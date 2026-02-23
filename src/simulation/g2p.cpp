// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"
#include <omp.h>


void Simulation::G2P(){

    unsigned int plastic_count = 0;

    #ifdef WARNINGS
        debug("G2P");
    #endif

    std::fill( particles.pic.begin(),  particles.pic.end(),  TV::Zero() );
    std::fill( particles.flip.begin(), particles.flip.end(), TV::Zero() );
    std::fill( particles.Bmat.begin(), particles.Bmat.end(), TM::Zero() );

    #pragma omp parallel num_threads(n_threads)
    {

        #pragma omp for nowait
        for(int p = 0; p < Np; p++){
            TV xp = particles.x[p];
            TV vp    = TV::Zero();
            TV flipp = TV::Zero();
            TM Bp    = TM::Zero();
            TM v_grad = TM::Zero();
            unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
            unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1;
        #ifdef THREEDIM
            unsigned int k_base = std::floor((xp(2)-grid.zc)*one_over_dx) - 1;
        #endif

            for(int i = i_base; i < i_base+4; i++){
                T xi = grid.x[i];
                T wi = N((xp(0)-xi)*one_over_dx);
                T wi_grad = dNdu((xp(0) - xi) * one_over_dx)  * one_over_dx;
                for(int j = j_base; j < j_base+4; j++){
                    T yi = grid.y[j];
                    T wj = N((xp(1) - yi)*one_over_dx);
                    T wj_grad = dNdu((xp(1) - yi) * one_over_dx)  * one_over_dx;
        #ifdef THREEDIM
                    for(int k = k_base; k < k_base+4; k++){
                        T zi = grid.z[k];
                        T wk = N((xp(2) - zi)*one_over_dx);
                        T wk_grad = dNdu((xp(2) - zi) * one_over_dx)  * one_over_dx;
                        
                        T weight = wi * wj * wk; //wip(xp(0), xp(1), xp(2), xi, yi, zi, one_over_dx);
                        TV weight_grad; 
                        weight_grad << wi_grad*wj*wk,
                                       wi*wj_grad*wk,
                                       wi*wj*wk_grad;

                        vp += grid.v[ind(i,j,k)] * weight;
                        v_grad += grid.v[ind(i,j,k)] * weight_grad.transpose();
                        if (flip_ratio < 0){ // APIC
                            TV posdiffvec = TV::Zero();
                            posdiffvec(0) = xi-xp(0);
                            posdiffvec(1) = yi-xp(1);
                            posdiffvec(2) = zi-xp(2);
                            Bp += grid.v[ind(i,j,k)] * posdiffvec.transpose() * weight;
                        }
                        if (flip_ratio >= -1){ // PIC-FLIP or AFLIP
                            flipp += grid.flip[ind(i,j,k)] * weight;
                        }
                    } // end loop k
        #else
                    T weight = wi * wj; //wip(xp(0), xp(1), xi, yi, one_over_dx);
                    TV weight_grad;
                    weight_grad << wi_grad*wj,
                                   wi*wj_grad;

                    vp += grid.v[ind(i,j)] * weight;
                    v_grad += grid.v[ind(i,j)] * weight_grad.transpose();
                    if (flip_ratio < 0){ // APIC
                        TV posdiffvec = TV::Zero();
                        posdiffvec(0) = xi-xp(0);
                        posdiffvec(1) = yi-xp(1);
                        Bp += grid.v[ind(i,j)] * posdiffvec.transpose() * weight;
                    }
                    if (flip_ratio >= -1){ // PIC-FLIP or AFLIP
                        flipp += grid.flip[ind(i,j)] * weight;
                    }
        #endif
                } // end loop j
            } // end loop i
            particles.pic[p] = vp;
            if (flip_ratio < 0){ // APIC
                particles.Bmat[p] = Bp;
            }
            if (flip_ratio >= -1){ // PIC-FLIP or AFLIP
                particles.flip[p] = flipp;
            }

            // Update material
            TM Fe_trial = particles.F[p];
            Fe_trial = (TM::Identity() + dt * v_grad )* Fe_trial;
            particles.F[p] = Fe_trial;

            plasticity(p, plastic_count, Fe_trial);


        } // end loop p

    } // end omp paralell

} // end G2P
