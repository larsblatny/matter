// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"
#include <omp.h>


void Simulation::G2P(){

    #ifdef WARNINGS
        debug("G2P");
    #endif

    std::fill( particles.pic.begin(),  particles.pic.end(),  TV::Zero() );
    std::fill( particles.flip.begin(), particles.flip.end(), TV::Zero() );
    std::fill( particles.Bmat.begin(), particles.Bmat.end(), TM::Zero() );

    #pragma omp parallel num_threads(n_threads)
    {
        std::vector<TV> particles_pic_local(Np);  std::fill( particles_pic_local.begin(),  particles_pic_local.end(),  TV::Zero() );
        std::vector<TV> particles_flip_local(Np); std::fill( particles_flip_local.begin(), particles_flip_local.end(), TV::Zero() );
        std::vector<TM> particles_Bmat_local(Np); std::fill( particles_Bmat_local.begin(), particles_Bmat_local.end(), TM::Zero() );

        #pragma omp for nowait
        for(int p = 0; p < Np; p++){
            TV xp = particles.x[p];
            TV vp    = TV::Zero();
            TV flipp = TV::Zero();
            TM Bp    = TM::Zero();
            unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
            unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1;
        #ifdef THREEDIM
            unsigned int k_base = std::floor((xp(2)-grid.zc)*one_over_dx) - 1;
        #endif

            for(int i = i_base; i < i_base+4; i++){
                T xi = grid.x[i];
                for(int j = j_base; j < j_base+4; j++){
                    T yi = grid.y[j];
        #ifdef THREEDIM
                    for(int k = k_base; k < k_base+4; k++){
                        T zi = grid.z[k];
                        T weight = wip(xp(0), xp(1), xp(2), xi, yi, zi, one_over_dx);
                        vp += grid.v[ind(i,j,k)] * weight;
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
                    T weight = wip(xp(0), xp(1), xi, yi, one_over_dx);
                    vp += grid.v[ind(i,j)] * weight;
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
            particles_pic_local[p] = vp;
            if (flip_ratio < 0){ // APIC
                particles_Bmat_local[p] = Bp;
            }
            if (flip_ratio >= -1){ // PIC-FLIP or AFLIP
                particles_flip_local[p] = flipp;
            }
        } // end loop p

        #pragma omp critical
        {
            for(int p = 0; p < Np; p++){
                particles.pic[p] += particles_pic_local[p];
                if (flip_ratio < 0){ // APIC
                    particles.Bmat[p] += particles_Bmat_local[p];
                }
                if (flip_ratio >= -1){ // PIC-FLIP or AFLIP
                    particles.flip[p] += particles_flip_local[p];
                }
            }
        } // end omp critical

    } // end omp paralell

} // end G2P
