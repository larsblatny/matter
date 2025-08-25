// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"
#include <omp.h>

void Simulation::MUSL(){

    #ifdef WARNINGS
        debug("MUSL");
    #endif


    #pragma omp parallel for schedule(static) num_threads(n_threads)
    for(int p=0; p<Np; p++){

        //// Velicity is updated
        if (flip_ratio < -1){ // APIC
            particles.v[p] = particles.pic[p];
        } else if (flip_ratio < 0){ // AFLIP
            particles.v[p] = (-flip_ratio) * ( particles.v[p] + particles.flip[p] ) + (1 - (-flip_ratio)) * particles.pic[p];
        } else{ // PIC-FLIP
            particles.v[p] =   flip_ratio  * ( particles.v[p] + particles.flip[p] ) + (1 -   flip_ratio)  * particles.pic[p];
        }

    } // end loop over particles


    std::fill( grid.v.begin(), grid.v.end(), TV::Zero() );

    #pragma omp parallel num_threads(n_threads)
    {
        std::vector<TV> grid_v_local(grid_nodes, TV::Zero() );

        #pragma omp for nowait
        for(int p = 0; p < Np; p++){
            TV xp = particles.x[p];
            unsigned int i_base = std::max(0, int(std::floor((xp(0)-grid.xc)*one_over_dx)) - 1); // i_base = std::min(i_base, Nx-4); // the subtraction of one is valid for both quadratic and cubic splines
            unsigned int j_base = std::max(0, int(std::floor((xp(1)-grid.yc)*one_over_dx)) - 1); // j_base = std::min(j_base, Ny-4);
        #ifdef THREEDIM
            unsigned int k_base = std::max(0, int(std::floor((xp(2)-grid.zc)*one_over_dx)) - 1); // k_base = std::min(k_base, Nz-4);
        #endif

            for(int i = i_base; i < i_base+4; i++){
                T xi = grid.x[i];
                for(int j = j_base; j < j_base+4; j++){
                    T yi = grid.y[j];
        #ifdef THREEDIM
                    for(int k = k_base; k < k_base+4; k++){
                        T zi = grid.z[k];
                        T weight = wip(xp(0), xp(1), xp(2), xi, yi, zi, one_over_dx);
                        if (weight > 1e-25){
                            grid_v_local[ind(i,j,k)] += particles.v[p] * weight;
                            if (flip_ratio < 0){ // APIC
                                TV posdiffvec = TV::Zero();
                                posdiffvec(0) = xi-xp(0);
                                posdiffvec(1) = yi-xp(1);
                                posdiffvec(2) = zi-xp(2);
                                grid_v_local[ind(i,j,k)] += particles.Bmat[p] * posdiffvec * apicDinverse * weight;
                            }
                        }
                    } // end for k
        #else
                    T weight = wip(xp(0), xp(1), xi, yi, one_over_dx);
                    if (weight > 1e-25){
                        grid_v_local[ind(i,j)]     += particles.v[p] * weight;
                        if (flip_ratio < 0){ // APIC
                            TV posdiffvec = TV::Zero();
                            posdiffvec(0) = xi-xp(0);
                            posdiffvec(1) = yi-xp(1);
                            grid_v_local[ind(i,j)] += particles.Bmat[p] * posdiffvec * apicDinverse * weight;
                        }
                    }
        #endif
                } // end for j
            } // end for i
        } // end for p

        #pragma omp critical
        {
            for (int l = 0; l<grid_nodes; l++){
                grid.v[l] += grid_v_local[l];
            } // end for l
        } // end omp critical

    } // end omp parallel


    for (int l = 0; l<grid_nodes; l++){
        T mi = grid.mass[l];
        if (mi > 0)
            grid.v[l] /= (mi/particle_mass);
        else
            grid.v[l].setZero();
    }

} // end P2G_MUSL
