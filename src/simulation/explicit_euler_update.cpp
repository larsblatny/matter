// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"
#include <omp.h>


void Simulation::explicitEulerUpdate(){

    #ifdef WARNINGS
        debug("explicitEulerUpdate");
    #endif

    // std::vector<TV> grid_force(grid_nodes, TV::Zero());

    // #pragma omp parallel num_threads(n_threads)
    // {
    //     std::vector<TV> grid_force_local(grid_nodes, TV::Zero());

    //     #pragma omp for nowait
    //     for(int p = 0; p < Np; p++){

    //         TM Fe = particles.F[p];

    //         TM dPsidF;
    //         if (elastic_model == ElasticModel::NeoHookean){
    //             dPsidF = NeoHookeanPiola(Fe);
    //         }
    //         else if (elastic_model == ElasticModel::Hencky){ // St Venant Kirchhoff with Hencky strain
    //             dPsidF = HenckyPiola(Fe);
    //         }
    //         else{
    //             debug("You specified an unvalid ELASTIC model!");
    //         }

    //         TM tau = dPsidF * Fe.transpose();

    //         TV xp = particles.x[p];
    //         unsigned int i_base = std::max(0, int(std::floor((xp(0)-grid.xc)*one_over_dx)) - 1); // i_base = std::min(i_base, Nx-4); // the subtraction of one is valid for both quadratic and cubic splines
    //         unsigned int j_base = std::max(0, int(std::floor((xp(1)-grid.yc)*one_over_dx)) - 1); // j_base = std::min(j_base, Ny-4);
    //     #ifdef THREEDIM
    //         unsigned int k_base = std::max(0, int(std::floor((xp(2)-grid.zc)*one_over_dx)) - 1); // k_base = std::min(k_base, Nz-4);
    //     #endif

    //         for(int i = i_base; i < i_base+4; i++){
    //             T xi = grid.x[i];
    //             for(int j = j_base; j < j_base+4; j++){
    //                 T yi = grid.y[j];
    //     #ifdef THREEDIM
    //                 for(int k = k_base; k < k_base+4; k++){
    //                     T zi = grid.z[k];
    //                     if ( grid.mass[ind(i,j,k)] > 0){
    //                         grid_force_local[ind(i,j,k)] += tau * grad_wip(xp(0), xp(1), xp(2), xi, yi, zi, one_over_dx);
    //                     } // end if non-zero grid mass
    //                 } // end for k
    //     #else
    //                 if ( grid.mass[ind(i,j)] > 0){
    //                     grid_force_local[ind(i,j)] += tau * grad_wip(xp(0), xp(1), xi, yi, one_over_dx);
    //                 } // end if non-zero grid mass
    //     #endif
    //              } // end for j
    //          } // end for i
    //     } // end for particles

    //     #pragma omp critical
    //     {
    //         for (int l = 0; l<grid_nodes; l++){
    //             grid_force[l] += grid_force_local[l];
    //         } // end for l
    //     } // end omp critical

    // } // end omp parallel

    //////////// if external grid gravity: //////////////////
    // std::pair<TMX, TMX> external_gravity_pair = createExternalGridGravity();
    ////////////////////////////////////////////////////////

    T dt_particle_volume = dt * particle_volume;
    TV dt_gravity = dt * gravity;

    #pragma omp parallel for collapse(DIMENSION) num_threads(n_threads)
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
#ifdef THREEDIM
            for(int k = 0; k < Nz; k++){
                unsigned int index = ind(i,j,k);
#else
                unsigned int index = ind(i,j);
#endif
                T mi = grid.mass[index];
                if (mi > 0){

                    TV velocity_increment = -dt_particle_volume * grid.force[index] / mi + dt_gravity;

                    //////////// if external grid gravity: //////////////////
                    // T external_gravity = external_gravity_pair.first(i,j);
                    // T external_gravity = external_gravity_pair.second(i,j);
                    // velocity_increment_x += dt * external_gravity(0);
                    // velocity_increment_y += dt * external_gravity(1);
                    ////////////////////////////////////////////////////////

                    TV old_vi = grid.v[index];
                    TV new_vi = old_vi + velocity_increment;

#ifdef THREEDIM
                    TV Xi(grid.x[i], grid.y[j], grid.z[k]);
#else
                    TV Xi(grid.x[i], grid.y[j]);
#endif

                    boundaryCollision(index, Xi, new_vi);

                    // TODO:
                    // boundaryCorrection(new_xi, new_yi, new_vxi, new_vyi);

                    // Only if impose velocity on certain grid nodes:
                    // overwriteGridVelocity(grid.x[i], grid.y[j], new_vi);

                    grid.v[index]    = new_vi;
                    grid.flip[index] = new_vi - old_vi;

                } // end if non-zero grid mass
#ifdef THREEDIM
            } // end for k
#endif
        } // end for j
    } // end for i

} // end explicitEulerUpdate
