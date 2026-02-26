// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"
#include <omp.h>


void Simulation::explicitEulerUpdate(){

    #ifdef WARNINGS
        debug("explicitEulerUpdate");
    #endif

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
