#include "simulation.hpp"
#include <omp.h>

void Simulation::P2G_nonlocal(){

    grid.delta_gamma.resize(Nx*Ny); std::fill( grid.delta_gamma.begin(), grid.delta_gamma.end(), 0.0 );

    #pragma omp parallel num_threads(n_threads)
    {
        std::vector<T> grid_delta_gamma_local(Nx*Ny);

        #pragma omp for
        for(int p = 0; p < Np; p++){
            TV xp = particles.x[p];
            unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - 1; // the subtraction of one is valid for both quadratic and cubic splines
            unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - 1;
            unsigned int k_base = std::floor((xp(2)-grid.zc)*one_over_dx) - 1;

            for(int i = i_base; i < i_base+4; i++){
                T xi = grid.x[i];
                for(int j = j_base; j < j_base+4; j++){
                    T yi = grid.y[j];
                    for(int k = k_base; k < k_base+4; k++){
                        T zi = grid.z[j];
                        T weight = wip(xp(0), xp(1), xp(2), xi, yi, zi, one_over_dx);
                        if (weight > 1e-25){
                            grid_delta_gamma_local[ind(i,j,k)] += particles.delta_gamma[p] * weight;
                        } // end if
                    } // end for k
                } // end for j
            } // end for i
        } // end for p


        #pragma omp critical
        {
            for (int l = 0; l<Nx*Ny*Nz; l++){
                grid.delta_gamma[l] += grid_delta_gamma_local[l];
            } // end for l
        } // end omp critical

    } // end omp parallel

    for (int l = 0; l<Nx*Ny*Ny; l++){
        T mi = grid.mass[l];
        if (mi > 0)
            grid.delta_gamma[l] /= (mi / particle_mass);
        else
            grid.delta_gamma[l] = 0;
        //grid.v[l] = (mi > 0) ? grid.v[l]/mi : TV::Zero(); // condition ? result_if_true : result_if_false
    }

} // end P2G_Optimized_Parallel
