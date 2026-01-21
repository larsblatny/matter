#include "simulation.hpp"
#include <omp.h>


void Simulation::nonlocalG2P(){

    std::fill( particles.delta_gamma_nonloc.begin(), particles.delta_gamma_nonloc.end(), 0.0 );
    std::vector<T> divisor; divisor.resize(Np); std::fill( divisor.begin(), divisor.end(), 0.0 );

    #pragma omp parallel num_threads(n_threads)
    {
        #pragma omp for
        for(int p = 0; p < Np; p++){
            TV xp = particles.x[p];
            unsigned int i_base = std::floor((xp(0)-grid.xc)*one_over_dx) - (nonlocal_support-1); // the subtraction of one is valid for both quadratic and cubic splines
            unsigned int j_base = std::floor((xp(1)-grid.yc)*one_over_dx) - (nonlocal_support-1);

            for(int i = i_base; i < i_base+(2*nonlocal_support); i++){
                T xi = grid.x[i];
                for(int j = j_base; j < j_base+(2*nonlocal_support); j++){
                    T yi = grid.y[j];
    #ifdef THREEDIM
                    for(int k = k_base; k < k_base+4; k++){
                        T zi = grid.z[k];
                        T dist_norm_sq = (xp(0)-xi)*(xp(0)-xi) + (xp(1)-yi)*(xp(1)-yi) + (xp(2)-zi)*(xp(2)-zi);
                        if (dist_norm_sq <= nonlocal_l_sq){
                            T kernel = std::exp(-4 * dist_norm_sq / nonlocal_l_sq);
                            T mi = grid.mass[ind(i,j,k)];
                            particles.delta_gamma_nonloc[p] += grid.delta_gamma[ind(i,j,k)] * kernel * mi;
                            divisor[p] += kernel * mi;
                        }
                    } // end for k
    #else    
                    T dist_norm_sq = (xp(0)-xi)*(xp(0)-xi) + (xp(1)-yi)*(xp(1)-yi);
                    if (dist_norm_sq <= nonlocal_l_sq){
                        T kernel = std::exp(-4 * dist_norm_sq / nonlocal_l_sq);
                        T mi = grid.mass[ind(i,j)];
                        particles.delta_gamma_nonloc[p] += grid.delta_gamma[ind(i,j)] * kernel * mi;
                        divisor[p] += kernel * mi;
                    }
    #endif
                } // end loop j
            } // end loop i
        } // end loop p
    } // end omp paralell

    for (int p = 0; p<Np; p++){
        T divisor_p = divisor[p];
        if (divisor_p > 0)
            particles.delta_gamma_nonloc[p] /= divisor_p;
        else
            particles.delta_gamma_nonloc[p] = 0;
    }

} // end nonlocalG2P
