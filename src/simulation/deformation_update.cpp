// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"
#include <omp.h>

// Deformation gradient is updated based on the NEW GRID VELOCITIES and the OLD PARTICLE POSITIONS
void Simulation::deformationUpdate(){

    #ifdef WARNINGS
        debug("deformationUpdate");
    #endif

    std::fill( particles.delta_gamma.begin(), particles.delta_gamma.end(), 0.0 );
    unsigned int plastic_count = 0;

    #pragma omp parallel for reduction(+:plastic_count) num_threads(n_threads)
    for(int p=0; p<Np; p++){

        TM v_grad = TM::Zero();
        TV xp = particles.x[p];
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
                    T weight = wi * wj * wk;
                    TV weight_grad; 
                    weight_grad << wi_grad*wj*wk,
                                    wi*wj_grad*wk,
                                    wi*wj*wk_grad;
                    v_grad += grid.v[ind(i,j,k)] * weight_grad.transpose();
                } // end loop k
    #else
                T weight = wi * wj;
                TV weight_grad;
                weight_grad << wi_grad*wj,
                               wi*wj_grad;
                v_grad += grid.v[ind(i,j)] * weight_grad.transpose();
    #endif
            } // end loop i
        } // end loop j

        TM Fe_trial = particles.F[p];
        Fe_trial = Fe_trial + dt * v_grad * Fe_trial;
        particles.F[p] = Fe_trial;

        plasticity(p, plastic_count, Fe_trial);

    } // end loop over particles

    if (!reduce_verbose)
        debug("               Proj: ", plastic_count, " / ", Np);

} // end deformationUpdate
