// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"
#include <omp.h>

void Simulation::P2G(){

    #ifdef WARNINGS
        debug("P2G");
    #endif

    #pragma omp parallel num_threads(n_threads)
    {
        std::vector<TV> grid_v_local(grid_nodes, TV::Zero() );
        std::vector<T> grid_mass_local(grid_nodes);
        std::vector<T> grid_friction_local(grid_nodes);
        std::vector<TV> grid_force_local(grid_nodes, TV::Zero() );

        #pragma omp for nowait
        for(int p = 0; p < Np; p++){
            TV xp = particles.x[p];
            unsigned int i_base = std::max(0, int(std::floor((xp(0)-grid.xc)*one_over_dx)) - 1); // i_base = std::min(i_base, Nx-4); // the subtraction of one is valid for both quadratic and cubic splines
            unsigned int j_base = std::max(0, int(std::floor((xp(1)-grid.yc)*one_over_dx)) - 1); // j_base = std::min(j_base, Ny-4);
        #ifdef THREEDIM
            unsigned int k_base = std::max(0, int(std::floor((xp(2)-grid.zc)*one_over_dx)) - 1); // k_base = std::min(k_base, Nz-4);
        #endif

            // ALT 1: Compute stress
            // TM Fe = particles.F[p];
            // TM dPsidF;
            // if (elastic_model == ElasticModel::NeoHookean){
            //     dPsidF = NeoHookeanPiola(Fe);
            // }
            // else if (elastic_model == ElasticModel::Hencky){ // St Venant Kirchhoff with Hencky strain
            //     dPsidF = HenckyPiola(Fe);
            // }
            // else{
            //     debug("You specified an unvalid ELASTIC model!");
            // }
            // TM tau = dPsidF * Fe.transpose();

            // ALT 2: Use the saved stress
            TM tau = particles.tau[p];

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

                        if (weight > 1e-25){
                            grid_mass_local[ind(i,j,k)]  += weight;
                            grid_v_local[ind(i,j,k)]     += particles.v[p] * weight;
                            grid_force_local[ind(i,j,k)] += tau * weight_grad;
                            if (flip_ratio < 0){ // APIC
                                TV posdiffvec = TV::Zero();
                                posdiffvec(0) = xi-xp(0);
                                posdiffvec(1) = yi-xp(1);
                                posdiffvec(2) = zi-xp(2);
                                grid_v_local[ind(i,j,k)] += particles.Bmat[p] * posdiffvec * apicDinverse * weight;
                            }
                            if (use_mibf)
                                grid_friction_local[ind(i,j,k)] += particles.muI[p] * weight;
                        }
                    } // end for k
        #else
                    T weight = wi * wj;
                    TV weight_grad;
                    weight_grad << wi_grad*wj,
                                   wi*wj_grad;

                    if (weight > 1e-25){
                        grid_mass_local[ind(i,j)]  += weight;
                        grid_v_local[ind(i,j)]     += particles.v[p] * weight;
                        grid_force_local[ind(i,j)] += tau * weight_grad;
                        if (flip_ratio < 0){ // APIC
                            TV posdiffvec = TV::Zero();
                            posdiffvec(0) = xi-xp(0);
                            posdiffvec(1) = yi-xp(1);
                            grid_v_local[ind(i,j)] += particles.Bmat[p] * posdiffvec * apicDinverse * weight;
                        }
                        if (use_mibf)
                            grid_friction_local[ind(i,j)] += particles.muI[p] * weight;
                    }
        #endif
                } // end for j
            } // end for i
        } // end for p

        #pragma omp critical
        {
            for (int l = 0; l<grid_nodes; l++){
                grid.mass[l]          += grid_mass_local[l];
                grid.v[l]             += grid_v_local[l];
                grid.force[l]         += grid_force_local[l];
                if (use_mibf)
                    grid.friction[l]  += grid_friction_local[l];
            } // end for l
        } // end omp critical

    } // end omp parallel

    ///////////////////////////////////////////////////////////
    // At this point in time grid.mass is equal to m_i / m_p //
    ///////////////////////////////////////////////////////////

    #pragma omp parallel for num_threads(n_threads)
    for (int l = 0; l<grid_nodes; l++){
        T mi = grid.mass[l];
        if (mi > 0)
            grid.v[l] /= mi;
        else
            grid.v[l].setZero();
        //grid.v[l] = (mi > 0) ? grid.v[l]/mi : TV::Zero(); // condition ? result_if_true : result_if_false
        if (use_mibf)
            grid.friction[l] /= mi;
    }

    for(auto&& m: grid.mass){
        m *= particle_mass;
    }

} // end P2G


void Simulation::P2G_sparse_scan() {
    #pragma omp parallel num_threads(n_threads)
    {
        // Thread local grid quantities
        std::vector<TV> grid_v_local(scan_sparse_grid.num_active_nodes(), TV::Zero() );
        std::vector<T> grid_mass_local(scan_sparse_grid.num_active_nodes(), T(0));
        std::vector<T> grid_friction_local(scan_sparse_grid.num_active_nodes(), T(0));
        std::vector<TV> grid_force_local(scan_sparse_grid.num_active_nodes(), TV::Zero() );

        // Parallel particles loop
        #pragma omp for nowait
        for(int p = 0; p < Np; p++){
            TV xp = particles.x[p];
            int i_base = std::floor(xp(0)*one_over_dx - 0.5); 
            int j_base = std::floor(xp(1)*one_over_dx - 0.5); 
            int k_base = std::floor(xp(2)*one_over_dx - 0.5);

            // Kirchhoff stress
            TM tau = particles.tau[p];

            // Grid loop for each particle
            for(int ix = i_base+0; ix<i_base+3; ix++){
                T xi = ix * dx;
                T wi = N((xp(0)-xi)*one_over_dx);
                T wi_grad = dNdu((xp(0) - xi) * one_over_dx)  * one_over_dx;
                for (int iy=j_base+0; iy<j_base+3; iy++) {
                    T yj = iy * dx;
                    T wj = N((xp(1) - yj)*one_over_dx);
                    T wj_grad = dNdu((xp(1) - yj) * one_over_dx)  * one_over_dx;
                    for (int iz=k_base+0; iz<k_base+3; iz++) {
                        T zk = iz * dx;
                        T wk = N((xp(2) - zk)*one_over_dx);
                        T wk_grad = dNdu((xp(2) - zk) * one_over_dx)  * one_over_dx;

                        // weight and gradient
                        T weight = wi * wj * wk; 
                        TV weight_grad; 
                        weight_grad << wi_grad*wj*wk,
                                       wi*wj_grad*wk,
                                       wi*wj*wk_grad;

                        // node (i, j, k) -> node flattened id
                        int gid = scan_sparse_grid.ijk_to_gid(ix, iy, iz);
                        assert(gid>=0);

                        // P2G
                        grid_mass_local[gid]  += weight; //weight * particle_mass;
                        grid_v_local[gid]     += particles.v[p] * weight; //particles.v[p] * weight * particle_mass;
                        grid_force_local[gid] += tau * weight_grad;
                        if (flip_ratio < 0) { // APIC
                            TV posdiffvec = TV::Zero();
                            posdiffvec(0) = xi-xp(0);
                            posdiffvec(1) = yj-xp(1);
                            posdiffvec(2) = zk-xp(2);
                            grid_v_local[gid] += particles.Bmat[p] * posdiffvec * apicDinverse * weight; //particles.Bmat[p] * posdiffvec * apicDinverse * weight * particle_mass;
                        }
                        if (use_mibf) {
                            grid_friction_local[gid] += particles.muI[p] * weight; //particles.muI[p] * weight * particle_mass;
                        }
                    } // end z
                } // end y
            } // end x
        } // end p

        // Local to global
        #pragma omp critical
        {
            for (int l = 0; l<scan_sparse_grid.num_active_nodes(); l++){
                grid.mass[l]          += grid_mass_local[l];
                grid.v[l]             += grid_v_local[l];
                grid.force[l]         += grid_force_local[l];
                if (use_mibf)
                    grid.friction[l]  += grid_friction_local[l];
            } // end for l
        } // end omp critical
    } // end opm parallel
}