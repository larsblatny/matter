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

    #pragma omp parallel num_threads(n_threads) reduction(+:plastic_count)
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
                        
                        T weight = wi * wj * wk;
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
                    T weight = wi * wj;
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

            if (!use_musl){
                TM Fe_trial = particles.F[p];
                Fe_trial = Fe_trial + dt * v_grad * Fe_trial;
                particles.F[p] = Fe_trial;

                plasticity(p, plastic_count, Fe_trial);
            }

        } // end loop p

    } // end omp paralell

    if (!reduce_verbose && !use_musl)
        debug("               Proj: ", plastic_count, " / ", Np);

} // end G2P



void Simulation::G2PSparseScan() {
    unsigned int plastic_count = 0;

    const int node_lo_offset = 0;
    const int node_hi_offset = 2;

    // min/max nodal id
    std::vector<int> tminx(n_threads, INT32_MAX);
    std::vector<int> tminy(n_threads, INT32_MAX);
    std::vector<int> tminz(n_threads, INT32_MAX);
    std::vector<int> tmaxx(n_threads, -INT32_MAX);
    std::vector<int> tmaxy(n_threads, -INT32_MAX);
    std::vector<int> tmaxz(n_threads, -INT32_MAX);

    #pragma omp parallel num_threads(n_threads) reduction(+:plastic_count)
    {
        int tid = omp_get_thread_num();
        int lminx = INT32_MAX;
        int lminy = INT32_MAX;
        int lminz = INT32_MAX;
        int lmaxx = -INT32_MAX;
        int lmaxy = -INT32_MAX;
        int lmaxz = -INT32_MAX;

        // Parallel particles loop
        #pragma omp for nowait
        for (int p=0; p<Np; p++) {
            TV xp = particles.x[p];
            TV vp    = TV::Zero();
            TV flipp = TV::Zero();
            TM Bp    = TM::Zero();
            TM v_grad = TM::Zero();
            int i_base = std::floor(xp(0)*one_over_dx - 0.5); 
            int j_base = std::floor(xp(1)*one_over_dx - 0.5); 
            int k_base = std::floor(xp(2)*one_over_dx - 0.5);

            // Get local min/max base
            int ilo = i_base + node_lo_offset;
            int ihi = i_base + node_hi_offset;
            int jlo = j_base + node_lo_offset;
            int jhi = j_base + node_hi_offset;
            int klo = k_base + node_lo_offset;
            int khi = k_base + node_hi_offset;

            lminx = std::min(lminx, ilo);
            lminy = std::min(lminy, jlo);
            lminz = std::min(lminz, klo);
            lmaxx = std::max(lmaxx, ihi);
            lmaxy = std::max(lmaxy, jhi);
            lmaxz = std::max(lmaxz, khi);

            // Local grid nodes loop
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

                        // node (i, j, k) -> flattened id
                        int gid = scan_sparse_grid.ijk_to_gid(ix, iy, iz);
                        assert(gid>=0);

                        // G2P
                        vp += grid.v[gid] * weight;
                        v_grad += grid.v[gid] * weight_grad.transpose();
                        if (flip_ratio < 0){ // APIC
                            TV posdiffvec = TV::Zero();
                            posdiffvec(0) = xi-xp(0);
                            posdiffvec(1) = yj-xp(1);
                            posdiffvec(2) = zk-xp(2);
                            Bp += grid.v[gid] * posdiffvec.transpose() * weight;
                        }
                        if (flip_ratio >= -1){ // PIC-FLIP or AFLIP
                            flipp += grid.flip[gid] * weight;
                        }
                    } // end k
                } // end j
            } // end i

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

            particles.x[p] = particles.x[p] + dt * particles.pic[p];
            // Velicity is updated
            if (flip_ratio < -1){ // APIC
                particles.v[p] = particles.pic[p];
            } else if (flip_ratio < 0){ // AFLIP
                particles.v[p] = (-flip_ratio) * ( particles.v[p] + particles.flip[p] ) + (1 - (-flip_ratio)) * particles.pic[p];
            } else{ // PIC-FLIP
                particles.v[p] =   flip_ratio  * ( particles.v[p] + particles.flip[p] ) + (1 -   flip_ratio)  * particles.pic[p];
            }
        } // end p

        tminx[tid] = lminx;
        tminy[tid] = lminy;
        tminz[tid] = lminz;
        tmaxx[tid] = lmaxx;
        tmaxy[tid] = lmaxy;
        tmaxz[tid] = lmaxz;
    } // end omp pragma

    particles.minx_id = *std::min_element(tminx.begin(), tminx.end());
    particles.miny_id = *std::min_element(tminy.begin(), tminy.end());
    particles.minz_id = *std::min_element(tminz.begin(), tminz.end());
    particles.maxx_id = *std::max_element(tmaxx.begin(), tmaxx.end());
    particles.maxy_id = *std::max_element(tmaxy.begin(), tmaxy.end());
    particles.maxz_id = *std::max_element(tmaxz.begin(), tmaxz.end());
}
