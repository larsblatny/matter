#include "simulation.hpp"
#include <omp.h>

void Simulation::nonlocalProjection(){

    nonlocalP2G();
    nonlocalG2P();

    int plastic_count = 0;

    #pragma omp parallel for reduction(+:plastic_count) num_threads(n_threads)
    for(int p=0; p<Np; p++){

        // if local approach
        // T delta_gamma_nonloc = particles.delta_gamma[p] * dt / d_prefac; // this is now eps_pl_dev_instant

        // if nonlocal approach
        T delta_gamma_nonloc = particles.delta_gamma_nonloc[p] * dt / d_prefac; // this is now eps_pl_dev_instant

        if (delta_gamma_nonloc > 0){ // If delta_gamma_nonloc is > 0, then delta_gamma > 0 also
            plastic_count++;

            T delta_gamma = particles.delta_gamma[p] * dt / d_prefac; // this is now eps_pl_dev_instant

            if (delta_gamma > 0){ // if originally PLASTIC particle, we perform the projection using the NONLOCAL delta gamma

                Eigen::JacobiSVD<TM> svd(particles.Fe_trial[p], Eigen::ComputeFullU | Eigen::ComputeFullV);
                TV hencky = svd.singularValues().array().abs().max(1e-4).log();
                T  hencky_trace = hencky.sum();
                TV hencky_deviatoric = hencky - (hencky_trace / dim) * TV::Ones();
                T  hencky_deviatoric_norm = hencky_deviatoric.norm();

                if (hencky_deviatoric_norm > 0)
                    hencky_deviatoric /= hencky_deviatoric_norm; // normalize the deviatoric vector so it gives a unit vector specifying the deviatoric direction

                hencky -= delta_gamma_nonloc * hencky_deviatoric;
                // hencky -= delta_gamma * hencky_deviatoric;
                particles.F[p] = svd.matrixU() * hencky.array().exp().matrix().asDiagonal() * svd.matrixV().transpose();
                particles.eps_pl_dev[p] += delta_gamma;

            }
            // In either case of an elastic or plastic particle, we record the plastic strain. This will change the yield stress in the next time step.
            particles.eps_pl_dev_nonloc[p] += delta_gamma_nonloc;
        }

    } // end for loop

    debug("               projected nonlocal = ", plastic_count, " / ", Np);

} // end nonlocalProjection
