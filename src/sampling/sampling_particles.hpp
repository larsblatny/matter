// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef SAMPLING_PARTICLES_HPP
#define SAMPLING_PARTICLES_HPP

#include "../tools.hpp"
#include "../data_structures.hpp"
#include "../../deps/tph_poisson-0.3/thinks/poisson_disk_sampling/poisson_disk_sampling.h"

#ifdef THREEDIM

    template <typename S>
    void sampleParticles(S& sim, T kRadius, T ppc = 8, unsigned int crop_to_shape = 0, std::uint32_t attempts = 30, std::uint32_t seed = 42) {
        const T Lx = sim.Lx;
        const T Ly = sim.Ly;
        const T Lz = sim.Lz;
        std::uint32_t kAttempts = attempts;
        std::uint32_t kSeed = seed;
        std::array<T, 3> kXMin = std::array<T, 3>{{0, 0, 0}};
        std::array<T, 3> kXMax = std::array<T, 3>{{Lx, Ly, Lz}};

        debug("Sampling particles...");
        std::vector<std::array<T, 3>> square_samples = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax, kAttempts, kSeed);
        std::vector<std::array<T, 3>> samples;

        debug("    Number of square samples: ", square_samples.size());

        sim.dx = std::cbrt(ppc / T(square_samples.size()) * Lx*Ly*Lz);
        sim.particle_volume = sim.dx * sim.dx * sim.dx / ppc; // = Lx*Ly*Lz / T(square_samples.size())
        sim.particle_mass = sim.rho * sim.particle_volume;

        debug("    dx set to ", sim.dx);

         /////// Cylinder
        if (crop_to_shape == 1){
            for(int p = 0; p < square_samples.size(); p++){
                if ( (square_samples[p][0]-Lx/2.0)*(square_samples[p][0]-Lx/2.0) + (square_samples[p][2]-Lz/2.0)*(square_samples[p][2]-Lz/2.0) < (Lx/2.0)*(Lx/2.0) ){
                    samples.push_back(square_samples[p]);
                }
            }
        }
        //////// Silo
        else if (crop_to_shape == 2){
            for(int p = 0; p < square_samples.size(); p++){
                T x = square_samples[p][0]-Lx/2.0;
                T y = square_samples[p][1];
                T z = square_samples[p][2]-Lz/2.0;
                T r_surface = std::tanh(y) + 1;
                T r_surface_sq = r_surface * r_surface;
                T r_point_sq = x*x + z*z;
                if (r_point_sq < r_surface_sq){
                    samples.push_back(square_samples[p]);
                }
            }
        }
        else{
            debug("    No shape specified, using just a square.");
            samples = square_samples;
        }

        sim.Np = samples.size();
        debug("    Number of particles samples: ", sim.Np);

        sim.particles = Particles(sim.Np);
        for(int p = 0; p < sim.Np; p++){
            for(int d = 0; d < 3; d++){
                sim.particles.x[p](d) = samples[p][d];
            }
        }

    } // end sampleParticles

#else // TWODIM

    template <typename S>
    void sampleParticles(S& sim, T kRadius, T ppc = 6, unsigned int crop_to_shape = 0, std::uint32_t attempts = 200, std::uint32_t seed = 42){
        const T Lx = sim.Lx;
        const T Ly = sim.Ly;
        std::uint32_t kAttempts = attempts;
        std::uint32_t kSeed = seed;
        std::array<T, 2> kXMin = std::array<T, 2>{{0, 0}};
        std::array<T, 2> kXMax = std::array<T, 2>{{Lx, Ly}};

        debug("Sampling particles...");
        std::vector<std::array<T, 2>> square_samples = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax, kAttempts, kSeed);
        std::vector<std::array<T, 2>> samples;

        debug("    Number of square samples: ", square_samples.size());

        sim.dx = std::sqrt(ppc / T(square_samples.size()) * Lx*Ly);
        sim.particle_volume = sim.dx * sim.dx / ppc; // = Lx*Ly / T(square_samples.size())
        sim.particle_mass = sim.rho * sim.particle_volume;

        debug("    dx set to ", sim.dx);

        /////// Quadratic Gate
        if (crop_to_shape == 1){
            T height = 0.05; // 0.016;
            for(int p = 0; p < square_samples.size(); p++){
                T xp = square_samples[p][0];
                T y_gate = height + 100*(xp-Lx)*(xp-Lx) - 0.5*sim.dx;
                if (square_samples[p][1] < y_gate){
                    samples.push_back(square_samples[p]);
                }
            }
        }
        else{
            debug("    No shape specified, using just a square.");
            samples = square_samples;
        }

        sim.Np = samples.size();
        debug("    Number of particles samples: ", sim.Np);

        sim.particles = Particles(sim.Np);
        for(int p = 0; p < sim.Np; p++){
            for(int d = 0; d < 2; d++){
                sim.particles.x[p](d) = samples[p][d];
            }
        }


    } // end sampleParticles

#endif 

#ifdef THREEDIM 

    template <typename S>
    void sampleParticles(S& sim, std::vector<TV> origins, std::vector<TV> sizes, T kRadius, T ppc = 8, std::uint32_t attempts = 30, std::uint32_t seed = 42) {
        std::uint32_t kAttempts = attempts;
        std::uint32_t kSeed     = seed;

        debug("Sampling particles from ", origins.size(), " boxes...");

        std::vector<std::array<T, 3>> samples;

        sim.sampling_start_idx.push_back(0);

        T total_volume = 0;   
        T total_square_samples = 0;   

        for (int i = 0; i < origins.size(); ++i) {

            TV mn = origins[i];
            TV mx = origins[i] + sizes[i];

            std::array<T, 3> kXMin = std::array<T, 3>{{mn(0), mn(1), mn(2)}};
            std::array<T, 3> kXMax = std::array<T, 3>{{mx(0), mx(1), mx(2)}};
            std::vector<std::array<T, 3>> box_samples = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax, kAttempts, kSeed);

            total_volume += sizes[i](0) * sizes[i](1) * sizes[i](2);
            total_square_samples += box_samples.size();

            for (int p = 0; p < box_samples.size(); ++p)
                samples.push_back(box_samples[p]);

            sim.sampling_start_idx.push_back(samples.size());

            debug("    box ", i, ": ", sim.sampling_start_idx[i + 1] - sim.sampling_start_idx[i], " particles");
        }

        sim.dx = std::cbrt(ppc / T(total_square_samples) * total_volume);
        sim.particle_volume = sim.dx * sim.dx * sim.dx / ppc;
        sim.particle_mass   = sim.rho * sim.particle_volume;

        debug("    dx set to ", sim.dx);

        TV dmin = origins[0];
        TV dmax = origins[0] + sizes[0];

        for (int i = 1; i < origins.size(); ++i) {

            dmin = dmin.cwiseMin(origins[i]);
            dmax = dmax.cwiseMax(origins[i] + sizes[i]);

        }

        sim.Lx = dmax(0) - dmin(0);
        sim.Ly = dmax(1) - dmin(1);
        sim.Lz = dmax(2) - dmin(2);

        debug("Lx Ly Lz = ", sim.Lx, " ", sim.Ly, " ", sim.Lz);

        sim.Np = samples.size();
        debug("    Number of particles samples: ", sim.Np);

        sim.particles = Particles(sim.Np);
        for (int p = 0; p < sim.Np; p++) {
            for (int d = 0; d < 3; d++) {
                sim.particles.x[p](d) = samples[p][d];
            }
        }

    } // end sampleParticles (multiple objects)

#else // TWODIM

    template <typename S>
    void sampleParticles(S& sim, std::vector<TV> origins, std::vector<TV> sizes,
                         T kRadius, T ppc = 6, std::uint32_t attempts = 200, std::uint32_t seed = 42) {
        std::uint32_t kAttempts = attempts;
        std::uint32_t kSeed     = seed;

        debug("Sampling particles from ", origins.size(), " boxes...");

        std::vector<std::array<T, 2>> samples;

        sim.sampling_start_idx.push_back(0);

        T total_volume = 0;  
        T total_square_samples = 0;   


        for (int i = 0; i < origins.size(); ++i) {

            TV mn = origins[i];
            TV mx = origins[i] + sizes[i];

            std::array<T, 2> kXMin = std::array<T, 2>{{mn(0), mn(1)}};
            std::array<T, 2> kXMax = std::array<T, 2>{{mx(0), mx(1)}};
            std::vector<std::array<T, 2>> box_samples = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax, kAttempts, kSeed);

            total_volume += sizes[i](0) * sizes[i](1);
            total_square_samples += box_samples.size();

            for (int p = 0; p < box_samples.size(); ++p)
                samples.push_back(box_samples[p]);  

            sim.sampling_start_idx.push_back(samples.size());

            debug("    box ", i, ": ", sim.sampling_start_idx[i + 1] - sim.sampling_start_idx[i], " particles");
        }

        sim.dx = std::sqrt(ppc / T(total_square_samples) * total_volume);
        sim.particle_volume = sim.dx * sim.dx / ppc;
        sim.particle_mass   = sim.rho * sim.particle_volume;

        debug("    dx set to ", sim.dx);

        TV dmin = origins[0];
        TV dmax = origins[0] + sizes[0];

        for (int i = 1; i < origins.size(); ++i) {
            dmin = dmin.cwiseMin(origins[i]);
            dmax = dmax.cwiseMax(origins[i] + sizes[i]);
        }

        sim.Lx = dmax(0) - dmin(0);
        sim.Ly = dmax(1) - dmin(1);

        sim.Np = samples.size();
        debug("    Number of particles samples: ", sim.Np);

        sim.particles = Particles(sim.Np);

        for (int p = 0; p < sim.Np; p++) {
            for (int d = 0; d < 2; d++) {
                sim.particles.x[p](d) = samples[p][d];
            }
        }

    } // end sampleParticles (multiple objects)

#endif

#endif  // SAMPLING_PARTICLES_HPP
