// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef REGULAR_SAMPLING_HPP
#define REGULAR_SAMPLING_HPP

#include "../tools.hpp"
#include "../data_structures.hpp"

template <typename S>
#ifdef THREEDIM
void RegularSampleParticles(S& sim, T ppc = 8, unsigned int crop_to_shape = 0)
#else
void RegularSampleParticles(S& sim, T ppc = 4, unsigned int crop_to_shape = 0)
#endif
{
    debug("Regular sampling of particles...");

    #ifdef THREEDIM

    const T Lx = sim.Lx;
    const T Ly = sim.Ly;
    const T Lz = sim.Lz;

    const T Lx0 = 0;
    const T Ly0 = 0;
    const T Lz0 = 0;

    const T pSpacing = sim.dx / std::cbrt(ppc);

    const int nx = int(Lx / pSpacing);
    const int ny = int(Ly / pSpacing);
    const int nz = int(Lz / pSpacing);

    std::vector<TV> square_samples;
    std::vector<TV> samples;

    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    for (int k = 0; k < nz; ++k)
    {
        TV point(
            Lx0 + (i + T(0.5)) * pSpacing,
            Ly0 + (j + T(0.5)) * pSpacing,
            Lz0 + (k + T(0.5)) * pSpacing
        );

        square_samples.push_back(point);
    }

    debug("    Number of square samples: ", square_samples.size());
    debug("    dx set to ", sim.dx);

    /////// Cylinder
    if (crop_to_shape == 1){
        for(int p = 0; p < square_samples.size(); p++){
            if ( (square_samples[p][0] - Lx / 2.0)*(square_samples[p][0] - Lx / 2.0) + (square_samples[p][2] - Lz / 2.0)*(square_samples[p][2] - Lz / 2.0) < (Lx / 2.0)*(Lx / 2.0) ){
                samples.push_back(square_samples[p]);
            }
        }
    }

    //////// Silo
    else if (crop_to_shape == 2){
        for(int p = 0; p < square_samples.size(); p++){
            T x = square_samples[p][0] - Lx / 2.0;
            T y = square_samples[p][1];
            T z = square_samples[p][2] - Lz / 2.0;
            T r_surface = std::tanh(y) + 1;
            T r_surface_sq = r_surface * r_surface;
            T r_point_sq = x * x + z * z;
            if (r_point_sq < r_surface_sq){
                samples.push_back(square_samples[p]);
            }
        }
    }

    else{
        debug("    No shape specified, using just a square.");
        samples = square_samples;
    }

    sim.particle_volume = pSpacing * pSpacing * pSpacing;
    sim.particle_mass   = sim.rho * sim.particle_volume;

    #else 

    const T Lx = sim.Lx;
    const T Ly = sim.Ly;

    const T Lx0 = 0;
    const T Ly0 = 0;

    const T pSpacing = sim.dx / std::sqrt(ppc);

    const int nx = int(Lx / pSpacing);
    const int ny = int(Ly / pSpacing);

    std::vector<TV> square_samples;
    std::vector<TV> samples;

    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    {
        TV point(
            Lx0 + (i + T(0.5)) * pSpacing,
            Ly0 + (j + T(0.5)) * pSpacing
        );

        samples.push_back(point);
    }

    sim.particle_volume = pSpacing * pSpacing;
    sim.particle_mass   = sim.rho * sim.particle_volume;

    debug("    Number of square samples: ", square_samples.size());
    debug("    dx set to ", sim.dx);

    /////// Quadratic Gate
    if (crop_to_shape == 1){
        T height = 0.05; // 0.016;
        for(int p = 0; p < square_samples.size(); p++){
            T xp = square_samples[p][0];
            T y_gate = height + 100 * (xp - Lx) * (xp - Lx) - 0.5 * sim.dx;
            if (square_samples[p][1] < y_gate){
                samples.push_back(square_samples[p]);
            }
        }
    }

    else{
        debug("    No shape specified, using just a square.");
        samples = square_samples;
    }

    sim.particle_volume = pSpacing * pSpacing;
    sim.particle_mass   = sim.rho * sim.particle_volume;

    #endif

    sim.Np = samples.size();
    debug("    Number of particle samples: ", sim.Np);
    debug("    Particle spacing: ", pSpacing);

    sim.particles = Particles(sim.Np);
    sim.particles.x = samples;

}

#endif // REGULAR_SAMPLING_HPP
