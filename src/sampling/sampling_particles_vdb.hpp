// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef SAMPLING_PARTICLES_VDB_HPP
#define SAMPLING_PARTICLES_VDB_HPP

#include "../tools.hpp"
#include "../data_structures.hpp"
#include "../../deps/tph_poisson-0.3/thinks/poisson_disk_sampling/poisson_disk_sampling.h"

#include "../objects/object_vdb.hpp"

template <typename S>
#ifdef THREEDIM
void sampleParticlesFromVdb(S& sim, ObjectVdb& obj, T kRadius, T ppc = 8)
#else // TWODIM
void sampleParticlesFromVdb(S& sim, ObjectVdb& obj, T kRadius, T ppc = 6)
#endif
{

    debug("Sampling particles from VDB...");

    std::uint32_t kAttempts = 30;
    std::uint32_t kSeed = 42;

    TV min_corner, max_corner;
    obj.bounds(min_corner, max_corner);
    TV L = max_corner - min_corner;

    #ifdef THREEDIM
        debug("    Min corner: ", min_corner(0), ", ", min_corner(1), ", ", min_corner(2));
        debug("    Max corner: ", max_corner(0), ", ", max_corner(1), ", ", max_corner(2));

        std::array<T, 3> kXMin = std::array<T, 3>{{min_corner(0), min_corner(1), min_corner(2)}};
        std::array<T, 3> kXMax = std::array<T, 3>{{max_corner(0), max_corner(1), max_corner(2)}};
        std::vector<std::array<T, 3>> square_samples = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax, kAttempts, kSeed);

        sim.dx = std::cbrt(ppc / T(square_samples.size()) * L(0)*L(1)*L(2));
        sim.particle_volume = sim.dx * sim.dx * sim.dx / ppc;
        sim.particle_mass = sim.rho * sim.particle_volume;

        sim.Lx = L(0);
        sim.Ly = L(1);
        sim.Lz = L(2);
    #else // TWODIM
        debug("    Min corner: ", min_corner(0), ", ", min_corner(1));
        debug("    Max corner: ", max_corner(0), ", ", max_corner(1));

        std::array<T, 2> kXMin = std::array<T, 2>{{min_corner(0), min_corner(1)}};
        std::array<T, 2> kXMax = std::array<T, 2>{{max_corner(0), max_corner(1)}};
        std::vector<std::array<T, 2>> square_samples = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax, kAttempts, kSeed);

        sim.dx = std::sqrt(ppc / T(square_samples.size()) * L(0)*L(1));
        sim.particle_volume = sim.dx * sim.dx / ppc;
        sim.particle_mass = sim.rho * sim.particle_volume;

        sim.Lx = L(0);
        sim.Ly = L(1);
    #endif

    debug("    Number of square samples: ", square_samples.size());
    debug("    dx set to ", sim.dx);

    std::vector<TV> samples;
    for(int p = 0; p < square_samples.size(); p++){

        #ifdef THREEDIM
            TV point(square_samples[p][0], square_samples[p][1], square_samples[p][2]);
        #else // TWODIM
            TV point(square_samples[p][0], square_samples[p][1]);
        #endif

        if ( obj.inside(point) ){
            samples.push_back(point);
        }
    }

    sim.Np = samples.size();
    debug("    Number of particles samples: ", sim.Np);

    sim.particles = Particles(sim.Np);
    sim.particles.x = samples;

} // end sampleParticlesFromVdb

template <typename S>
#ifdef THREEDIM
void sampleParticlesFromVdb(S& sim, std::vector<ObjectVdb> objects, T kRadius, T ppc = 8)
#else // TWODIM
void sampleParticlesFromVdb(S& sim, std::vector<ObjectVdb> objects, T kRadius, T ppc = 6)
#endif
{
    debug("Sampling particles from ", objects.size(), " VDB objects...");

    std::uint32_t kAttempts = 5;

    std::vector<TV> samples;  

    sim.sampling_start_idx.push_back(0);

    T total_volume = 0; 
    T total_square_samples = 0; 

    for (int i = 0; i < objects.size(); ++i) {
        ObjectVdb&    obj   = objects[i];
        std::uint32_t kSeed = 42;    

        TV min_corner, max_corner;
        obj.bounds(min_corner, max_corner);
        TV L = max_corner - min_corner;

        debug("    Min corner: ", min_corner(0), ", ", min_corner(1), ", ", min_corner(2));
        debug("    Max corner: ", max_corner(0), ", ", max_corner(1), ", ", max_corner(2));

    #ifdef THREEDIM
        std::array<T, 3> kXMin = std::array<T, 3>{{min_corner(0), min_corner(1), min_corner(2)}};
        std::array<T, 3> kXMax = std::array<T, 3>{{max_corner(0), max_corner(1), max_corner(2)}};
        std::vector<std::array<T, 3>> square_samples = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax, kAttempts, kSeed);
        total_volume += L(0) * L(1) * L(2);
    #else // TWODIM
        std::array<T, 2> kXMin = std::array<T, 2>{{min_corner(0), min_corner(1)}};
        std::array<T, 2> kXMax = std::array<T, 2>{{max_corner(0), max_corner(1)}};
        std::vector<std::array<T, 2>> square_samples = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax, kAttempts, kSeed);
        total_volume += L(0) * L(1);
    #endif
        total_square_samples += square_samples.size();

        for (int p = 0; p < square_samples.size(); ++p) {
    #ifdef THREEDIM
        TV point(square_samples[p][0], square_samples[p][1], square_samples[p][2]);
    #else // TWODIM
        TV point(square_samples[p][0], square_samples[p][1]);
    #endif
        if (obj.inside(point)) {
            samples.push_back(point);
        }
    }

        sim.sampling_start_idx.push_back(samples.size());        

        debug("    object ", i, ": ", sim.sampling_start_idx[i + 1] - sim.sampling_start_idx[i], " particles");

    }

    #ifdef THREEDIM
        sim.dx              = std::cbrt(ppc / T(total_square_samples) * total_volume);
        sim.particle_volume = sim.dx * sim.dx * sim.dx / ppc;
    #else // TWODIM
        sim.dx              = std::sqrt(ppc / T(total_square_samples) * total_volume);
        sim.particle_volume = sim.dx * sim.dx / ppc;
    #endif

    sim.particle_mass = sim.rho * sim.particle_volume;

    TV union_min;
    TV union_max;

    for (int i = 0; i < objects.size(); ++i){
        TV min;
        TV max;
        objects[i].bounds(min, max);
        union_min = union_min.cwiseMin(min);
        union_max = union_max.cwiseMax(max);
    }
    TV LU = union_max - union_min;
    sim.Lx = LU(0);
    sim.Ly = LU(1);
    #ifdef THREEDIM
        sim.Lz = LU(2);
    #endif

    debug("    dx set to ", sim.dx);
    sim.Np = samples.size();
    debug("    Number of particles sampled: ", sim.Np);

    sim.particles     = Particles(sim.Np);
    sim.particles.x   = samples;

} // end sampleParticlesFromVdb

#endif  // SAMPLING_PARTICLES_VDB_HPP
