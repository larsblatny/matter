// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef REGULAR_SAMPLING_VDB_HPP
#define REGULAR_SAMPLING_VDB_HPP

#include "../tools.hpp"
#include "../data_structures.hpp"
#include "../objects/object_vdb.hpp"

template <typename S>
#ifdef THREEDIM
void regularSampleParticlesFromVdb(S& sim, ObjectVdb& obj, T ppc = 8)
#else
void regularSampleParticlesFromVdb(S& sim, ObjectVdb& obj, T ppc = 4)
#endif
{

    debug("Regular sampling of particles from VDB...");

    TV min_corner, max_corner;
    obj.bounds(min_corner, max_corner);
    TV L = max_corner - min_corner;

    #ifdef THREEDIM
    debug("    Min corner: ", min_corner(0), ", ", min_corner(1), ", ", min_corner(2));
    debug("    Max corner: ", max_corner(0), ", ", max_corner(1), ", ", max_corner(2));

    const T pSpacing = sim.dx / std::cbrt(ppc);

    const int nx = int(L(0) / pSpacing);
    const int ny = int(L(1) / pSpacing);
    const int nz = int(L(2) / pSpacing);

    std::vector<TV> samples;

    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    for (int k = 0; k < nz; ++k)
    {
        TV point(
            min_corner(0) + (i + T(0.5)) * pSpacing,
            min_corner(1) + (j + T(0.5)) * pSpacing,
            min_corner(2) + (k + T(0.5)) * pSpacing
        );

        if (obj.inside(point))
            samples.push_back(point);
    }

    sim.particle_volume = pSpacing * pSpacing * pSpacing;
    sim.particle_mass   = sim.rho * sim.particle_volume;

    sim.Lx = L(0);
    sim.Ly = L(1);
    sim.Lz = L(2);

    #else 

    debug("    Min corner: ", min_corner(0), ", ", min_corner(1));
    debug("    Max corner: ", max_corner(0), ", ", max_corner(1));

    const T pSpacing = sim.dx / std::sqrt(ppc);

    const int nx = int(L(0) / pSpacing);
    const int ny = int(L(1) / pSpacing);

    std::vector<TV> samples;

    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    {
        TV point(
            min_corner(0) + (i + T(0.5)) * pSpacing,
            min_corner(1) + (j + T(0.5)) * pSpacing
        );

        if (obj.inside(point))
            samples.push_back(point);
    }

    sim.particle_volume = pSpacing * pSpacing;
    sim.particle_mass   = sim.rho * sim.particle_volume;

    sim.Lx = L(0);
    sim.Ly = L(1);

    #endif

    sim.Np = samples.size();
    debug("    Number of particle samples: ", sim.Np);
    debug("    Particle spacing: ", pSpacing);

    sim.particles = Particles(sim.Np);
    sim.particles.x = samples;
}

#endif // REGULAR_SAMPLING_VDB_HPP
