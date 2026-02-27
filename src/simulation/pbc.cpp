// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"

void Simulation::PBCAddParticles1D(){

    TV incr = TV::Zero();
    incr(0) = dx;

    num_add_pbc_particles = 0;
    for(int p = 0; p < Np; p++){

        TV part_x = particles.x[p]; // copy

        particles.x.push_back(part_x+incr);
        particles.v.push_back(particles.v[p]);
        particles.F.push_back(particles.F[p]);
        particles.Bmat.push_back(particles.Bmat[p]);
        particles.tau.push_back(particles.tau[p]);
        num_add_pbc_particles++;

        particles.x.push_back(part_x-incr);
        particles.v.push_back(particles.v[p]);
        particles.F.push_back(particles.F[p]);
        particles.Bmat.push_back(particles.Bmat[p]);
        particles.tau.push_back(particles.tau[p]);
        num_add_pbc_particles++;

    } // end loop over p

    Np = particles.x.size();

#ifdef WARNINGS
    debug("Added particles:                    ", num_add_pbc_particles);
    debug("Number of particles after addition: ", particles.x.size()   );
#endif
} // end PBCAddParticles1D()

void Simulation::PBCAddParticles(unsigned int safety_factor){

    num_add_pbc_particles = 0;

    #pragma omp parallel for reduction(+:num_add_pbc_particles) num_threads(n_threads)
    for(int p = 0; p < Np; p++){

        TV part_x = particles.x[p]; // copy
        T diff;

        diff = Lx - part_x(0);
        if ( diff < safety_factor*dx && diff > 0){
            part_x(0) = -diff;
            #pragma omp critical
            {
                particles.x.push_back(part_x);
                particles.v.push_back(particles.v[p]);
                particles.F.push_back(particles.F[p]);
                particles.Bmat.push_back(particles.Bmat[p]);
                particles.tau.push_back(particles.tau[p]);
            }
            num_add_pbc_particles++;
            continue; // go directly to next p
        }

        diff = part_x(0);
        if ( diff < safety_factor*dx && diff >= 0 ){
            part_x(0) = Lx + diff;
            #pragma omp critical
            {
                particles.x.push_back(part_x);
                particles.v.push_back(particles.v[p]);
                particles.F.push_back(particles.F[p]);
                particles.Bmat.push_back(particles.Bmat[p]);
                particles.tau.push_back(particles.tau[p]);
            }
            num_add_pbc_particles++;
            continue; // go directly to next p
        }

    } // end loop over p

    Np = particles.x.size();

#ifdef WARNINGS
    debug("Added particles:                    ", num_add_pbc_particles);
    debug("Number of particles after addition: ", particles.x.size()   );
#endif
} // end PBCAddParticles()

void Simulation::PBCDelParticles(){

    #ifdef WARNINGS
        debug("PBCDelParticles");
    #endif

    particles.x.erase(   particles.x.end()   -num_add_pbc_particles, particles.x.end()   );
    particles.v.erase(   particles.v.end()   -num_add_pbc_particles, particles.v.end()   );
    particles.F.erase(   particles.F.end()   -num_add_pbc_particles, particles.F.end()   );
    particles.Bmat.erase(particles.Bmat.end()-num_add_pbc_particles, particles.Bmat.end());
    particles.tau.erase( particles.tau.end() -num_add_pbc_particles, particles.tau.end() );
    Np = particles.x.size();

#ifdef WARNINGS
    debug("Number of particles after delition: ", particles.x.size());
    debug("Np:                                 ", Np                );
#endif
} // end PBCDelParticles()
