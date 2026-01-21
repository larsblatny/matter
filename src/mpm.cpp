// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "tools.hpp"
#include "simulation/simulation.hpp"
#include "sampling/sampling_particles.hpp"

#include "objects/object_bump.hpp"
#include "objects/object_gate.hpp"
#include "objects/object_ramp.hpp"
#include "objects/object_plate.hpp"

// Comment if not compiling with OpenVDB:
// #include "objects/object_vdb.hpp"
// #include "sampling/sampling_particles_vdb.hpp"


int main(){
    // openvdb::initialize(); // Comment if not using openvdb

    Simulation sim;

    sim.initialize(/*save to file*/ true, /*path*/ "output/", /*name*/ "nonlocal_k001_dgproj");

    sim.save_grid = true;
    sim.end_frame = 100;     // last frame to simulate
    sim.fps = 7500;           // frames per second
    sim.n_threads = 8;      // number of threads in parallel
    sim.cfl = 0.5;          // CFL constant, typically around 0.5
    sim.flip_ratio = 0.99; // (A)PIC-(A)FLIP ratio in [-1,1].


    // INITILIZE ELASTICITY
    sim.elastic_model = ElasticModel::Hencky;
    sim.E = 69e9;     // Young's modulus (Pa)
    sim.nu = 0.33;   // Poisson's ratio (-)
    sim.rho = 1000; // Density (kg/m3)

    ////// GRAVITY ANGLE [default: gravity is 0]
    T theta_deg = 0; // angle in degrees of gravity vector
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero(); //
    // sim.gravity[0] = +9.81 * std::sin(theta);
    // sim.gravity[1] = -9.81 * std::cos(theta);

    ////// INITIAL PARTICLE POSITIONS
    sim.Lx = 4;
    sim.Ly = 5;
    #ifdef THREEDIM
        sim.Lz = 5;
    #endif

    T k_rad = 0.01;
    sampleParticles(sim, k_rad);

    sim.grid_reference_point = TV::Zero();

    // //////////////////////////////////////////////////////////////////////////////////
    // // int Npx = 52+1;
    // // int Npy = 65+1;

    // // int Npx = 104+1;
    // // int Npy = 130+1;

    // int Npx = 208+1;
    // int Npy = 260+1;

    // // int Npx = 416+1;
    // // int Npy = 520+1;

    // T dxp = sim.Lx / (Npx-1.0);

    // sim.dx = 2.0 * dxp;

    // const unsigned int ppc = 4;
    // sim.Np = Npx * Npy;
    // sim.particles = Particles(sim.Np);

    // sim.particle_volume = sim.dx * sim.dx * sim.dx / ppc; // = Lx*Ly*Lz / T(square_samples.size())
    // sim.particle_mass = sim.rho * sim.particle_volume;

    // int p = -1;
    // for(int i = 0; i < Npx; i++){
    //     for(int j = 0; j < Npy; j++){
    //             p++;

    //             // T px = (i+disp_i[d])*sim.dx;
    //             // T py = (j+disp_j[d])*sim.dx;
    //             T px = (i)*dxp;
    //             T py = (j)*dxp;

    //             sim.particles.x[p](0) = px;
    //             sim.particles.x[p](1) = py;
    //             sim.particles.v[p](0) = 0;
    //             sim.particles.v[p](1) = 0;

    //     } // end for i
    // } // end for j
    // debug("Added particles = ", p);
    // if ((p+1) != sim.Np){
    //     debug("Particle number mismatch!!!");
    //     return 0;
    // }
    //////////////////////////////////////////////////////////////////////////////////

    ////// OBJECTS AND TERRAINS
    sim.plates.push_back(std::make_unique<ObjectPlate>(sim.Ly, PlateType::top, BC::NoSlip, 0, -1e15, 1e15, 0, -7, 0)); 
    sim.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::SlipStick, 0)); 



    ////// PLASTICITY
    sim.plastic_model = PlasticModel::VM; 
    sim.q_max = 400e6;
    sim.xi = 2e9;
    sim.xi_nonloc = -4e9;
    sim.nonlocal_l = 0.5;

    sim.q_prefac = std::sqrt(3.0/2.0);

    for(int p = 0; p < sim.Np; p++){
        if (sim.particles.x[p](0) < 1.0 && sim.particles.x[p](1) < 1.0) {
            sim.particles.q_max[p] = 0.9375 * sim.q_max;
        }
        else{
            sim.particles.q_max[p] = sim.q_max;
        }
    }


    sim.simulate();

	return 0;
}
