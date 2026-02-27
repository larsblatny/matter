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

    sim.initialize(/*save to file*/ true, /*path*/ "output/", /*name*/ "collapse");

    sim.save_grid = true;
    sim.end_frame = 20;     // last frame to simulate
    sim.fps = 10;           // frames per second
    sim.n_threads = 8;      // number of threads in parallel
    sim.cfl = 0.5;          // CFL constant, typically around 0.5
    sim.flip_ratio = -0.95; // (A)PIC-(A)FLIP ratio in [-1,1].

    // INITILIZE ELASTICITY
    sim.elastic_model = ElasticModel::Hencky;
    sim.E = 1e6;     // Young's modulus (Pa)
    sim.nu = 0.3;   // Poisson's ratio (-)
    sim.rho = 1000; // Density (kg/m3)

    ////// GRAVITY ANGLE [default: gravity is 0]
    T theta_deg = 0; // angle in degrees of gravity vector
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero(); //
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    ////// INITIAL PARTICLE POSITIONS
    sim.Lx = 1;
    sim.Ly = 1;
    T k_rad = 0.01;
    #ifdef THREEDIM
        sim.Lz = 5;
    #endif
    sampleParticles(sim, k_rad);

    ////// OPTIONAL: CHANGE INITIAL PARTICLE POSITIONS
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) += 0.5*sim.dx;
    }
    sim.grid_reference_point = TV::Zero();

    ////// OPTIONAL: INITIAL PARTICLE VELOCITIES
    // sim.particles.v = ...

    ////// OBJECTS AND TERRAINS
    sim.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::NoSlip)); 

    /////// Here are some examples how to use the objects derived from ObjectGeneral:
    // T friction = 0.2; 
    // sim.objects.push_back(std::make_unique<ObjectBump>(BC::SlipFree, friction));
    // sim.objects.push_back(std::make_unique<ObjectGate>(BC::SlipFree, friction));

    /////// Here is an example how to use ObjectVdb (uncomment includes and openvdb::initialize() above):
    // sim.objects.push_back(std::make_unique<ObjectVdb>("../levelsets/vdb_file_name.vdb", BC::NoSlip, friction));

    ////// PLASTICITY
    sim.plastic_model = PlasticModel::DPVisc; // Viscous model with Drucker_Prager yield surface

    sim.use_pradhana = true; // Supress unwanted volume expansion in Drucker-Prager models
    sim.q_prefac = 1.0 / std::sqrt(2.0); // [default: sqrt(1/2)] Prefactor in def. of q, here q = sqrt(1/2 * s:s)

    sim.M = std::tan(30*M_PI/180.0); // Internal friction
    sim.q_cohesion = 0; // Yield surface's intercection of q-axis (in Pa), 0 is the cohesionless case
    sim.visc_exponent = 1; // Exponent in viscous models
    sim.visc_time = 0; // Viscous time parameter

    sim.simulate();

	return 0;
}
