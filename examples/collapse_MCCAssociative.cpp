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

    sim.initialize(/*save to file*/ true, /*path*/ "output/", /*name*/ "collapse_MCCAssociative");

    sim.save_grid = true;
    sim.end_frame = 40;     // last frame to simulate
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
    T theta_deg = 20; // angle in degrees of gravity vector
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero(); //
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    ////// INITIAL PARTICLE POSITIONS
    sim.Lx = 1;
    sim.Ly = 1;
    T k_rad = 0.015;
    #ifdef THREEDIM
        sim.Lz = 5;
    #endif
    sampleParticles(sim, k_rad);

    sim.grid_reference_point = TV::Zero();

    ////// OBJECTS AND TERRAINS
    T friction = 0.2;
    sim.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::left, BC::SlipFree, friction)); 

    ////// PLASTICITY
    sim.plastic_model = PlasticModel::MCCAssociative; // MCC yield surface with associative viscoplastic flow rule
    sim.hardening_law = HardeningLaw::ExpoExpl;

    sim.M = std::tan(25*M_PI/180.0); // Internal friction
    sim.p0 = 5.0e3;
    sim.beta = 0.0;
    sim.xi = 0.02;
    sim.eta = 2.0; // Comment to remove viscosity and compare

    sim.simulate();

	return 0;
}
