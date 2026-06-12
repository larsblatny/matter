// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "tools.hpp"
#include "simulation/simulation.hpp"
#include "sampling/sampling_particles.hpp"

#include "objects/object_bump.hpp"
#include "objects/object_gate.hpp"
#include "objects/object_ramp.hpp"
#include "objects/object_plate.hpp"

// Comment if not compiling with OpenVDB:
#include "objects/object_vdb.hpp"
#include "sampling/sampling_particles_vdb.hpp"


int main(){
    openvdb::initialize(); // Comment if not using openvdb

    Simulation sim;

    sim.initialize(/*save to file*/ true, /*path*/ "output/", /*name*/ "mountain");

    sim.save_grid = false; //true;
    sim.end_frame = 100; //500;     // last frame to simulate
    sim.fps = 5;           // frames per second
    sim.n_threads = 8;      // number of threads in parallel
    sim.cfl = 0.5;          // CFL constant, typically around 0.5
    sim.flip_ratio = -0.95; // (A)PIC-(A)FLIP ratio in [-1,1].
    sim.reduce_verbose = true;

    sim.use_sparse = false;

    // INITILIZE ELASTICITY
    sim.elastic_model = ElasticModel::Hencky;
    sim.E = 5e5;     // Young's modulus (Pa)
    sim.nu = 0.3;   // Poisson's ratio (-)
    sim.rho = 2000; // Density (kg/m3)

    ////// GRAVITY ANGLE [default: gravity is 0]
    T theta_deg = 30.; // angle in degrees of gravity vector
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero(); //
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    ////// INITIAL PARTICLE POSITIONS    
    std::string release_path = std::string("../data/mountain/mountain_release_shifted.vdb");
    ObjectVdb release = ObjectVdb(release_path);
    sampleParticlesFromVdb(sim, release, 0.18); // h = 0.431492478463 m
    // sampleParticlesFromVdb(sim, release, 0.1668627); // h = 0.400046583858 m
    // sampleParticlesFromVdb(sim, release, 0.06663595279); // h = 0.159996934438 m

   
    sim.grid_reference_point = TV::Zero();

    ////// OPTIONAL: INITIAL PARTICLE VELOCITIES
    // sim.particles.v = ...

    ////// OBJECTS AND TERRAINS
    // sim.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::NoSlip)); 
    // sim.plates.push_back(std::make_unique<ObjectPlate>(2., PlateType::bottom, BC::NoSlip)); 

    std::string terrain_path = std::string("../data/mountain/mountain_terrain_shifted.vdb");
    std::string name;
    name = "Terrain"; 
    sim.objects.push_back(std::make_unique<ObjectVdb>(terrain_path, BC::SlipFree, 0.35, name));

    /////// Here are some examples how to use the objects derived from ObjectGeneral:
    // T friction = 0.2; 
    // sim.objects.push_back(std::make_unique<ObjectBump>(BC::SlipFree, friction));
    // sim.objects.push_back(std::make_unique<ObjectGate>(BC::SlipFree, friction));

    /////// Here is an example how to use ObjectVdb (uncomment includes and openvdb::initialize() above):
    // sim.objects.push_back(std::make_unique<ObjectVdb>("../levelsets/vdb_file_name.vdb", BC::NoSlip, friction));

    ////// PLASTICITY
    // sim.plastic_model = PlasticModel::DPVisc; // Perzyna model with Drucker_Prager yield surface

    // sim.use_pradhana = true; // Supress unwanted volume expansion in Drucker-Prager models
    // sim.q_prefac = 1.0 / std::sqrt(2.0); // [default: sqrt(1/2)] Prefactor in def. of q, here q = sqrt(1/2 * s:s)

    // sim.M = std::tan(20*M_PI/180.0); // Internal friction
    // sim.q_cohesion = 0; // Yield surface's intercection of q-axis (in Pa), 0 is the cohesionless case
    // sim.perzyna_exp = 1; // Exponent in Perzyna models
    // sim.perzyna_visc = 0.12; // Viscous time parameter is Perzyna models




    sim.plastic_model = PlasticModel::DP; // Perzyna model with Drucker_Prager yield surface

    sim.q_prefac = 1.0 / std::sqrt(2.0); // [default: sqrt(1/2)] Prefactor in def. of q, here q = sqrt(1/2 * s:s)

    sim.M = std::tan(30*M_PI/180.0); // Internal friction
    sim.q_cohesion = 0; // Yield surface's intercection of q-axis (in Pa), 0 is the cohesionless case


    sim.simulate();

	return 0;
}
