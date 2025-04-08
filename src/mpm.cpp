// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "tools.hpp"
#include "simulation/simulation.hpp"
#include "sampling/sampling_particles.hpp"
#include "objects/object_plate.hpp"



int main(){

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

    ////// PLASTICITY
    sim.plastic_model = PlasticModel::DPVisc; // Perzyna model with Drucker_Prager yield surface

    sim.use_pradhana = true; // Supress unwanted volume expansion in Drucker-Prager models
    sim.q_prefac = 1.0 / std::sqrt(2.0); // [default: sqrt(1/2)] Prefactor in def. of q, here q = sqrt(1/2 * s:s)

    sim.M = std::tan(30*M_PI/180.0); // Internal friction
    sim.q_cohesion = 0; // Yield surface's intercection of q-axis (in Pa), 0 is the cohesionless case
    sim.perzyna_exp = 1; // Exponent in Perzyna models
    sim.perzyna_visc = 0; // Viscous time parameter is Perzyna models

    ////// INITIAL PARTICLE POSITIONS
    sim.Lx = 1;
    sim.Ly = 1;
    T k_rad = 0.01;
    #ifdef THREEDIM
        sim.Lz = 0.1;
    #endif
    sampleParticles(sim, k_rad);

    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        #ifdef THREEDIM
            sim.particles.x[p](2) -= 0.5*sim.Lz;
        #endif
    }

    sim.grid_reference_point = TV::Zero();

    ////// OBJECTS AND TERRAINS
    sim.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::NoSlip));
    #ifdef THREEDIM
        sim.plates.push_back(std::make_unique<ObjectPlate>(-0.5*sim.Lz, PlateType::back, BC::SlipFree));
        sim.plates.push_back(std::make_unique<ObjectPlate>(0.5*sim.Lz, PlateType::front, BC::SlipFree));
    #endif

    ////// FIXED GRID
    #ifdef THREEDIM
        TV Lmin(-3, 0,-0.05);
        TV Lmax(3,  1, 0.05);
    #else
        TV Lmin(-3, 0);
        TV Lmax(3,  1);
    #endif
    sim.remeshFixed(1, Lmin, Lmax); 


    /////////////////////////////////
    /////////// DIMMMATER ///////////
    /////////////////////////////////

    T drag_constant = 2;
    T material_diameter = 0.02;
    sim.dim_drag_constant = (1.0/8.0) * drag_constant * M_PI * std::pow(material_diameter,2) * sim.rho;
    sim.dim_vels_filepath = sim.directory + sim.sim_name + "/nodevels.txt";
    sim.dim_drag_filepath = sim.directory + sim.sim_name + "/nodedrag.txt";
    sim.dim_inds_filepath = sim.directory + sim.sim_name + "/nodeinds.txt";


    ///// CREATE LIST OF COUPLING INDICES
    std::vector<int> coupling_indices_temp;
    for (int i = 0; i<sim.Nx; i++){
        T xi = sim.grid.xc + i*sim.dx;
        for (int j = 0; j<sim.Ny; j++){
            T yi = sim.grid.yc + j*sim.dx;
            #ifdef THREEDIM
            for (int k = 0; k<sim.Nz; k++){
                if (xi > 1 && yi < 0.25)
                    coupling_indices_temp.push_back((i*sim.Ny + j) * sim.Nz + k);
            }
            #else
                if (xi > 1 && yi < 0.25)
                    coupling_indices_temp.push_back(i*sim.Ny+j);
            #endif
                
        }
    }

    ////// SAVE COUPLING INDICES TO FILE
    std::ofstream out_file(sim.dim_inds_filepath);
        if (!out_file.is_open()) {
            std::cerr << "Unable to save coupling inds" << std::endl;
            return 0;
        }
        bool firstline = true;
        for (const int& index : coupling_indices_temp){ 
            if (!firstline){
                out_file << "\n";
            }
            out_file << index;
            firstline = false;
        }
    out_file.close();


    // READ COUPLING INDICES FROM FILE (UNNECESSARY, EXAMPLE ONLY)
    std::ifstream in_file(sim.dim_inds_filepath);
    if (!in_file.is_open()) {
        std::cerr << "Unable to open coupling indices" << std::endl;
    }
    int value;
    while (in_file >> value) {
        sim.coupling_indices.push_back(value);
    }
    in_file.close();

    debug("Num of coupling inds = ", sim.coupling_indices.size());

    sim.simulate();

	return 0;
}
