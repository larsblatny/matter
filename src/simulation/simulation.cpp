// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"

Simulation::Simulation(){
    current_time_step = 0;
    time  = 0;
    frame = 0;
    exit  = 0;

    is_initialized = false;

    runtime_p2g = 0;
    runtime_g2p = 0;
    runtime_euler = 0;
    runtime_defgrad = 0;
    runtime_total = 0;

    // create grid
    grid = Grid();
}


void Simulation::initialize(bool save, std::string dir, std::string name){

    std::cout << "-----------------------------------------------------------------------------------" << std::endl;
    std::cout << "    88b           d88                                                              " << std::endl;
    std::cout << "    888b         d888                  aa          aa                              " << std::endl;
    std::cout << "    88 8b       d8 88                  88          88                              " << std::endl;
    std::cout << "    88  8b     d8  88   adPPYYba   aaaa88aaaa  aaaa88aaaa    adPPYba   8b dPPYba   " << std::endl;
    std::cout << "    88   8b   d8   88  aa      Y8  aaaa88aaaa  8888888888  a8P     88  88P     Y8  " << std::endl;
    std::cout << "    88    8b d8    88   adPPPPP88      88          88      adPPPPP88   88          " << std::endl;
    std::cout << "    88     888     88  88      88      aa          aa      a8b         88          " << std::endl;
    std::cout << "    88      8      88    adPPYba                            adPPYba    88          " << std::endl;
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;

    save_sim = save;
    directory = dir;
    sim_name = name;

    if (save_sim)
        createDirectory();
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;

    is_initialized = true;
}


void Simulation::createDirectory(){

    if (mkdir(directory.c_str(), 0777) == -1)
        std::cerr << "Directory " << directory << " has already been created:  " << strerror(errno) << std::endl;
    else
        std::cout << "Directory " << directory << " was created now"<< std::endl;
    if (mkdir((directory + sim_name).c_str(), 0777) == -1)
        std::cerr << "Simulation " << sim_name << " has already been created: " << strerror(errno) << std::endl;
    else
        std::cout << "Simulation " << sim_name << " was created now" << std::endl;

    std::string file_in  = std::string(SRC_DIR) + "/mpm.cpp";
    std::string file_out = directory + sim_name + "/initial_setup.cpp";

    std::ifstream in(file_in, std::ios::binary);
    std::ofstream out(file_out, std::ios::binary);
    out << in.rdbuf();
    if (!(in && out)){
        std::cerr << "Initial setup " << file_in << " was NOT successfully copied to " << file_out << std::endl;
        exit = 1;
    }

}

void Simulation::simulate(){

    if (!is_initialized){
        debug("Simulation not initialized. Call the initialize(...) function in your input file.");
        return;
    }

    if (elastic_model != ElasticModel::Hencky && plastic_model != PlasticModel::NoPlasticity){
        debug("This plastic model is only compatible with Hencky's elasticity model");
        debug("Please use: elastic_model = Hencky");
        return;
    }

    if (dim == 3){
        debug("This is a 3D simulation.");
    }
    else if (dim == 2){
        debug("This is a 2D simulation.");
    }
    else{
        debug("Unsupported spline degree");
        return;
    }

    #if SPLINEDEG == 3
      apicDinverse = 3.0/(dx*dx);
      debug("Using cubic splines.");
    #elif SPLINEDEG == 2
      apicDinverse = 4.0/(dx*dx);
      debug("Using quadratic splines.");
    #elif SPLINEDEG == 1
        apicDinverse = 0; // NB not implemented
        debug("Using linear hat functions.");
    #else
        #error Unsupported spline degree
    #endif

    if (exit)
        return;

    lambda = nu * E / ( (1.0 + nu) * (1.0 - 2.0*nu) ); // first Lame parameter
    mu = E / (2.0*(1.0+nu)); // shear modulus
    K = calculateBulkModulus(); // bulk modulus
    wave_speed = std::sqrt(E/rho); // elastic wave speed

    dt_max = cfl_elastic * dx / wave_speed;

    frame_dt = 1.0 / fps;

    gravity_final = gravity;

    one_over_dx = 1.0 / dx;
    one_over_dx_square = one_over_dx * one_over_dx;

    d_prefac = 1 / q_prefac;
    e_mu_prefac = 2*q_prefac            * mu;
    f_mu_prefac = 2*q_prefac * q_prefac * mu;

    fac_Q = I_ref / (grain_diameter*std::sqrt(rho_s));

    if (use_mibf){
        if (plastic_model == PlasticModel::DPMui || plastic_model == PlasticModel::MCCMui){
            std::fill(particles.muI.begin(), particles.muI.end(), mu_1);
        } 
        else {// e.g., DPVisc
            std::fill(particles.muI.begin(), particles.muI.end(), M);
        }
    }

    debug("Number of particles: ", Np);
    debug("Grid spacing dx:     ", dx);
    debug("Elastic wave speed:  ", wave_speed);
    debug("Maximum dt:          ", dt_max);
    debug("Particle volume:     ", particle_volume);
    debug("Particle mass:       ", particle_mass);

    time = 0;
    frame = 0;
    final_time = end_frame * frame_dt;

    T time_tolerance = 1e-15;

    if (save_sim){
        saveInfo();
        saveParticleData();
    }

    // Total runtime of simulation
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    while (frame < end_frame){
        if (!reduce_verbose){
            std::cout << "Frame: "               << frame  << " / "    << end_frame  << std::endl;
            std::cout << "               Name: " << sim_name           << std::endl;
            std::cout << "               Step: " << current_time_step  << std::endl;
            std::cout << "               Time: " << time   << " -> "   << (frame+1)*frame_dt << std::endl;
        }
        advanceStep();
        if (exit == 1)
            return;
        time += dt;
        current_time_step++;
        if( frame_dt*(frame+1) - time < time_tolerance ){
            frame++;
            std::cout << "End of frame " << frame << std::endl;
            if (save_sim){
                saveParticleData();
                if (save_grid)
                    saveGridData();
                // saveAvgData();
            }

        }
        if (std::abs(final_time-time) < time_tolerance || final_time < time){
            std::cout << "The simulation ended at time = " << time << std::endl;
            if (save_sim){
                saveParticleData();
                if (save_grid)
                    saveGridData();
                // saveAvgData();
            }
            break;
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    runtime_total = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "Simulation took " << runtime_total << " milliseconds" << std::endl;
    debug("Runtime P2G     = ", runtime_p2g     * 1000.0, " milliseconds");
    debug("Runtime G2P     = ", runtime_g2p     * 1000.0, " milliseconds");
    debug("Runtime Euler   = ", runtime_euler   * 1000.0, " milliseconds");
    debug("Runtime DefGrad = ", runtime_defgrad * 1000.0, " milliseconds");

    if (save_sim)
        saveTiming();
}





void Simulation::advanceStep(){
    updateDt();

    if (pbc){
        PBCAddParticles(4);
    }

    resizeGrid();

    timer t_p2g; t_p2g.start();
    P2G();
    t_p2g.stop(); runtime_p2g += t_p2g.get_timing();



    #ifdef DIMMATTER

        // Write velocities at coupling indices to file
        std::ofstream out_file(dim_vels_filepath);
        if (!out_file.is_open()) {
            std::cerr << "Unable to save velocities" << std::endl;
            exit = 1;
            return;
        }
        bool firstline = true;
        for (const int& index : coupling_indices){ 
            if (!firstline){
                out_file << "\n";
            }
            TV v = grid.v[index];
            #ifdef THREEDIM
                out_file << v(0) << " " << v(1) << " " << v(2);
            #else
                out_file << v(0) << " " << v(1);
            #endif
            firstline = false;
        }
        out_file.close();

        // At this point, DIM will be calculating the drag force
        // code should pause here while this is done
        // before the drag force is read by Matter in explicitEulerUpdate

    #endif // DIMMATTER



    // checkMassConservation();
    // checkMomentumConservation();

    timer t_euler; t_euler.start();
    explicitEulerUpdate();
    t_euler.stop(); runtime_euler += t_euler.get_timing();

    // addExternalParticleGravity();

    if (pbc){
        PBCDelParticles();
    }


    timer t_g2p; t_g2p.start();
    G2P();
    t_g2p.stop(); runtime_g2p += t_g2p.get_timing();

    if (use_musl == true)
        MUSL();

    timer t_defgrad; t_defgrad.start();
    deformationUpdate();
    t_defgrad.stop(); runtime_defgrad += t_defgrad.get_timing();

    positionUpdate();

    moveObjects();

} // end advanceStep



void Simulation::moveObjects(){
    for (auto& obj : plates) {
        obj->move(dt, frame_dt, time);
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// EXTRA HELPER FUNCTIONS ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

T Simulation::calculateBulkModulus(){
    return lambda + 2.0 * mu / dim;
}

void Simulation::checkMomentumConservation(){
    TV momentum_grid = TV::Zero();
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
        #ifdef THREEDIM
            for(int k=0; k<Nz; k++)
                momentum_grid += grid.mass[ind(i,j,k)] * grid.v[ind(i,j,k)];
        #else
            momentum_grid += grid.mass[ind(i,j)] * grid.v[ind(i,j)];
        #endif
        }
    }

    TV momentum_particle = TV::Zero();
    for(int p=0; p<Np; p++){
        momentum_particle += particle_mass * particles.v[p];
    }
    debug("               Total part momentum = ", momentum_grid.norm());
    debug("               Total grid momentum = ", momentum_particle.norm());

    if ( (momentum_grid-momentum_particle).norm() > 1e-10 * momentum_particle.norm() ){
        debug("MOMENTUM NOT CONSERVED!!!");
        exit = 1;
        return;
    }
}

void Simulation::checkMassConservation(){
    T particle_mass_total = particle_mass*Np;
    T grid_mass_total = 0;
    for(auto&& m: grid.mass)
        grid_mass_total += m;

    debug("               Total grid mass = ", grid_mass_total    );
    debug("               Total part mass = ", particle_mass_total);

    if ( std::abs(grid_mass_total-particle_mass_total) > 1e-10 * particle_mass_total ){
        debug("MASS NOT CONSERVED!!!");
        exit = 1;
        return;
    }
}


// This function is to be used in explicitEulerUpdate after boundaryCollision
// It must be hard-coded to choice
void Simulation::overwriteGridVelocity(TV Xi, TV& vi){
    T y_start = Ly - 0.25*dx;
    T width = 2*dx;
    T v_imp = 0.1; // positive value means tension
    if (Xi(1) > y_start - width + v_imp * time)
        vi(1) = v_imp;
    if (Xi(1) < 0       + width - v_imp * time)
        vi(1) = -v_imp;
}


// This function must be hard-coded to choice
void Simulation::addExternalParticleGravity(){
    // // 1. Transfer grid velocity to particles
    // G2P();
    // // 2. Apply gravity on particle velocity
    // for(int p=0; p<Np; p++)
    //     particles.v[p] += dt * (-2.0*amplitude*particles.x0[p]);
    // // 3. Transfer particle velocity back to grid
    // P2G();
} // end addExternalParticleGravity


// // TODO:
// void Simulation::boundaryCorrection(T xi, T yi, T& vxi, T& vyi){
//
//     // trial step
//     T x_next = xi + vxi * dt;
//     T y_next = yi + vyi * dt;
//     moveObjects();
//
//     for (ObjectPlate &obj : objects_plate) {
//         bool colliding = obj.inside(x_next, y_next);
//         if (colliding) {
//             if (obj.plate_type == upper ){
//                 debug(obj.name, " dist = ", obj.distance(x_next, y_next));
//                 vyi += obj.distance(x_next, y_next) / dt; // distance is negative since grid point is inside object
//             }
//             else if (obj.plate_type == lower){
//                 debug(obj.name, " dist = ", obj.distance(x_next, y_next));
//                 vyi -= obj.distance(x_next, y_next) / dt;
//             }
//         } // end if colliding
//
//     } // end iterator over objects
//
//     // Correct back
//     moveObjects(-dt);
//
// } // end boundaryCorrection
