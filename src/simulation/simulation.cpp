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

    setup_file = std::string(SRC_DIR) + "/mpm.cpp"; // overwritten by Python
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

    if (setup_file.empty()){ // e.g. an interactive Python session, no script to copy
        std::cout << "No setup file to copy" << std::endl;
        return;
    }

    std::string file_in = setup_file;

    // Keep the extension of the input file, so a Python run gets initial_setup.py
    std::string extension = "";
    size_t dot   = file_in.find_last_of('.');
    size_t slash = file_in.find_last_of("/\\");
    if (dot != std::string::npos && (slash == std::string::npos || dot > slash))
        extension = file_in.substr(dot);

    std::string file_out = directory + sim_name + "/initial_setup" + extension;

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

    if (use_mibf && use_basal_friction_field){
        debug("use_mibf and use_basal_friction_field cannot both be true, please set one of them to false.");
        return;
    }

    if (use_basal_friction_field && !basal_friction_field.isSet()){
        debug("use_basal_friction_field is true but no friction field has been set.");
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

    if (max_dt < 0) {
        max_dt = cfl_elastic * dx / wave_speed;
    }
    
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

    // Initialize grid
    if (use_sparse){
        debug("Using sparse grid.");

        if (pbc){
            debug("ERROR: Sparse grid (use_sparse) is not yet compatible with periodic boundary conditions (pbc), please set one of them to false.");
            return;
        }
        if (use_musl){
            debug("ERROR: Sparse grid (use_sparse) is not yet compatible with MUSL (use_musl), please set one of them to false.");
            return;
        }
        if (dim == 2){
            debug("ERROR: Sparse grid (use_sparse) is not yet compatible with 2D, please set use_sparse = false, or use 3D.");
            return;
        }

        int B = 4; // Block size
        BlockScanGrid::I3 bmin{0, 0, 0}; // Needs to be reset before calling scan_sparse_grid.init()
        BlockScanGrid::I3 bmax{8, 8, 8}; // Needs to be reset before calling scan_sparse_grid.init()
        scan_sparse_grid.init(B, bmin, bmax, int(n_threads));
    }
    else{
        grid = Grid();
    }

    debug("Number of particles: ", Np);
    debug("Grid spacing dx:     ", dx);
    debug("Elastic wave speed:  ", wave_speed);
    debug("Maximum dt:          ", max_dt);
    debug("Particle volume:     ", particle_volume);
    debug("Particle mass:       ", particle_mass);

    current_time_step = 0;
    time = 0;
    frame = 0;
    final_time = end_frame * frame_dt;

    if (save_sim){
        saveInfo();
        saveParticleData();
        saveForces();
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
        if (interrupt_check && interrupt_check())
            return;
        time += dt;
        current_time_step++;
        if(frame_dt*(frame+1) - time < min_dt*1.1){
            frame++;
            std::cout << "End of frame " << frame << std::endl;
            if (save_sim){
                saveParticleData();
                saveForces();
                if (save_grid && !use_sparse) // TODO: save sparse grid data has not been implemented
                    saveGridData();
                if (save_avg)
                    saveAvgData();
            }
        }
        if (std::abs(final_time-time) < min_dt*1.1 || final_time < time){
            std::cout << "The simulation ended at time = " << time << std::endl;
            break;
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    runtime_total = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "Simulation took " << runtime_total << " milliseconds" << std::endl;
    debug("Runtime P2G     = ", runtime_p2g     * 1000.0, " milliseconds");
    debug("Runtime G2P     = ", runtime_g2p     * 1000.0, " milliseconds");
    debug("Runtime Euler   = ", runtime_euler   * 1000.0, " milliseconds");
    if (use_musl)
        debug("Runtime DefGrad = ", runtime_defgrad * 1000.0, " milliseconds");

    if (save_sim)
        saveTiming();
}




void Simulation::advanceStep(){

    for (auto& obj: objects) obj->force.setZero();
    for (auto& obj: plates) obj->force.setZero();

    updateDt();

    if (pbc){ // TODO: not yet supported by sparse grid
        if (current_time_step == 0)
            remeshFixed(4);
    }
    else{
        if (!use_sparse) {
            if (current_time_step == 0) {
                remeshFixedInit(2,2,2);
            } else {
                remeshFixedCont();
            }
        } 
        else {
            if (current_time_step==0) {
                getParticlesMinMax();
            }

            // Reset the dense scan domain
            int B = scan_sparse_grid.get_B();
            BlockScanGrid::I3 new_bmin{BlockScanGrid::floor_div(particles.minx_id, B)-1, BlockScanGrid::floor_div(particles.miny_id, B)-1, BlockScanGrid::floor_div(particles.minz_id, B)-1};
            BlockScanGrid::I3 new_bmax{BlockScanGrid::floor_div(particles.maxx_id, B)+2, BlockScanGrid::floor_div(particles.maxy_id, B)+2, BlockScanGrid::floor_div(particles.maxz_id, B)+2};
            scan_sparse_grid.init(B, new_bmin, new_bmax, n_threads);

            // scan-based approach to mark active blocks and construct sparse grid
            markActiveBlocksScan();
            resetSparseGridScan();
        }
    }

    if (pbc){ // TODO: not yet supported by sparse grid
        PBCAddParticles(4);
    }

    if (!use_sparse) {
        resizeGrid();
    }

    timer t_p2g; t_p2g.start();
    if (!use_sparse) {
        P2G();
    } 
    else {
        P2GSparseScan();
    }
    t_p2g.stop(); runtime_p2g += t_p2g.get_timing();

    timer t_euler; t_euler.start();
    if (!use_sparse) {
        explicitEulerUpdate();
    } 
    else {
        explicitEulerUpdateSparseScan();
    }
    t_euler.stop(); runtime_euler += t_euler.get_timing();

    if (pbc){
        PBCDelParticles();
    }

    timer t_g2p; t_g2p.start();
    if (!use_sparse) {
        G2P();
    } 
    else {
        G2PSparseScan();
    }
    t_g2p.stop(); runtime_g2p += t_g2p.get_timing();

    if (use_musl){ // TODO: not yet supported by sparse grid
        MUSL();

        timer t_defgrad; t_defgrad.start();
        deformationUpdate();
        t_defgrad.stop(); runtime_defgrad += t_defgrad.get_timing();
    }

    if (!use_sparse) {
        positionUpdate();
    }

    moveObjects();

} // end advanceStep



void Simulation::moveObjects(){
    if ( !gravity_special || (gravity_special && time >= gravity_time) ){
        for (auto& obj : plates) {
            obj->move(dt, frame_dt, time);
        }
        for (auto& obj : objects) {
            obj->move(time);
        }
    } // endif
} // end moveObjects


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
