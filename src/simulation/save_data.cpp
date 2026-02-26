// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"
#include "../../deps/tinyply.h"
#include "../tools.hpp"

void Simulation::saveParticleData(std::string extra){

    std::vector<T> pressure_vec(Np);
    std::vector<T> devstress_vec(Np);
    std::vector<T> Je_vec(Np);

    #pragma omp parallel for num_threads(n_threads)
    for(int p = 0; p < Np; p++){

        TM Fe = particles.F[p];

        // ALT 1: Compute stress
        // TM tau;
        // if (elastic_model == ElasticModel::NeoHookean)
        //     tau = NeoHookeanPiola(Fe) * Fe.transpose();
        // else if (elastic_model == ElasticModel::Hencky)
        //     tau = HenckyPiola(Fe) * Fe.transpose();

        // ALT 2: Use the saved stress
        TM tau = particles.tau[p];

        T Je = Fe.determinant();

        T pressure  = -tau.trace() / dim;
        TM tau_dev = tau + pressure * TM::Identity();
        T devstress = q_prefac * std::sqrt(selfDoubleDot(tau_dev));

        pressure_vec[p] = pressure;
        devstress_vec[p] = devstress;
        Je_vec[p] = Je;
    }


    std::string filename = directory + sim_name + "/particles_f" + extra + std::to_string(frame) + ".ply";
    std::ofstream out;
    out.open(filename, std::ios::out | std::ios::binary);
    tinyply::PlyFile file;
    auto type = std::is_same<T, float>::value ? tinyply::Type::FLOAT32 : tinyply::Type::FLOAT64;

    #ifdef THREEDIM
    file.add_properties_to_element(
        "vertex",
        { "x", "y", "z"},
        type,
        particles.x.size(),
        reinterpret_cast<uint8_t*>(particles.x.data()),
        tinyply::Type::INVALID,
        0);
    file.add_properties_to_element(
        "vertex",
        { "vx", "vy", "vz" },
        type,
        particles.v.size(),
        reinterpret_cast<uint8_t*>(particles.v.data()),
        tinyply::Type::INVALID,
        0);
    #else
    file.add_properties_to_element(
        "vertex",
        { "x", "y"},
        type,
        particles.x.size(),
        reinterpret_cast<uint8_t*>(particles.x.data()),
        tinyply::Type::INVALID,
        0);

    file.add_properties_to_element(
        "vertex",
        { "vx", "vy"},
        type,
        particles.v.size(),
        reinterpret_cast<uint8_t*>(particles.v.data()),
        tinyply::Type::INVALID,
        0);
    #endif

        file.add_properties_to_element(
        "vertex",
        { "p" },
        type,
        pressure_vec.size(),
        reinterpret_cast<uint8_t*>(pressure_vec.data()),
        tinyply::Type::INVALID,
        0);

    file.add_properties_to_element(
        "vertex",
        { "q" },
        type,
        devstress_vec.size(),
        reinterpret_cast<uint8_t*>(devstress_vec.data()),
        tinyply::Type::INVALID,
        0);

    file.add_properties_to_element(
        "vertex",
        { "Je" },
        type,
        Je_vec.size(),
        reinterpret_cast<uint8_t*>(Je_vec.data()),
        tinyply::Type::INVALID,
        0);

    if (plastic_model != PlasticModel::NoPlasticity){

        file.add_properties_to_element(
            "vertex",
            { "eps_pl_vol" },
            type,
            particles.eps_pl_vol.size(),
            reinterpret_cast<uint8_t*>(particles.eps_pl_vol.data()),
            tinyply::Type::INVALID,
            0);

        file.add_properties_to_element(
            "vertex",
            { "eps_pl_dev" },
            type,
            particles.eps_pl_dev.size(),
            reinterpret_cast<uint8_t*>(particles.eps_pl_dev.data()),
            tinyply::Type::INVALID,
            0);

        file.add_properties_to_element(
            "vertex",
            { "delta_gamma" },
            type,
            particles.delta_gamma.size(),
            reinterpret_cast<uint8_t*>(particles.delta_gamma.data()),
            tinyply::Type::INVALID,
            0);

        if ( (plastic_model == PlasticModel::DPVisc && use_mibf) || (plastic_model == PlasticModel::MCCVisc && use_mibf) || plastic_model == PlasticModel::DPMui || plastic_model == PlasticModel::MCCMui){
            file.add_properties_to_element(
                "vertex",
                { "muI" },
                type,
                particles.muI.size(),
                reinterpret_cast<uint8_t*>(particles.muI.data()),
                tinyply::Type::INVALID,
                0);
        } // end if MIBF models

        if (plastic_model == PlasticModel::DPMui || plastic_model == PlasticModel::MCCMui){

            file.add_properties_to_element(
                "vertex",
                { "viscosity" },
                type,
                particles.viscosity.size(),
                reinterpret_cast<uint8_t*>(particles.viscosity.data()),
                tinyply::Type::INVALID,
                0);
        } // end if MuI models

    } // end if plasticity

    file.write(out, true);

    std::ofstream outFile(directory + sim_name + "/last_written.txt");
    outFile << std::to_string(frame) << "\n";
    outFile.close();

} // end saveParticleData()

void Simulation::saveGridData(std::string extra){

    std::string filename = directory + sim_name + "/grid_f" + extra + std::to_string(frame) + ".ply";
    std::ofstream out;
    out.open(filename, std::ios::out | std::ios::binary);
    tinyply::PlyFile file;
    auto type = std::is_same<T, float>::value ? tinyply::Type::FLOAT32 : tinyply::Type::FLOAT64;

    int counter = 0;
    std::vector<TV> grid_x_save;
    grid_x_save.resize(grid_nodes);
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            #ifdef THREEDIM
            for(int k = 0; k < Nz; k++){
                TV grid_x_node(grid.x[i],grid.y[j],grid.z[k]);
            #else
                TV grid_x_node(grid.x[i],grid.y[j]);
            #endif
                grid_x_save[counter] = grid_x_node;
                counter++;
            #ifdef THREEDIM
            }
            #endif
        }
    }

    #ifdef THREEDIM
    file.add_properties_to_element(
        "vertex",
        { "x", "y", "z"},
        type,
        grid_x_save.size(),
        reinterpret_cast<uint8_t*>(grid_x_save.data()),
        tinyply::Type::INVALID,
        0);
    file.add_properties_to_element(
        "vertex",
        { "vx", "vy", "vz" },
        type,
        grid.v.size(),
        reinterpret_cast<uint8_t*>(grid.v.data()),
        tinyply::Type::INVALID,
        0);
    #else
    file.add_properties_to_element(
        "vertex",
        { "x", "y"},
        type,
        grid_x_save.size(),
        reinterpret_cast<uint8_t*>(grid_x_save.data()),
        tinyply::Type::INVALID,
        0);

    file.add_properties_to_element(
        "vertex",
        { "vx", "vy"},
        type,
        grid.v.size(),
        reinterpret_cast<uint8_t*>(grid.v.data()),
        tinyply::Type::INVALID,
        0);
    #endif

    file.add_properties_to_element(
        "vertex",
        { "mass" },
        type,
        grid.mass.size(),
        reinterpret_cast<uint8_t*>(grid.mass.data()),
        tinyply::Type::INVALID,
        0);

    if (use_mibf){
        file.add_properties_to_element(
            "vertex",
            { "friction" },
            type,
            grid.friction.size(),
            reinterpret_cast<uint8_t*>(grid.friction.data()),
            tinyply::Type::INVALID,
            0);
    }

    file.write(out, true);

} // end saveGridData()




void Simulation::computeAvgData(TM& volavg_cauchy, TM& volavg_kirchh, T& Javg){

    volavg_cauchy = TM::Zero();
    volavg_kirchh = TM::Zero();

    T Jsum = 0;
    for(int p = 0; p < Np; p++){

        TM Fe = particles.F[p];

        // ALT 1: Compute stress
        // TM tau;
        // if (elastic_model == ElasticModel::NeoHookean)
        //     tau = NeoHookeanPiola(Fe) * Fe.transpose();
        // else if (elastic_model == ElasticModel::Hencky)
        //     tau = HenckyPiola(Fe) * Fe.transpose();

        // ALT 2: Use the saved stress
        TM tau = particles.tau[p];

        T Je = Fe.determinant();
        T J = Je * std::exp( particles.eps_pl_vol[p] );
        Jsum += J;

        volavg_cauchy += tau;
        volavg_kirchh += tau * J;
    }

    volavg_cauchy /= Jsum;
    volavg_kirchh /= Jsum;
    Javg = Jsum / Np;

}


void Simulation::saveAvgData(){

    TM volavg_cauchy = TM::Zero();
    TM volavg_kirchh = TM::Zero();
    T Javg;

    computeAvgData(volavg_cauchy, volavg_kirchh, Javg);

    std::ofstream outFile1(directory + sim_name + "/avg_cauchy_frame_" + std::to_string(frame) + ".csv");
    #ifdef THREEDIM
        outFile1 << volavg_cauchy(0,0)    << ","
                 << volavg_cauchy(0,1)    << ","
                 << volavg_cauchy(0,2)    << ","
                 << volavg_cauchy(1,0)    << ","
                 << volavg_cauchy(1,1)    << ","
                 << volavg_cauchy(1,2)    << ","
                 << volavg_cauchy(2,0)    << ","
                 << volavg_cauchy(2,1)    << ","
                 << volavg_cauchy(2,2)    << "\n";
    #else
        outFile1 << volavg_cauchy(0,0)    << ","
                 << volavg_cauchy(0,1)    << ","
                 << volavg_cauchy(1,0)    << ","
                 << volavg_cauchy(1,1)    << "\n";
    #endif
    outFile1.close();

    std::ofstream outFile2(directory + sim_name + "/avg_kirchh_frame_" + std::to_string(frame) + ".csv");
    #ifdef THREEDIM
        outFile2 << volavg_kirchh(0,0)    << ","
                 << volavg_kirchh(0,1)    << ","
                 << volavg_kirchh(0,2)    << ","
                 << volavg_kirchh(1,0)    << ","
                 << volavg_kirchh(1,1)    << ","
                 << volavg_kirchh(1,2)    << ","
                 << volavg_kirchh(2,0)    << ","
                 << volavg_kirchh(2,1)    << ","
                 << volavg_kirchh(2,2)    << "\n";
    #else
        outFile2 << volavg_kirchh(0,0)    << ","
                 << volavg_kirchh(0,1)    << ","
                 << volavg_kirchh(1,0)    << ","
                 << volavg_kirchh(1,1)    << "\n";
    #endif
    outFile2.close();

    std::ofstream outFile3(directory + sim_name + "/avg_J_frame_" + std::to_string(frame) + ".csv");
    outFile3 << Javg  << "\n";
    outFile3.close();

    std::ofstream outFile4(directory + sim_name + "/last_saved_frame.txt");
    outFile4 << std::to_string(frame) << "\n";
    outFile4.close();

}


void Simulation::saveInfo(){
    
    std::ofstream infoFile(directory + sim_name + "/info.txt");
    infoFile << end_frame           << "\n"   // 0
             << fps                 << "\n"   // 1
             << dx                  << "\n"   // 2
             << Np                  << "\n"   // 3
             << particle_volume     << "\n";  // 4
    infoFile.close();
}

void Simulation::saveTiming(){

    std::ofstream timeFile(directory + sim_name + "/info_timing.txt");
    timeFile << current_time_step          << "\n"   // 0
             << runtime_total / ((T)1000)  << "\n"   // 1
             << runtime_p2g                << "\n"   // 2
             << runtime_g2p                << "\n"   // 3
             << runtime_euler              << "\n"   // 4
             << runtime_defgrad            << "\n";  // 5
    timeFile.close();
}

