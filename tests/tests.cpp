// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include <gtest/gtest.h>

#include "../src/tools.hpp"
#include "../src/simulation/simulation.hpp"
#include "../src/sampling/sampling_particles.hpp"

#include "../src/objects/object_curve.hpp"
#include "../src/objects/object_ground.hpp"
#include "../src/objects/object_ground_rotated.hpp"

#include <random>


TEST(BoundaryTest, AnalyticSlipFree) {

    Simulation sim;
    sim.initialize(false);
    sim.reduce_verbose = true;

    T half_rounds = 1;
    int num_frames_in_half_round = 100;
    sim.end_frame = half_rounds*num_frames_in_half_round;
    sim.fps = num_frames_in_half_round / 0.5949238427;

    sim.n_threads = 1;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1000;

    sim.gravity[1] = -9.81;

    sim.dx = 0.001;
    sim.particle_volume = std::pow(sim.dx, sim.dim);
    sim.particle_mass = sim.rho * sim.particle_volume;
    sim.Np = 1;
    sim.particles = Particles(sim.Np);

    sim.particles.x[0](0) = -1;
    sim.particles.x[0](1) = 1;

    sim.objects.push_back(std::make_unique<ObjectCurve>(BC::SlipFree, 0.0)); 

    sim.simulate();

    T v_sim = sim.particles.v[0](0);
    T v_true = std::sqrt(2*9.81);
    T diff = std::abs(v_sim-v_true)/v_true;
    debug("diff: ", diff);
    ASSERT_NEAR(diff, 0.0, 1e-3);
}


TEST(BoundaryTest, AnalyticSlipStick) {

    Simulation sim;
    sim.initialize(false);
    sim.reduce_verbose = true;

    T half_rounds = 1;
    int num_frames_in_half_round = 100;
    sim.end_frame = half_rounds*num_frames_in_half_round;
    sim.fps = num_frames_in_half_round / 0.5949238427;

    sim.n_threads = 1;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1000;

    sim.gravity[1] = -9.81;

    sim.dx = 0.001;
    sim.particle_volume = std::pow(sim.dx, sim.dim);
    sim.particle_mass = sim.rho * sim.particle_volume;
    sim.Np = 1;
    sim.particles = Particles(sim.Np);

    sim.particles.x[0](0) = -1;
    sim.particles.x[0](1) = 1;

    sim.objects.push_back(std::make_unique<ObjectCurve>(BC::SlipStick, 0.0)); 

    sim.simulate();

    T v_sim = sim.particles.v[0](0);
    T v_true = std::sqrt(2*9.81);
    T diff = std::abs(v_sim-v_true)/v_true;
    debug("diff: ", diff);
    ASSERT_NEAR(diff, 0.0, 1e-3);
}

TEST(CoulombFrictionTest, PlateSlipFree) {

    Simulation sim;
    sim.initialize(false);
    sim.reduce_verbose = true;

    sim.end_frame = 5;
    sim.fps = 1;
    sim.n_threads = 1;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    T theta_deg = 24;
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero();
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;
    sim.E = 1e6; 
    sim.nu = 0.3; 
    sim.rho = 1000; 

    sim.Lx = 0.1;
    sim.Ly = 0.05;
    #ifdef THREEDIM
        sim.Lz = 0.05;
    #endif
    T k_rad = 0.005;
    sampleParticles(sim, k_rad);
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) += 0.5*sim.dx;
    }

    T friction = std::tan(15.0 * M_PI / 180.0);

    sim.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::SlipFree, friction)); 

    sim.simulate();

    T mean_sim_x;
    for(int p = 0; p < sim.Np; p++)
        mean_sim_x += sim.particles.x[p](0);
    mean_sim_x /= (T)sim.Np;

    T final_time = (T)sim.end_frame / sim.fps;
    T mean_true_x = 0.5*9.81*final_time*final_time*(std::sin(theta) - friction*std::cos(theta));
    T diff = std::abs(mean_sim_x-mean_true_x) / mean_true_x;
    debug("diff: ", diff);
    ASSERT_NEAR(diff, 0.0, 1e-3);
}

TEST(CoulombFrictionTest, PlateSlipStick) {

    Simulation sim;
    sim.initialize(false);
    sim.reduce_verbose = true;

    sim.end_frame = 5;
    sim.fps = 1;
    sim.n_threads = 1;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    T theta_deg = 24;
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero();
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;
    sim.E = 1e6; 
    sim.nu = 0.3; 
    sim.rho = 1000; 

    sim.Lx = 0.1;
    sim.Ly = 0.05;
    #ifdef THREEDIM
        sim.Lz = 0.05;
    #endif
    T k_rad = 0.005;
    sampleParticles(sim, k_rad);
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) += 0.5*sim.dx;
    }

    T friction = std::tan(15.0 * M_PI / 180.0);

    sim.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::SlipStick, friction)); 

    sim.simulate();

    T mean_sim_x;
    for(int p = 0; p < sim.Np; p++)
        mean_sim_x += sim.particles.x[p](0);
    mean_sim_x /= (T)sim.Np;

    T final_time = (T)sim.end_frame / sim.fps;
    T mean_true_x = 0.5*9.81*final_time*final_time*(std::sin(theta) - friction*std::cos(theta));
    T diff = std::abs(mean_sim_x-mean_true_x) / mean_true_x;
    debug("diff: ", diff);
    ASSERT_NEAR(diff, 0.0, 1e-3);
}

TEST(CoulombFrictionTest, GeneralSlipFree) {

    Simulation sim;
    sim.initialize(false);
    sim.reduce_verbose = true;

    sim.end_frame = 5;
    sim.fps = 1;
    sim.n_threads = 1;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    T theta_deg = 24;
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero();
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;
    sim.E = 1e6; 
    sim.nu = 0.3; 
    sim.rho = 1000; 

    sim.Lx = 0.1;
    sim.Ly = 0.05;
    #ifdef THREEDIM
        sim.Lz = 0.05;
    #endif
    T k_rad = 0.005;
    sampleParticles(sim, k_rad);
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) += 0.5*sim.dx;
    }

    T friction = std::tan(15.0 * M_PI / 180.0);

    sim.objects.push_back(std::make_unique<ObjectGround>(BC::SlipFree, friction)); 

    sim.simulate();

    T mean_sim_x;
    for(int p = 0; p < sim.Np; p++)
        mean_sim_x += sim.particles.x[p](0);
    mean_sim_x /= (T)sim.Np;

    T final_time = (T)sim.end_frame / sim.fps;
    T mean_true_x = 0.5*9.81*final_time*final_time*(std::sin(theta) - friction*std::cos(theta));
    T diff = std::abs(mean_sim_x-mean_true_x) / mean_true_x;
    debug("diff: ", diff);
    ASSERT_NEAR(diff, 0.0, 1e-3);
}

TEST(CoulombFrictionTest, GeneralSlipStick) {

    Simulation sim;
    sim.initialize(false);
    sim.reduce_verbose = true;

    sim.end_frame = 5;
    sim.fps = 1;
    sim.n_threads = 1;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    T theta_deg = 24;
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero();
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;
    sim.E = 1e6; 
    sim.nu = 0.3; 
    sim.rho = 1000; 

    sim.Lx = 0.1;
    sim.Ly = 0.05;
    #ifdef THREEDIM
        sim.Lz = 0.05;
    #endif
    T k_rad = 0.005;
    sampleParticles(sim, k_rad);
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) += 0.5*sim.dx;
    }

    T friction = std::tan(15.0 * M_PI / 180.0);

    sim.objects.push_back(std::make_unique<ObjectGround>(BC::SlipStick, friction)); 

    sim.simulate();

    T mean_sim_x;
    for(int p = 0; p < sim.Np; p++)
        mean_sim_x += sim.particles.x[p](0);
    mean_sim_x /= (T)sim.Np;

    T final_time = (T)sim.end_frame / sim.fps;
    T mean_true_x = 0.5*9.81*final_time*final_time*(std::sin(theta) - friction*std::cos(theta));
    T diff = std::abs(mean_sim_x-mean_true_x) / mean_true_x;
    debug("diff: ", diff);
    ASSERT_NEAR(diff, 0.0, 1e-3);
}

TEST(BoundaryTest, MIBF) {

    T friction = std::tan(30.0 * M_PI / 180.0);
    T theta = 32 * M_PI / 180;
        
    Simulation sim_one;
    sim_one.initialize(false);

    sim_one.save_grid = true;
    sim_one.end_frame = 1;
    sim_one.fps = 2;
    sim_one.n_threads = 8;   
    sim_one.cfl = 0.5;     
    sim_one.flip_ratio = -0.95; 

    sim_one.gravity = TV::Zero();
    sim_one.gravity[0] = +9.81 * std::sin(theta);
    sim_one.gravity[1] = -9.81 * std::cos(theta);

    sim_one.elastic_model = ElasticModel::Hencky;
    sim_one.E = 1e5;     
    sim_one.nu = 0.3;   
    sim_one.rho = 1000; 

    sim_one.Lx = 0.1;
    sim_one.Ly = 0.05;
    #ifdef THREEDIM
        sim_one.Lz = 0.05;
    #endif
    sampleParticles(sim_one, 0.001);
    for(int p = 0; p < sim_one.Np; p++){
        sim_one.particles.x[p](0) -= 0.5*sim_one.Lx;
        sim_one.particles.x[p](1) += 0.5*sim_one.dx;
    }

    sim_one.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::SlipFree, friction)); 

    sim_one.plastic_model = PlasticModel::DPVisc; 

    sim_one.use_pradhana = false; 
    sim_one.use_mibf = false;       // NB

    sim_one.M = friction;
    sim_one.q_cohesion = 0;
    sim_one.visc_exponent = 1;
    sim_one.visc_time = 0;

    sim_one.simulate();

    auto max_x_it_1 = std::max_element( sim_one.particles.x.begin(), sim_one.particles.x.end(), [](const TV &x1, const TV &x2){return x1(0) < x2(0);} );
    T max_x_1 = (*max_x_it_1)(0);




    Simulation sim_two;
    sim_two.initialize(false);

    sim_two.save_grid = true;
    sim_two.end_frame = 1;
    sim_two.fps = 2;
    sim_two.n_threads = 8;   
    sim_two.cfl = 0.5;     
    sim_two.flip_ratio = -0.95; 

    sim_two.gravity = TV::Zero();
    sim_two.gravity[0] = +9.81 * std::sin(theta);
    sim_two.gravity[1] = -9.81 * std::cos(theta);

    sim_two.elastic_model = ElasticModel::Hencky;
    sim_two.E = 1e5;     
    sim_two.nu = 0.3;   
    sim_two.rho = 1000; 

    sim_two.Lx = 0.1;
    sim_two.Ly = 0.05;
    #ifdef THREEDIM
        sim_two.Lz = 0.05;
    #endif
    sampleParticles(sim_two, 0.001);
    for(int p = 0; p < sim_two.Np; p++){
        sim_two.particles.x[p](0) -= 0.5*sim_two.Lx;
        sim_two.particles.x[p](1) += 0.5*sim_two.dx;
    }

    sim_two.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::SlipFree, friction)); 

    sim_two.plastic_model = PlasticModel::DPVisc; 

    sim_two.use_pradhana = false; 
    sim_two.use_mibf = true;       // NB       

    sim_two.M = friction;
    sim_two.q_cohesion = 0;
    sim_two.visc_exponent = 1;
    sim_two.visc_time = 0;

    sim_two.simulate();

    auto max_x_it_2 = std::max_element( sim_two.particles.x.begin(), sim_two.particles.x.end(), [](const TV &x1, const TV &x2){return x1(0) < x2(0);} );
    T max_x_2 = (*max_x_it_2)(0);

    T diff = std::abs(max_x_1 - max_x_2);
    debug(diff);
    ASSERT_NEAR(diff, 0.0, 1e-12);




    Simulation sim_three;
    sim_three.initialize(false);

    sim_three.save_grid = true;
    sim_three.end_frame = 1;
    sim_three.fps = 2;
    sim_three.n_threads = 8;   
    sim_three.cfl = 0.5;     
    sim_three.flip_ratio = -0.95; 

    sim_three.gravity = TV::Zero();
    sim_three.gravity[0] = +9.81 * std::sin(theta);
    sim_three.gravity[1] = -9.81 * std::cos(theta);

    sim_three.elastic_model = ElasticModel::Hencky;
    sim_three.E = 1e5;     
    sim_three.nu = 0.3;   
    sim_three.rho = 1000; 

    sim_three.Lx = 0.1;
    sim_three.Ly = 0.05;
    #ifdef THREEDIM
        sim_three.Lz = 0.05;
    #endif
    sampleParticles(sim_three, 0.001);
    for(int p = 0; p < sim_three.Np; p++){
        sim_three.particles.x[p](0) -= 0.5*sim_three.Lx;
        sim_three.particles.x[p](1) += 0.5*sim_three.dx;
    }

    sim_three.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::SlipFree, friction)); 

    sim_three.plastic_model = PlasticModel::DPVisc; 

    sim_three.use_pradhana = false; 
    sim_three.use_mibf = false;                // NB
    sim_three.q_prefac = std::sqrt(3.0/2.0);   // NB


    sim_three.M = std::sqrt(3.0)*friction;     // NB
    sim_three.q_cohesion = 0;
    sim_three.visc_exponent = 1;
    sim_three.visc_time = 0;

    sim_three.simulate();

    auto max_x_it_3 = std::max_element( sim_three.particles.x.begin(), sim_three.particles.x.end(), [](const TV &x1, const TV &x2){return x1(0) < x2(0);} );
    T max_x_3 = (*max_x_it_3)(0);

    T difff = std::abs(max_x_1 - max_x_3);
    debug(difff);
    ASSERT_NEAR(difff, 0.0, 1e-12);

#ifdef THREEDIM
    Simulation sim_four;
    sim_four.initialize(false);

    sim_four.use_sparse = true; // NB

    sim_four.save_grid = true;
    sim_four.end_frame = 1;
    sim_four.fps = 2;
    sim_four.n_threads = 8;   
    sim_four.cfl = 0.5;     
    sim_four.flip_ratio = -0.95; 

    sim_four.gravity = TV::Zero();
    sim_four.gravity[0] = +9.81 * std::sin(theta);
    sim_four.gravity[1] = -9.81 * std::cos(theta);

    sim_four.elastic_model = ElasticModel::Hencky;
    sim_four.E = 1e5;     
    sim_four.nu = 0.3;   
    sim_four.rho = 1000; 

    sim_four.Lx = 0.1;
    sim_four.Ly = 0.05;
    
        sim_four.Lz = 0.05;
    
    sampleParticles(sim_four, 0.001);
    for(int p = 0; p < sim_four.Np; p++){
        sim_four.particles.x[p](0) -= 0.5*sim_four.Lx;
        sim_four.particles.x[p](1) += 0.5*sim_four.dx;
    }

    sim_four.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::SlipFree, friction)); 

    sim_four.plastic_model = PlasticModel::DPVisc; 

    sim_four.use_pradhana = false; 
    sim_four.use_mibf = false;             

    sim_four.M = friction;
    sim_four.q_cohesion = 0;
    sim_four.visc_exponent = 1;
    sim_four.visc_time = 0;

    sim_four.simulate();

    auto max_x_it_4 = std::max_element( sim_four.particles.x.begin(), sim_four.particles.x.end(), [](const TV &x1, const TV &x2){return x1(0) < x2(0);} );
    T max_x_4 = (*max_x_it_4)(0);

    T diffff = std::abs(max_x_1 - max_x_4);
    debug(diffff);
    ASSERT_NEAR(diffff, 0.0, 1e-12);
#endif

}


TEST(ForceTest, NoSlipGround) {

    Simulation sim;
    sim.initialize(false);

    sim.end_frame = 10;
    sim.fps = 10;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = 0;

    sim.gravity[1] = -9.81;

    sim.Lx = 1;
    sim.Ly = 1;
    #ifdef THREEDIM
        sim.Lz = 0.1;
    #endif

    T friction = 0;
    T y_ground = 0;
    std::string name = "ground";

    sim.objects.push_back(std::make_unique<ObjectGround>(BC::NoSlip, friction, name, true, y_ground));

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;

    sim.E = 1e10;
    sim.nu = 0.3;
    sim.rho = 400;

    sampleParticles(sim, 0.2);

    sim.simulate();

    TV obj_force = TV::Zero();
    for (auto& obj: sim.objects){
        obj_force = obj->force;
    }

    #ifdef THREEDIM
        T exp_force = sim.rho * sim.Lx * sim.Ly * sim.Lz * sim.gravity[1];
    #else
        T exp_force = sim.rho * sim.Lx * sim.Ly * sim.gravity[1];
    #endif

    T com_force = obj_force(1);
    T diff = std::abs(exp_force - com_force)/exp_force;
    debug("diff: ", diff);

    ASSERT_NEAR(diff, 0.0, 1e-3);
}

TEST(ForceTest, SlipFreeRotatedGround) {

    Simulation sim;
    sim.initialize(false);

    sim.end_frame = 10;
    sim.fps = 10;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = 0;

    T theta_deg = 10.0;
    T theta = theta_deg * (M_PI / 180.0);
    T cos_theta = std::cos(theta);
    T sin_theta = std::sin(theta);

    T theta_r_deg = 0;
    T theta_r = theta_r_deg * (M_PI / 180.0);

    T friction = std::tan(theta_r);

    sim.gravity[1] = -9.81;

    sim.Lx = 1;
    sim.Ly = 1;
    #ifdef THREEDIM
        sim.Lz = 0.1;
    #endif

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;

    sim.E = 1e10;
    sim.rho = 400;
    sim.nu = 0.3;

    sampleParticles(sim, 0.1);

    for(int p = 0; p < sim.Np; p++){

        T x_new = sim.particles.x[p](1) * cos_theta - sim.particles.x[p](0) * sin_theta;
        T z_new = sim.particles.x[p](1) * sin_theta + sim.particles.x[p](0) * cos_theta;

        sim.particles.x[p](1) = x_new;
        sim.particles.x[p](0) = z_new;

    }

    sim.objects.push_back(std::make_unique<ObjectGroundRotated>(BC::SlipFree, friction, "ground_rotated", true, theta));

    sim.simulate();

    TV obj_force = TV::Zero();
    for (auto& obj: sim.objects){
        obj_force = obj->force;
    }

    #ifdef THREEDIM
        T exp_force_1 = (sim.rho * sim.Lx * sim.Ly * sim.Lz * sim.gravity[1] * cos_theta) * (sin_theta - friction * cos_theta);
        T exp_force_2 = (sim.rho * sim.Lx * sim.Ly * sim.Lz * sim.gravity[1] * cos_theta) * (cos_theta + friction * sin_theta);
    #else
        T exp_force_1 = (sim.rho * sim.Lx * sim.Ly * sim.gravity[1] * cos_theta) * (sin_theta - friction * cos_theta);
        T exp_force_2 = (sim.rho * sim.Lx * sim.Ly * sim.gravity[1] * cos_theta) * (cos_theta + friction * sin_theta);
    #endif

    T com_force_1 = obj_force(0);
    T com_force_2 = obj_force(1);

    T diff_1 = std::abs(exp_force_1 - com_force_1) / exp_force_1;
    T diff_2 = std::abs(exp_force_2 - com_force_2) / exp_force_2;

    debug("diff_1: ", diff_1);
    debug("diff_2: ", diff_2);
    ASSERT_NEAR(diff_1, 0.0, 1.2e-3);
    ASSERT_NEAR(diff_2, 0.0, 1.2e-3);
}

TEST(ForceTest, NoSlipPlate) {

    Simulation sim;
    sim.initialize(false);

    sim.end_frame = 10;
    sim.fps = 10;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = 0;

    sim.gravity[1] = -9.81;

    sim.Lx = 1;
    sim.Ly = 1;
    #ifdef THREEDIM
        sim.Lz = 0.1;
    #endif

    T friction = 0;

    #ifdef THREEDIM
        sim.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::NoSlip, friction, -1, 2, 0.0, 0.0, 0.0, 0.0, 0.0, "bottom_plate", true));
    #else
        sim.plates.push_back(std::make_unique<ObjectPlate>(0, PlateType::bottom, BC::NoSlip, friction, -1, 2, 0.0, 0.0, 0.0, 0.0, "bottom_plate", true));
    #endif

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;

    sim.E = 1e10;
    sim.nu = 0.3;
    sim.rho = 400;

    sampleParticles(sim, 0.2);

    sim.simulate();

    TV plate_force = TV::Zero();
    for (auto& obj: sim.plates){
        plate_force = obj->force;
    }

    #ifdef THREEDIM
        T exp_force = sim.rho * sim.Lx * sim.Ly * sim.Lz * sim.gravity[1];
    #else
        T exp_force = sim.rho * sim.Lx * sim.Ly * sim.gravity[1];
    #endif
    T com_force = plate_force(1);
    T diff = std::abs(exp_force - com_force) / exp_force;
    debug("diff: ", diff);

    ASSERT_NEAR(diff, 0.0, 1e-3);
}


TEST(ElasticityTest, BulkModulus) {

    Simulation sim;
    sim.initialize(false);
    sim.reduce_verbose = true;
    sim.end_frame = 100;
    sim.fps = 1;
    sim.n_threads = 8;
    sim.cfl = 0.6;
    sim.flip_ratio = -0.95;
    sim.gravity = TV::Zero();
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 10000;

    sim.Lx = 1;
    sim.Ly = 1;
    #ifdef THREEDIM
    sim.Lz = 0.2;
        sampleParticles(sim, 0.02, 8);
    #else
        sampleParticles(sim, 0.02, 4);
    #endif

    T vel = 0.001;

    T vmin_factor = 10;
    T load_factor = 1;

    #ifdef THREEDIM
        sim.plates.push_back(std::make_unique<ObjectPlate>(0-0.5*sim.dx,       PlateType::bottom, BC::NoSlip, 0, -1e15, 1e15,   0,  vel, 0, vmin_factor, load_factor));
        sim.plates.push_back(std::make_unique<ObjectPlate>(sim.Ly+0.5*sim.dx,  PlateType::top,    BC::NoSlip, 0, -1e15, 1e15,   0, -vel, 0, vmin_factor, load_factor));
    #else
        sim.plates.push_back(std::make_unique<ObjectPlate>(0-0.5*sim.dx,       PlateType::bottom, BC::NoSlip, 0, -1e15, 1e15,   0,  vel,    vmin_factor, load_factor)); 
        sim.plates.push_back(std::make_unique<ObjectPlate>(sim.Ly+0.5*sim.dx,  PlateType::top,    BC::NoSlip, 0, -1e15, 1e15,   0, -vel,    vmin_factor, load_factor)); 
    #endif

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;

    sim.simulate();

    TM volavg_cauchy = TM::Zero();
    TM volavg_kirchh = TM::Zero();
    T Javg;
    sim.computeAvgData(volavg_cauchy, volavg_kirchh, Javg);

    #ifdef THREEDIM
        T volavg_p = -1.0 * (volavg_kirchh(0,0) + volavg_kirchh(1,1) + volavg_kirchh(2,2)) / 3;
    #else
        T volavg_p = -1.0 * (volavg_kirchh(0,0) + volavg_kirchh(1,1)) / 2;
    #endif

    T volavg_epsv = std::log(Javg);

    T measured_K = volavg_p / (-volavg_epsv);
    T true_K = sim.calculateBulkModulus();

    T rel_diff = std::abs(measured_K - true_K) / true_K;

    debug("true_K:     ", true_K);
    debug("measured_K: ", measured_K);
    debug("rel_diff:   ", rel_diff);

    ASSERT_NEAR(rel_diff, 0.0, 0.03);
}

TEST(ElasticityTest, BulkModulusMUSL) {

    Simulation sim;
    sim.initialize(false);
    sim.reduce_verbose = true;
    sim.end_frame = 100;
    sim.fps = 1;
    sim.n_threads = 8;
    sim.cfl = 0.6;
    sim.flip_ratio = -0.95;
    sim.gravity = TV::Zero();
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 10000;

    sim.use_musl = true; // NB

    sim.Lx = 1;
    sim.Ly = 1;
    #ifdef THREEDIM
    sim.Lz = 0.2;
        sampleParticles(sim, 0.02, 8);
    #else
        sampleParticles(sim, 0.02, 4);
    #endif

    T vel = 0.001;

    T vmin_factor = 10;
    T load_factor = 1;

    #ifdef THREEDIM
        sim.plates.push_back(std::make_unique<ObjectPlate>(0-0.5*sim.dx,       PlateType::bottom, BC::NoSlip, 0, -1e15, 1e15,   0,  vel, 0, vmin_factor, load_factor));
        sim.plates.push_back(std::make_unique<ObjectPlate>(sim.Ly+0.5*sim.dx,  PlateType::top,    BC::NoSlip, 0, -1e15, 1e15,   0, -vel, 0, vmin_factor, load_factor));
    #else
        sim.plates.push_back(std::make_unique<ObjectPlate>(0-0.5*sim.dx,       PlateType::bottom, BC::NoSlip, 0, -1e15, 1e15,   0,  vel,    vmin_factor, load_factor)); 
        sim.plates.push_back(std::make_unique<ObjectPlate>(sim.Ly+0.5*sim.dx,  PlateType::top,    BC::NoSlip, 0, -1e15, 1e15,   0, -vel,    vmin_factor, load_factor)); 
    #endif

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;

    sim.simulate();

    TM volavg_cauchy = TM::Zero();
    TM volavg_kirchh = TM::Zero();
    T Javg;
    sim.computeAvgData(volavg_cauchy, volavg_kirchh, Javg);

    #ifdef THREEDIM
        T volavg_p = -1.0 * (volavg_kirchh(0,0) + volavg_kirchh(1,1) + volavg_kirchh(2,2)) / 3;
    #else
        T volavg_p = -1.0 * (volavg_kirchh(0,0) + volavg_kirchh(1,1)) / 2;
    #endif

    T volavg_epsv = std::log(Javg);

    T measured_K = volavg_p / (-volavg_epsv);
    T true_K = sim.calculateBulkModulus();

    T rel_diff = std::abs(measured_K - true_K) / true_K;

    debug("true_K:     ", true_K);
    debug("measured_K: ", measured_K);
    debug("rel_diff:   ", rel_diff);

    ASSERT_NEAR(rel_diff, 0.0, 0.03);
}

TEST(EnergyTest, Rotation) {

    Simulation sim;
    sim.initialize(false);
    sim.reduce_verbose = true;
    sim.end_frame = 20;
    sim.fps = 1;
    sim.gravity = TV::Zero();
    sim.cfl = 0.5;
    sim.flip_ratio = -1;
    sim.n_threads = 8;
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1550;

    T h_gate, l_gate;
    sim.Lx = 1;
    sim.Ly = 1;
    #ifdef THREEDIM
        sim.Lz = 0.05;
    #endif
    sampleParticles(sim, 0.01);

    T total_energy_init = 0;
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) -= 0.5*sim.Ly;

        T vx = -1.0*sim.particles.x[p](1) + 0.5;
        T vy =  1.0*sim.particles.x[p](0) + 0.5;
        sim.particles.v[p](0) = vx;
        sim.particles.v[p](1) = vy;

        total_energy_init += 0.5*(vx*vx + vy*vy); // per unit mass
    }

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;

    sim.simulate();

    T total_energy_last = 0;
    for(int p = 0; p < sim.Np; p++){
        T vx = sim.particles.v[p](0);
        T vy = sim.particles.v[p](1);
        #ifdef THREEDIM
            T vz = sim.particles.v[p](2);
            total_energy_last += 0.5*(vx*vx + vy*vy + vz*vz); // per unit mass
        #else
            total_energy_last += 0.5*(vx*vx + vy*vy); // per unit mass
        #endif
    }

    T rel_diff = (total_energy_init - total_energy_last) / total_energy_init;

    debug("rel_diff:   ", rel_diff);

    EXPECT_TRUE((rel_diff >= 0) && (rel_diff <= 1e-3));
    // Can not have energy increase!
    // Energy decrease within a relative tolerance
}

TEST(EnergyTest, RotationMUSL) {

    Simulation sim;
    sim.initialize(false);
    sim.reduce_verbose = true;
    sim.end_frame = 20;
    sim.fps = 1;
    sim.gravity = TV::Zero();
    sim.cfl = 0.5;
    sim.flip_ratio = -1;
    sim.n_threads = 8;
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1550;

    sim.use_musl = true; // NB

    T h_gate, l_gate;
    sim.Lx = 1;
    sim.Ly = 1;
    #ifdef THREEDIM
        sim.Lz = 0.05;
    #endif
    sampleParticles(sim, 0.01);

    T total_energy_init = 0;
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) -= 0.5*sim.Ly;

        T vx = -1.0*sim.particles.x[p](1) + 0.5;
        T vy =  1.0*sim.particles.x[p](0) + 0.5;
        sim.particles.v[p](0) = vx;
        sim.particles.v[p](1) = vy;

        total_energy_init += 0.5*(vx*vx + vy*vy); // per unit mass
    }

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::NoPlasticity;

    sim.simulate();

    T total_energy_last = 0;
    for(int p = 0; p < sim.Np; p++){
        T vx = sim.particles.v[p](0);
        T vy = sim.particles.v[p](1);
        #ifdef THREEDIM
            T vz = sim.particles.v[p](2);
            total_energy_last += 0.5*(vx*vx + vy*vy + vz*vz); // per unit mass
        #else
            total_energy_last += 0.5*(vx*vx + vy*vy); // per unit mass
        #endif
    }

    T rel_diff = (total_energy_init - total_energy_last) / total_energy_init;

    debug("rel_diff:   ", rel_diff);

    EXPECT_TRUE((rel_diff >= 0) && (rel_diff <= 1e-3));
    // Can not have energy increase!
    // Energy decrease within a relative tolerance
}

TEST(CollapseTest, DruckerPragerOne) {

    Simulation sim;
    sim.initialize(false);

    sim.reduce_verbose = true;
    sim.end_frame = 70;
    sim.fps = 50;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    T theta_deg = 10;
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero();
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1000;

    sim.Lx = 0.20;
    sim.Ly = 0.15;
    T k_rad = 0.0015;
    #ifdef THREEDIM
        sim.Lz = 0.10;
        sampleParticles(sim, k_rad, 8);
    #else
        sampleParticles(sim, k_rad, 4);
    #endif
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](1) += 0.5*sim.dx;
    }

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::DPSoft;

    sim.use_pradhana = true;

    sim.xi = 0;

    sim.q_cohesion = 0;
    sim.M = std::tan(30.0 * M_PI / 180.0);

    #ifdef THREEDIM
        sim.plates.push_back(std::make_unique<ObjectPlate>(0,       PlateType::back,   BC::SlipFree));
        sim.plates.push_back(std::make_unique<ObjectPlate>(sim.Lz,  PlateType::front,  BC::SlipFree));
        sim.plates.push_back(std::make_unique<ObjectPlate>(0,       PlateType::left,   BC::SlipFree));
        sim.plates.push_back(std::make_unique<ObjectPlate>(0,       PlateType::bottom, BC::NoSlip));
    #else
        sim.plates.push_back(std::make_unique<ObjectPlate>(0,  PlateType::left,   BC::SlipFree));
        sim.plates.push_back(std::make_unique<ObjectPlate>(0,  PlateType::bottom, BC::NoSlip));
    #endif

    sim.simulate();

    auto max_x_it = std::max_element( sim.particles.x.begin(), sim.particles.x.end(), [](const TV &x1, const TV &x2){return x1(0) < x2(0);} );
    T max_x = (*max_x_it)(0);
    T diff = std::abs(max_x - 0.56);
    debug("diff:   ", diff);
    ASSERT_NEAR(diff, 0.0, 0.011);
}

TEST(CollapseTest, DruckerPragerTwo) {

    Simulation sim;
    sim.initialize(false);

    sim.reduce_verbose = true;
    sim.end_frame = 70;
    sim.fps = 50;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    T theta_deg = 10;
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero();
    sim.gravity[0] = +9.81 * std::sin(theta);
    sim.gravity[1] = -9.81 * std::cos(theta);

    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1000;

    sim.Lx = 0.20;
    sim.Ly = 0.15;
    T k_rad = 0.0015;
    #ifdef THREEDIM
        sim.Lz = 0.10;
        sampleParticles(sim, k_rad, 8);
    #else
        sampleParticles(sim, k_rad, 4);
    #endif
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](1) += 0.5*sim.dx;
    }

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::DPVisc;

    sim.use_pradhana = true;

    sim.visc_time = 0;
    sim.visc_exponent = 1;

    sim.q_cohesion = 0;
    sim.M = std::tan(30.0 * M_PI / 180.0);

    #ifdef THREEDIM
        sim.plates.push_back(std::make_unique<ObjectPlate>(0,       PlateType::back,   BC::SlipFree));
        sim.plates.push_back(std::make_unique<ObjectPlate>(sim.Lz,  PlateType::front,  BC::SlipFree));
        sim.plates.push_back(std::make_unique<ObjectPlate>(0,       PlateType::left,   BC::SlipFree));
        sim.plates.push_back(std::make_unique<ObjectPlate>(0,       PlateType::bottom, BC::NoSlip));
    #else
        sim.plates.push_back(std::make_unique<ObjectPlate>(0,  PlateType::left,   BC::SlipFree));
        sim.plates.push_back(std::make_unique<ObjectPlate>(0,  PlateType::bottom, BC::NoSlip));
    #endif

    sim.simulate();

    auto max_x_it = std::max_element( sim.particles.x.begin(), sim.particles.x.end(), [](const TV &x1, const TV &x2){return x1(0) < x2(0);} );
    T max_x = (*max_x_it)(0);
    T diff = std::abs(max_x - 0.56);
    debug("diff:   ", diff);
    ASSERT_NEAR(diff, 0.0, 0.011);
}


template <int D>
using MatT = Eigen::Matrix<T, D, D>;
template <int D>
using VecT = Eigen::Matrix<T, D, 1>;

template <int D>
static void expectSVDInvariants(const MatT<D>& F, const std::string& label) {
    FastSVDImpl<D> svd(F);
    const MatT<D> U = svd.matrixU();
    const MatT<D> V = svd.matrixV();
    const VecT<D> s = svd.singularValues();
    const std::string at = " [D=" + std::to_string(D) + " " + label + "]";

    const T maxabs = F.cwiseAbs().maxCoeff();
    const T nrm = (maxabs > 0 && std::isfinite(maxabs)) ? maxabs : T(1);

    EXPECT_LT((U * s.asDiagonal() * V.transpose() - F).cwiseAbs().maxCoeff() / nrm, 1e-12)
        << "reconstruction" << at;
    EXPECT_LT((U.transpose() * U - MatT<D>::Identity()).cwiseAbs().maxCoeff(), 1e-12)
        << "U orthogonality" << at;
    EXPECT_LT((V.transpose() * V - MatT<D>::Identity()).cwiseAbs().maxCoeff(), 1e-12)
        << "V orthogonality" << at;

    for (int i = 0; i < D; i++)
        EXPECT_GE(s(i), T(0)) << "singular value " << i << " negative" << at;
    for (int i = 0; i + 1 < D; i++)
        EXPECT_GE(s(i), s(i + 1)) << "singular values not decreasing" << at;

    const T detF = F.determinant();
    if (std::isfinite(detF) && std::abs(detF) > T(1e-6) * std::pow(nrm, D))
        EXPECT_NEAR(U.determinant() * V.determinant(), detF > 0 ? T(1) : T(-1), 1e-12)
            << "det(U)det(V) inconsistent with sign(det F)" << at;
}

template <int D>
static MatT<D> randMat(std::mt19937& rng, T amp) {
    std::normal_distribution<T> gauss(0, 1);
    MatT<D> M;
    for (int i = 0; i < D; i++)
        for (int j = 0; j < D; j++) M(i, j) = amp * gauss(rng);
    return M;
}

template <int D>
static MatT<D> randRot(std::mt19937& rng) {
    Eigen::JacobiSVD<MatT<D>> sv(randMat<D>(rng, 1), Eigen::ComputeFullU | Eigen::ComputeFullV);
    return MatT<D>(sv.matrixU() * sv.matrixV().transpose());
}

template <int D>
static void runMatchesEigen() {
    std::mt19937 rng(20240517);
    // Amplitudes spanning near-rigid to heavily distorted deformation gradients.
    for (T amp : {T(0.01), T(0.1), T(0.5), T(1.0)}) {
        for (int trial = 0; trial < 3000; trial++) {
            const MatT<D> F = MatT<D>::Identity() + randMat<D>(rng, amp);
            expectSVDInvariants<D>(F, "random");

            Eigen::JacobiSVD<MatT<D>> ref(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
            FastSVDImpl<D> svd(F);
            const VecT<D> se = ref.singularValues();
            if (se(D - 1) <= T(1e-6) * se(0))
                continue; // not well determined; invariants above still apply

            EXPECT_LT((svd.singularValues() - se).cwiseAbs().maxCoeff() / se(0), 1e-11)
                << "singular values disagree with Eigen [D=" << D << "]";

            // The two quantities plasticity.cpp actually builds from the SVD.
            const VecT<D> h  = svd.singularValues().array().abs().max(1e-4).log();
            const VecT<D> he = se.array().abs().max(1e-4).log();
            EXPECT_LT((svd.matrixU() * h.asDiagonal() * svd.matrixU().transpose()
                       - ref.matrixU() * he.asDiagonal() * ref.matrixU().transpose())
                          .cwiseAbs().maxCoeff(), 1e-9)
                << "Kirchhoff stress disagrees [D=" << D << "]";
            EXPECT_LT((svd.matrixU() * h.array().exp().matrix().asDiagonal() * svd.matrixV().transpose()
                       - ref.matrixU() * he.array().exp().matrix().asDiagonal() * ref.matrixV().transpose())
                          .cwiseAbs().maxCoeff() / se(0), 1e-10)
                << "projected F disagrees [D=" << D << "]";
        }
    }
}

TEST(FastSVDTest, MatchesEigenOnDeformationGradients) {
    runMatchesEigen<2>();
    runMatchesEigen<3>();
}

template <int D>
static void runDegenerate() {
    std::mt19937 rng(11);

    expectSVDInvariants<D>(MatT<D>::Zero(), "zero");
    expectSVDInvariants<D>(MatT<D>::Identity(), "identity");
    expectSVDInvariants<D>(MatT<D>(-MatT<D>::Identity()), "negative identity");

    for (int trial = 0; trial < 1500; trial++) {
        // Pure rotations and reflections: all singular values equal, the worst case
        // for any method relying on eigenvector separation.
        const MatT<D> R = randRot<D>(rng);
        expectSVDInvariants<D>(R, "rotation");
        MatT<D> refl = R;
        refl.col(0) *= -1;
        expectSVDInvariants<D>(refl, "reflection");

        // Rank deficient by construction.
        MatT<D> rd = randMat<D>(rng, 1);
        rd.col(1) = rd.col(0);
        expectSVDInvariants<D>(rd, "rank deficient");

        // Prescribed spectra: repeated, collapsed, and extreme condition numbers.
        VecT<D> sp;
        if constexpr (D == 3) {
            sp = VecT<D>(1, 1, 1);
            expectSVDInvariants<D>(randRot<D>(rng) * sp.asDiagonal() * randRot<D>(rng), "sigma 1,1,1");
            sp = VecT<D>(2, 2, 1);
            expectSVDInvariants<D>(randRot<D>(rng) * sp.asDiagonal() * randRot<D>(rng), "sigma 2,2,1");
            sp = VecT<D>(1, 1, 0);
            expectSVDInvariants<D>(randRot<D>(rng) * sp.asDiagonal() * randRot<D>(rng), "sigma 1,1,0");
            sp = VecT<D>(1, 0, 0);
            expectSVDInvariants<D>(randRot<D>(rng) * sp.asDiagonal() * randRot<D>(rng), "sigma 1,0,0");
            sp = VecT<D>(1e8, 1, 1e-8);
            expectSVDInvariants<D>(randRot<D>(rng) * sp.asDiagonal() * randRot<D>(rng), "ill conditioned");
            sp = VecT<D>(1e-30, 1e-30, 1e-30);
            expectSVDInvariants<D>(randRot<D>(rng) * sp.asDiagonal() * randRot<D>(rng), "tiny");
        } else {
            sp = VecT<D>(1, 1);
            expectSVDInvariants<D>(randRot<D>(rng) * sp.asDiagonal() * randRot<D>(rng), "sigma 1,1");
            sp = VecT<D>(1, 0);
            expectSVDInvariants<D>(randRot<D>(rng) * sp.asDiagonal() * randRot<D>(rng), "sigma 1,0");
            sp = VecT<D>(1e8, 1e-8);
            expectSVDInvariants<D>(randRot<D>(rng) * sp.asDiagonal() * randRot<D>(rng), "ill conditioned");
            sp = VecT<D>(1e-30, 1e-30);
            expectSVDInvariants<D>(randRot<D>(rng) * sp.asDiagonal() * randRot<D>(rng), "tiny");
        }
        // Extreme magnitudes, which the power-of-two prescaling exists to handle.
        expectSVDInvariants<D>(randMat<D>(rng, 1e150),  "huge");
        expectSVDInvariants<D>(randMat<D>(rng, 1e-150), "small");
    }
}

TEST(FastSVDTest, DegenerateAndSingularInputs) {
    runDegenerate<2>();
    runDegenerate<3>();
}

template <int D>
static void runExtremeMagnitudes() {
    // frexp reports exponents outside the representable power-of-two range for
    // subnormal inputs and for inputs near the largest finite double. Before the
    // prescaling exponent was clamped, the first produced NaN (0 * inf) and the
    // second Inf (the undo step overflowing), where Eigen returns finite results.
    // Reconstruction is not asserted: a subnormal input carries only a few
    // significant bits, so its own representation granularity dominates.
    std::mt19937 rng(5);
    std::uniform_real_distribution<T> unif(-1, 1);
    const T mags[] = {T(1e-320), T(1e-310), std::numeric_limits<T>::denorm_min(),
                      std::numeric_limits<T>::min(), T(1e300), T(1e307)};

    for (T mag : mags) {
        for (int trial = 0; trial < 200; trial++) {
            MatT<D> F;
            for (int i = 0; i < D; i++)
                for (int j = 0; j < D; j++) F(i, j) = mag * unif(rng);

            FastSVDImpl<D> svd(F);
            ASSERT_TRUE(svd.singularValues().allFinite()) << "magnitude " << mag << " D=" << D;
            ASSERT_TRUE(svd.matrixU().allFinite())        << "magnitude " << mag << " D=" << D;
            ASSERT_TRUE(svd.matrixV().allFinite())        << "magnitude " << mag << " D=" << D;

            EXPECT_LT((svd.matrixU().transpose() * svd.matrixU() - MatT<D>::Identity())
                          .cwiseAbs().maxCoeff(), 1e-12) << "magnitude " << mag << " D=" << D;
            EXPECT_LT((svd.matrixV().transpose() * svd.matrixV() - MatT<D>::Identity())
                          .cwiseAbs().maxCoeff(), 1e-12) << "magnitude " << mag << " D=" << D;

            Eigen::JacobiSVD<MatT<D>> ref(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
            const VecT<D> se = ref.singularValues();
            if (se(0) > 0)
                EXPECT_LT((svd.singularValues() - se).cwiseAbs().maxCoeff() / se(0), 1e-10)
                    << "magnitude " << mag << " D=" << D;
        }
    }
}

TEST(FastSVDTest, SubnormalAndNearOverflowInputs) {
    runExtremeMagnitudes<2>();
    runExtremeMagnitudes<3>();
}

// Checks the singular values against EXACT ALGEBRA rather than another solver.
// For an integer F, C = F^T F is an exact integer matrix, so the coefficients of
// its characteristic polynomial are exact integers, and the singular values must
// satisfy the elementary symmetric identities
//     sum sigma_i^2 = tr C,   sum_{i<j} sigma_i^2 sigma_j^2 = M,   prod sigma_i^2 = det C
// where M is the sum of the principal 2x2 minors. Those three determine the
// eigenvalue multiset uniquely, so no root finding or reference SVD is involved.
template <int D>
static void runExactAlgebra() {
    std::mt19937 rng(777);
    std::uniform_int_distribution<int> coeff(-10, 10);

    for (int trial = 0; trial < 4000; trial++) {
        MatT<D> F;
        long long Fi[D][D];
        for (int i = 0; i < D; i++)
            for (int j = 0; j < D; j++) { Fi[i][j] = coeff(rng); F(i, j) = T(Fi[i][j]); }

        long long C[D][D];
        for (int i = 0; i < D; i++)
            for (int j = 0; j < D; j++) {
                long long acc = 0;
                for (int k = 0; k < D; k++) acc += Fi[k][i] * Fi[k][j];
                C[i][j] = acc;
            }

        long long tr = 0, minors = 0, det = 0;
        for (int i = 0; i < D; i++) tr += C[i][i];
        if constexpr (D == 2) {
            minors = C[0][0] * C[1][1] - C[0][1] * C[1][0];
        } else {
            minors = (C[0][0]*C[1][1] - C[0][1]*C[1][0]) + (C[0][0]*C[2][2] - C[0][2]*C[2][0])
                   + (C[1][1]*C[2][2] - C[1][2]*C[2][1]);
            det = C[0][0]*(C[1][1]*C[2][2] - C[1][2]*C[2][1])
                - C[0][1]*(C[1][0]*C[2][2] - C[1][2]*C[2][0])
                + C[0][2]*(C[1][0]*C[2][1] - C[1][1]*C[2][0]);
        }

        FastSVDImpl<D> svd(F);
        const VecT<D> s = svd.singularValues();
        long double l[D];
        for (int i = 0; i < D; i++) l[i] = (long double)s(i) * (long double)s(i);

        // Normalize by the natural scale of each identity so tolerances are relative.
        const long double sc = std::max<long double>(1.0L, (long double)tr);
        long double e1 = 0, e2 = 0, e3 = 0;
        for (int i = 0; i < D; i++) e1 += l[i];
        e1 -= (long double)tr;
        if constexpr (D == 2) {
            e2 = l[0]*l[1] - (long double)minors;
        } else {
            e2 = l[0]*l[1] + l[0]*l[2] + l[1]*l[2] - (long double)minors;
            e3 = l[0]*l[1]*l[2] - (long double)det;
        }
        EXPECT_LT(std::abs(e1) / sc,          1e-12) << "trace identity, D=" << D;
        EXPECT_LT(std::abs(e2) / (sc*sc),     1e-12) << "minor identity, D=" << D;
        if constexpr (D == 3)
            EXPECT_LT(std::abs(e3) / (sc*sc*sc), 1e-12) << "determinant identity, D=" << D;
    }
}

TEST(FastSVDTest, ExactAlgebraIdentitiesOnIntegerMatrices) {
    runExactAlgebra<2>();
    runExactAlgebra<3>();
}

// Matrices whose SVD is known in closed form.
TEST(FastSVDTest, ClosedFormCases) {
    const T phi = (T(1) + std::sqrt(T(5))) / T(2);   // golden ratio

    {   // [[1,1],[0,1]] has singular values phi and 1/phi exactly
        Eigen::Matrix<T,2,2> F; F << 1, 1, 0, 1;
        FastSVDImpl<2> s(F);
        EXPECT_NEAR(s.singularValues()(0), phi,        1e-15);
        EXPECT_NEAR(s.singularValues()(1), T(1) / phi, 1e-15);
    }
    {   // 5 x (3-4-5 rotation): both singular values exactly 5
        Eigen::Matrix<T,2,2> F; F << 3, -4, 4, 3;
        FastSVDImpl<2> s(F);
        EXPECT_NEAR(s.singularValues()(0), T(5), 1e-14);
        EXPECT_NEAR(s.singularValues()(1), T(5), 1e-14);
    }
    {   // antidiagonal: singular values 3 and 2
        Eigen::Matrix<T,2,2> F; F << 0, 2, 3, 0;
        FastSVDImpl<2> s(F);
        EXPECT_NEAR(s.singularValues()(0), T(3), 1e-14);
        EXPECT_NEAR(s.singularValues()(1), T(2), 1e-14);
    }
    {   // 5 x (rotation about z): all three singular values exactly 5
        Eigen::Matrix<T,3,3> F; F << 3, -4, 0, 4, 3, 0, 0, 0, 5;
        FastSVDImpl<3> s(F);
        for (int i = 0; i < 3; i++) EXPECT_NEAR(s.singularValues()(i), T(5), 1e-14);
    }
    {   // signed diagonal: singular values are the sorted magnitudes
        Eigen::Matrix<T,3,3> F; F << 2, 0, 0, 0, -6, 0, 0, 0, 3;
        FastSVDImpl<3> s(F);
        EXPECT_NEAR(s.singularValues()(0), T(6), 1e-14);
        EXPECT_NEAR(s.singularValues()(1), T(3), 1e-14);
        EXPECT_NEAR(s.singularValues()(2), T(2), 1e-14);
    }
    {   // golden-ratio block plus a unit direction: phi, 1, 1/phi
        Eigen::Matrix<T,3,3> F; F << 1, 1, 0, 0, 1, 0, 0, 0, 1;
        FastSVDImpl<3> s(F);
        EXPECT_NEAR(s.singularValues()(0), phi,        1e-15);
        EXPECT_NEAR(s.singularValues()(1), T(1),       1e-15);
        EXPECT_NEAR(s.singularValues()(2), T(1) / phi, 1e-15);
    }
    {   // 7 x permutation: all singular values exactly 7
        Eigen::Matrix<T,3,3> F; F << 0, 0, 7, 0, 7, 0, 7, 0, 0;
        FastSVDImpl<3> s(F);
        for (int i = 0; i < 3; i++) EXPECT_NEAR(s.singularValues()(i), T(7), 1e-14);
    }
}
