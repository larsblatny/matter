// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include <gtest/gtest.h>

#include "../src/tools.hpp"
#include "../src/simulation/simulation.hpp"
#include "../src/sampling/sampling_particles.hpp"

#include "../src/objects/object_curve.hpp"
#include "../src/objects/object_ground.hpp"


TEST(BoundaryTest, AnalyticSLIPFREE) {

    Simulation sim;
    sim.initialize(false);

    T half_rounds = 1;
    int num_frames_in_half_round = 100;
    sim.end_frame = half_rounds*num_frames_in_half_round;
    sim.fps = num_frames_in_half_round / 0.5949238427;

    sim.n_threads = 1;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    sim.elastic_model = Hencky;
    sim.plastic_model = NoPlasticity;
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1000;

    sim.gravity[1] = -9.81;

    sim.dx = 0.001;
    sim.particle_volume = std::pow(sim.dx, sim.dim);
    #ifdef MULTIMATERIAL
        for(int p = 0; p < sim.Np; p++)
            sim.particles.mass[p] = sim.rho;
    #else
        sim.particle_mass = sim.rho * sim.particle_volume;
    #endif
    sim.Np = 1;
    sim.particles = Particles(sim.Np);

    sim.particles.x[0](0) = -1;
    sim.particles.x[0](1) = 1;

    ObjectCurve curve = ObjectCurve(SLIPFREE, 0.0);  sim.objects.push_back(&curve);

    sim.simulate();

    T v_sim = sim.particles.v[0](0);
    T v_true = std::sqrt(2*9.81);
    T diff = std::abs(v_sim-v_true)/v_true;
    ASSERT_NEAR(diff, 0.0, 1e-4);
}


TEST(BoundaryTest, AnalyticSLIPSTICK) {

    Simulation sim;
    sim.initialize(false);

    T half_rounds = 1;
    int num_frames_in_half_round = 100;
    sim.end_frame = half_rounds*num_frames_in_half_round;
    sim.fps = num_frames_in_half_round / 0.5949238427;

    sim.n_threads = 1;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    sim.elastic_model = Hencky;
    sim.plastic_model = NoPlasticity;
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1000;

    sim.gravity[1] = -9.81;

    sim.dx = 0.001;
    sim.particle_volume = std::pow(sim.dx, sim.dim);
    #ifdef MULTIMATERIAL
        for(int p = 0; p < sim.Np; p++)
            sim.particles.mass[p] = sim.rho;
    #else
        sim.particle_mass = sim.rho * sim.particle_volume;
    #endif
    sim.Np = 1;
    sim.particles = Particles(sim.Np);

    sim.particles.x[0](0) = -1;
    sim.particles.x[0](1) = 0.9999;

    ObjectCurve curve = ObjectCurve(SLIPSTICK, 0.0);  sim.objects.push_back(&curve);

    sim.simulate();

    T v_sim = sim.particles.v[0](0);
    T v_true = std::sqrt(2*9.81);
    T diff = std::abs(v_sim-v_true)/v_true;
    ASSERT_NEAR(diff, 0.0, 1e-3);
}

TEST(CoulombFrictionTest, PlateSLIPFREE) {

    Simulation sim;
    sim.initialize(false);

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

    sim.elastic_model = Hencky;
    sim.plastic_model = NoPlasticity;
    sim.E = 1e6; 
    sim.nu = 0.3; 
    sim.rho = 1000; 

    sim.Lx = 0.1;
    sim.Ly = 0.05;
    T k_rad = 0.005;
    sampleParticles(sim, k_rad);
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) += 0.5*sim.dx;
    }
    sim.grid_reference_point = TV::Zero();

    T friction = std::tan(15.0 * M_PI / 180.0);

    ObjectPlate ground = ObjectPlate(0, bottom, SLIPFREE, friction); 
    sim.plates.push_back(ground);

    sim.simulate();

    T mean_sim_x;
    for(int p = 0; p < sim.Np; p++)
        mean_sim_x += sim.particles.x[p](0);
    mean_sim_x /= (T)sim.Np;

    T final_time = (T)sim.end_frame / sim.fps;
    T mean_true_x = 0.5*9.81*final_time*final_time*(std::sin(theta) - friction*std::cos(theta));
    T diff = std::abs(mean_sim_x-mean_true_x) / mean_true_x;
    ASSERT_NEAR(diff, 0.0, 1e-3);
}

TEST(CoulombFrictionTest, PlateSLIPSTICK) {

    Simulation sim;
    sim.initialize(false);

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

    sim.elastic_model = Hencky;
    sim.plastic_model = NoPlasticity;
    sim.E = 1e6; 
    sim.nu = 0.3; 
    sim.rho = 1000; 

    sim.Lx = 0.1;
    sim.Ly = 0.05;
    T k_rad = 0.005;
    sampleParticles(sim, k_rad);
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) += 0.5*sim.dx;
    }
    sim.grid_reference_point = TV::Zero();

    T friction = std::tan(15.0 * M_PI / 180.0);

    ObjectPlate ground = ObjectPlate(0, bottom, SLIPSTICK, friction); 
    sim.plates.push_back(ground);

    sim.simulate();

    T mean_sim_x;
    for(int p = 0; p < sim.Np; p++)
        mean_sim_x += sim.particles.x[p](0);
    mean_sim_x /= (T)sim.Np;

    T final_time = (T)sim.end_frame / sim.fps;
    T mean_true_x = 0.5*9.81*final_time*final_time*(std::sin(theta) - friction*std::cos(theta));
    T diff = std::abs(mean_sim_x-mean_true_x) / mean_true_x;
    ASSERT_NEAR(diff, 0.0, 1e-3);
}



TEST(CoulombFrictionTest, GeneralSLIPFREE) {

    Simulation sim;
    sim.initialize(false);

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

    sim.elastic_model = Hencky;
    sim.plastic_model = NoPlasticity;
    sim.E = 1e6; 
    sim.nu = 0.3; 
    sim.rho = 1000; 

    sim.Lx = 0.1;
    sim.Ly = 0.05;
    T k_rad = 0.005;
    sampleParticles(sim, k_rad);
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) += 0.5*sim.dx;
    }
    sim.grid_reference_point = TV::Zero();

    T friction = std::tan(15.0 * M_PI / 180.0);

    ObjectGround ground = ObjectGround(SLIPFREE, friction); 
    sim.objects.push_back(&ground);

    sim.simulate();

    T mean_sim_x;
    for(int p = 0; p < sim.Np; p++)
        mean_sim_x += sim.particles.x[p](0);
    mean_sim_x /= (T)sim.Np;

    T final_time = (T)sim.end_frame / sim.fps;
    T mean_true_x = 0.5*9.81*final_time*final_time*(std::sin(theta) - friction*std::cos(theta));
    T diff = std::abs(mean_sim_x-mean_true_x) / mean_true_x;
    ASSERT_NEAR(diff, 0.0, 1e-3);
}

TEST(CoulombFrictionTest, GeneralSLIPSTICK) {

    Simulation sim;
    sim.initialize(false);

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

    sim.elastic_model = Hencky;
    sim.plastic_model = NoPlasticity;
    sim.E = 1e6; 
    sim.nu = 0.3; 
    sim.rho = 1000; 

    sim.Lx = 0.1;
    sim.Ly = 0.05;
    T k_rad = 0.005;
    sampleParticles(sim, k_rad);
    for(int p = 0; p < sim.Np; p++){
        sim.particles.x[p](0) -= 0.5*sim.Lx;
        sim.particles.x[p](1) += 0.5*sim.dx;
    }
    sim.grid_reference_point = TV::Zero();

    T friction = std::tan(15.0 * M_PI / 180.0);

    ObjectGround ground = ObjectGround(SLIPSTICK, friction); 
    sim.objects.push_back(&ground);

    sim.simulate();

    T mean_sim_x;
    for(int p = 0; p < sim.Np; p++)
        mean_sim_x += sim.particles.x[p](0);
    mean_sim_x /= (T)sim.Np;

    T final_time = (T)sim.end_frame / sim.fps;
    T mean_true_x = 0.5*9.81*final_time*final_time*(std::sin(theta) - friction*std::cos(theta));
    T diff = std::abs(mean_sim_x-mean_true_x) / mean_true_x;
    ASSERT_NEAR(diff, 0.0, 1e-3);
}




TEST(BoundaryTest, MIBF) {

    T friction = std::tan(30.0 * M_PI / 180.0);
    T theta = 32 * M_PI / 180;
    
    ObjectPlate ground = ObjectPlate(0,  bottom, SLIPFREE, friction);  
    
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

    sim_one.elastic_model = Hencky;
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
    sim_one.grid_reference_point = TV::Zero();

    sim_one.plates.push_back(ground);

    sim_one.plastic_model = DPVisc; 

    sim_one.use_pradhana = false; 
    sim_one.use_mises_q = false;
    sim_one.use_mibf = false;          // NB

    sim_one.M = friction;
    sim_one.q_cohesion = 0;
    sim_one.perzyna_exp = 1;
    sim_one.perzyna_visc = 0;

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

    sim_two.elastic_model = Hencky;
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
    sim_two.grid_reference_point = TV::Zero();

    sim_two.plates.push_back(ground);

    sim_two.plastic_model = DPVisc; 

    sim_two.use_pradhana = false; 
    sim_two.use_mises_q = false;
    sim_two.use_mibf = true;          // NB       

    sim_two.M = friction;
    sim_two.q_cohesion = 0;
    sim_two.perzyna_exp = 1;
    sim_two.perzyna_visc = 0;

    sim_two.simulate();

    auto max_x_it_2 = std::max_element( sim_two.particles.x.begin(), sim_two.particles.x.end(), [](const TV &x1, const TV &x2){return x1(0) < x2(0);} );
    T max_x_2 = (*max_x_it_2)(0);

    T diff = std::abs(max_x_1 - max_x_2);
    debug(diff);
    ASSERT_NEAR(diff, 0.0, 1e-13);


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

    sim_three.elastic_model = Hencky;
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
    sim_three.grid_reference_point = TV::Zero();

    sim_three.plates.push_back(ground);

    sim_three.plastic_model = DPVisc; 

    sim_three.use_pradhana = false; 
    sim_three.use_mises_q = true;        // NB
    sim_three.use_mibf = false;          // NB

    sim_three.M = std::sqrt(3.0)*friction; // NB
    sim_three.q_cohesion = 0;
    sim_three.perzyna_exp = 1;
    sim_three.perzyna_visc = 0;

    sim_three.simulate();

    auto max_x_it_3 = std::max_element( sim_three.particles.x.begin(), sim_three.particles.x.end(), [](const TV &x1, const TV &x2){return x1(0) < x2(0);} );
    T max_x_3 = (*max_x_it_3)(0);

    T difff = std::abs(max_x_1 - max_x_3);
    debug(difff);
    ASSERT_NEAR(difff, 0.0, 1e-13);

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
        ObjectPlate ground = ObjectPlate(0-0.5*sim.dx,       bottom, NOSLIP, 0, -1e15, 1e15,   0,  vel, 0, vmin_factor, load_factor);  sim.plates.push_back(ground);
        ObjectPlate compre = ObjectPlate(sim.Ly+0.5*sim.dx,  top,    NOSLIP, 0, -1e15, 1e15,   0, -vel, 0, vmin_factor, load_factor);  sim.plates.push_back(compre);
    #else
        ObjectPlate ground = ObjectPlate(0-0.5*sim.dx,       bottom, NOSLIP, 0, -1e15, 1e15,   0,  vel,    vmin_factor, load_factor);  sim.plates.push_back(ground);
        ObjectPlate compre = ObjectPlate(sim.Ly+0.5*sim.dx,  top,    NOSLIP, 0, -1e15, 1e15,   0, -vel,    vmin_factor, load_factor);  sim.plates.push_back(compre);
    #endif

    sim.elastic_model = Hencky;
    sim.plastic_model = NoPlasticity;

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

    sim.elastic_model = Hencky;
    sim.plastic_model = NoPlasticity;

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
    sim.grid_reference_point = TV::Zero();

    sim.elastic_model = Hencky;
    sim.plastic_model = DPSoft;

    sim.use_mises_q = false;
    sim.use_pradhana = true;

    sim.xi = 0;

    sim.q_cohesion = 0;
    sim.M = std::tan(30.0 * M_PI / 180.0);

    #ifdef THREEDIM
    ObjectPlate sideback   = ObjectPlate(0,       back,   SLIPFREE);  sim.plates.push_back(sideback);
    ObjectPlate sidefront  = ObjectPlate(sim.Lz,  front,  SLIPFREE);  sim.plates.push_back(sidefront);
    ObjectPlate sideleft   = ObjectPlate(0,       left,   SLIPFREE);  sim.plates.push_back(sideleft);
    ObjectPlate ground     = ObjectPlate(0,       bottom, NOSLIP);    sim.plates.push_back(ground);
    #else
    ObjectPlate sideleft   = ObjectPlate(0,  left,   SLIPFREE);  sim.plates.push_back(sideleft);
    ObjectPlate ground     = ObjectPlate(0,  bottom, NOSLIP);    sim.plates.push_back(ground);
    #endif

    sim.simulate();

    auto max_x_it = std::max_element( sim.particles.x.begin(), sim.particles.x.end(), [](const TV &x1, const TV &x2){return x1(0) < x2(0);} );
    T max_x = (*max_x_it)(0);
    T diff = std::abs(max_x - 0.56);
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
    sim.grid_reference_point = TV::Zero();

    sim.elastic_model = Hencky;
    sim.plastic_model = DPVisc;

    sim.use_mises_q = false;
    sim.use_pradhana = true;

    sim.perzyna_visc = 0;
    sim.perzyna_exp = 1;

    sim.q_cohesion = 0;
    sim.M = std::tan(30.0 * M_PI / 180.0);

    #ifdef THREEDIM
    ObjectPlate sideback   = ObjectPlate(0,       back,   SLIPFREE);  sim.plates.push_back(sideback);
    ObjectPlate sidefront  = ObjectPlate(sim.Lz,  front,  SLIPFREE);  sim.plates.push_back(sidefront);
    ObjectPlate sideleft   = ObjectPlate(0,       left,   SLIPFREE);  sim.plates.push_back(sideleft);
    ObjectPlate ground     = ObjectPlate(0,       bottom, NOSLIP);    sim.plates.push_back(ground);
    #else
    ObjectPlate sideleft   = ObjectPlate(0,  left,   SLIPFREE);  sim.plates.push_back(sideleft);
    ObjectPlate ground     = ObjectPlate(0,  bottom, NOSLIP);    sim.plates.push_back(ground);
    #endif

    sim.simulate();

    auto max_x_it = std::max_element( sim.particles.x.begin(), sim.particles.x.end(), [](const TV &x1, const TV &x2){return x1(0) < x2(0);} );
    T max_x = (*max_x_it)(0);
    T diff = std::abs(max_x - 0.56);
    ASSERT_NEAR(diff, 0.0, 0.011);
}
