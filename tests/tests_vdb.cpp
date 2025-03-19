// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include <gtest/gtest.h>

#include "../src/tools.hpp"
#include "../src/simulation/simulation.hpp"

#include "../src/objects/object_curve.hpp"
#include "../src/objects/object_vdb.hpp"

TEST(FileIOTest, LoadFile) {
    std::string file_path = std::string(INCLUDE_DIR) + "/curve.vdb";
    std::ifstream file(file_path);
    ASSERT_TRUE(file.is_open());
}


TEST(BoundaryTest, CalcNormal) {

    openvdb::initialize(); 

    ObjectVdb obj1 = ObjectVdb(std::string(INCLUDE_DIR) + "/curve.vdb"); 
    ObjectCurve obj2 = ObjectCurve();

    T xpos = -0.5;
    #ifdef THREEDIM
        TV pos(xpos, xpos*xpos, 0.0);
    #else
        TV pos(xpos, xpos*xpos);
    #endif
    
    T diff1 = std::abs(obj1.normal(pos)[0]-obj2.normal(pos)[0]);
    T diff2 = std::abs(obj1.normal(pos)[1]-obj2.normal(pos)[1]);
    T diff = diff1 + diff2;

    EXPECT_TRUE(obj1.inside(pos) == obj2.inside(pos));

    ASSERT_NEAR(diff, 0.0, 1e-4);

}


TEST(BoundaryTest, VDB) {

    openvdb::initialize();

    Simulation sim;
    sim.initialize(false);

    T half_rounds = 1;
    int num_frames_in_half_round = 100;
    sim.end_frame = half_rounds*num_frames_in_half_round;
    sim.fps = num_frames_in_half_round / 0.5949238427;

    sim.n_threads = 8;
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

    std::string file_path = std::string(INCLUDE_DIR) + "/curve.vdb";
    ObjectVdb curve  = ObjectVdb(file_path, SLIPFREE, 0.0); sim.objects.push_back(&curve);


    sim.simulate();

    T v_sim = sim.particles.v[0](0);
    T v_true = std::sqrt(2*9.81);
    T diff = std::abs(v_sim-v_true)/v_true;
    ASSERT_NEAR(diff, 0.0, 1e-3);
}
