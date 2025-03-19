// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <cmath>
#include <assert.h>
#include <vector>
#include <chrono>

#include "../tools.hpp"
#include "../data_structures.hpp"
#include "../timer.hpp"

#include "../objects/object_plate.hpp"
#include "../objects/object_general.hpp"
#include "../objects/object_bump.hpp"
#include "../objects/object_gate.hpp"
#include "../objects/object_ramp.hpp"

class Simulation{
public:
  Simulation();
  ~Simulation(){};

  int exit = 0;

  unsigned int n_threads = 1;
  unsigned int end_frame = 1;
  T fps = 1;

  bool is_initialized = false;
  bool save_sim = true;
  bool reduce_verbose = false;
  bool pbc = false;
  bool change_particle_positions = false;
  bool gravity_special = false;
  bool save_grid = false;
  bool use_mibf = false;
  bool use_musl = false;

  T cfl = 0.5;
  T cfl_elastic = 0.5;
  T flip_ratio = -0.95;

  T rho = 1000;

  TV gravity = TV::Zero();
  T gravity_time = 0;
  // bool no_liftoff = true;

  T Lx = 1;
  T Ly = 1;
#ifdef THREEDIM
  T Lz = 1;
#endif

  TV grid_reference_point = 2e10 * TV::Ones();

  // Particle data
  Particles particles;
  unsigned int Np;
#ifndef MULTIMATERIAL
  T particle_mass;
#endif
  T particle_volume; // initial particle volume
  T dx;

  // Elastoplasticity
  ElasticModel elastic_model = Hencky;
  PlasticModel plastic_model = NoPlasticity;
  HardeningLaw hardening_law = ExpoImpl;

  T E = 1e6; // Young's modulus (3D)
  T nu = 0.3; // Poisson's ratio (3D)

  bool use_pradhana = true;
  bool use_mises_q = false;

  // Von Mises:
  T q_max = 100;
  T q_min = 100;
  T p_min = -1e20;
  T xi = 0;

  // Drucker Prager
  T M = 1;
  T q_cohesion = 0;

  // Perzyna
  T perzyna_exp = 1;
  T perzyna_visc = 0;

  // MCC
  T beta = 0;
  T p0 = 1000;

  // mu(I) rheology
  T rho_s = 2500;
  T grain_diameter = 1e-3;
  T I_ref = 0.279;
  T mu_1 = std::tan(20.9*M_PI/180.0);
  T mu_2 = std::tan(32.8*M_PI/180.0);;

  T stress_tolerance = 1e-5;
  
  // Making wave_speed public to declare it in mpm.cpp when MULTIMATERIAL
  #ifdef MULTIMATERIAL
      T wave_speed;
      T dt_max;
  #endif

  // Objects
  std::vector<ObjectPlate> plates;
  std::vector<ObjectGeneral*> objects;

  // Functions
  void initialize(bool save = true, std::string dir = "output/", std::string name = "dummy");
  void simulate();
  void saveInfo();
  void saveTiming();
  void saveAvgData();
  void computeAvgData(TM& volavg_cauchy, TM& volavg_kirchh, T& Javg);
  void saveParticleData(std::string extra = "");
  void saveGridData(std::string extra = "");

  void createDirectory();

  void advanceStep();
  void updateDt();

  void resizeGrid();
  void remeshFixed(unsigned int extra_nodes);
  void remeshFixedInit(unsigned int sfx, unsigned int sfy, unsigned int sfz);
  void remeshFixedCont();

  void P2G();
  void explicitEulerUpdate();
  void G2P();
  void deformationUpdate();

  void MUSL();

  void positionUpdate();
  void PBCAddParticles1D();
  void PBCAddParticles(unsigned int safety_factor);
  void PBCDelParticles();
  unsigned int num_add_pbc_particles;

  TM NeoHookeanPiola(TM & Fe);
  TM HenckyPiola(TM & Fe);
  void plasticity(unsigned int p, unsigned int & plastic_count, TM & Fe_trial);

  void moveObjects();

  void boundaryCollision(int index, TV Xi, TV& vi);
  void overwriteGridVelocity(TV Xi, TV& vi);

  T calculateBulkModulus();
  void checkMomentumConservation();
  void checkMassConservation();

  void addExternalParticleGravity();
  std::pair<TMX, TMX> createExternalGridGravity();


  #ifdef THREEDIM
    const unsigned int dim = 3;
  #else
    const unsigned int dim = 2;
  #endif

private:

  unsigned int current_time_step = 0;
  unsigned int frame = 0;
  T time = 0;

  T runtime_p2g = 0;
  T runtime_g2p = 0;
  T runtime_euler = 0;
  T runtime_defgrad = 0;
  T runtime_total = 0;

  std::string sim_name;
  std::string directory;

  T final_time;
  T frame_dt;
  T dt;
  #ifndef MULTIMATERIAL
      T wave_speed;
      T dt_max;
  #endif

  TV gravity_final;

  T mu; // shear modulus
  T lambda; // first Lame parameter
  T K; // bulk modulus

  T fac_Q; // for mu(I) rheology

  // Prefactors for plasticity models
  T q_prefac;    // q        = factor * ||dev(tau)||
  T d_prefac;    // gamma    = factor * ||dev(eps)||
  T e_mu_prefac; // q        = factor * ||dev(eps)||
  T f_mu_prefac; // q^tr - q = factor * dt * gamma_dot
  T rma_prefac;

  // Precomputations
  T one_over_dx;
  T one_over_dx_square;
  T apicDinverse;

  // Grid handling and remeshing
  Grid grid;
  unsigned int Nx, Ny;
#ifdef THREEDIM
  unsigned int Nz;
#endif
  unsigned int grid_nodes;

#ifdef THREEDIM
    inline unsigned int ind(unsigned int i, unsigned int j, unsigned int k){
      return (i*Ny + j) * Nz + k; // 3D
    }
#else
    inline unsigned int ind(unsigned int i, unsigned int j){
        return (i*Ny + j); // 3D
    }
#endif

  T min_x_init;
  T max_x_init;
  T Nx_init;
  T low_x_init;
  T high_x_init;

  T min_y_init;
  T max_y_init;
  T Ny_init;
  T low_y_init;
  T high_y_init;

#ifdef THREEDIM
  T min_z_init;
  T max_z_init;
  T Nz_init;
  T low_z_init;
  T high_z_init;
#endif


}; // END Simulation class



inline TM Simulation::NeoHookeanPiola(TM & Fe){
    return mu * (Fe - Fe.transpose().inverse()) + lambda * std::log(Fe.determinant()) * Fe.transpose().inverse();
} // end NeoHookeanPiola

inline TM Simulation::HenckyPiola(TM & Fe){
    Eigen::JacobiSVD<TM> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
    TA sigma = svd.singularValues().array(); // abs() for inverse also??
    TM logSigma = sigma.abs().log().matrix().asDiagonal();
    TM invSigma = sigma.inverse().matrix().asDiagonal();
    TM dPsidF = svd.matrixU() * ( 2*mu*invSigma*logSigma + lambda*logSigma.trace()*invSigma ) * svd.matrixV().transpose();
    return dPsidF;
} // end HenckyPiola



#endif  // SIMULATION_HPP
