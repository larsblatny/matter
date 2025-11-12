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
#include <memory>

#include "../tools.hpp"
#include "../data_structures.hpp"
#include "../timer.hpp"

#include "../objects/object_general.hpp"
#include "../objects/object_plate.hpp"

class Simulation{
public:
  Simulation();
  ~Simulation(){};

  int exit = 0;

#ifdef THREEDIM
  const unsigned int dim = 3;
#else
  const unsigned int dim = 2;
#endif

  unsigned int n_threads = 1; // number of OMP threads
  unsigned int end_frame = 1; // last frame in the simulation

  bool is_initialized = false;
  bool save_sim = true;
  bool reduce_verbose = false;
  bool pbc = false;
  bool change_particle_positions = false;
  bool gravity_special = false;
  bool save_grid = false;
  bool use_mibf = false;
  bool use_musl = false;

  TV grid_reference_point = 2e10 * TV::Ones();
  TV gravity = TV::Zero();

  T min_dt = 1e-14; // minimum dt, also used to check  end of frame, use with caution
  T fps = 1; // frames per second
  T cfl = 0.5; // classical CFL coefficient
  T cfl_elastic = 0.5; // CFL-like coffefficient for elastic wave speed
  T flip_ratio = -0.95; // [0,1]: PIC/FLIP where 1 is FLIP and 0 is PIC. [-1,0): APIC/AFLIP where -1 is AFLIP.
  T rho = 1000; // density
  T gravity_time = 0; // used if gravity_special = true
  // bool no_liftoff = true; // can be used if gravity_special = true
  T Lx = 1;
  T Ly = 1;
#ifdef THREEDIM
  T Lz = 1;
#endif

  // Particle data
  Particles particles;
  unsigned int Np; // number of particles
  unsigned int num_add_pbc_particles; // used for periodic boundary conditions
  T particle_mass;   // constant particle mass
  T particle_volume; // initial particle volume
  T dx; // grid cell size

  // Elastoplasticity
  ElasticModel elastic_model = ElasticModel::Hencky;
  PlasticModel plastic_model = PlasticModel::NoPlasticity;
  HardeningLaw hardening_law = HardeningLaw::ExpoImpl;
  bool use_pradhana = true; // volume correction for Drucker-Prager-based models
  T E = 1e6; // Young's modulus (3D)
  T nu = 0.3; // Poisson's ratio (3D)
  T stress_tolerance = 1e-5; // only used in some models
  T xi = 0; // softening/hardening parameter

  // Von Mises:
  T q_max = 100;
  T q_min = 100;
  T p_min = -1e20;

  // Drucker Prager
  T M = 1;
  T q_cohesion = 0;

  // Perzyna
  T perzyna_exp = 1;
  T perzyna_visc = 0;
  
  // Associative viscoplastic rule
  T eta = 1.0; // 1.0 equivalent to no viscosity

  // MCC
  T beta = 0;
  T p0 = 1000;

  // mu(I) rheology
  T rho_s = 2500;
  T grain_diameter = 1e-3;
  T I_ref = 0.279;
  T mu_1 = std::tan(20.9*M_PI/180.0);
  T mu_2 = std::tan(32.8*M_PI/180.0);;

  // Prefactor for q in plasticity models
  T q_prefac  = 1.0 / std::sqrt(2.0); // q = factor * ||dev(tau)||

  // Objects
  std::vector<std::unique_ptr<ObjectPlate>> plates;
  std::vector<std::unique_ptr<ObjectGeneral>> objects;

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
  void plasticity(unsigned int p, unsigned int & plastic_count, TM & Fe_trial);
  void moveObjects();
  void boundaryCollision(int index, TV Xi, TV& vi);
  void overwriteGridVelocity(TV Xi, TV& vi);
  void checkMomentumConservation();
  void checkMassConservation();
  void addExternalParticleGravity();
  std::pair<TMX, TMX> createExternalGridGravity();
  T calculateBulkModulus();
  TM NeoHookeanPiola(TM & Fe);
  TM HenckyPiola(TM & Fe);

private:

  unsigned int current_time_step = 0;
  unsigned int frame = 0;

  T time = 0;
  T runtime_p2g = 0;
  T runtime_g2p = 0;
  T runtime_euler = 0;
  T runtime_defgrad = 0;
  T runtime_total = 0;

  T final_time;
  T frame_dt;
  T dt;
  T dt_max;

  std::string sim_name;
  std::string directory;
 
  T wave_speed; // elastic wave speed
  T mu; // shear modulus
  T lambda; // first Lame parameter
  T K; // bulk modulus
  T fac_Q; // for mu(I) rheology

  TV gravity_final; // used if gravity_special = true

  // Prefactors for plasticity models
  T d_prefac;    // gamma    = factor * ||dev(eps)||
  T e_mu_prefac; // q        = factor * ||dev(eps)||
  T f_mu_prefac; // q^tr - q = factor * dt * gamma_dot

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
    return (i*Ny + j) * Nz + k; 
  }
#else
  inline unsigned int ind(unsigned int i, unsigned int j){
      return (i*Ny + j); 
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
}; // end Simulation class

inline TM Simulation::NeoHookeanPiola(TM & Fe){
    return mu * (Fe - Fe.transpose().inverse()) + lambda * std::log(Fe.determinant()) * Fe.transpose().inverse();
} // end NeoHookeanPiola

inline TM Simulation::HenckyPiola(TM & Fe){
    Eigen::JacobiSVD<TM> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
    TA sigma = svd.singularValues().array(); 
    TM logSigma = sigma.abs().log().matrix().asDiagonal();
    TM invSigma = sigma.inverse().matrix().asDiagonal();
    TM dPsidF = svd.matrixU() * ( 2*mu*invSigma*logSigma + lambda*logSigma.trace()*invSigma ) * svd.matrixV().transpose();
    return dPsidF;
} // end HenckyPiola

#endif  // SIMULATION_HPP
