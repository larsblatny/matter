// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef DATA_STRUCTURES_HPP
#define DATA_STRUCTURES_HPP

#include "tools.hpp"

class Particles{
public:
  Particles(unsigned int Np = 1){
      x.resize(Np); std::fill( x.begin(), x.end(), TV::Zero() );
      v.resize(Np); std::fill( v.begin(), v.end(), TV::Zero() );
      pic.resize(Np); std::fill( pic.begin(), pic.end(), TV::Zero() );
      flip.resize(Np); std::fill( flip.begin(), flip.end(), TV::Zero() );

      eps_pl_dev.resize(Np); std::fill( eps_pl_dev.begin(), eps_pl_dev.end(), 0.0 );
      eps_pl_vol.resize(Np); std::fill( eps_pl_vol.begin(), eps_pl_vol.end(), 0.0 );
      eps_pl_vol_pradhana.resize(Np);std::fill( eps_pl_vol_pradhana.begin(),eps_pl_vol_pradhana.end(), 0.0 );
      delta_gamma.resize(Np); std::fill( delta_gamma.begin(), delta_gamma.end(), 0.0 );
      viscosity.resize(Np); std::fill( viscosity.begin(), viscosity.end(), 0.0 );
      muI.resize(Np); std::fill( muI.begin(), muI.end(), 0.0 );

      F.resize(Np); std::fill( F.begin(), F.end(), TM::Identity() );
      Bmat.resize(Np); std::fill( Bmat.begin(), Bmat.end(), TM::Zero() );
      tau.resize(Np); std::fill( tau.begin(), tau.end(), TM::Zero() );
  }

  std::vector<TV> x;
  std::vector<TV> v;
  std::vector<TV> pic;
  std::vector<TV> flip;

  std::vector<T> eps_pl_dev;
  std::vector<T> eps_pl_vol;
  std::vector<T> eps_pl_vol_pradhana;
  std::vector<T> delta_gamma;
  std::vector<T> viscosity;
  std::vector<T> muI;

  std::vector<TM> F;
  std::vector<TM> Bmat;
  std::vector<TM> tau;


};

class Grid{
public:
    Grid(){}
    std::vector<T> x;
    std::vector<T> y;
    #ifdef THREEDIM
    std::vector<T> z;
    #endif
    std::vector<TV> v;
    std::vector<TV> flip;
    std::vector<TV> force;
    std::vector<T> mass;
    std::vector<T> friction;
    T xc;
    T yc;
    #ifdef THREEDIM
    T zc;
    #endif
};

#endif // DATA_STRUCTURES_HPP
