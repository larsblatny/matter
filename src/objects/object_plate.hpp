// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef OBJECTPLATE_HPP
#define OBJECTPLATE_HPP

#include "../tools.hpp"

class ObjectPlate{
public:

#ifdef THREEDIM

    ObjectPlate(T pos_object, PlateType plate_type, BC bc = BC::NoSlip, T friction = 0, T pos_lower = -1e15, T pos_upper = 1e15, T vx_object = 0.0, T vy_object = 0.0, T vz_object = 0.0, T vmin_factor = 1.0, T load_factor = 0.0, std::string name = "plate", bool force_calc = false) :
              pos_object(pos_object),
              plate_type(plate_type),
              bc(bc),
              friction(friction),
              pos_lower(pos_lower),
              pos_upper(pos_upper),
              vx_object(vx_object),
              vy_object(vy_object),
              vz_object(vz_object),
              vx_object_original(vx_object),
              vy_object_original(vy_object),
              vz_object_original(vz_object),
              vmin_factor(vmin_factor),
              load_factor(load_factor),
              name(name),
              force_calc(force_calc) {}

      bool inside(const TV& X_in) const {
          if (plate_type == PlateType::left){
              if (X_in(1) < pos_upper && X_in(1) > pos_lower && (X_in(0) - pos_object) <= 0){
                  return true;
              } else{
                  return false;
              }
          }
          else if (plate_type == PlateType::right){
              if (X_in(1) < pos_upper && X_in(1) > pos_lower && (pos_object - X_in(0)) <= 0){
                  return true;
              } else{
                  return false;
              }
          }
          else if (plate_type == PlateType::bottom){
              if (X_in(0) < pos_upper && X_in(0) > pos_lower && (X_in(1) - pos_object) <= 0){
                  return true;
              } else{
                  return false;
              }
          }
          else if (plate_type == PlateType::top){
              if (X_in(0) < pos_upper && X_in(0) > pos_lower && (pos_object - X_in(1)) <= 0){
                  return true;
              } else{
                  return false;
              }
          }
          else if (plate_type == PlateType::back){
              if (X_in(2) - pos_object <= 0){
                  return true;
              } else{
                  return false;
              }
          }
          else if (plate_type == PlateType::front){
              if (pos_object - X_in(2) <= 0){
                  return true;
              } else{
                  return false;
              }
          }
          else{
            debug("FATAL: invalid plate type");
            return false;
          }

        } // end inside(x,y,z)


#else // TWODIM

    ObjectPlate(T pos_object, PlateType plate_type, BC bc = BC::NoSlip, T friction = 0, T pos_lower = -1e15, T pos_upper = 1e15, T vx_object = 0.0, T vy_object = 0.0, T vmin_factor = 1.0, T load_factor = 0.0, std::string name = "plate", bool force_calc = false) :
              pos_object(pos_object),
              plate_type(plate_type),
              bc(bc),
              friction(friction),
              pos_lower(pos_lower),
              pos_upper(pos_upper),
              vx_object(vx_object),
              vy_object(vy_object),
              vx_object_original(vx_object),
              vy_object_original(vy_object),
              vmin_factor(vmin_factor),
              load_factor(load_factor),
              name(name),
              force_calc(force_calc) {}

    bool inside(const TV& X_in) const {
        if (plate_type == PlateType::left){
            if (X_in(1) < pos_upper && X_in(1) > pos_lower && (X_in(0) - pos_object) <= 0){
                return true;
            } else{
                return false;
            }
        }
        else if (plate_type == PlateType::right){
            if (X_in(1) < pos_upper && X_in(1) > pos_lower && (pos_object - X_in(0)) <= 0){
                return true;
            } else{
                return false;
            }
        }
        else if (plate_type == PlateType::bottom){
            if (X_in(0) < pos_upper && X_in(0) > pos_lower && (X_in(1) - pos_object) <= 0){
                return true;
            } else{
                return false;
            }
        }
        else if (plate_type == PlateType::top){
            if (X_in(0) < pos_upper && X_in(0) > pos_lower && (pos_object - X_in(1)) <= 0){
                return true;
            } else{
                return false;
            }
        }
        else{
            debug("FATAL: invalid plate type");
            return false;
        }
    }

#endif


  void move(T dt, T frame_dt, T time){

      T factor = (time < load_factor * frame_dt) ? (1.0 / vmin_factor) : 1.0;
      vx_object = vx_object_original * factor;
      vy_object = vy_object_original * factor;
    #ifdef THREEDIM
      vz_object = vz_object_original * factor;
    #endif

      if (plate_type == PlateType::left || plate_type == PlateType::right){
          pos_object += dt * vx_object;
          pos_upper  += dt * vy_object;
          pos_lower  += dt * vy_object;
      }
      else if (plate_type == PlateType::bottom || plate_type == PlateType::top){
          pos_object += dt * vy_object;
          pos_upper  += dt * vx_object;
          pos_lower  += dt * vx_object;
      }
#ifdef THREEDIM
      else if (plate_type == PlateType::back || plate_type == PlateType::front){
          pos_object += dt * vz_object;
      }
#endif


  } // end move(dt, frame_dt, time)

  T pos_object;

  PlateType plate_type;
  BC bc;
  T friction;

  T pos_upper;
  T pos_lower;

  T vx_object;
  T vy_object;
  #ifdef THREEDIM
    T vz_object;
  #endif

  T vx_object_original;
  T vy_object_original;
  #ifdef THREEDIM
    T vz_object_original;
  #endif

  T vmin_factor;
  T load_factor;

  std::string name;

  bool force_calc;

  TV force = TV::Zero();

};


#endif  // OBJECTPLATE_HPP
