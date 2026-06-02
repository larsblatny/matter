// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef OBJECTSILO_HPP
#define OBJECTSILO_HPP

#include "object_general.hpp"

class ObjectSilo : public ObjectGeneral {
public:
    T cut; // must be a negative number

    ~ObjectSilo(){}

    ObjectSilo(BC bc_in, T friction_in, std::string name_in = "silo", bool force_calc_in = false, T cut_in = -1) : ObjectGeneral(bc_in, friction_in, name_in, force_calc_in), cut(cut_in)  {}

    bool inside(const TV& X_in) const override {

        T x = X_in(0);
        T y = X_in(1);
        T z = X_in(2);

        if (y<cut){
            return false;
        }

        T r_surface = std::tanh(y) + 1;
        T r_surface_sq = r_surface * r_surface;
        T r_point_sq = x*x + z*z;

        if (r_point_sq < r_surface_sq)
            return false; 
        else 
            return true; 

    } // end inside

    TV normal(const TV& X_in) const override {

        TV n;

        T x = X_in(0);
        T y = X_in(1);
        T z = X_in(2);

        T theta = std::atan2(z, x);
        T sech_y = 1.0 / std::cosh(y); 

        n(0) = -std::cos(theta);
        n(1) =  sech_y * sech_y; // always >= 0
        n(2) = -std::sin(theta);

        return n.normalized();

    } // end normal

};

#endif  // OBJECTSILO_HPP
