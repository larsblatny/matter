// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef OBJECTGROUNDROTATED_HPP
#define OBJECTGROUNDROTATED_HPP

#include "object_general.hpp"

class ObjectGroundRotated : public ObjectGeneral {
public:
    T theta;

    ~ObjectGroundRotated(){}

    ObjectGroundRotated(BC bc_in = BC::NoSlip, T friction_in = 0.0, std::string name_in = "", bool force_calc_in = false, T theta_in = 0) : ObjectGeneral(bc_in, friction_in, name_in, force_calc_in), theta(theta_in) {}

    bool inside(const TV& X_in) const override {
        
        T x = X_in(0);
        T y = X_in(1);

        if (y * std::cos(theta) + x * std::sin(theta) > 0) 
            return false;
        else
            return true;

    }

    TV normal(const TV& X_in) const override {

        TV n = TV::Zero();

        n(0) = std::sin(theta);
        n(1) = std::cos(theta);

        return n;
    }

};

#endif  // OBJECTGROUNDROTATED_HPP
