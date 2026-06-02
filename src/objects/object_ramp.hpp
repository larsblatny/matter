// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef OBJECTRAMP_HPP
#define OBJECTRAMP_HPP

#include "object_general.hpp"

class ObjectRamp : public ObjectGeneral {
public:

    ~ObjectRamp(){}

    ObjectRamp(BC bc_in, T friction_in, std::string name_in = "ramp", bool force_calc_in = false) : ObjectGeneral(bc_in, friction_in, name_in, force_calc_in) {}

    bool inside(const TV& X_in) const override {

        T x = X_in(0);
        T y = X_in(1);

        T y_limit = 0.1 * std::tanh(-100*x);
        if (y < y_limit)
            return true;
        else
            return false;

    }

    TV normal(const TV& X_in) const override {

        T x = X_in(0);
        TV n = TV::Zero();

        T tmp = std::cosh(-100*x);
        T b_der = 0.1*(-100)/(tmp*tmp);
        n(0) = -b_der;
        n(1) = 1;

        return n.normalized();
    }

};

#endif  // OBJECTRAMP_HPP
