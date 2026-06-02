// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef OBJECTGATE_HPP
#define OBJECTGATE_HPP

#include "object_general.hpp"

class ObjectGate : public ObjectGeneral {
public:
    T height;

    ~ObjectGate(){}

    ObjectGate(BC bc_in, T friction_in, std::string name_in = "gate", bool force_calc_in = false, T height_in = 0.016) : ObjectGeneral(bc_in, friction_in, name_in, force_calc_in), height(height_in) {}

    bool inside(const TV& X_in) const override {

        T x = X_in(0);
        T y = X_in(1);

        T y_limit = height + 100 * x*x;

        if (y > y_limit)
            return true;
        else
            return false;

    }

    TV normal(const TV& X_in) const override {

        T x = X_in(0);
        TV n = TV::Zero();

        T b_der = 2.0 * 100 * x;
        n(0) = -b_der;
        n(1) = 1;
        n *= -1;

        return n.normalized();
    }

};

#endif  // OBJECTGATE_HPP
