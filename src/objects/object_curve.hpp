// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef OBJECTCURVE_HPP
#define OBJECTCURVE_HPP

#include "object_general.hpp"

class ObjectCurve : public ObjectGeneral {
public:

    ~ObjectCurve(){}

    ObjectCurve(BoundaryCondition bc_in = NOSLIP, T friction_in = 0.0, std::string name_in = "") : ObjectGeneral(bc_in, friction_in, name_in) {}

    bool inside(const TV& X_in) override {

        T x = X_in(0);
        T y = X_in(1);

        if (y > x*x)
            return false;
        else
            return true;

    }

    TV normal(const TV& X_in) override {

        T x = X_in(0);
        TV n;

        n(0) = -2*x;
        n(1) = 1;

        return n.normalized();

    }

};

#endif  // OBJECTCURVE_HPP
