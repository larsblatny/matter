// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef OBJECTGROUND_HPP
#define OBJECTGROUND_HPP

#include "object_general.hpp"

class ObjectGround : public ObjectGeneral {
public:
    T y_ground;

    ~ObjectGround(){}

    ObjectGround(BoundaryCondition bc_in = NOSLIP, T friction_in = 0.0, std::string name_in = "", T y_ground_in = 0) : ObjectGeneral(bc_in, friction_in, name_in), y_ground(y_ground_in) {}

    bool inside(const TV& X_in) override {

        T y = X_in(1);

        if (y > y_ground)
            return false;
        else
            return true;

    }

    TV normal(const TV& X_in) override {

        TV n;
        
        n(0) = 0;
        n(1) = 1;

        return n;

    }

};

#endif  // OBJECTGROUND_HPP
