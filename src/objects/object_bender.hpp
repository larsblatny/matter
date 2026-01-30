// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef OBJECTBENDER_HPP
#define OBJECTBENDER_HPP

#include "object_general.hpp"

class ObjectBender : public ObjectGeneral {
public:
    T r;
    TV c;

    T ri;
    T rf;
    T tf;

    ~ObjectBender(){}

    ObjectBender(BC bc_in = BC::NoSlip, T friction_in = 0.0, std::string name_in = "", T ri_in = 1, T rf_in = 1, T tf_in = 1) : ObjectGeneral(bc_in, friction_in, name_in), ri(ri_in), rf(rf_in), tf(tf_in), r(ri_in), c(TV(0,-ri_in)) {}

    bool inside(const TV& X_in) const override {

        TV d = X_in - c;

        return d.squaredNorm() <= (r*r);
    }

    TV normal(const TV& X_in) const override {

        TV n = X_in - c; // center to point gives OUTWARD normal

        double len = n.norm();
        if (len > 0.0)
            n /= len;

        return n;
    }

    void move(T time) override {

        // update radius
        T stime = time / tf;
        r = 1.0 / ( (1-stime)/ri + stime/rf );

        // update center
        c(1) = -r;
    }

    TV v_object(T time, const TV& X_in) const override {

        T stime = time / tf;
        T rad = 1.0 / ( (1-stime)/ri + stime/rf );

        // dr/dt (always negative)
        T drds = -rad * rad * (1.0 / rf - 1.0 / ri);
        T drdt = drds / tf;

        TV ey = TV::Zero();
        ey(1) = 1;

        TV vel = drdt * (ey - normal(X_in));
        vel(0) = 0;

        // debug("Position: \n", X_in, "\n");
        // debug("Velocity: \n", vel, "\n");

        return vel;
    }



};

#endif  // OBJECTBENDER_HPP
