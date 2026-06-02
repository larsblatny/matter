// Copyright (C) 2026 Lars Blatny. Released under GPL-3.0 license.

#ifndef OBJECTBOXROTATED_HPP
#define OBJECTBOXROTATED_HPP

#include "object_general.hpp"

class ObjectBoxRotated : public ObjectGeneral {
public:
    TV L;
    TV c;
    TV c0;
    TV v;
    T ct; // cos theta
    T st; // sin theta

    ~ObjectBoxRotated(){}

    ObjectBoxRotated(BC bc_in = BC::NoSlip, T friction_in = 0.0, std::string name_in = "box_rotated", bool force_calc_in = false, TV L_in = TV::Ones(), TV c_in = TV::Zero(), TV v_in = TV::Zero(), T theta_in = 0) : ObjectGeneral(bc_in, friction_in, name_in, force_calc_in), L(L_in), c(c_in), c0(c_in), v(v_in), ct(std::cos(theta_in)), st(std::sin(theta_in)) {}

    bool inside(const TV& X_in) const override {

        TV d = X_in - c;

        TV local_d;
        local_d[0] =  ct * d[0] + st * d[1];
        local_d[1] = -st * d[0] + ct * d[1];

        TV h = (T)0.5 * L;

        for (int i = 0; i < TV::RowsAtCompileTime; ++i) {
            if (std::abs(local_d[i]) > h[i]) {
                return false;
            }
        }

        return true;
    }

    TV normal(const TV& X_in) const override {

        // World -> local
        TV d = X_in - c;

        TV local_d;
        local_d[0] =  ct * d[0] + st * d[1];
        local_d[1] = -st * d[0] + ct * d[1];

        TV h = (T)0.5 * L;

        // Local-space normal
        TV n_local = TV::Zero();

        int axis = 0;
        T max_dist = std::abs(local_d[0] / h[0]);

        for (int i = 1; i < TV::RowsAtCompileTime; ++i) {

            T dist = std::abs(local_d[i] / h[i]);

            if (dist > max_dist) {
                max_dist = dist;
                axis = i;
            }
        }

        n_local[axis] = (local_d[axis] >= 0) ? (T)1 : (T)-1;

        // Local -> world
        TV n_world;

        n_world[0] = ct * n_local[0] - st * n_local[1];
        n_world[1] = st * n_local[0] + ct * n_local[1];

        return n_world;
    }

    void move(T time) override {
        c = c0 + v * time;
    }

    TV v_object(T time, const TV& X_in) const override {
        return v; // return zero if no velocity of object. Same as not defining the function.
    }

};

#endif // OBJECTBOXROTATED_HPP
