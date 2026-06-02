// Copyright (C) 2026 Lars Blatny. Released under GPL-3.0 license.

#ifndef OBJECTBOX_HPP
#define OBJECTBOX_HPP

#include "object_general.hpp"

class ObjectBox : public ObjectGeneral {
public:
    TV L;
    TV c;
    TV c0;
    TV v;

    ~ObjectBox(){}

    ObjectBox(BC bc_in = BC::NoSlip, T friction_in = 0.0, std::string name_in = "box", bool force_calc_in = false, TV L_in = TV::Ones(), TV c_in = TV::Zero(), TV v_in = TV::Zero()) : ObjectGeneral(bc_in, friction_in, name_in, force_calc_in), L(L_in), c(c_in), c0(c_in), v(v_in) {}

    bool inside(const TV& X_in) const override {

        TV d = X_in - c;
        TV h = (T)0.5 * L; // half extension

        for (int i = 0; i < TV::RowsAtCompileTime; ++i) {
            if (std::abs(d[i]) > h[i]) {
                return false;
            }
        }
        return true;
    }

    TV normal(const TV& X_in) const override {

        TV d = X_in - c;
        TV h = (T)0.5 * L; // half extension

        TV n = TV::Zero();

        // Find the closest face
        int axis = 0;
        T max_dist = std::abs(d[0] / h[0]);

        for (int i = 1; i < TV::RowsAtCompileTime; ++i) {
            T dist = std::abs(d[i] / h[i]);

            if (dist > max_dist) {
                max_dist = dist;
                axis = i;
            }
        }

        // Outward normal of the closest face
        n[axis] = (d[axis] >= 0) ? (T)1 : (T)-1;

        return n;

    }

    void move(T time) override {
        c = c0 + v * time;
    }

    TV v_object(T time, const TV& X_in) const override {
        return v; // return zero if no velocity of object. Same as not defining the function.
    }



};

#endif  // OBJECTBOX_HPP
