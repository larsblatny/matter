// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef BASAL_FRICTION_FIELD_HPP
#define BASAL_FRICTION_FIELD_HPP

#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include "tools.hpp"

// Assumes in 3D the horizontal plane is (x,z)
class BasalFrictionField{
public:
    BasalFrictionField(){}

    bool isSet() const { return initialized; }

#ifdef THREEDIM
    // values must be given row-major with z varying slowest, i.e.
    // values[k*nx + i] is the friction coefficient at grid point
    // (x0 + i*d, z0 + k*d). The grid spacing d is the same along both axes.
    void set(const std::vector<T>& values_in, unsigned int nx_in, unsigned int nz_in,
             T x0_in, T z0_in, T d_in){
        assert(nx_in >= 2 && nz_in >= 2);
        assert(d_in > 0);
        assert(values_in.size() == static_cast<size_t>(nx_in) * static_cast<size_t>(nz_in));

        values = values_in;
        nx = nx_in;
        nz = nz_in;
        x0 = x0_in;
        z0 = z0_in;
        d = d_in;
        inv_d = T(1) / d_in;
        initialized = true;
    } // end set

    // Bilinear interpolation at (x,z). Points outside clamped to the nearest edge
    T interpolate(T x, T z) const{
        assert(initialized);

        // Fractional grid coordinates, clamped to the data bounds (constant
        // extrapolation outside the grid).
        T fx = std::min(std::max((x - x0) * inv_d, T(0)), T(nx - 1));
        T fz = std::min(std::max((z - z0) * inv_d, T(0)), T(nz - 1));

        unsigned int i0 = static_cast<unsigned int>(std::floor(fx));
        unsigned int k0 = static_cast<unsigned int>(std::floor(fz));
        unsigned int i1 = std::min(i0 + 1, nx - 1);
        unsigned int k1 = std::min(k0 + 1, nz - 1);

        T tx = fx - static_cast<T>(i0);
        T tz = fz - static_cast<T>(k0);

        T f00 = values[k0 * nx + i0];
        T f10 = values[k0 * nx + i1];
        T f01 = values[k1 * nx + i0];
        T f11 = values[k1 * nx + i1];

        T f0 = f00 * (1 - tx) + f10 * tx; // interpolate along x at z0
        T f1 = f01 * (1 - tx) + f11 * tx; // interpolate along x at z1

        return f0 * (1 - tz) + f1 * tz; // interpolate along z
    } // end interpolate

    T interpolate(const TV& X) const { return interpolate(X(0), X(2)); }

#else // TWODIM

    // values must be given along a regular line in x, i.e. values[i] is the
    // friction coefficient at grid point (x0 + i*d).
    void set(const std::vector<T>& values_in, unsigned int nx_in, T x0_in, T d_in){
        assert(nx_in >= 2);
        assert(d_in > 0);
        assert(values_in.size() == static_cast<size_t>(nx_in));

        values = values_in;
        nx = nx_in;
        x0 = x0_in;
        d = d_in;
        inv_d = T(1) / d_in;
        initialized = true;
    } // end set

    // Linear interpolation at x. Points outside the data bounds are clamped to the nearest edge
    T interpolate(T x) const{
        assert(initialized);

        T fx = std::min(std::max((x - x0) * inv_d, T(0)), T(nx - 1));

        unsigned int i0 = static_cast<unsigned int>(std::floor(fx));
        unsigned int i1 = std::min(i0 + 1, nx - 1);

        T tx = fx - static_cast<T>(i0);

        return values[i0] * (1 - tx) + values[i1] * tx;
    } // end interpolate

    T interpolate(const TV& X) const { return interpolate(X(0)); }

#endif

private:
    std::vector<T> values;
    unsigned int nx = 0;
#ifdef THREEDIM
    unsigned int nz = 0;
    T z0 = 0;
#endif
    T x0 = 0, d = 1, inv_d = 1;
    bool initialized = false;
};

#endif // BASAL_FRICTION_FIELD_HPP
