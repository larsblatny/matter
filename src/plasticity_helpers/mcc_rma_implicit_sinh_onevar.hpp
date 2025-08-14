// Copyright (C) 2024 Tobias Verheijen, Lars Blatny. Released under GPL-3.0 license.

#ifndef MCCRMAIMPLICITSINHONEVAR_HPP
#define MCCRMAIMPLICITSINHONEVAR_HPP

#include "../tools.hpp"
#include "limited_search_sinh.hpp"

bool MCCRMAImplicitSinhOnevar(T& p, T& q, int& exit, T M, T p00, T beta, T mu, T K, T xi, T epv)
/**
 **** WRITTEN BY TOBIAS VERHEIJEN ****
 *
 * Implicit Return Mapping Algorithm:
 * This method uses a 1-d Newton-Raphson scheme iterating over *epv* values from the following set of equations and the
 * sinh hardening law for p0:
 *
 * p = (epv - epv^{n})K + p_trial
 * q = q_trial / (1 + 2 * mu * (1 + 2 * beta) * delta_gamma)
 * delta_gamma = (p_trial - p) / (2 * K * M^2 * p + M^2 * (beta - 1) * K * p0)
 * p0 = K * sinh(-xi*epv + sinh^{-1}(p00/K))
 * y = M^2 * (p - p0) * (p + beta * p0) + (1 + 2 * beta) * q^2
 *
 * @param p: Trial p: -K * Tr(epsilon)
 * @param q: Trial q: 2 * G * sqrt(dev(eps):dev(eps))
 * @param exit: Exit flag for program, if error encountered, exit set to 1
 * @param M: Slope of the critical state line
 * @param p00: Initial consolidation pressure value
 * @param beta: Cohesion parameter
 * @param mu:  G = Lame parameter
 * @param K: Bulk modulus, (lambda + 2*mu/dim)
 * @param xi: xi value for simulation
 * @param epv: Plastic Volumetric Hencky Strain
 * @return true if deformation is plastic, false if deformation is not plastic
 */
{
    T p0 = std::max(T(0.01), K * std::sinh(-xi * epv + std::asinh(p00 / K)));
    T Msq = M * M;
    // T damping_factor = 1. + 2. * beta;
    T damping_factor = 1;
    T y = Msq * (p*p + (beta - 1.) * p0*p - beta * p0*p0) + damping_factor * q*q;

    // If y <= 0 indicates that the p, q position falls within the yield surface
    if (y <= 0) {
        return false;
    }

    // store initial values and rescale
    T p0_t = p0;
    T scale_factor = 1./p0_t;
    T p0_min = 0.01 * scale_factor; p00 *= scale_factor; p0 = 1.; p *= scale_factor; q *= scale_factor; mu *= scale_factor; K *= scale_factor;

    T p_t = p;  T q_t = q;  T epv_t = epv;

    // Determine if trial point is right or left of apex of yield surface
    int right = 0;
    if (p_t > 0.5 * (p0 - beta*p0)) {
        right = 1;
    }

    // Newton-Raphson iterative root finding
    int max_iter = 40;
    for (int iter = 0;  iter < max_iter;  iter++) {
        T midpoint = (p0 - beta * p0) * 0.5;  // Find p-center of yield surface to avoid div-by-0 error
        if (std::abs(p - midpoint) < T(1.e-4)) {
            // if p_t == midpoint then we map q to top of yield surface since no plastic change
            q = M * std::sqrt((p0 - p) * (beta * p0 + p) / damping_factor);

            // rescale values
            p *= p0_t;
            q *= p0_t;
            return true;
        }

        p = (epv - epv_t) * K + p_t;  // compute current p value

        // Compute p0 value
        T p0_helper = K * std::sinh(-xi * epv + std::asinh(p00/K));
        p0 = std::max(p0_min, p0_helper);

        // define delta_gamma helper
        T denom = 2. * K * Msq * p + Msq * (beta - 1.) * K * p0;
        T dg_step = (p_t - p) / denom;  // compute current delta_gamma value

        // define q helper
        T q_denom = 1. + 2. * mu * damping_factor * dg_step;
        q = q_t / q_denom;  // compute current q value

        // compute yield surface
        y = Msq * (p - p0) * (p + beta * p0) + damping_factor * q * q;

        // Check for convergence, iteration > 3 check is done to ensure some displacement is measured.
        if (iter > 5 && std::abs(y) < 1e-4) {
            break;
        }

        if (iter == max_iter - 1) {  // did not break loop
            if (p0 > p0_min) {
                debug("RMA: FATAL did not exit loop at iteration = ", iter, ", iteration = ", iter);
                debug("y    = ", y);
                debug("p0   = ", p0);
                debug("K    = ", K);
                debug("G    = ", mu);
                debug("pt   = ", p_t);
                debug("qt   = ", q_t);
                debug("epvt = ", epv_t);
                debug("epv  = ", epv);
                debug("p    = ", p);
                debug("q    = ", q);
                debug("p0_t = ", p0_t);
                debug("scale factor = ", scale_factor); 
                exit = 1;
            }
        }

        // Compute derivative helpers
        T dp_depv = K;  // compute dp/depv

        // Compute dp0/depv, the if statement helps to make a continuous max() approximation
        T dp0_depv;
        T dp0_depv_helper = -xi * K * std::cosh(-xi * epv + std::asinh(p00 / K));
        if (p0_helper > p0_min){
            dp0_depv = 0.5 * (dp0_depv_helper + dp0_depv_helper - p0_min);
        } else {
            dp0_depv = 0.5 * (dp0_depv_helper + p0_min - dp0_depv_helper);
        }

        T ddelta_gamma_depv = (-dp_depv / denom) - (p_t - p) * (2. * K * Msq * dp_depv + Msq * (beta - 1.) * K * dp0_depv) \
 / (denom * denom);  // compute ddg/depv
        T dq_depv = -q_t * 2. * mu * damping_factor * ddelta_gamma_depv / (q_denom * q_denom);  // compute dq/depv

        // Get direction for epv update
        T dy_depv = Msq * (dp_depv - dp0_depv) * (p + beta * p0) + \
                     Msq * (p - p0) * (beta * dp0_depv + dp_depv) + \
                     2. * damping_factor * q * dq_depv;  // compute dy/depv

        // Take scaled-admissible Newton-Raphson step: x^{n+1} = x^{n} - scale * f(x^{n}) / f'(x^{n})
        epv = limited_search_sinh(epv, epv_t, y / dy_depv, xi, K, p_t, p00, p0_min, beta, right);

        if (epv == std::numeric_limits<double>::infinity()) {
            debug("RMA: FATAL no valid epv value found = ", iter, ", iteration = ", iter);
            debug(iter, ":  r1   = ", y);
            debug(iter, ":  y    = ", y);
            debug(iter, ":  p0   = ", p0);
            debug(iter, ":  p0_t   = ", p0_t);
            debug(iter, ":  pt   = ", p_t);
            debug(iter, ":  qt   = ", q_t);
            debug(iter, ":  epvt   = ", epv_t);
            debug(iter, ":  epv   = ", epv);
            debug(iter, ":  p    = ", p);
            debug(iter, ":  q    = ", q);
            debug(iter, ":  K    = ", K);
            debug(iter, ":  mu    = ", mu);
            exit = 1;
            break;
        }
    }

    // ensure p is within valid range
    p = std::max(p, -beta * p0);
    p = std::min(p, p0);

    // directly solve for q
    q = M * std::sqrt((p0 - p) * (beta*p0 + p)/damping_factor);

    // Rescale values to original scale
    p *= p0_t;
    q *= p0_t;

    return true;  // This indicates plasticity
}

#endif // MCCRMAIMPLICITSINHONEVAR_HPP
