// Copyright (C) 2024 Tobias Verheijen, Lars Blatny. Released under GPL-3.0 license.

#ifndef MCCRMAEXPLICITONEVAR_HPP
#define MCCRMAEXPLICITONEVAR_HPP

#include "../tools.hpp"

bool MCCRMAExplicitOnevar(T& p, T& q, int& exit, T M, T p0, T beta, T mu, T K)
/**
 **** WRITTEN BY TOBIAS VERHEIJEN ****
 *
 * @param p: p_stress: -K * Tr(epsilon)
 * @param q: q_stress: 2 * G * sqrt(dev(eps):dev(eps))
 * @param exit: Exit flag for program, if error encountered, exit set to 1
 * @param M: Slope of the critical state line
 * @param p0: Consolidation pressure
 * @param beta: Cohesion parameter
 * @param mu: G = Lame parameter
 * @param K: Bulk modulus, (lambda + 2*mu/dim)
 * @return true if deformation is plastic, false if deformation is not plastic
 */
{
    T Msq = M * M;
    // T damping_factor = 1. + 2.*beta;
    T damping_factor = 1;
    T y = Msq * (p - p0) * (p + beta * p0) + damping_factor * (q * q);

    // If y <= 0 indicates that the p, q position falls within the yield surface
    if (y <= 0) {
        return false;
    }

    T delta_gamma = 0;
    T p0_t = p0;
    T scale_factor = 1./p0_t;
    p0 = 1.; p *= scale_factor; q *= scale_factor; mu *= scale_factor; K *= scale_factor;
    T pt = p;  T qt = q;  T bm1 = beta - 1.;

    int max_iter = 100;
    for (int iter = 0;  iter < max_iter;  iter++) {
        // Compute new p and q values
        p = (pt - Msq * K * delta_gamma * bm1 * p0) / (1. + 2. * Msq * K * delta_gamma);
        q = qt / (1. + 2. * damping_factor * mu * delta_gamma);

        // Recompute y to see if new (p, q) on yield surface
        y = Msq * (p - p0) * (p + beta * p0) + damping_factor * (q * q);

        if (iter > 3 && y < 1.e-6) {  // Check for convergence
            break;
        }

        if (iter == max_iter - 1){ // did not break loop
            debug("RMA: FATAL did not exit loop at iter = ", iter);
            debug(iter, ":  y    = ", y);
            debug(iter, ":  p0   = ", p0);
            debug(iter, ":  pt   = ", pt);
            debug(iter, ":  qt   = ", qt);
            debug(iter, ":  p    = ", p);
            debug(iter, ":  q    = ", q);
            exit = 1;
            break;
        }

        // Compute derivative helpers
        T dp_ddelta_gamma = (-2. * Msq * K * pt - Msq * K * bm1 * p0) / ((1. + 2. * Msq * K * delta_gamma) * (1. + 2. * Msq * K * delta_gamma));
        T dq_ddelta_gamma = -2. * damping_factor * mu * qt / ((1. + 2. * damping_factor * mu * delta_gamma) * (1. + 2. * damping_factor * mu * delta_gamma));

        // Take Newton-Raphson Step
        T dyddelta_gamma = Msq * (2. * p * dp_ddelta_gamma + bm1 * p0 * dp_ddelta_gamma) + 2. * damping_factor * q * dq_ddelta_gamma;
        delta_gamma -= y / dyddelta_gamma;
    }

    // ensure p is within valid range
    p = std::max(p, -beta * p0);
    p = std::min(p, p0);

    // directly solve for q
    q = M * std::sqrt((p0 - p) * (beta*p0 + p) / damping_factor);

    // Rescale values to original scale
    p *= p0_t;
    q *= p0_t;

    return true;
}


#endif // MCCRMAEXPLICITONEVAR_HPP
