// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <assert.h>
#include <iomanip>
#include <string>


/////////////// USER PARAMETERS (GLOBAL) //////////
typedef double T; // float or double
#define THREEDIM // Uncomment for 2D
#define SPLINEDEG 2 // Quadratic B-spline
// #define WARNINGS // Write more debug info to screen
#define DIMMATTER
///////////////////////////////////////////////////



// Needed for tinyply (do NOT uncomment):
#define TINYPLY_IMPLEMENTATION

// Needed for OMP collapse (do NOT uncomment):
#ifdef THREEDIM
    #define DIMENSION 3
#else
    #define DIMENSION 2
#endif


#ifdef THREEDIM
    typedef Eigen::Matrix<T, 3, 3> TM;
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TMX;
    typedef Eigen::Matrix<T, 3, 1> TV;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> TVX;
    typedef Eigen::Array<T,3,1> TA;
#else
    typedef Eigen::Matrix<T, 2, 2> TM;
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TMX;
    typedef Eigen::Matrix<T, 2, 1> TV;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> TVX;
    typedef Eigen::Array<T,2,1> TA;
#endif
////////////////////////

enum class PlateType { top, bottom, left, right, front, back };
enum class ElasticModel { Hencky, NeoHookean };
enum class PlasticModel { NoPlasticity, VM, DP, DPSoft, MCC, VMVisc, DPVisc, MCCVisc, DPMui, MCCMui};
enum class HardeningLaw { NoHard, ExpoExpl, ExpoImpl, SinhExpl, SinhImpl };
enum class BC { NoSlip, SlipStick, SlipFree };

///////////////////// TOOLS ////////////////////////

template <typename T>
void debug(T in){
  std::cout << std::setprecision(12) << in << std::endl;
}
template <typename T, typename U>
void debug(T in1, U in2){
  std::cout << in1 << in2 << std::endl;
}
template <typename T, typename U, typename V>
void debug(T in1, U in2, V in3){
  std::cout << in1 << in2 << in3 << std::endl;
}
template <typename T, typename U, typename V, typename W>
void debug(T in1, U in2, V in3, W in4){
  std::cout << in1 << in2 << in3 << in4 << std::endl;
}
template <typename T, typename U, typename V, typename W, typename X>
void debug(T in1, U in2, V in3, W in4, X in5){
  std::cout << in1 << in2 << in3 << in4 << in5 << std::endl;
}
template <typename T, typename U, typename V, typename W, typename X, typename Y>
void debug(T in1, U in2, V in3, W in4, X in5, Y in6){
  std::cout << in1 << in2 << in3 << in4 << in5 << in6 << std::endl;
}
template <typename T, typename U, typename V, typename W, typename X, typename Y, typename Z>
void debug(T in1, U in2, V in3, W in4, X in5, Y in6, Z in7){
  std::cout << in1 << in2 << in3 << in4 << in5 << in6 << in7 << std::endl;
}

inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

inline T selfDoubleDot(TM& A){

    #ifdef THREEDIM
    T out = A(0,0)*A(0,0) + A(0,1)*A(0,1) + A(0,2)*A(0,2)
          + A(1,0)*A(1,0) + A(1,1)*A(1,1) + A(1,2)*A(1,2)
          + A(2,0)*A(2,0) + A(2,1)*A(2,1) + A(2,2)*A(2,2);
    #else
    T out = A(0,0)*A(0,0) + A(0,1)*A(0,1)
          + A(1,0)*A(1,0) + A(1,1)*A(1,1);
    #endif

    return out;
}

#if SPLINEDEG == 3

    inline T N(T u){
        T uabs = std::abs(u);
        if (uabs < 1.0){
            return 0.5 * uabs*uabs*uabs - uabs*uabs + 0.6666666666666666666666666666666666666666666666666;
        }
            else if (uabs < 2.0){
            return 0.1666666666666666666666666666666666666666666 * (2.0 - uabs) * (2.0 - uabs) * (2.0 - uabs);
        }
            else {
            return 0;
        }
    }

    inline T dNdu(T u){
        T uabs = std::abs(u);
        if (uabs < 1.0){
            return u * (1.5 * uabs - 2.0);
        }
        else if (uabs < 2.0){
            return -0.5 * sgn(u) * (2.0 - uabs) * (2.0 - uabs);
        }
        else {
            return 0;
        }
    }

    inline T d2Ndu2(T u){
        T uabs = std::abs(u);
        if (uabs < 1.0){
            return (3.0 * uabs - 2.0);
        }
        else if (uabs < 2.0){
            return (2.0 - uabs);
        }
        else {
            return 0;
        }
    }

#elif SPLINEDEG == 2 // QUADRATIC SPLINES

    inline T N(T u){
        T uabs = std::abs(u);
        if (uabs < 0.5){
            return 0.75 - uabs * uabs;
        }
    		else if (uabs < 1.5){
            return 0.5 * (1.5 - uabs) * (1.5 - uabs);
        }
    		else {
            return 0;
    	}
    }

    inline T dNdu(T u){
        T uabs = std::abs(u);
        if (uabs < 0.5){
            return (-2*u);
        }
    	else if (uabs < 1.5){
            return (u - 1.5*sgn(u));
        }
    	else {
            return 0;
    	}
    }

    inline T d2Ndu2(T u){
        T uabs = std::abs(u);
        if (uabs < 0.5){
            return -2.0;
        }
    	else if (uabs < 1.5){
            return 1;
        }
    	else {
            return 0;
    	}
    }

#elif SPLINEDEG == 1

inline T N(T u){
    T uabs = std::abs(u);
    return std::max(T(0), 1.0-uabs);
}

inline T dNdu(T u){
    T uabs = std::abs(u);
    if (uabs < 1.0){
        if (u > 0){
            return -1;
        } else{
            return 1;
        }
    }
    else {
        return 0;
    }
}

inline T d2Ndu2(T u){
    return 0; // NB: Not implemented
}

#else
 #error Unsupported spline degree
#endif


#ifdef THREEDIM

    inline T wip(T xp, T yp, T zp, T xi, T yi, T zi, T one_over_h){
        return N( (xp - xi) * one_over_h ) * N( (yp - yi) * one_over_h ) * N( (zp - zi) * one_over_h );
    }

    inline TV grad_wip(T xp, T yp, T zp, T xi, T yi, T zi, T one_over_h){
        TV out;
        out << dNdu((xp - xi) * one_over_h) * N((yp - yi) * one_over_h) * N((zp - zi) * one_over_h) * one_over_h,
               dNdu((yp - yi) * one_over_h) * N((xp - xi) * one_over_h) * N((zp - zi) * one_over_h) * one_over_h,
               dNdu((zp - zi) * one_over_h) * N((xp - xi) * one_over_h) * N((yp - yi) * one_over_h) * one_over_h;
        return out;
    }

    inline T laplace_wip(T xp, T yp, T zp, T xi, T yi, T zi, T one_over_h, T one_over_h_square){
        T term1 = d2Ndu2((xp - xi) * one_over_h) * N((yp - yi) * one_over_h) *  N((zp - zi) * one_over_h);
        T term2 = d2Ndu2((yp - yi) * one_over_h) * N((xp - xi) * one_over_h) *  N((zp - zi) * one_over_h);
        T term3 = d2Ndu2((zp - zi) * one_over_h) * N((xp - xi) * one_over_h) *  N((yp - yi) * one_over_h);
        return ( term1 + term2 + term3 ) * one_over_h_square;
    }

#else // TWODIM

    inline T wip(T xp, T yp, T xi, T yi, T one_over_h){
        return N( (xp - xi) * one_over_h ) * N( (yp - yi) * one_over_h );
    }

    inline TV grad_wip(T xp, T yp, T xi, T yi, T one_over_h){
        TV out;
        out << dNdu((xp - xi) * one_over_h) * N((yp - yi) * one_over_h) * one_over_h,
               dNdu((yp - yi) * one_over_h) * N((xp - xi) * one_over_h) * one_over_h;
        return out;
    }

    // inline T wip(T xp, T yp, T xi, T yi, T one_over_h){
    //     T Lx = 0.1;
    //     return N( fmod(xp - xi, Lx) * one_over_h ) * N( (yp - yi) * one_over_h );
    // }
    // inline TV grad_wip(T xp, T yp, T xi, T yi, T one_over_h){
    //     T Lx = 0.1;
    //     T xpminxi = xp-xi;
    //     if (xpminxi <= -Lx || xpminxi >= Lx)
    //         xpminxi = -fmod(xpminxi, Lx);
    //     TV out;
    //     out << dNdu( xpminxi  * one_over_h) * N((yp - yi) * one_over_h) * one_over_h,
    //            dNdu((yp - yi) * one_over_h) * N( xpminxi  * one_over_h) * one_over_h;
    //     return out;
    // }

    inline T laplace_wip(T xp, T yp, T xi, T yi, T one_over_h, T one_over_h_square){
        T term1 = d2Ndu2((xp - xi) * one_over_h) * N((yp - yi) * one_over_h);
        T term2 = d2Ndu2((yp - yi) * one_over_h) * N((xp - xi) * one_over_h);
        return ( term1 + term2 ) * one_over_h_square;
    }

#endif

// Taken from: https://gist.github.com/lorenzoriano/5414671
// linspace(a,b,N) return std::vector of length N in the closed interval [a,b], i.e., including b
inline std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

// Taken from: https://stackoverflow.com/questions/21216909/these-python-functions-in-c
// Works like numpy.arange, does NOT include stop value
inline std::vector<T> arange(T start, T stop, T step) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}


#endif  // TOOLS_HPP
