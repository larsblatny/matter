// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"

void Simulation::updateDt(){

    auto max_velocity_it = std::max_element( particles.v.begin(), particles.v.end(),
                                             []( const TV &v1, const TV &v2 )
                                             {
                                                 return v1.squaredNorm() < v2.squaredNorm();
                                             } );
    T max_speed = (*max_velocity_it).norm();

    if (max_speed >= wave_speed){
        debug("               FYI the particle speed ", max_speed, " is larger than elastic wave speed ", wave_speed);
    }

#ifdef WARNINGS
    debug("               max_dt = ", max_dt);
#endif

    if (std::abs(max_speed) > 1e-10){
        T cfl_dt = cfl * dx / max_speed;
#ifdef WARNINGS
        debug("               cfl_dt = ", cfl_dt);
#endif
        dt = std::min(cfl_dt, max_dt);
    } else {
        dt = max_dt;
#ifdef WARNINGS
        debug("               cfl_dt = not computed, max_speed too low");
#endif
    }

    dt = std::min(dt, frame_dt*(frame+1) - time);
    dt = std::min(dt, final_time         - time);
    dt = std::max(dt, min_dt);

#ifdef WARNINGS
    debug("               dt     = ", dt);
#endif


    // Here one may hard-code a special gravity evolution with time
    // Example here: linear gravity increase until "gravity_time"
    if (gravity_special){

        if (time < gravity_time){
            gravity = gravity_final * time/gravity_time;
        }
        else{
            gravity = gravity_final;
            // if (no_liftoff){
            //     for(int p=0; p<Np; p++){
            //         if (particles.x[p](0) > 0.5*Lx){
            //             particles.x[p](0) -= 0.5*Lx;
            //             particles.x[p](1) -= Ly+10*dx;
            //         }
            //     }
            //     no_liftoff = false;
            }

    } // end gravity_special



} // end updateDt
