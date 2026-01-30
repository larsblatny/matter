// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "simulation.hpp"

void Simulation::boundaryCollision(int index, TV Xi, TV& vi){

    // Make a copy
    TV vi_orig = vi;

    // New positions
    Xi += dt * vi; // recall Xi was passed by value

    for(auto& obj : objects) {
        bool colliding = obj->inside(Xi);
        if (colliding) {
            TV v_object = obj->v_object(time, Xi);
            TV v_rel = vi_orig - v_object;

            if (obj->bc == BC::NoSlip) {
                v_rel.setZero();
            } // end BC::NoSlip

            else if (obj->bc == BC::SlipFree) {
                TV n = obj->normal(Xi);
                T dot = v_rel.dot(n);
                if (dot < 0){ // if moving towards object

                    T friction = obj->friction;
                    if (use_mibf)
                        friction = grid.friction[index];

                    TV v_tang = v_rel - dot * n;
                    if (friction > 0){
                        if( -dot * friction < v_tang.norm() )
                            v_rel = v_tang + v_tang.normalized() * dot * friction;
                        else
                            v_rel.setZero();
                    } else{
                        v_rel = v_tang;
                    } // end non-zero friction
                } // end moving towards object

            } // end BC::SlipFree

            else if (obj->bc == BC::SlipStick) {
                TV n = obj->normal(Xi);
                T dot = v_rel.dot(n);

                v_rel = v_rel - dot * n; //  v_rel is now v_tang

                if (dot < 0){ // if moving towards object

                    T friction = obj->friction;
                    if (use_mibf)
                        friction = grid.friction[index];

                    if (friction > 0){
                        if( -dot * friction < v_rel.norm() )
                            v_rel = v_rel + v_rel.normalized() * dot * friction;
                        else
                            v_rel.setZero();
                    } // end non-zero friction
                } // end moving towards object

            } // end BC::SlipStick

            else {
                debug("INVALID BOUNDARY CONDITION!!!");
                exit = 1;
                return;
            }

            vi = v_rel + v_object;

            // update velocity copy before next iteration
            vi_orig = vi; // Comment this line to enforce ordering of objects (i.e., use only last object in list)

        } // end if colliding

    } // end iterator over general objects


#ifdef THREEDIM

    for (auto& obj : plates) {
        bool colliding = obj->inside(Xi);
        if (colliding) {
            T vx_rel = vi_orig(0) - obj->vx_object;
            T vy_rel = vi_orig(1) - obj->vy_object;
            T vz_rel = vi_orig(2) - obj->vz_object;

            if (obj->bc == BC::NoSlip) {
                vx_rel = 0;
                vy_rel = 0;
                vz_rel = 0;
            } // end BC::NoSlip

            else if (obj->bc == BC::SlipStick) {

                T friction = obj->friction;
                    if (use_mibf)
                        friction = grid.friction[index];

                if (obj->plate_type == PlateType::top || obj->plate_type == PlateType::bottom){
                    // tangential velocity is the (x,z) components

                    T vel_t      = std::sqrt(vx_rel*vx_rel + vz_rel*vz_rel);
                    T fric_vel_n = friction * std::abs(vy_rel);

                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( vel_t <= fric_vel_n ){
                        vx_rel = 0; // tangential component also set to zero
                        vz_rel = 0;
                    }
                    else { // just reduce tangential component
                        vx_rel = sgn(vx_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vz_rel*vz_rel/(vx_rel*vx_rel));
                        vz_rel = sgn(vz_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vx_rel*vx_rel/(vz_rel*vz_rel));
                    }

                    // normal component (y) must be set to zero
                    vy_rel = 0;

                } // end PlateType::top and PlateType::bottom plate
                else if (obj->plate_type == PlateType::left || obj->plate_type == PlateType::right){
                    // tangential velocity is the (y,z) components

                    T vel_t      = std::sqrt(vy_rel*vy_rel + vz_rel*vz_rel);
                    T fric_vel_n = friction * std::abs(vx_rel);

                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( vel_t <= fric_vel_n ){
                        vy_rel = 0; // tangential component also set to zero
                        vz_rel = 0;
                    }
                    else { // just reduce tangential component
                        vy_rel = sgn(vy_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vz_rel*vz_rel/(vy_rel*vy_rel));
                        vz_rel = sgn(vz_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vy_rel*vy_rel/(vz_rel*vz_rel));
                    }

                    // normal component (x) must be set to zero
                    vx_rel = 0;

                } // end PlateType::left or PlateType::right plate
                else if (obj->plate_type == PlateType::front || obj->plate_type == PlateType::back){
                    // tangential velocity is the (x,y) components

                    T vel_t      = std::sqrt(vx_rel*vx_rel + vy_rel*vy_rel);
                    T fric_vel_n = friction * std::abs(vz_rel);

                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( vel_t <= fric_vel_n ){
                        vx_rel = 0; // tangential component also set to zero
                        vy_rel = 0;
                    }
                    else { // just reduce tangential component
                        vx_rel = sgn(vx_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vy_rel*vy_rel/(vx_rel*vx_rel));
                        vy_rel = sgn(vy_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vx_rel*vx_rel/(vy_rel*vy_rel));
                    }

                    // normal component (z) must be set to zero
                    vz_rel = 0;

                } // end PlateType::front or PlateType::back plate
            } // end BC::SlipStick



            else if (obj->bc == BC::SlipFree) {

                T friction = obj->friction;
                    if (use_mibf)
                        friction = grid.friction[index];

                if ((obj->plate_type == PlateType::top && vy_rel > 0) || (obj->plate_type == PlateType::bottom && vy_rel < 0)){
                    // tangential velocity is the (x,z) components
                    T vel_t      = std::sqrt(vx_rel*vx_rel + vz_rel*vz_rel);
                    T fric_vel_n = friction * std::abs(vy_rel);

                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( vel_t <= fric_vel_n ){
                        vx_rel = 0; // tangential component also set to zero
                        vz_rel = 0;
                    }
                    else { // just reduce tangential component
                        vx_rel = sgn(vx_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vz_rel*vz_rel/(vx_rel*vx_rel));
                        vz_rel = sgn(vz_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vx_rel*vx_rel/(vz_rel*vz_rel));
                    }

                    // normal component (y) must be set to zero
                    vy_rel = 0;

                } // end PlateType::top and PlateType::bottom plate
                else if ((obj->plate_type == PlateType::left && vx_rel < 0) || (obj->plate_type == PlateType::right && vx_rel > 0)){
                    // tangential velocity is the (y,z) components
                    T vel_t      = std::sqrt(vy_rel*vy_rel + vz_rel*vz_rel);
                    T fric_vel_n = friction * std::abs(vx_rel);

                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( vel_t <= fric_vel_n ){
                        vy_rel = 0; // tangential component also set to zero
                        vz_rel = 0;
                    }
                    else { // just reduce tangential component
                        vy_rel = sgn(vy_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vz_rel*vz_rel/(vy_rel*vy_rel));
                        vz_rel = sgn(vz_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vy_rel*vy_rel/(vz_rel*vz_rel));
                    }

                    // normal component (x) must be set to zero
                    vx_rel = 0;

                } // end PlateType::left or PlateType::right plate
                else if ((obj->plate_type == PlateType::front && vz_rel > 0) || (obj->plate_type == PlateType::back && vz_rel < 0 )){
                    // tangential velocity is the (x,y) components
                    T vel_t      = std::sqrt(vx_rel*vx_rel + vy_rel*vy_rel);
                    T fric_vel_n = friction * std::abs(vz_rel);

                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( vel_t <= fric_vel_n ){
                        vx_rel = 0; // tangential component also set to zero
                        vy_rel = 0;
                    }
                    else { // just reduce tangential component
                        vx_rel = sgn(vx_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vy_rel*vy_rel/(vx_rel*vx_rel));
                        vy_rel = sgn(vy_rel) * (vel_t - fric_vel_n) / std::sqrt(1 + vx_rel*vx_rel/(vy_rel*vy_rel));
                    }

                    // normal component (z) must be set to zero
                    vz_rel = 0;

                } // end PlateType::front or PlateType::back plate
            } // end BC::SlipFree



            else {
                debug("INVALID BOUNDARY CONDITION!!!");
                exit = 1;
                return;
            }

            vi(0) = vx_rel + obj->vx_object;
            vi(1) = vy_rel + obj->vy_object;
            vi(2) = vz_rel + obj->vz_object;

            // update velocity copy before next iteration
            vi_orig = vi; // Comment this line to enforce ordering of objects (i.e., use only last object in list)

        } // end if colliding

    } // end iterator over 3D plate objects


#else // TWODIM


    for (auto& obj : plates) {
        bool colliding = obj->inside(Xi);
        if (colliding) {
            T vx_rel = vi_orig(0) - obj->vx_object;
            T vy_rel = vi_orig(1) - obj->vy_object;

            if (obj->bc == BC::NoSlip) {
                vx_rel = 0;
                vy_rel = 0;
            } // end BC::NoSlip

            else if (obj->bc == BC::SlipStick) {

                T friction = obj->friction;
                    if (use_mibf)
                        friction = grid.friction[index];

                if (obj->plate_type == PlateType::top || obj->plate_type == PlateType::bottom){
                    // tangential velocity is the (x) component and must be changed
                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( std::abs(vx_rel) < friction * std::abs(vy_rel) ){
                        vx_rel = 0; // tangential component also set to zero
                    }
                    else { // just reduce tangential component
                        if (vx_rel > 0)
                            vx_rel -= friction * std::abs(vy_rel);
                        else
                            vx_rel += friction * std::abs(vy_rel);
                    }

                    // normal component (y) must be set to zero
                    vy_rel = 0;

                } // end PlateType::top and PlateType::bottom plate
                else if (obj->plate_type == PlateType::left || obj->plate_type == PlateType::right){
                    // tangential velocity is the (y) component and must be changed

                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( std::abs(vy_rel) < friction * std::abs(vx_rel) ){
                        vy_rel = 0; // tangential component also set to zero
                    }
                    else { // just reduce tangential component
                        if (vy_rel > 0)
                            vy_rel -= friction * std::abs(vx_rel);
                        else
                            vy_rel += friction * std::abs(vx_rel);
                    }

                    // normal component (x) must be set to zero
                    vx_rel = 0;

                } // end PlateType::left or PlateType::right plate

            } // end BC::SlipStick



            else if (obj->bc == BC::SlipFree) {

                T friction = obj->friction;
                    if (use_mibf)
                        friction = grid.friction[index];

                if ((obj->plate_type == PlateType::top && vy_rel > 0) || (obj->plate_type == PlateType::bottom && vy_rel < 0)){
                    // tangential velocity is the (x) component and must be changed
                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( std::abs(vx_rel) < friction * std::abs(vy_rel) ){
                        vx_rel = 0; // tangential component also set to zero
                    }
                    else { // just reduce tangential component
                        if (vx_rel > 0)
                            vx_rel -= friction * std::abs(vy_rel);
                        else
                            vx_rel += friction * std::abs(vy_rel);
                    }

                    // normal component (y) must be set to zero
                    vy_rel = 0;

                } // end PlateType::top and PlateType::bottom plate
                else if ((obj->plate_type == PlateType::left && vx_rel < 0) || (obj->plate_type == PlateType::right && vx_rel > 0)){
                    // tangential velocity is the (y) component and must be changed
                    if (friction == 0){
                        // Do nothing
                    }
                    else if ( std::abs(vy_rel) < friction * std::abs(vx_rel) ){
                        vy_rel = 0; // tangential component also set to zero
                    }
                    else { // just reduce tangential component
                        if (vy_rel > 0)
                            vy_rel -= friction * std::abs(vx_rel);
                        else
                            vy_rel += friction * std::abs(vx_rel);
                    }

                    // normal component (x) must be set to zero
                    vx_rel = 0;

                } // end PlateType::left or PlateType::right plate
            } // end BC::SlipFree



            else {
                debug("INVALID BOUNDARY CONDITION!!!");
                exit = 1;
                return;
            }

            vi(0) = vx_rel + obj->vx_object;
            vi(1) = vy_rel + obj->vy_object;

            // update velocity copy before next iteration 
            vi_orig = vi; // Comment this line to enforce ordering of objects (i.e., use only last object in list)

        } // end if colliding

    } // end iterator over 2D plate objects

#endif

} // end boundaryCollision
