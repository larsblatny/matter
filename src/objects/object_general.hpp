// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef OBJECTGENERAL_HPP
#define OBJECTGENERAL_HPP

#include "../tools.hpp"

class ObjectGeneral{
public:

    BC bc;
    T friction;
    std::string name;

    ObjectGeneral(BC bc, T friction, std::string name, TV v_object = TV::Zero()) : bc(bc), friction(friction), name(name) {}

    virtual ~ObjectGeneral(){}

    virtual bool inside(const TV& X_in) const = 0;

    virtual TV normal(const TV& X_in) const = 0;

    virtual void move(T time) {
        // do nothing by default
    }

    virtual TV v_object(T time, const TV& X_in) const {
        return TV::Zero(); // return zero by default
    }

};

#endif  // OBJECTGENERAL_HPP
