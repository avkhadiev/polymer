// 2017 Artur Avkhadiev
/*! \file atom.cpp
*/
#include <utility>      /* std::pair, std::make_pair */
#include <iostream>
#include <string>
#include "../include/vector.h"
#include <../include/atom.h>
Atom initialize_atom(double mass, Vector r, Vector v, double t){
    Vector f = {.x = 0.0, .y = 0.0, .z = 0.0};
    std::pair<Vector, double> position = std::make_pair(r, t);
    std::pair<Vector, double> velocity = std::make_pair(v, t);
    std::pair<Vector, double> force = std::make_pair(f, t);
    Atom atom = {.mass = mass,
        .position = position,
        .velocity = velocity,
        .force = force};
    return atom;
}
std::string atom_to_string(Atom a) {
    std::string a_str_header = "ATOM:\n";
    // mass string
    std::string a_str_m_value = "m = " + std::to_string(a.mass);
    std::string a_str_m = a_str_m_value + "\n";
    // position string
    std::string a_str_position_value = vector_to_string(a.position.first);
    std::string a_str_position_time = std::to_string(a.position.second);
    std::string a_str_position = "r(t): t = " +
        a_str_position_time + " r = "
        + a_str_position_value + "\n";
    // velocity string
    std::string a_str_velocity_value = vector_to_string(a.velocity.first);
    std::string a_str_velocity_time = std::to_string(a.velocity.second);
    std::string a_str_velocity = "v(t): t = "
        + a_str_velocity_time +
        " r = " + a_str_velocity_value + "\n";
    // force string
    std::string a_str_force_value = vector_to_string(a.force.first);
    std::string a_str_force_time = std::to_string(a.force.second);
    std::string a_str_force = "f(t): t = "
        + a_str_force_time
        + " r = " + a_str_force_value;
    std::string a_str = a_str_header
        + a_str_m
        + a_str_position
        + a_str_velocity
        + a_str_force;
    return a_str;
};
::std::ostream& operator<<(::std::ostream& os, const Atom& a) {
    return os << atom_to_string(a).c_str();
};
