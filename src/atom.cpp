// 2017 Artur Avkhadiev
/*! \file atom.cpp
*/
#include <utility>      /* std::pair, std::make_pair */
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "../include/vector.h"
#include <../include/atom.h>
Atom initialize_atom(double mass, Vector r, Vector v, double t){
    Vector f = {.x = 0.0, .y = 0.0, .z = 0.0};
    std::pair<Vector, double> position = std::make_pair(r, t);
    std::pair<Vector, double> velocity = std::make_pair(v, t);
    // force always has to be recalculated; time is not defined until then
    std::pair<Vector, double> force = std::make_pair(f, -1.0);
    Atom atom = {.mass = mass,
        .position = position,
        .velocity = velocity,
        .force = force};
    return atom;
}
bool is_time_consistent(Atom atom, double time) {
    double position_time = atom.position.second;
    double velocity_time = atom.velocity.second;
    bool consistent = false;
    if (position_time == velocity_time) {
        if (time == -1) {
            consistent = true;
        }
        else {
            if (time == position_time){
                consistent = true;
            }
        }
    }
    return consistent;
}
void set_time(Atom *atom, double time) {
    atom->position.second = time;
    atom->velocity.second = time;
}
std::string atom_to_string(Atom a, bool verbose) {
    std::string a_str_header = "ATOM:\n";
    // mass string
    std::string a_str_m_header = "m = ";
    std::string a_str_m_value = std::to_string(a.mass);
    // position string
    std::string a_str_position_header = "r(t):";
    std::string a_str_position_value_header = "r = ";
    std::string a_str_position_value = vector_to_string(a.position.first);
    std::string a_str_position_time_header = "t = ";
    std::string a_str_position_time = std::to_string(a.position.second);
    // velocity string
    std::string a_str_velocity_header = "v(t):";
    std::string a_str_velocity_value_header = "v = ";
    std::string a_str_velocity_value = vector_to_string(a.velocity.first);
    std::string a_str_velocity_time_header = "t = ";
    std::string a_str_velocity_time = std::to_string(a.velocity.second);
    // force string
    std::string a_str_force_header = "f(t):";
    std::string a_str_force_value_header = "f = ";
    std::string a_str_force_value = vector_to_string(a.force.first);
    std::string a_str_force_time_header = "t = ";;
    std::string a_str_force_time = std::to_string(a.force.second);
    // make a string
    std::string a_str;
    if (verbose){
        // include all headers
        a_str = a_str_header
            + a_str_m_header + a_str_m_value + "\n"
            + a_str_position_header + " "
                + a_str_position_time_header + a_str_position_time + " "
                + a_str_position_value_header + a_str_position_value + "\n"
            + a_str_velocity_header + " "
                + a_str_velocity_time_header + a_str_velocity_time + " "
                + a_str_velocity_value_header + a_str_velocity_value + "\n"
            + a_str_force_header + " "
                + a_str_force_time_header + a_str_force_time + " "
                + a_str_force_value_header + a_str_force_value;
    }
    else {
        // do not include headers nor times nor forces
        a_str = a_str_m_value + " "
            + a_str_position_value + " "
            + a_str_velocity_value;
    }
    return a_str;
};
::std::ostream& operator<<(::std::ostream& os, const Atom& a) {
    return os << atom_to_string(a).c_str();
};
Atom string_to_atom(std::string nonverbose_str){
    // in non-verbose lines, the information about an atom is one line
    std::istringstream ss(nonverbose_str.c_str());
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> words(begin, end);
    // convert strings to data
    double mass = atof(words.at(0).c_str());
    double rx = atof(words.at(1).c_str());
    double ry = atof(words.at(2).c_str());
    double rz = atof(words.at(3).c_str());
    double vx = atof(words.at(4).c_str());
    double vy = atof(words.at(5).c_str());
    double vz = atof(words.at(6).c_str());
    Vector r = {.x = rx, .y = ry, .z = rz};
    Vector v = {.x = vx, .y = vy, .z = vz};
    // force will be calculated later
    Vector f = {.x = 0.0, .y = 0.0, .z = 0.0};
    // time will be assigned later
    Atom atom = initialize_atom(mass, r, v);
    return atom;
};
