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
    std::pair<Vector, double> force = std::make_pair(f, t);
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
        a_str = a_str_m_value + "\n"
            + a_str_position_value + "\n"
            + a_str_velocity_value;
    }
    return a_str;
};
::std::ostream& operator<<(::std::ostream& os, const Atom& a) {
    return os << atom_to_string(a).c_str();
};
Atom string_to_atoms(std::string nonverbose_str){
    std::stringstream ss(nonverbose_str.c_str());
    std::string line;
    std::vector<std::string> line_vector;
    while(std::getline(ss, line, '\n')){
      // do something with the line
      line_vector.push_back(line);
    }
    // retrieve information line-by-line;
    // see this
    //https://stackoverflow.com/questions/5607589/right-way-to-split-an-stdstring-into-a-vectorstring
};
