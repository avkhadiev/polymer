// 2017 Artur Avkhadiev
/*! \file simple_atom.cpp
*/
#include <string>
#include <sstream>
#include <vector>

#include "../include/vector.h"
#include "../include/simple_atom.h"
namespace simple{
    Atom::Atom(Vector position, Vector velocity, Vector force) :
        position (position),
        velocity (velocity),
        force (force){}
    Atom::~Atom(){}
    Atom& Atom::operator=(const Atom &other){
        position = other.position;
        velocity = other.velocity;
        force = other.force;
        return *this;
    }
    bool Atom::operator==(const Atom &other) const {
        return (position == other.position && velocity == other.velocity);
    }
    bool Atom::operator!=(const simple::Atom &other) const {
        return !(*this == other);
    }
    std::string Atom::to_string(bool verbose) const{
        std::string a_str_header = "SIMPLE ATOM:";
        // headers
        std::string a_str_position_header = "r:";
        std::string a_str_velocity_header = "v:";
        std::string a_str_force_header = "f:";
        // values
        std::string a_str_position_value = vector_to_string(position);
        std::string a_str_velocity_value = vector_to_string(velocity);
        std::string a_str_force_value = vector_to_string(force);
        // make a string
        std::string a_str;
        if (verbose){
            // include all headers
            a_str = a_str_header + "\n"
                + a_str_position_header + " " + a_str_position_value + "\n"
                + a_str_velocity_header + " " + a_str_velocity_value + "\n"
                + a_str_force_header + " " + a_str_force_value;
        }
        else {
            // do not include headers or forces
            a_str = a_str_position_value;// + " " + a_str_velocity_value;
        }
        return a_str;
    }
    Atom string_to_atom(std::string nonverbose_str){
        // in non-verbose lines, the information about an atom is one line
        std::istringstream ss(nonverbose_str.c_str());
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> words(begin, end);
        double rx, ry, rz;
        double vx, vy, vz;
        // convert strings to data
        rx = atof(words.at(0).c_str());
        ry = atof(words.at(1).c_str());
        rz = atof(words.at(2).c_str());
        if (words.size() == 6){
            vx = atof(words.at(3).c_str());
            vy = atof(words.at(4).c_str());
            vz = atof(words.at(5).c_str());
        }
        else{
            vx = vy = vz = 0.0;
        }
        Vector r = vector(rx, ry, rz);
        Vector v = vector(vx, vy, vz);
        return Atom(r, v);
    }
    ::std::ostream& operator<<(::std::ostream& os, const simple::Atom& a){
        // output at default verbosity
        return os << a.to_string().c_str();
    }
}   // namespace simple
