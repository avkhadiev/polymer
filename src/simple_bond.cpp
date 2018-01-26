// 2017 Artur Avkhadiev
/*! \file simple_bond.cpp
*/
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include "../include/vector.h"
#include "../include/simple_atom.h"
#include "../include/simple_bond.h"
namespace simple {
    Bond::Bond(Vector new_position, Vector new_velocity) :
        position (divide(new_position, norm(new_position))),
        velocity (new_velocity) {}
    Bond::Bond(Atom start, Atom end)
    {
        Vector new_position = subtract(end.position, start.position);
        Vector new_velocity = subtract(end.velocity, start.velocity);
        double bondlength = norm(new_position);
        position = divide(new_position, bondlength);
        velocity = divide(new_velocity, bondlength);
    }
    Bond::~Bond(){}
    Bond& Bond::operator=(const Bond &other){
        position = other.position;
        velocity = other.velocity;
        return *this;
    }
    bool Bond::operator==(const Bond &other) const {
        return (position == other.position && velocity == other.velocity);
    }
    bool Bond::operator!=(const simple::Bond &other) const {
        return !(*this == other);
    }
    std::string Bond::to_string(bool verbose) const{
        std::string b_str_header = "SIMPLE BOND:";
        // headers
        std::string b_str_position_header = "omega_r:";
        std::string b_str_velocity_header = "omega_v:";
        // values
        std::string b_str_position_value = vector_to_string(position);
        std::string b_str_velocity_value = vector_to_string(velocity);
        // make a string
        std::string b_str;
        if (verbose){
            // include all headers
            b_str = b_str_header + "\n"
                + b_str_position_header + " " + b_str_position_value + "\n"
                + b_str_velocity_header + " " + b_str_velocity_value;
        }
        else {
            b_str = b_str_position_value;// + " " + b_str_velocity_value;
        }
        return b_str;
    }
    Bond string_to_bond(std::string nonverbose_str){
        // in non-verbose lines, the information about a bond is one line
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
        // if velocities are recorded
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
        return Bond(r, v);      // position will be normalized
    }
    ::std::ostream& operator<<(::std::ostream& os, const simple::Bond& b){
        // output at default verbosity
        return os << b.to_string().c_str();
    }
} // namespace simple
