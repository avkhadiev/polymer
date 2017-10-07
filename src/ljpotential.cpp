// 2017 Artur Avkhadiev
/*! \file ljpotential.cpp
*/
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "../include/parsing.h"
#include "../include/ljpotential.h"
#include "../include/vector.h"
LJPotential::LJPotential() :
    _epsilon (1.0),
    _sigma   (1.0) {};
LJPotential::~LJPotential() {
};
double LJPotential::get_epsilon() const {
    return _epsilon;
};
double LJPotential::get_sigma() const {
    return _sigma;
};
double LJPotential::calculate_pair_potential(double inv_rijsq) {
    double pair_potential;
    double frac = pow(inv_rijsq, 3.0);
    pair_potential = 4.0 * (pow(frac, 2.0) - frac);
    return pair_potential;
};
double LJPotential::calculate_neg_pair_virial(double inv_rijsq) {
    double neg_pair_virial;
    double frac = pow(inv_rijsq, 3.0);
    neg_pair_virial = - 24.0 * (2 * pow(frac, 2.0) - frac);
    return neg_pair_virial;
};
double LJPotential::calculate_fstrength_over_r(double inv_rijsq) {
    // Allen & Tildesley Eq. 5.3
    double neg_pair_virial = calculate_neg_pair_virial(inv_rijsq);
    double fstrength_over_r = -inv_rijsq * neg_pair_virial;
    return fstrength_over_r;
};
std::string LJPotential::get_str() const {
    std::string potential_str;
    std::string epsilon_str = "\\epsilon " + std::to_string(_epsilon);
    std::string sigma_str = "\\sigma " + std::to_string(_sigma);
    potential_str = epsilon_str + "\n" + sigma_str + "\n";
    return potential_str;
};
::std::ostream& operator<<(::std::ostream& os, const LJPotential& potential){
    return os << potential.get_str().c_str();
}
void LJPotential::writeout_parameters_to_file(std::string outdir,
    std::string sim_name) {
    // delimeter for csv file
    std::string delim = ",";
    // stores path to output file
    std::string fout;
    fout = outdir + parse_string(sim_name) + "_potential.csv";
    std::ofstream writeout;
    // open writeout for output operations and s
    // try to open in truncate mode
    writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
    // now file has to be opened
    if (writeout.is_open()) {
        // if file could be opened...
        // writeout potential's string representation
        // writeout the header
        writeout << "\\epsilon" << delim;
        writeout << "\\sigma" << std::endl;
        // writeout the parameters
        writeout << get_epsilon() << delim;
        writeout << get_sigma() << std::endl;
    }
    else {
        // if file could not be opened
        std::string err_msg = "LJPotential::writeout_parameters_to_file: unable to open file at";
        fprintf(stderr, "%s %s\n", err_msg.c_str(), fout.c_str());
        perror("open");
    }
    writeout.close();
};
