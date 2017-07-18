// 2017 Artur Avkhadiev
/*! \file ljpotential.cpp
* does NOT measure energy, distance in units of epsilon, sigma
* (they are included in the calculations)
*/
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "../include/parsing.h"
#include "../include/ljpotential_not_unity.h"
#include "../include/vector.h"
LJPotential::LJPotential() :
    _epsilon (1.0),
    _epsilon4 (4.0 * _epsilon),
    _epsilon24 (24.0 * _epsilon),
    _sigma   (1.0),
    _sigma6  (pow(_sigma, 6.0)),
    _sigma12 (pow(_sigma6, 2.0)),
    _calculate_prefactors (true) {};
LJPotential::LJPotential(double epsilon,
    double sigma,
    bool calculate_prefactors) :
    _epsilon (epsilon),
    _epsilon4 (4.0 * _epsilon),
    _epsilon24 (24.0 * _epsilon),
    _sigma   (sigma),
    _sigma6  (pow(_sigma, 6.0)),
    _sigma12 (pow(_sigma6, 2.0)),
    _calculate_prefactors (calculate_prefactors) {};
LJPotential::~LJPotential() {
};
double LJPotential::get_epsilon() const {
    return _epsilon;
};
double LJPotential::get_sigma() const {
    return _sigma;
};
bool LJPotential::get_calculate_prefactors(){
    return _calculate_prefactors;
};
void LJPotential::set_epsilon(double epsilon){
    _epsilon = epsilon;
    _epsilon4 = 4.0 * epsilon;
    _epsilon24 = 24.0 * epsilon;
};
void LJPotential::set_sigma(double sigma){
    _sigma = sigma;
    _sigma6  = pow(_sigma, 6.0);
    _sigma12 = pow(_sigma6, 2.0);
};
void LJPotential::set_calculate_prefactors(bool calculate_prefactors){
    _calculate_prefactors = calculate_prefactors;
};
double LJPotential::calculate_pair_potential(double inv_rijsq) {
    double pair_potential;
    double frac = _sigma6 * pow(inv_rijsq, 3.0);
    // Allen & Tildesley Eq. 1.6
    if (_calculate_prefactors){
        pair_potential = _epsilon4 * (pow(frac, 2.0) - frac);
    }
    else {
        pair_potential = (pow(frac, 2.0) - frac);
    }
    return pair_potential;
};
double LJPotential::calculate_neg_pair_virial(double inv_rijsq) {
    double pair_virial;
    double frac = _sigma6 * pow(inv_rijsq, 3.0);
    // Allen & Tildesley Eq. 2.59 - 2.63
    if (_calculate_prefactors){
        pair_virial = _epsilon24 * (2 * pow(frac, 2.0) - frac);
    }
    else {
        pair_virial = (2 * pow(frac, 2.0) - frac);
    }
    return pair_virial;
};
double LJPotential::calculate_fstrength_over_r(double inv_rijsq) {
    // Allen & Tildesley Eq. 5.3
    double neg_pair_virial = calculate_neg_pair_virial(inv_rijsq);
    double fstrength_over_r = inv_rijsq * neg_pair_virial;
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
void LJPotential::read_parameters_from_file(std::string indir,
    std::string sim_name) {
    // number of paramaters = 2
    int nparams = 2;
    double epsilon;
    double sigma;
    // stores path to input file
    std::string fin;
    fin = indir + sim_name + "_potential.cfg";
    std::ifstream readout;
    // stores parameter lines
    std::string parameter_line;
    // for storing (epsilon, sigma)
    std::vector<double> parameters;
    // try to open the file
    readout.open(fin, std::ifstream::in);
    if (!readout.is_open()) {
        // if the file could not be opened
        std::string err_msg = "read_state_to_file: unable to open file at";
        fprintf(stderr, "%s %s\n", err_msg.c_str(), fin.c_str());
        perror("open");
    }
    else {
        // if the file was opened...
        // get parameters
        for (int i = 0; i < nparams; ++i) {
            std::getline(readout, parameter_line);
            // break parameter line into words (name value)
            std::istringstream ss(parameter_line.c_str());
            std::istream_iterator<std::string> begin(ss);
            std::istream_iterator<std::string> end;
            std::vector<std::string> words(begin, end);
            parameters.push_back(atof(words.at(1).c_str()));
        }
        epsilon = parameters.at(0);
        sigma = parameters.at(1);
        // set parameters
        set_epsilon(epsilon);
        set_sigma(sigma);
    }
};
