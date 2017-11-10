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
    Potential(),
    _epsilon (1.0),
    _sigma   (1.0),
    _sigmasq (pow(_sigma, 2.0)){};
LJPotential::LJPotential(double epsilon, double sigma) :
    Potential(),          
    _epsilon (epsilon),
    _sigma (sigma),
    _sigmasq (pow(_sigma, 2.0)){};
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
    double frac = pow(_sigmasq * inv_rijsq, 3.0);
    pair_potential = 4.0 * (pow(frac, 2.0) - frac);
    return pair_potential;
};
double LJPotential::calculate_neg_pair_virial(double inv_rijsq) {
    double neg_pair_virial;
    double frac = pow(_sigmasq * inv_rijsq, 3.0);
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
/** from Stoddard and Ford, 1973. The advantage of this form is that
  * one does not need to compute a square root of r in adjusting the
  * potential
  */
double AdjustedLJPotential::_calculate_corr1(){
    double e = LJPotential::_epsilon;
    double ssq = LJPotential::_sigmasq;
    double frac = pow(ssq / _rcsq, 3.0);
    double frac2 = pow(frac, 2.0);
    double corr1 = - 4.0 * e * (7.0 * frac2 - 4.0 * frac);
    return corr1;
};
double AdjustedLJPotential::_calculate_corr2(){
    double e = LJPotential::_epsilon;
    double ssq = LJPotential::_sigmasq;
    double frac = pow(ssq / _rcsq, 3.0);
    double frac2 = pow(frac, 2.0);
    double corr2 = 4.0 * e * (6.0 * frac2 - 3.0 * frac);
    return corr2 ;
};
AdjustedLJPotential::AdjustedLJPotential(double cutoff) :
    LJPotential(),
    _rc(cutoff),
    _rcsq(pow(cutoff, 2.0)),
    _corr1(_calculate_corr1()),
    _corr2(_calculate_corr2()){
};
AdjustedLJPotential::AdjustedLJPotential(double epsilon,
    double sigma, double cutoff) :
    LJPotential(epsilon, sigma),
    _rc(cutoff),
    _rcsq(pow(cutoff, 2.0)),
    _corr1(_calculate_corr1()),
    _corr2(_calculate_corr2()){
};
AdjustedLJPotential::~AdjustedLJPotential() {
};
std::string AdjustedLJPotential::get_str() const {
    std::string cutoff_str = "cutoff" + std::to_string(_rc);
    std::string potential_str = LJPotential::get_str() + cutoff_str + "\n";
    return potential_str;
};
double AdjustedLJPotential::calculate_pair_potential(double inv_rijsq) {
    double unadj = LJPotential::calculate_pair_potential(inv_rijsq);
    double adjst = _corr1 + _corr2 * (1.0 / (inv_rijsq * _rcsq));
    double pair_potential = unadj + adjst;
    return pair_potential;
};
double AdjustedLJPotential::calculate_neg_pair_virial(double inv_rijsq) {
    double unadj = LJPotential::calculate_neg_pair_virial(inv_rijsq);
    double adjst = 2.0 * _corr2 * (1.0 / (inv_rijsq * _rcsq));
    double neg_pair_virial = unadj + adjst;
    return neg_pair_virial;
}
double AdjustedLJPotential::calculate_fstrength_over_r(double inv_rijsq){
    double neg_pair_virial = calculate_neg_pair_virial(inv_rijsq);
    double fstrength_over_r = -inv_rijsq * neg_pair_virial;
    return fstrength_over_r;
}
void AdjustedLJPotential::writeout_parameters_to_file(std::string outdir,
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
        writeout << "\\sigma" << delim;
        writeout << "cutoff" << std::endl;
        // writeout the parameters
        writeout << get_epsilon() << delim;
        writeout << get_sigma() << delim;
        writeout << get_rc() << std::endl;
    }
    else {
        // if file could not be opened
        std::string err_msg = "AdjustedLJPotential::writeout_parameters_to_file: unable to open file at";
        fprintf(stderr, "%s %s\n", err_msg.c_str(), fout.c_str());
        perror("open");
    }
    writeout.close();
};
