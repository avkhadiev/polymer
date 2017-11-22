// 2017 Artur Avkhadiev
/*! \file ljpotential.cpp
*/
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>                                   /* round */
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
// TODO
double LJPotential::rijsq(Vector rij){
    double rijsq = normsq(rij);
    if (rijsq == 0){
        std::string err_msg = "potential: distance between atoms is 0.\n";
        throw std::invalid_argument(err_msg);
    }
    return rijsq;
};
double LJPotential::vij(Vector ri, Vector rj) {
    double pair_potential;
    double frac = pow(_sigmasq / rijsq(rij(ri, rj)), 3.0);
    pair_potential = 4.0 * _epsilon * (pow(frac, 2.0) - frac);
    return pair_potential;
};
double LJPotential::wij(Vector ri, Vector rj) {
    double frac = pow(_sigmasq * rijsq(rij(ri, rj)), 3.0);
    double wij = 24.0 * _epsilon * (2 * pow(frac, 2.0) - frac);
    return wij;
};
Vector LJPotential::fij(Vector ri, Vector rj) {
    // Allen & Tildesley Eq. 5.3
    Vector r = rij(ri, rj);
    Vector fij = multiply(r, wij(ri, rj) / rijsq(r));
    return fij;
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
AdjustedLJPotential::AdjustedLJPotential(double cutoff, double box) :
    LJPotential(),
    _rc(cutoff),
    _rcsq(pow(cutoff, 2.0)),
    _box(box),
    _corr1(_calculate_corr1()),
    _corr2(_calculate_corr2()){
        if (cutoff > 0.5 * box){
            fprintf(stderr, "%s: %s (%f) > %s (%f), %s.\n", "ljpotential", "cutoff distance", cutoff, "half box", box/2.0, "will set cutoff to half box length");
            _rc = 0.5 * box;
            _rcsq = pow(_rc, 2.0);
        }
};
AdjustedLJPotential::AdjustedLJPotential(double epsilon,
    double sigma, double cutoff, double box) :
    LJPotential(epsilon, sigma),
    _rc(cutoff),
    _rcsq(pow(cutoff, 2.0)),
    _box(box),
    _corr1(_calculate_corr1()),
    _corr2(_calculate_corr2()){
        if (cutoff > 0.5 * box){
            fprintf(stderr, "%s: %s (%f) > %s (%f), %s.\n", "ljpotential", "cutoff distance", cutoff, "half box", box/2.0, "will set cutoff to half box length");
            _rc = 0.5 * box;
            _rcsq = pow(_rc, 2.0);
        }
};
AdjustedLJPotential::~AdjustedLJPotential() {
};
std::string AdjustedLJPotential::get_str() const {
    std::string cutoff_str = "cutoff " + std::to_string(_rc);
    std::string potential_str = LJPotential::get_str() + cutoff_str + "\n";
    return potential_str;
};
Vector AdjustedLJPotential::rij(Vector ri, Vector rj){
    Vector rij = LJPotential::rij(ri, rj);
    rij.x = rij.x - _box * round( rij.x / _box );
    rij.y = rij.y - _box * round( rij.y / _box );
    rij.z = rij.z - _box * round( rij.z / _box );
    return rij;
}
double AdjustedLJPotential::rijsq(Vector rij){
    return LJPotential::rijsq(rij);
};
double AdjustedLJPotential::vij(Vector ri, Vector rj) {
    double vij = 0.0;
    double rsq = rijsq(rij(ri, rj));
    if (rsq < _rcsq){
        double frac = pow(_sigmasq / rsq, 3.0);
        double unadj = 4.0 * _epsilon * (pow(frac, 2.0) - frac);
        double adjst = _corr1 + _corr2 * (rsq / _rcsq);
        vij = unadj + adjst;
    }
    return vij;
};
double AdjustedLJPotential::wij(Vector ri, Vector rj) {
    double wij = 0.0;
    double rsq = rijsq(rij(ri, rj));
    if (rsq < _rcsq){
        // if within cutoff distance
        double frac = pow(_sigmasq / rsq, 3.0);
        double unadj = 24.0 * _epsilon * (2 * pow(frac, 2.0) - frac);
        double adjst = - 2.0 * _corr2 * (rsq / _rcsq);
        wij = unadj + adjst;
    }
    return wij;
}
Vector AdjustedLJPotential::fij(Vector ri, Vector rj){
    Vector r = rij(ri, rj);
    Vector fij = multiply(r, wij(ri, rj) / rijsq(r));
    return fij;
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
