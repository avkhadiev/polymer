// 2017 Artur Avkhadiev
/*! \file ljpotential.cpp
*/
#include <cmath>
#include "../include/ljpotential.h"
#include "../include/vector.h"
#include "../include/atom.h"
LJPotential::LJPotential() :
    _epsilon (1.0),
    _sigma   (1.0),
    _sigma6  (pow(_sigma, 6.0)),
    _sigma12 (pow(_sigma6, 2.0)) {};
LJPotential::LJPotential(double epsilon, double sigma) :
    _epsilon (epsilon),
    _sigma   (sigma),
    _sigma6  (pow(_sigma, 6.0)),
    _sigma12 (pow(_sigma6, 2.0)) {};
LJPotential::~LJPotential() {
};
double LJPotential::get_epsilon(){
    return _epsilon;
};
double LJPotential::get_sigma(){
    return _sigma;
};
void LJPotential::set_epsilon(double epsilon){
    _epsilon = epsilon;
};
void LJPotential::set_sigma(double sigma){
    _sigma = sigma;
    _sigma6  = pow(_sigma, 6.0);
    _sigma12 = pow(_sigma6, 2.0);
};
Vector LJPotential::calculate_force(Atom &on_atom, Atom &from_atom){
    Vector r1 = on_atom.position;
    Vector r2 = from_atom.position;
    Vector r12 = subtract(r1, r2);
    double rrsq = normsq(divide(1, r12));
    double frac = _sigma6 * pow(rrsq, 3.0);
    double force_strength = 24 * _epsilon * rrsq * (2 * pow(frac, 2.0) - frac);
    Vector force = multiply(r12, force_strength);
    return force;
};
