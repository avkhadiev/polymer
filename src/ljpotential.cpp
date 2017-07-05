// 2017 Artur Avkhadiev
/*! \file ljpotential.cpp
*/
#include <cmath>
#include <stdexcept>
#include "../include/ljpotential.h"
#include "../include/vector.h"
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
Vector LJPotential::calculate_force(Vector r1, Vector r2){
    Vector r12 = subtract(r1, r2);
    double normsq_r12 = normsq(r12);
    Vector f12;
    if (normsq_r12 == 0){
        std::string err_message = "calculate_force: distance is zero";
        throw std::invalid_argument(err_message);
    }
    else {
        double rrsq = 1 / normsq_r12;
        double frac = _sigma6 * pow(rrsq, 3.0);
        double force_strength = 24 * _epsilon * rrsq * (2 * pow(frac, 2.0) - frac);
        f12 = multiply(r12, force_strength);
    }
    return f12;
};
