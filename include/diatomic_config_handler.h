// 2017 Artur Avkhadiev
/*! \file config_handler.h
*/
#ifndef POLYMER_DIATOMIC_CONFIG_HANDLER_H
#define POLYMER_DIATOMIC_CONFIG_HANDLER_H
#include <cassert>
#include <cmath>                            /**> trig functions */
#include <functional>                       /**> std::bind      */
#include <chrono>                           /**> seed generation for RNG */
#include <random>                           /**> Mersenne Twister RNG */
#include "config_handler.h"
class DiatomicConfigHandler :
    public ConfigHandler {
protected:
    Vector _ini_bp_axis;                    /**> initial bond position axis */
    Vector _ini_bv_axis;                    /**> initial bond velocity axis */
    // RNG
    unsigned seed();                        /**> generate seed */
    typedef struct rng_t {
        unsigned seed;
        std::mt19937 generator;             /**> Standard Mersenne Twister */
        std::uniform_real_distribution<double> phi_dist;
        std::uniform_real_distribution<double> theta_dist;
    } RNG;
    void setup_rng();
    void generate_bond_position(double phi, double theta);
    void generate_bond_position();          /**> random orientation */
    void generate_bond_velocity();          /**> random normal to position */
    // stores current bond velocity & position as initial axes --- used to
    // set up observables.  e.g. projection of angular momentum onto r x v has
    // to be conserved.
    void save_initial_axes();
    RNG rng;
public:
    DiatomicConfigHandler();
    Vector ini_bp_axis() const {return _ini_bp_axis;};
    Vector ini_bv_axis() const {return _ini_bv_axis;};
    /**
    * phi \in [0, 2] in units of pi
    * theta \in [0, 1] in units of pi
    */
    DiatomicConfigHandler(double phi , double theta);
    ~DiatomicConfigHandler(){};
};
#endif
