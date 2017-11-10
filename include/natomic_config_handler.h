// 2017 Artur Avkhadiev
/*! \file natomic_config_handler.h
*/
#ifndef POLYMER_NATOMIC_CONFIG_HANDLER_H
#define POLYMER_NATOMIC_CONFIG_HANDLER_H
#include <cassert>
#include <cmath>                            /**> trig functions */
#include <functional>                       /**> std::bind      */
#include <chrono>                           /**> seed generation for RNG */
#include <random>                           /**> Mersenne Twister RNG */
#include "config_handler.h"
class NAtomicConfigHandler :
    public ConfigHandler {
protected:
    // Generate the molecule in a plane?
    struct gen_configuration_mode {
        bool is_planar;
        double phi;
    } mode;
    // RNG
    unsigned seed();                        /**> generate seed */
    typedef struct rng_t {
        unsigned seed;
        std::mt19937 generator;             /**> Standard Mersenne Twister */
        const double phi_min = 0.0;
        const double phi_max = 2.0;
        std::uniform_real_distribution<double> phi_dist;
        const double theta_min = 0.0;
        const double theta_max = 1.0;
        std::uniform_real_distribution<double> theta_dist;
        const double speed_mu = 0.0;
        const double speed_sigma = 0.25;
        std::lognormal_distribution<double> speed_dist;
    } RNG;
    void setup_rng();
    void generate_bond_position(simple::Bond& b, double phi, double theta);
    //random orientation
    void generate_bond_position(simple::Bond& b);
    // random normal to position, with lognormal
    void generate_bond_velocity(simple::Bond& b);
    void scale_velocities(double target_kinetic_energy);
    RNG rng;
public:
    NAtomicConfigHandler(double target_energy = 0.0, bool is_planar = false);
    ~NAtomicConfigHandler(){};
};
#endif
