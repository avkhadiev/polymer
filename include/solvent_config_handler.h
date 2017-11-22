// 2017 Artur Avkhadiev
/*! \file solvent_config_handler.h
*       Generates solvent molecules in a cubic box of dimension L
*       Do not (!) create this config handler if you do not plan to have any
*       solvents in the simulation.
*       If you create a solvent config handler with N = 0, NC = 0, it will
*       set NC to a default value and calculate N = 4 NC^3
*/
#ifndef POLYMER_SOLVENT_CONFIG_HANDLER_H
#define POLYMER_SOLVENT_CONFIG_HANDLER_H
#include <cassert>
#include <cmath>                            /**> trig functions */
#include <functional>                       /**> std::bind      */
#include <chrono>                           /**> seed generation for RNG */
#include <random>                           /**> Mersenne Twister RNG */
#include "vector.h"
#include "config_handler.h"
class SolventConfigHandler :
    public ConfigHandler {
protected:
    double _solvent_sigma;                  /**> radius of solvent molecules  */
    double _density;                        /**> density of solvent molecules */
    int _nc;                                /**> n of lattice cells per length*/
    double _box;                            /**> determined from density */
    // RNG
    unsigned seed();                        /**> generate seed */
    typedef struct rng_t {
        unsigned seed;
        std::mt19937 generator;             /**> Standard Mersenne Twister */
        std::uniform_real_distribution<double> coordinate;
        const double min = 0.0;
        const double max = 1.0;
    } RNG;
    void setup_rng();
    Vector uniform_on_sphere();             /**> unit vector uniform on sphere*/
    void subtract_momentum();               /**> remove CM momentum */
    RNG rng;
public:
    int nc() const{return _nc;};
    double box() const{return _box;};
    void fcc_positions();                   /**> generate solvents on lattice */
    double solvent_kinetic_energy();        /**> kinetic energy of solvents */
    double temperature();
    /**> temperature is in units of k_B x K */
    void rescale_velocities(double temperature);
    void ran_velocities(double temperature);/**> generate random velocities */
    virtual std::string get_info_str();
    /** number of cells takes precedence unless n is explicitly specified (non-zero) **/
    SolventConfigHandler(double solvent_sigma, double density = 0.0, int nc = 0, int n = 0.0);
    ~SolventConfigHandler(){};
};
#endif
