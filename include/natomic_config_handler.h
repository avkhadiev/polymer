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
#include "../include/polymer_observables.h"
#include "config_handler.h"
class NAtomicConfigHandler :
    public ConfigHandler {
public:
    NAtomicConfigHandler(ForceUpdater* fupd = NULL,
        bool remove_p = true, bool remove_l = false, bool check_overlap = true,
        bool init_planar = false, bool randomize_phi = true,
        double new_phi = 0.0);
    NAtomicConfigHandler(simple::BondState& state, ForceUpdater* fupd = NULL);
    NAtomicConfigHandler(simple::AtomState& state, ForceUpdater* fupd = NULL);
    ~NAtomicConfigHandler(){};
    bool should_remove_linear_momentum;
    bool should_remove_angular_momentum;
    bool ovp;                                       // check overlap?
    size_t ndof_polymer() const {return temp_kin_polymer.ndof();};
    std::string print_observables();
    void initialize2D(double temperature = 5.0);
    void initialize3D(double temperature = 5.0);
    void initialize2D_given_tot_energy(double target_energy_per_atom = 3.0);
    void initialize3D_given_tot_energy(double target_energy_per_atom = 3.0);
    void scale_velocities_to_temperature(double temperature);
    void scale_velocities_to_energy(double energy);
    // helpful quantities
    double polymer_kinetic_energy();
    double polymer_potential_energy();
    double polymer_kin_temperature();
    Vector polymer_linear_momentum();
    Vector polymer_angular_momentum();
    Vector polymer_rcm();
protected:
    // Observables
    polymer::KE ke_polymer;
    polymer::PE pe_polymer;
    polymer::KineticTemperature temp_kin_polymer;
    polymer::LinMomComponent px_polymer;
    polymer::LinMomComponent py_polymer;
    polymer::LinMomComponent pz_polymer;
    polymer::AngMomComponent lx_polymer;
    polymer::AngMomComponent ly_polymer;
    polymer::AngMomComponent lz_polymer;
    void zero_observables();
    // Generate the molecule in a plane?
    struct gen_configuration_mode {
        bool is_planar;
        bool is_phi_random;
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
        const double speedsq_mean = 1.0;
        std::exponential_distribution<double> speedsq_dist;
    } RNG;
    void setup_rng();
    void generate_bond_position(simple::Bond& b, double phi, double theta);
    void generate_bond_position(simple::Bond& b);
    void generate_bond_velocity(simple::Bond& b);
    // initialization methods
    void random_positions();
    void random_velocities();
    const size_t maxiter_overlap = 10000;
    bool is_overlap(Vector ri, Vector rj, double sigmasq);
    // given a polymer atom index, checks for a potential overlap
    // of that atom with any previous (< index) atom in the molecule
    bool check_overlaps(size_t polymer_atom_index, double sigmasq);
    void remove_linear_momentum();
    // TODO
    void remove_angular_momentum();
    RNG rng;
private:
    const size_t _maxiter_init_given_tot_energy = 10000;
};
#endif
