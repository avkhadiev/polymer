// 2017 Artur Avkhadiev
/*! \file diatomic_config_handler.h
*/
#include "../include/diatomic_config_handler.h"
#include "../include/default_macros.h"
unsigned DiatomicConfigHandler::seed(){
    typedef std::chrono::high_resolution_clock clock;
    clock::time_point tp = clock::now();
    clock::duration d = tp.time_since_epoch();
    return d.count();
}
void DiatomicConfigHandler::setup_rng(){
    rng.seed = seed();
    rng.generator = std::mt19937(rng.seed);
    rng.phi_dist = std::uniform_real_distribution<double>(0.0, 2.0);
    rng.theta_dist = std::uniform_real_distribution<double>(0.0, 1.0);
}
void DiatomicConfigHandler::generate_bond_position(double phi, double theta)
{
    assert(phi <= 2.0 && phi >= 0.0);
    assert(theta <= 1.0 &&  theta > 0.0);
    phi = phi * PI;
    theta = theta * PI;
    Vector &r = _bond_state.polymers.at(0).bonds.at(0).position;
    r = vector(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
    // normalize (just in case)
    r = divide(r, norm(r));
}
void DiatomicConfigHandler::generate_bond_position()
{
    double phi = rng.phi_dist(rng.generator);
    double theta = rng.theta_dist(rng.generator);
    generate_bond_position(phi, theta);
}
void DiatomicConfigHandler::generate_bond_velocity()
{
    // orient vector randomly
    double phi = rng.phi_dist(rng.generator);
    double theta = rng.theta_dist(rng.generator);
    phi = phi * PI;
    theta = theta * PI;
    Vector &v = _bond_state.polymers.at(0).bonds.at(0).velocity;
    v = vector(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
    // subtract projection onto bond
    Vector &r = _bond_state.polymers.at(0).bonds.at(0).position;
    v = subtract(v, multiply(r, dot(v, r)));
    // normalize remainder
    v = divide(v, norm(v));
}
void DiatomicConfigHandler::save_initial_axes(){
    _ini_bp_axis = _bond_state.polymers.at(0).bonds.at(0).position;
    _ini_bv_axis = _bond_state.polymers.at(0).bonds.at(0).velocity;
}
DiatomicConfigHandler::DiatomicConfigHandler():
    ConfigHandler()
{
    assert(polymer_nb() == 1);
    assert(polymer_d() == 1.0);
    setup_rng();
    generate_bond_position();
    generate_bond_velocity();
    _atom_state.update(_bond_state);
    save_initial_axes();
}
DiatomicConfigHandler::DiatomicConfigHandler(double phi, double theta):
    ConfigHandler() {
    assert(polymer_nb() == 1);
    assert(polymer_d() == 1.0);
    setup_rng();
    generate_bond_position(phi, theta);
    generate_bond_velocity();
    _atom_state.update(_bond_state);
    save_initial_axes();
}
