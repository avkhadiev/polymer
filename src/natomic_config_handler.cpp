// 2017 Artur Avkhadiev
/*! \file natomic_config_handler.h
*/
#include "../include/natomic_config_handler.h"
#include "../include/default_macros.h"
unsigned NAtomicConfigHandler::seed(){
    typedef std::chrono::high_resolution_clock clock;
    clock::time_point tp = clock::now();
    clock::duration d = tp.time_since_epoch();
    return d.count();
}
void NAtomicConfigHandler::setup_rng(){
    rng.seed = seed();
    rng.generator = std::mt19937(rng.seed);
    rng.phi_dist
        = std::uniform_real_distribution<double>(rng.phi_min, rng.phi_max);
    rng.theta_dist
        = std::uniform_real_distribution<double>(rng.theta_min, rng.theta_max);
    rng.speed_dist
        = std::lognormal_distribution<double>(rng.speed_mu, rng.speed_sigma);
}
void NAtomicConfigHandler::generate_bond_position(simple::Bond& b,
    double phi, double theta)
{
    assert(phi <= rng.phi_max && phi >= rng.phi_min);
    assert(theta <= rng.theta_max &&  theta > rng.theta_min);
    phi = phi * PI;
    theta = theta * PI;
    Vector& r = b.position;
    r = vector(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
    r = divide(r, norm(r));  // normalize just in case
}
void NAtomicConfigHandler::generate_bond_position(simple::Bond& b)
{
    double phi;
    if (mode.is_planar){
        // generating molecule in a fixed plane
        phi = mode.phi;
    }
    else {
        phi = rng.phi_dist(rng.generator);
    }
    double theta = rng.theta_dist(rng.generator);
    generate_bond_position(b, phi, theta);
}
void NAtomicConfigHandler::generate_bond_velocity(simple::Bond& b)
{
    double phi;
    if (mode.is_planar){
        // generating molecule in a fixed plane
        phi = mode.phi;
    }
    else {
        phi = rng.phi_dist(rng.generator);
    }
    double theta = rng.theta_dist(rng.generator);
    double speed = rng.speed_dist(rng.generator);
    phi = phi * PI;
    theta = theta * PI;
    Vector &r = b.position;
    Vector &v = b.velocity;
    // orient velocity randomly and subtract projection onto bond,
    // then normalize and scale to random magnitude
    v = vector(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
    v = subtract(v, multiply(r, dot(v, r)));
    v = divide(v, norm(v));
    v = multiply(v, speed);
    //check that v is perpendicular to r
    //fprintf(stderr, "%s * %s = %1.5f\n",
    //    vector_to_string(v).c_str(),
    //    vector_to_string(r).c_str(),
    //    dot(v, r));
}
void NAtomicConfigHandler::scale_velocities(double target_kenergy){
    // so far, only work with a sigle polymer whose CM is at the origin
    assert(nmolecules() == 1);
    double K = 0.0;
    double m = polymer_m();
    for(const simple::Atom& atom : atom_state().polymers.at(0).atoms){
        K += m * normsq(atom.velocity) / 2;
    }
    assert(K != 0.0);
    double scale_factor = sqrt(target_kenergy / K);
    for(const simple::Atom& atom : atom_state().polymers.at(0).atoms){
        multiply(atom.velocity, scale_factor);
    }
}
NAtomicConfigHandler::NAtomicConfigHandler(double target_energy, bool is_planar):
    ConfigHandler(),
    mode({is_planar, 0.0})
{
    // so far, only work with a sigle polymer whose CM is at the origin
    assert(nmolecules() == 1);
    setup_rng();
    if (mode.is_planar) {
        mode.phi = rng.phi_dist(rng.generator);
    }
    for (simple::Bond& b : _bond_state.polymers.at(0).bonds){
        generate_bond_position(b);
        generate_bond_velocity(b);
    }
    _atom_state.update(_bond_state);
    // calculate potential energy, substract from target energy
    // set up a force updater and feed it the atomic state to get potential
    double target_kinetic_energy = target_energy; // FIXME
    if (target_kinetic_energy > 0) {
        scale_velocities(target_kinetic_energy);
    }
}
