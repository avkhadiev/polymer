// 2017 Artur Avkhadiev
/*! \file natomic_config_handler.h
*/
#include <math.h>
#include "../include/natomic_config_handler.h"
#include "../include/default_macros.h"
#include <exception>
NAtomicConfigHandler::NAtomicConfigHandler(ForceUpdater* fupd,
    bool remove_p, bool remove_l, bool overlap,
    bool init_planar, bool randomize_phi,
    double new_phi
):
    ConfigHandler(fupd),
    should_remove_linear_momentum(remove_p),
    should_remove_angular_momentum(remove_l),
    ovp(overlap),
    ke_polymer( cf_observables::calc_mean,
                cf_observables::calc_err,
                cf_observables::print_val),
    pe_polymer( cf_observables::calc_mean,
                cf_observables::calc_err,
                cf_observables::print_val),
    temp_kin_polymer(
                should_remove_linear_momentum,
                should_remove_angular_momentum,
                cf_observables::calc_mean,
                cf_observables::calc_err,
                cf_observables::print_val),
    px_polymer( observable::x_vec,
                cf_observables::calc_mean,
                cf_observables::calc_err,
                cf_observables::print_val),
    py_polymer( observable::y_vec,
                cf_observables::calc_mean,
                cf_observables::calc_err,
                cf_observables::print_val),
    pz_polymer( observable::z_vec,
                cf_observables::calc_mean,
                cf_observables::calc_err,
                cf_observables::print_val),
    lx_polymer( observable::x_vec,
                cf_observables::calc_mean,
                cf_observables::calc_err,
                cf_observables::print_val),
    ly_polymer( observable::y_vec,
                cf_observables::calc_mean,
                cf_observables::calc_err,
                cf_observables::print_val),
    lz_polymer( observable::z_vec,
                cf_observables::calc_mean,
                cf_observables::calc_err,
                cf_observables::print_val),
    mode({.is_planar = init_planar,
        .is_phi_random = randomize_phi,
        .phi = new_phi})
{
    // so far, only work with a sigle polymer whose CM is at the origin
    assert(nmolecules() == 1);
    setup_rng();
    zero_observables();
}
NAtomicConfigHandler::NAtomicConfigHandler(simple::BondState& state, ForceUpdater* fupd) : NAtomicConfigHandler() {
    _bond_state = state;
    _atom_state = simple::AtomState(_bond_state);
    // so far, only work with a sigle polymer whose CM is at the origin
    assert(nmolecules() == 1);
    setup_rng();
    zero_observables();
}

NAtomicConfigHandler::NAtomicConfigHandler(simple::AtomState& state, ForceUpdater* fupd) : NAtomicConfigHandler() {
    _atom_state = state;
    _bond_state = simple::BondState(_atom_state);
    // so far, only work with a sigle polymer whose CM is at the origin
    assert(nmolecules() == 1);
    setup_rng();
    zero_observables();
}

void NAtomicConfigHandler::zero_observables(){
    pe_polymer.value = 0.0;
    ke_polymer.value = 0.0;
    temp_kin_polymer.value = 0.0;
    px_polymer.value = 0.0;
    py_polymer.value = 0.0;
    pz_polymer.value = 0.0;
    lx_polymer.value = 0.0;
    ly_polymer.value = 0.0;
    lz_polymer.value = 0.0;
}
std::string NAtomicConfigHandler::print_observables(){
    std::string str = "";
    str += "NAtomicConfigHandler:\n";
    str += "PE " + std::to_string(polymer_potential_energy()) + "\n";
    str += "KE " + std::to_string(polymer_kinetic_energy()) + "\n";
    str += "Tkin " + std::to_string(polymer_kin_temperature()) + "\n";
    str += "Rcm " + vector_to_string(polymer_rcm());
    return str;
}
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
    rng.speedsq_dist
        = std::exponential_distribution<double>(rng.speedsq_mean);
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
    double speedsq = rng.speedsq_dist(rng.generator);
    phi = phi * PI;
    theta = theta * PI;
    Vector &r = b.position;
    Vector &v = b.velocity;
    // orient velocity randomly and subtract projection onto bond,
    // then normalize and scale to random magnitude according to distribution
    v = vector(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
    v = subtract(v, multiply(r, dot(v, r)));
    v = divide(v, norm(v));
    v = multiply(v, pow(speedsq, 0.5));
    //check that v is perpendicular to r
    //fprintf(stderr, "%s * %s = %1.5f\n",
    //    vector_to_string(v).c_str(),
    //    vector_to_string(r).c_str(),
    //    dot(v, r));
}
double NAtomicConfigHandler::polymer_kinetic_energy(){
    ke_polymer.value = 0.0;
    ke_polymer.update(atom_state());
    return ke_polymer.value;
}
double NAtomicConfigHandler::polymer_potential_energy(){
    pe_polymer.value = 0.0;
    if (_fupd != NULL){
        pe_polymer.value
            = _fupd->polymer_potential_energy(_atom_state.polymers.at(0));
    }
    else {
        fprintf(stderr, "%s\n", "NAtomicConfigHandler: potential energy requested, but ForceUpdater is not set. Returning V = 0");
    }
    return pe_polymer.value;
}
double NAtomicConfigHandler::polymer_kin_temperature(){
    temp_kin_polymer.value = 0;
    temp_kin_polymer.update(atom_state());
    return temp_kin_polymer.value;
}
Vector NAtomicConfigHandler::polymer_linear_momentum(){
    px_polymer.value = 0.0; px_polymer.update(_atom_state.polymers.at(0));
    py_polymer.value = 0.0; py_polymer.update(_atom_state.polymers.at(0));
    pz_polymer.value = 0.0; pz_polymer.update(_atom_state.polymers.at(0));
    Vector p = vector(px_polymer.value, py_polymer.value, pz_polymer.value);
    return p;
}
Vector NAtomicConfigHandler::polymer_angular_momentum(){
    lx_polymer.value = 0.0; lx_polymer.update(_atom_state.polymers.at(0));
    ly_polymer.value = 0.0; ly_polymer.update(_atom_state.polymers.at(0));
    lz_polymer.value = 0.0; lz_polymer.update(_atom_state.polymers.at(0));
    Vector l = vector(lx_polymer.value, ly_polymer.value, lz_polymer.value);
    return l;
}
Vector NAtomicConfigHandler::polymer_rcm(){
    Vector rcm = vector(0.0, 0.0, 0.0);
    Vector ri;
    size_t natoms = simple::BasePolymer::nb() + 1;
    for(size_t ia = 0; ia < natoms; ++ia){
        ri = atom_state().polymers.at(0).atoms.at(ia).position;
        add(rcm, ri);
    }
    divide(rcm, natoms);
    return rcm;
}
void NAtomicConfigHandler::remove_linear_momentum(){
    Vector pcm = polymer_linear_momentum();
    Vector vcm = divide(pcm, simple::BasePolymer::m());
    for(const simple::Atom& atom : atom_state().polymers.at(0).atoms){
        subtract(atom.velocity, vcm);
    }
}
void NAtomicConfigHandler::scale_velocities_to_temperature(double temperature){
    temperature = abs(temperature);
    temp_kin_polymer.value = 0; temp_kin_polymer.update(atom_state());
    double scale_factor = sqrt(temperature/polymer_kin_temperature());
    for(simple::Atom& atom : atom_state().polymers.at(0).atoms){
        atom.velocity = multiply(atom.velocity, scale_factor);
    }
    if (DEBUG) {
        fprintf(stdout, "%s\n", print_observables().c_str());
        fprintf(stderr, "%s %f\n", "target_temperature", temperature);
        fprintf(stderr, "%s %f\n", "scale_factor", scale_factor);
    }
}
void NAtomicConfigHandler::scale_velocities_to_energy(double energy){
    energy = abs(energy);
    double PE = polymer_potential_energy(); assert(PE != 0.0);
    double KE = polymer_kinetic_energy(); assert(KE != 0.0);
    double target_KE = energy - PE;
    size_t tries = 1;
    if (target_KE < 0.0){
        size_t maxiter = _maxiter_init_given_tot_energy;
        size_t iter = 0;
        while ((target_KE < 0.0) && (iter < maxiter)){
            random_positions();
            PE = polymer_potential_energy();
            KE = polymer_kinetic_energy();
            target_KE = energy - PE;
            ++iter;
        }
        if (iter >= maxiter){
            fprintf(stderr, "%s\n", "NAtomicHandler: could not initialize polymer at required energy (microcanonical)");
            EXIT_FAILURE;
        }
        tries += iter;
    }
    if (DEBUG) {
        fprintf(stdout, "%s: %s %zu %s\n",
            "NAtomicHandler scale_velocities_to_energy", "It took", tries, "tries to configure the polymer at required energy");
    }
    double scale_factor = sqrt(target_KE/KE);
    for(simple::Atom& atom : atom_state().polymers.at(0).atoms){
        atom.velocity = multiply(atom.velocity, scale_factor);
    }
    if (DEBUG) {
        fprintf(stdout, "%s\n", print_observables().c_str());
        fprintf(stderr, "%s %f\n", "target energy", energy);
        fprintf(stderr, "%s %f\n", "target_KE", target_KE);
        fprintf(stderr, "%s %f\n", "scale_factor", scale_factor);
    }
}
void NAtomicConfigHandler::random_positions(){
    if (_fupd != NULL){
        simple::Bond b;
        simple::BondPolymer& bpolymer = bond_state().polymers.at(0);
        double sigmasq
            = pow(simple::BasePolymer::d() * _fupd->get_polymer_potential()->get_sigma(), 2.0);
        size_t maxiter = maxiter_overlap;
        size_t iter;
        size_t tries = 0;
        for (size_t ib = 0; ib < simple::BasePolymer::nb(); ++ib){
            iter = 0;
            // generate next director
            ++tries;
            generate_bond_position(bpolymer.bonds.at(ib));
            // update the atomic representation
            _atom_state.update(_bond_state);
            // if overlap check required, then repeat trying to generate
            // the next director until the latest atom does not overlap
            // with any of the previous atoms in a polymer, or until
            // the maximum number of iterations is reached
            if (ovp){
                while (iter < maxiter && (check_overlaps(ib + 1, sigmasq)))
                {
                    ++iter;
                    ++tries;
                    generate_bond_position(bpolymer.bonds.at(ib));
                    _atom_state.update(_bond_state);
                }
            }
            if (iter >= maxiter)
            {
                std::string err = "NAtomicConfigHandler: while initializing random positions, number of iterations exceeded the allotted max because of overlap checks\n";
                throw std::out_of_range(err);
            }
        }
        if (DEBUG) {
            fprintf(stdout, "%s: %s %zu %s\n",
                "NAtomicHandler random_positions", "It took", tries, "tries to orient the polymer without overlaps");
        }
        _atom_state.update(_bond_state);
    }
    else {
        fprintf(stderr, "%s\n", "NAtomicConfigHandler: cannot initialize state without the force updater specified");
        return;
    }
}
void NAtomicConfigHandler::random_velocities(){
    // ensure COM is at origin
    double epsilon = pow(10, -6.0);
    atom_state().update(bond_state());
    if (normsq(polymer_rcm()) > epsilon){
        fprintf(stderr, "%s\n", "NAtomicConfigHandler: COM of polymer is not at origin.");
        EXIT_FAILURE;
    }
    else{
        for (simple::Bond& b : _bond_state.polymers.at(0).bonds){
            generate_bond_velocity(b);
        }
    }
    _atom_state.update(_bond_state);
}
bool NAtomicConfigHandler::is_overlap(Vector ri, Vector rj, double sigmasq){
    bool overlap = false;
    if (_fupd != NULL){
        double rijsq = normsq(_fupd->get_polymer_potential()->rij(ri, rj));
        if (rijsq < sigmasq){
            overlap = true;
        }
    }
    else {
        fprintf(stderr, "%s\n", "NAtomicConfigHandler: overlap check requested, but ForceUpdater is not set. Returning overlap = false");
    }
    return overlap;
}
bool NAtomicConfigHandler::check_overlaps(size_t polymer_atom_index,
    double sigmasq){
    bool overlap = false;
    if (polymer_atom_index > 0){
        Vector ri
            = atom_state().polymers.at(0).atoms.at(polymer_atom_index).position;
        Vector rj;
        for(size_t ia = 0; ia < polymer_atom_index - 1; ++ia){
            rj = atom_state().polymers.at(0).atoms.at(ia).position;
            if (is_overlap(ri, rj, sigmasq)) {
                overlap = true;
                break;
            }
        }
    }
    return overlap;
}
void NAtomicConfigHandler::initialize2D(double temperature){
    mode.is_planar = true;
    if (mode.is_phi_random){
        mode.phi = rng.phi_dist(rng.generator);
    }
    random_positions();
    random_velocities();
    scale_velocities_to_temperature(temperature);
}
void NAtomicConfigHandler::initialize3D(double temperature){
    mode.is_planar = false;
    random_positions();
    random_velocities();
    scale_velocities_to_temperature(temperature);
}
void NAtomicConfigHandler::initialize2D_given_tot_energy(double target_energy){
    mode.is_planar = true;
    if (mode.is_phi_random){
        mode.phi = rng.phi_dist(rng.generator);
    }
    random_positions();
    random_velocities();
    scale_velocities_to_energy(target_energy);
}
void NAtomicConfigHandler::initialize3D_given_tot_energy(double target_energy){
    mode.is_planar = false;
    random_positions();
    random_velocities();
    scale_velocities_to_energy(target_energy);
}
