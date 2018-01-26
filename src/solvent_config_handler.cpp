// 2017 Artur Avkhadievi
/*! \file solvent_config_handler.h
*/
#include <cmath>                            /**> cubic root */
#include <stdexcept>
#include "../include/solvent_config_handler.h"
#include "../include/default_macros.h"
SolventConfigHandler::SolventConfigHandler(double solvent_sigma,
    double density,
    int nc, int n) :
    ConfigHandler(),
    _solvent_sigma(solvent_sigma),
    _density (density),
    _nc (nc) {
        if (density <= 0.0) density = DEFAULT_SOLVENT_DENSITY;
        simple::BaseState::set_density(density);
        if ((n < 0) or ((n > 0) and (n != 4 * pow(_nc, 3.0)))) {
            std::string err_msg = "solvent config initialization: incorrect number of cells / solvents specified";
            throw std::invalid_argument(err_msg);
        }
        else {
            if (n == 0) {
                if (_nc < 0) _nc = DEFAULT_NC;
                n = 4 * pow(_nc, 3);
            }
            _box = cbrt(n / density);
            simple::BaseState::set_nsolvents(n);
        }
    }
unsigned SolventConfigHandler::seed(){
    typedef std::chrono::high_resolution_clock clock;
    clock::time_point tp = clock::now();
    clock::duration d = tp.time_since_epoch();
    return d.count();
}
void SolventConfigHandler::setup_rng(){
    rng.seed = seed();
    rng.generator = std::mt19937(rng.seed);
    rng.coordinate
        = std::uniform_real_distribution<double>(rng.min, rng.max);
}
Vector SolventConfigHandler::uniform_on_sphere(){
    double x = rng.coordinate(rng.generator);
    double y = rng.coordinate(rng.generator);
    double z = rng.coordinate(rng.generator);
    Vector random = vector(x, y, z);
    Vector unit = divide(random, norm(random));
    return unit;
}
void SolventConfigHandler::subtract_momentum(){
    // calculate CM velocity
    Vector vcm = vector(0.0, 0.0, 0.0);
    for(const simple::Solvent& solvent : atom_state().solvents){
        vcm += solvent.v();
    }
    vcm = divide(vcm, nsolvents());
    // subtract CM velocity
    for(simple::Solvent& solvent : atom_state().solvents){
        solvent.set_v(subtract(solvent.v(), vcm));
    }
}
void SolventConfigHandler::fcc_positions() {
    double box2 = _box / 2.0;
    double cell = _box / _nc;
    /* positions in unit cell are worked out as follows:
    *  there are 4 fcc points per cell: one at the vertex and three and the
    *  centers of the cubic faces adjacent to the vertex.
    * set the vertex lattice point at 0.25, 0.25, 0.25
    */
    /* coordinates along rows, lattice points along columns */
    double r_fcc[3][4] =
        {
        //    1st   2nd   3rd   4th     lattice pts per cell
            {0.25, 0.25, 0.75, 0.75},   // x-coordinates
            {0.25, 0.75, 0.75, 0.25},   // y-coordinates
            {0.25, 0.75, 0.25, 0.75}    // z-coordinates
        };
    // loop over x, y, z directions for each cubic cell
    int isolvent = 0;
    Vector r;
    for (size_t iz = 0; iz < _nc; ++iz){
        for (size_t iy = 0; iy < _nc; ++iy){
            for (size_t ix = 0; ix < _nc; ++ix){
                for (size_t pt = 0; pt < 4; ++pt){  // 4 lattice pts per cell
                    // shift position by (ix, iy, iz),
                    // multiply by cell length
                    // then shift by box2 to center everything at the origin
                    // TODO check for overlap
                    r.x = r_fcc[0][pt] + ix;
                    r.y = r_fcc[1][pt] + iy;
                    r.z = r_fcc[2][pt] + iz;
                    r = multiply(r, cell);
                    r = subtract(r, box2);
                    _atom_state.solvents.at(isolvent).set_r(r);
                    ++isolvent;
                }
            }
        }
    }
    _bond_state.update(_atom_state);
}
double SolventConfigHandler::solvent_kinetic_energy(){
    double K = 0.0;
    for(const simple::Solvent& solvent : atom_state().solvents){
        K += normsq(solvent.v());
    }
    K = solvent_m() * K / 2.0;
    return K;
}
double SolventConfigHandler::temperature(){
    int n = nsolvents();
    // K = (3N - 3) kT / 2
    double temperature = 2.0 * solvent_kinetic_energy() / (3.0 * n - 3);
    return temperature;
}
void SolventConfigHandler::rescale_velocities(double target_temperature){
    // calculate scaling factor
    double factor = sqrt(target_temperature / temperature());
    // rescale
    for(simple::Solvent& solvent : atom_state().solvents){
        solvent.set_v( multiply(solvent.v(), factor) );
    }
}
void SolventConfigHandler::ran_velocities(double temperature){
    for(simple::Solvent& solvent : atom_state().solvents){
        solvent.set_v( uniform_on_sphere() );
    }
    subtract_momentum();
    rescale_velocities(temperature);
    _bond_state.update(_atom_state);
}
std::string SolventConfigHandler::get_info_str(){
    std::string basic_info = ConfigHandler::get_info_str();
    std::string additional_info;
    std::string density
        = "density of solvents (reduced units): "
        + std::to_string(SolventConfigHandler::_density);
    std::string box
        = "box length (reduced units): "
        + std::to_string(SolventConfigHandler::_box);
    std::string temperature
        = "current temperature (reduced units): "
        + std::to_string(SolventConfigHandler::temperature());
    additional_info = density + "\n" + box + "\n" + temperature;
    std::string info_str = basic_info + "\n" + additional_info;
    return info_str;
}
