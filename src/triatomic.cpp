// 2017 Artur Avkhadiev
/*! \file triatomic.cpp
*/
#include "../include/simple_simulation.h"
namespace parameters{
    int NM = 1;
    int NB = 1;
    double M = 1.0;
    double D = 3.0;
}
namespace settings{
    std::string sim_name = "triatomic";
    std::string cndir = "/Users/Arthur/stratt/polymer/test/";
    std::string tpdir = "/Users/Arthur/stratt/polymer/test/";
    std::string dtdir = "/Users/Arthur/stratt/polymer/test/";
    double dt = 0.001;
    size_t icalc = 10;
    size_t iprint = 10;
    size_t isave = 0;
    size_t idata = 0;
    size_t itape = 0;
}
namespace rattle{
    double tol = pow(10, -3);
    double tiny = pow(10, -7);
    int maxiter = pow(10, 3);
}
int main(int argc, char **argv){
    // state setup
    simple::BaseState::set_nm(parameters::NM);
    simple::BasePolymer::set_nb(parameters::NB);
    simple::BasePolymer::set_m(parameters::M);
    simple::BasePolymer::set_d(parameters::D);
    ConfigHandler config = ConfigHandler();
    ConfigHandler& cfg = config;
    // integrator setup
    ForceUpdater force_loop = ForceUpdater(LJPotential());
    VerletIntegrator verlet = VerletIntegrator(force_loop);
    RattleIntegrator rattle = RattleIntegrator(force_loop,
        rattle::tol, rattle::tiny, rattle::maxiter);
    Integrator& integrator = rattle;
    // observable setup
    ObservableContainer container = ObservableContainer();
    ObservableContainer& obs = container;
    // simulation setup
    simple::Simulation sim
        = simple::Simulation(settings::sim_name,
            settings::cndir,
            settings::tpdir,
            settings::dtdir,
            cfg,
            integrator,
            obs,
            settings::dt,
            settings::icalc,
            settings::iprint,
            settings::isave,
            settings::idata,
            settings::itape);
    sim.evolve(100.0);
}
