// 2017 Artur Avkhadiev
/*! \file run_simulation.cpp
*/
#include<string>
#include<vector>
#include "../include/simple_simulation.h"
#include "../include/diatomic_observables.h"
#include "../include/dynamic_observables.h"
namespace parameters{
    int NM = 1;
    int NB = 1;
    double M = 1.0;
    double D = 4.0;
}
namespace settings{
    std::string sim_name = "test";
    std::string cndir = "/Users/Arthur/stratt/polymer/test/";
    std::string tpdir = "/Users/Arthur/stratt/polymer/test/";
    std::string dtdir = "/Users/Arthur/stratt/polymer/test/";
    std::string infile = "/Users/Arthur/stratt/polymer/test/readin_cn.txt";
    double runtime;
    double dt = 0.001;
    size_t icalc = 10;
    size_t iprint = 10;
    size_t isave = 100;
    size_t idata = 100;
    size_t itape = 10;
}
namespace rattle{
    double tol = pow(10, -7);
    double tiny = pow(10, -7);
    int maxiter = pow(10, 7);
}
namespace observables{
    // declare observables
    std::vector<container::Unit> vec
        = {};
}
namespace messages{
    std::string usage = "Usage: run_simulation <name> <runtime/tau> <tol>.";
    std::string defaults = "default bond length = "
        + std::to_string(parameters::D)
        + ", default nbonds = "
        + std::to_string(parameters::NB);
}
int main(int argc, char **argv){
    if (argc < 3){
        fprintf(stderr, "%s\n", messages::usage.c_str());
        fprintf(stderr, "%s\n", messages::defaults.c_str());
    }
    else {
        settings::sim_name = argv[1];
        settings::runtime = atof(argv[2]);      /**>  run time */
        rattle::tol = atof(argv[3]);            /**> tolerance */
        // include parameters in name
        //settings::sim_name
        //    += "_T_" + std::to_string(settings::runtime)
        //    + "_D_" + std::to_string(parameters::D)
        //    + "_NB_" + std::to_string(parameters::NB);
        //settings::sim_name = parse_string(settings::sim_name);
        // state setup
        simple::BaseState::set_nm(parameters::NM);
        simple::BasePolymer::set_nb(parameters::NB);
        simple::BasePolymer::set_m(parameters::M);
        simple::BasePolymer::set_d(parameters::D);
        ConfigHandler config = ConfigHandler();
        ConfigHandler& cfg = config;
        // observable setup
        ObservableContainer container = ObservableContainer(observables::vec);
        ObservableContainer& obs = container;
        // integrator setup
        ForceUpdater force_loop = ForceUpdater(LJPotential());
        VerletIntegrator verlet = VerletIntegrator(force_loop);
        RattleIntegrator rattle = RattleIntegrator(force_loop,
            rattle::tol, rattle::tiny, rattle::maxiter);
        Integrator& integrator = rattle;
        // simulation setup
        simple::Simulation sim
            = simple::Simulation(settings::sim_name,
                settings::cndir,
                settings::tpdir,
                settings::dtdir,
                settings::infile,
                cfg,
                integrator,
                obs,
                settings::dt,
                settings::icalc,
                settings::iprint,
                settings::isave,
                settings::idata,
                settings::itape);
        sim.evolve(settings::runtime);
    }
}
