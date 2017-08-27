// 2017 Artur Avkhadiev
/*! \file run_simulation.cpp
* Usage: run_simulation
*    <name> <infile>
*    <optional zero_observables>
*    <optional nbonds bond/sigma>
*    <optional runtime/tau dt/tau rattle_tol>
*    <optional icalc iprint isave idata itape>
*    <optional cndir dtdir tpdir>
*/
#include<string>
#include<vector>
#include "../include/simple_simulation.h"
#include "../include/diatomic_observables.h"
#include "../include/dynamic_observables.h"
namespace parameters{
    int NM = 1;                         /**> single-molecule simulation */
    // default parameters
    // may be changed via command line
    int NB = 2;                         /**> is a triatomic             */
    double M = 1.0;                     /**> specifies unit of mass     */
    double D = 4.0;                     /**> (-.-)--(-.-)--(-.-)--(-.-) */
}
namespace settings{
    // mandatory arguments, suggested values
    std::string sim_name = "bond" + std::to_string(parameters::NB);
    std::string infile = "/Users/Arthur/stratt/polymer/test/readin_cn.txt";
    std::string infile_default = "default";
    // optional arguments, default values
    int zero_observables = 0;           /**> don't reset observables */
    std::string cndir = "/Users/Arthur/stratt/polymer/test/";
    std::string dtdir = "/Users/Arthur/stratt/polymer/test/";
    std::string tpdir = "/Users/Arthur/stratt/polymer/test/";
    double runtime = 4.0;
    double dt = 0.001;
    // intervals, in number of time steps
    size_t icalc = 10;                  /**> observable update      */
    size_t iprint = 10;                 /**> status printout        */
    size_t isave = 100;                 /**> configuration update   */
    size_t idata = 100;                 /**> observable writeout    */
    size_t itape = 10;                  /**> tapefile update        */
}
namespace rattle{
    double tol = pow(settings::dt, 2.0);
    double tiny = pow(10, -7);
    int maxiter = pow(10, 7);
}
namespace observables{
    dynamic::IntKE ik = dynamic::IntKE(parameters::M);
    dynamic::AvgIntKE avg_ik = dynamic::AvgIntKE(ik);
    dynamic::LNorm ln = dynamic::LNorm(parameters::M);
    dynamic::LProj lp = dynamic::LProj(parameters::M);
    dynamic::V v = dynamic::V();
    dynamic::AvgV avg_v = dynamic::AvgV(v);
    dynamic::NegW w = dynamic::NegW();
    dynamic::AvgNegW avg_w = dynamic::AvgNegW(w);
    dynamic::NegWC wc = dynamic::NegWC();
    dynamic::AvgNegWC avg_wc = dynamic::AvgNegWC(wc);
    container::Unit ik_unit = {.obs = ik,
        .status = true, .timelog = true, .average = false,
        .eval_time = container::EvalTime::integrator};
    container::Unit avg_ik_unit = {.obs = avg_ik,
        .status = true, .timelog = true, .average = true,
        .eval_time = container::EvalTime::sim_loop};
    container::Unit ln_unit = {.obs = ln,
        .status = true, .timelog = true, .average = false,
        .eval_time = container::EvalTime::sim_loop};
    container::Unit lp_unit = {.obs = lp,
        .status = true, .timelog = true, .average = false,
        .eval_time = container::EvalTime::sim_loop};
    container::Unit v_unit = {.obs = v,
        .status = true, .timelog = true, .average = false,
        .eval_time = container::EvalTime::force_loop};
    container::Unit avg_v_unit = {.obs = avg_v,
        .status = true, .timelog = true, .average = true,
        .eval_time = container::EvalTime::sim_loop};
    container::Unit w_unit = {.obs = w,
        .status = true, .timelog = true, .average = false,
        .eval_time = container::EvalTime::force_loop};
    container::Unit avg_w_unit = {.obs = avg_w,
        .status = true, .timelog = true, .average = true,
        .eval_time = container::EvalTime::sim_loop};
    container::Unit wc_unit = {.obs = wc,
        .status = true, .timelog = true, .average = false,
        .eval_time = container::EvalTime::integrator};
    container::Unit avg_wc_unit = {.obs = avg_wc,
        .status = true, .timelog = true, .average = true,
        .eval_time = container::EvalTime::sim_loop};
    std::vector<container::Unit> vec
        = {ik_unit, avg_ik_unit, v_unit, avg_v_unit,
        w_unit, avg_w_unit, wc_unit, avg_wc_unit,
        ln_unit, lp_unit};
    void setup(ConfigHandler& cfg, bool zero_observables){
        if (zero_observables) {
            for (int i = 0; i < vec.size(); ++i) {
                vec.at(i).obs.zero();
            }
        }
        double m = parameters::M;
        Vector laxis = lp.L(cfg.atom_state());// outputs angular momentum vector
        laxis = divide(laxis, norm(laxis));
        ik = dynamic::IntKE(m);
        ln = dynamic::LNorm(m);
        lp = dynamic::LProj(m, laxis);
    }
}
namespace messages{
    std::string usage = "Usage: run_simulation <name> <infile>\n<optional zero_observables>\n<optional nbonds bond/sigma>\n<optional runtime/tau dt/tau rattle_tol>\n<optional icalc iprint isave idata itape>\n<optional cndir dtdir tpdir>";
    std::string defaults =
        "Defaults:\nzero_observables "
            + std::to_string(settings::zero_observables) + "\n"
        "nbonds " + std::to_string(parameters::NB) + "\n"
        + "bond/sigma " + std::to_string(parameters::D) + "\n"
        + "runtime/tau " + std::to_string(settings::runtime) + "\n"
        + "dt " + std::to_string(settings::dt) + "\n"
        + "tol " + std::to_string(rattle::tol) + "\n"
        + "icalc " + std::to_string(settings::icalc) + "\n"
        + "iprint " + std::to_string(settings::iprint) + "\n"
        + "isave " + std::to_string(settings::isave) + "\n"
        + "idata " + std::to_string(settings::idata) + "\n"
        + "itape " + std::to_string(settings::itape) + "\n"
        + "cndir " + settings::cndir + "\n"
        + "dtdir " + settings::dtdir + "\n"
        + "tpdir " + settings::tpdir;
    std::string infile = "If infile =='" + settings::infile_default + "', the state will be initialized with the default initialization method and requires arguments <zero_observables> <nbonds> and <bond/sigma>,\nwhere zero_observables != 0 erases the memory of all average observables; otherwise, the infile should specify a path to the configuration file,\nand arguments <zero_observables>, <nbonds>, and <bond/sigma> will be ignored.";
}

bool read_arguments(int argc, char **argv){
    bool success = true;
    switch (argc) {
        case 17: settings::tpdir = argv[16];
        case 16: settings::dtdir = argv[15];
        case 15: settings::cndir = argv[14];
        case 14: settings::itape = atoi(argv[13]);
        case 13: settings::idata = atoi(argv[12]);
        case 12: settings::isave = atoi(argv[11]);
        case 11: settings::iprint = atoi(argv[10]);
        case 10: settings::icalc = atoi(argv[9]);
        case 9: rattle::tol = atof(argv[8]);
        case 8: settings::dt = atof(argv[7]);
        case 7: settings::runtime = atof(argv[6]);
        case 6: parameters::D = atof(argv[5]);
        case 5: parameters::NB = atoi(argv[4]);
        case 4: settings::zero_observables = atoi(argv[3]);
        case 3:
            settings::infile = argv[2];
            settings::sim_name = argv[1];
            break;
        default:
            fprintf(stderr, "%s\n", messages::usage.c_str());
            fprintf(stderr, "%s\n", messages::defaults.c_str());
            fprintf(stderr, "%s\n", messages::infile.c_str());
            success = false;
    }
    return success;
}

int main(int argc, char **argv){
    int res = 0;

    if (!read_arguments(argc, argv)) res = 1;
    else {
        if (settings::infile == settings::infile_default){
            settings::infile = "";
        }
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
        // integrator setup
        ForceUpdater force_loop
            = ForceUpdater(LJPotential(), &observables::v, &observables::w);
        RattleIntegrator rattle = RattleIntegrator(force_loop,
            rattle::tol, rattle::tiny, rattle::maxiter,
            &observables::ik, &observables::wc);
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
        if (!sim.is_input_given()){
        /**
        * The state is now set up.
        * If it was NOT set up from input file,
        *  observables may need the initial-state information. For example,
        *  the orientation of the angular momentum vector initially needs
        *  to be given to the LProj observable.
        * If it was set up from the input file, this information need not
        *  be re-written.
        */
            observables::setup(config, settings::zero_observables);
        }
        sim.evolve(settings::runtime);
    }
    return res;
}
