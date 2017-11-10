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
#include "../include/natomic_config_handler.h"
#include "../include/diatomic_observables.h"
#include "../include/dynamic_observables.h"
#include "../include/diatomic_config_handler.h"
namespace parameters{
    // default parameters
    double EPP = 1.0;                   /**> epsilon polymer-polymer      */
    double ESS = 1.0;                   /**> epsilon solvent-solvent      */
    double EPS = 1.0;                   /**> epsilon polymer-solvent      */
    double SP  = 1.0;                   /**> sigma polymer                */
    double SS  = 1.0;                   /**> sigma solvent                */
    double RC  = 5.0;                   /**> potential cutoff for solvents*/
    double MP = 1.0;                    /**> polymer mass                 */
    double MS = 1.0;                    /**> solvent mass                 */
    int NM = 1;                         /**> number of polymers           */
    int NB = 2;                         /**> number of bonds              */
    int NC = 4;                         /**> nsolvents = 4 * nc^3         */
    double D = 0.86;                    /**> reduced in units of SP!      */
    double RHOS = 0.7;                  /**> density in reduced units     */
    double T = 1.0;                     /**> temperature in reduced units */
}
namespace settings{
    // mandatory arguments, suggested values
    std::string sim_name = "sim";
    std::string infile = "/Users/Arthur/stratt/polymer/test/readin_cn.txt";
    std::string infile_default = "default";
    std::string infile_planar = "planar";
    // optional arguments, default values
    int zero_observables = 0;           /**> don't reset observables */
    std::string polymer_potential = "polymer potential";
    std::string solvent_potential = "solvent potential";
    std::string inter_potential = "polymer-solvent potential";
    std::string potential_lj = "lj";
    std::string potential_lj_sf = "lj-sf";
    std::string potential_none = "none";
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
    dynamic::IntKE ik = dynamic::IntKE(parameters::MP);
    dynamic::AvgIntKE avg_ik = dynamic::AvgIntKE(ik);
    dynamic::LNorm ln = dynamic::LNorm(parameters::MP);
    dynamic::LProj lp = dynamic::LProj(parameters::MP);
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
        double m = parameters::MP;
        Vector laxis = lp.L(cfg.atom_state());// outputs angular momentum vector
        laxis = divide(laxis, norm(laxis));
        ik = dynamic::IntKE(m); // ik.update(cfg.atom_state());
        ln = dynamic::LNorm(m); // ln.update(cfg.atom_state());
        lp = dynamic::LProj(m, laxis);
    }
}
namespace messages{
    std::string usage = "Usage: run_simulation <name> <infile>\n<optional zero_observables target_kenergy>\n<optional potential>\n<optional nbonds bond/sigma>\n<optional runtime/tau dt/tau rattle_tol>\n<optional icalc iprint isave idata itape>\n<optional cndir dtdir tpdir>";
    std::string defaults =
        "Defaults:\nzero_observables "
            + std::to_string(settings::zero_observables) + "\n"
        + "target_kenergy " + std::to_string(settings::target_kenergy) + "\n"
        + "potential " + settings::potential_default + "\n"
        + "nbonds " + std::to_string(parameters::NB) + "\n"
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
    std::string potential = "If potential =='" + settings::potential_default + "', the simulation will use a default LJ potential with epsilon, sigma = 1.0; if potential == '" + settings::potential_none + "', there will be no potential at all. Units of mass and length will be determined by atomic mass and bond length, respectively. Time unit will be equal to 2 pi * d / rms(bond velocity)";
    std::string infile = "If infile =='" + settings::infile_default + "' (or '" + settings::infile_planar + "'), the state will be initialized with the default initialization method (or with molecules in planar configurations) and requires <nbonds> and <bond/sigma>; otherwise, the infile should specify a path to the configuration file, and arguments <nbonds>, and <bond/sigma> will be ignored. argument zero_observables != 0 erases the memory of all average observables and target_kenergy = 0 will not scale velocities.";
}

bool read_arguments(int argc, char **argv){
    bool success = true;
    switch (argc) {
        case 19: settings::tpdir = argv[18];
        case 18: settings::dtdir = argv[17];
        case 17: settings::cndir = argv[16];
        case 16: settings::itape = atoi(argv[15]);
        case 15: settings::idata = atoi(argv[14]);
        case 14: settings::isave = atoi(argv[13]);
        case 13: settings::iprint = atoi(argv[12]);
        case 12: settings::icalc = atoi(argv[11]);
        case 11: rattle::tol = atof(argv[10]);
        case 10: settings::dt = atof(argv[9]);
        case 9: settings::runtime = atof(argv[8]);
        case 8: parameters::D = atof(argv[7]);
        case 7: parameters::NB = atoi(argv[6]);
        case 6: settings::potential = argv[5];
        case 5: settings::target_kenergy = atof(argv[4]);
        case 4: settings::zero_observables = atoi(argv[3]);
        case 3:
            settings::infile = argv[2];
            settings::sim_name = argv[1];
            break;
        default:
            fprintf(stderr, "%s\n", messages::usage.c_str());
            fprintf(stderr, "%s\n", messages::defaults.c_str());
            fprintf(stderr, "%s\n", messages::potential.c_str());
            fprintf(stderr, "%s\n", messages::infile.c_str());
            success = false;
    }
    return success;
}

int main(int argc, char **argv){
    int res = 0;
    if (!read_arguments(argc, argv)) res = 1;
    else {
        // state setup
        simple::BaseState::set_nm(parameters::NM);
        simple::BasePolymer::set_nb(parameters::NB);
        simple::BasePolymer::set_m(parameters::M);
        simple::BasePolymer::set_d(parameters::D);
        // observable setup
        ObservableContainer container = ObservableContainer(observables::vec);
        ObservableContainer& obs = container;
        // potential and force updater setup
        LJPotential lj = LJPotential(1.0, 1.0);
        Potential empty = Potential();
        Potential *potential = &lj;
        ForceUpdater force_loop
            = ForceUpdater(potential, &observables::v, &observables::w);
        // state set up: FIXME link force updater to the configuration handler
        bool is_planar = false;
        if (settings::infile == settings::infile_default){
            settings::infile = "";
        }
        else if (settings::infile == settings::infile_planar){
            settings::infile = "";
            is_planar = true;
        }
        NAtomicConfigHandler config =
            NAtomicConfigHandler(settings::target_kenergy, is_planar);
        ConfigHandler& cfg = config;
        fprintf(stdout, "%s\n%s",
            "State initialized by the Configuration Handler:",
            cfg.bond_state().to_string(true, true).c_str());
        if (settings::potential != settings::potential_default){
            if (settings::potential == settings::potential_none){
                force_loop.set_potential(&empty);
                /*
                * having no LJ potential requires a change in units
                */
                parameters::D = 1.0;
                parameters::M = 1.0;
                // positions of polymer atoms in the atomic representation
                // have to be recalculated once bondlength is reset
                cfg.atom_state().update(cfg.bond_state());
                // find rms bond velocity
                double acc = 0.0;
                for (int i = 0; i < cfg.nmolecules(); ++i){
                    for (int j = 0; j < cfg.polymer_nb(); ++j){
                        acc += norm(cfg.bond_state().polymers.at(i).bonds.at(j).velocity);
                    }
                }
                double rms_bondv = acc / (cfg.nmolecules() * cfg.polymer_nb());
                double scalef = 2 * PI / rms_bondv;
                settings::dt = settings::dt * scalef;
                settings::runtime = settings::runtime * scalef;
            }
            else {
                fprintf(stderr, "potential '%s' is not recognized\n", settings::potential.c_str());
                exit(1);
            }
        }
        fprintf(stderr, "will use potential:\n%s", force_loop.get_potential()->get_str().c_str());
        // integrator setup
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
