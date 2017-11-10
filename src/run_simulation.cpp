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
#include "../include/settings_parser.h"
// particular potentials and integrators
#include "../include/ljpotential.h"
#include "../include/verlet_integrator.h"
#include "../include/rattle_integrator.h"
// state configuration handlers
#include "../include/natomic_config_handler.h"
#include "../include/diatomic_config_handler.h"
#include "../include/solvent_config_handler.h"
// observables handlers
#include "../include/diatomic_observables.h"
#include "../include/dynamic_observables.h"
namespace observables{
//    dynamic::IntKE ik = dynamic::IntKE(parameters::MP);
//    dynamic::AvgIntKE avg_ik = dynamic::AvgIntKE(ik);
//    dynamic::LNorm ln = dynamic::LNorm(parameters::MP);
//    dynamic::LProj lp = dynamic::LProj(parameters::MP);
//    dynamic::V v = dynamic::V();
//    dynamic::AvgV avg_v = dynamic::AvgV(v);
//    dynamic::NegW w = dynamic::NegW();
//    dynamic::AvgNegW avg_w = dynamic::AvgNegW(w);
//    dynamic::NegWC wc = dynamic::NegWC();
//    dynamic::AvgNegWC avg_wc = dynamic::AvgNegWC(wc);
//    container::Unit ik_unit = {.obs = ik,
//        .status = true, .timelog = true, .average = false,
//        .eval_time = container::EvalTime::integrator};
//    container::Unit avg_ik_unit = {.obs = avg_ik,
//        .status = true, .timelog = true, .average = true,
//        .eval_time = container::EvalTime::sim_loop};
//    container::Unit ln_unit = {.obs = ln,
//        .status = true, .timelog = true, .average = false,
//        .eval_time = container::EvalTime::sim_loop};
//    container::Unit lp_unit = {.obs = lp,
//        .status = true, .timelog = true, .average = false,
//        .eval_time = container::EvalTime::sim_loop};
//    container::Unit v_unit = {.obs = v,
//        .status = true, .timelog = true, .average = false,
//        .eval_time = container::EvalTime::force_loop};
//    container::Unit avg_v_unit = {.obs = avg_v,
//        .status = true, .timelog = true, .average = true,
//        .eval_time = container::EvalTime::sim_loop};
//    container::Unit w_unit = {.obs = w,
//        .status = true, .timelog = true, .average = false,
//        .eval_time = container::EvalTime::force_loop};
//    container::Unit avg_w_unit = {.obs = avg_w,
//        .status = true, .timelog = true, .average = true,
//        .eval_time = container::EvalTime::sim_loop};
//    container::Unit wc_unit = {.obs = wc,
//        .status = true, .timelog = true, .average = false,
//        .eval_time = container::EvalTime::integrator};
//    container::Unit avg_wc_unit = {.obs = avg_wc,
//        .status = true, .timelog = true, .average = true,
//        .eval_time = container::EvalTime::sim_loop};
std::vector<container::Unit> vec = {};
//    void setup(ConfigHandler& cfg, bool zero_observables){
//        if (zero_observables) {
//            for (int i = 0; i < vec.size(); ++i) {
//                vec.at(i).obs.zero();
//            }
//        }
//        double m = parameters::MP;
//        Vector laxis = lp.L(cfg.atom_state());// outputs angular momentum vector
//        laxis = divide(laxis, norm(laxis));
//        ik = dynamic::IntKE(m); // ik.update(cfg.atom_state());
//        ln = dynamic::LNorm(m); // ln.update(cfg.atom_state());
//        lp = dynamic::LProj(m, laxis);
//    }
}   // namespace observables
namespace args{
    std::string initial_state = "";
    std::string config_file;
    std::string sim_name;
    bool is_input_given = false;
}
namespace messages{
    std::string usage = "Usage: run_simulation <name> <config_file> <optional initial_state_file>";
}
bool read_arguments(int argc, char **argv){
    bool success = true;
    switch (argc) {
        case 4:
            args::initial_state = argv[3];
            args::is_input_given = true;
        case 3:
            args::config_file = argv[2];
            args::sim_name = argv[1];
            break;
        default:
            fprintf(stderr, "%s\n", messages::usage.c_str());
            success = false;
    }
    return success;
}

int main(int argc, char **argv){
    int res = 0;
    if (!read_arguments(argc, argv)) res = 1;
    else {
        SettingsParser settings = SettingsParser(args::config_file);
        // potential setup
        LJPotential pp = LJPotential(settings.epp, settings.sp);
        AdjustedLJPotential ss = AdjustedLJPotential(settings.ess, settings.ss, settings.rc_ss);
        AdjustedLJPotential ps = AdjustedLJPotential(settings.eps, 0.5 * (settings.ss + settings.sp), settings.rc_ps);
        // state setup
        simple::BaseState::set_nm(settings.np);
        simple::BaseState::set_nsolvents(4 * pow(settings.nc, 3.0));
        simple::BasePolymer::set_nb(settings.nb);
        simple::BasePolymer::set_m(settings.mp);
        simple::BasePolymer::set_d(settings.d);
        ObservableContainer container = ObservableContainer(observables::vec);
        ObservableContainer& obs = container;
        // state set up: FIXME link force updater to the configuration handler
        bool is_planar;
        if (settings.polymer_planar_conformation){
            is_planar = true;
        }
        else {
            is_planar = false;
        }
        // prepare state configuration
        SolventConfigHandler config = SolventConfigHandler(settings.rho_s, settings.nc);
        ConfigHandler& cfg = config;
        if (args::is_input_given){
            fprintf(stderr, "%s\n",
            "Don't yet know how to work with when initial state configuration is given");
            exit(EXIT_FAILURE);
        }
        else {
            config.fcc_positions();
            config.ran_velocities(settings.temperature);
        }
        fprintf(stdout, "%s\n%s",
            "State initialized by the Configuration Handler:",
            cfg.bond_state().to_string(true, true).c_str());
        //if (settings::potential != settings::potential_default){
            //if (settings::potential == settings::potential_none){
                //force_loop.set_potential(&empty);
                ///*
                //* having no LJ potential requires a change in units
                //*/
                //parameters::D = 1.0;
                //parameters::M = 1.0;
                //// positions of polymer atoms in the atomic representation
                //// have to be recalculated once bondlength is reset
                //cfg.atom_state().update(cfg.bond_state());
                //// find rms bond velocity
                //double acc = 0.0;
                //for (int i = 0; i < cfg.nmolecules(); ++i){
                    //for (int j = 0; j < cfg.polymer_nb(); ++j){
                        //acc += norm(cfg.bond_state().polymers.at(i).bonds.at(j).velocity);
                    //}
                //}
                //double rms_bondv = acc / (cfg.nmolecules() * cfg.polymer_nb());
                //double scalef = 2 * PI / rms_bondv;
                //settings::dt = settings::dt * scalef;
                //settings::runtime = settings::runtime * scalef;
            //}
            //else {
                //fprintf(stderr, "potential '%s' is not recognized\n", settings::potential.c_str());
                //exit(1);
            //}
        //}
        ForceUpdater force_loop = ForceUpdater(&pp, &ss, &ps);
        RattleIntegrator rattle = RattleIntegrator(force_loop,
            settings.tol, settings.tiny, settings.maxiter);
        VerletIntegrator verlet
            = VerletIntegrator(force_loop, NULL, config.box());
        Integrator& integrator = verlet;
        // simulation setup
        simple::Simulation sim
            = simple::Simulation(args::sim_name,
                settings.cndir,
                settings.tpdir,
                settings.dtdir,
                args::initial_state,
                cfg,
                integrator,
                obs,
                settings.dt,
                settings.icalc,
                settings.iprint,
                settings.isave,
                settings.idata,
                settings.itape);
        //if (!sim.is_input_given()){
        /**
        * The state is now set up.
        * If it was NOT set up from input file,
        *  observables may need the initial-state information. For example,
        *  the orientation of the angular momentum vector initially needs
        *  to be given to the LProj observable.
        * If it was set up from the input file, this information need not
        *  be re-written.
        */
        //    observables::setup(config, settings::zero_observables);
        //}
        sim.evolve(settings.runtime);
    }
    return res;
}
