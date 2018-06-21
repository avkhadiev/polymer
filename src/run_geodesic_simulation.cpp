// 2018 Artur Avkhadiev
/*! \file run_geodesic_simulation.cpp
* Usage: run_geodesic_simulation <name> <config_file> <optional initial_state_file>
*/
#include<string>
#include<vector>
#include "../include/default_macros.h"
#include "../include/settings_parser.h"
// particular potentials and integrators
#include "../include/potential.h"
#include "../include/ljpotential.h"
// state configuration handlers
#include "../include/natomic_config_handler.h"
#include "../include/diatomic_config_handler.h"
#include "../include/solvent_config_handler.h"
// observables handlers
#include "../include/general_observables.h"
#include "../include/polymer_observables.h"
//#include "../include/solvent_observables.h"
// simulation class
#include "../include/simple_simulation.h"
namespace args{
    std::string config_file;
    std::string sim_name;
}
namespace messages{
    std::string usage = "Usage: run_md_simulation <name> <config_file>";
}
bool read_arguments(int argc, char **argv){
    bool success = true;
    switch (argc) {
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
        // state setup
        simple::BaseState::set_nm(settings.np);
        simple::BaseState::set_nsolvents(4 * pow(settings.nc, 3.0));
        simple::BasePolymer::set_nb(settings.nb);
        // normalize tolerance for entire chain
        settings.epsilon = settings.epsilon / sqrt(settings.nb + 1);
        fprintf(stderr, "%s %f\n", "normalized epsilon (per atom)", settings.epsilon);
        // change in units if potential is zero
        bool no_potential = false;
        if ((settings.epp == 0.0)
            && (settings.ess == 0.0)
            && (settings.eps == 0.0)){
            /*
            * having no LJ potential requires a change in units
            */
            settings.np = 1.0;
            settings.mp = 1.0;
            settings.d = 1.0;
            no_potential = true;
        }
        simple::BasePolymer::set_m(settings.mp);
        simple::BasePolymer::set_d(settings.d);
        // setup observables
        bool calc_mean = true;
        bool calc_err = true;
        bool print_val = true;
        polymer::PE pe_polymer = polymer::PE(!calc_mean, !calc_err, !print_val);
        polymer::Virial w_polymer = polymer::Virial(!calc_mean, !calc_err, !print_val);
        int a1 = 0;
        int a2 = settings.nb / 3;
        int a3 = settings.nb / 2;
        polymer::ConstraintVirial wc_polymer
            = polymer::ConstraintVirial(!calc_mean, !calc_err, !print_val);
        polymer::PolAtomPosComponent r1x
            = polymer::PolAtomPosComponent(vector(1.0, 0.0, 0.0), 0, a1);
        polymer::PolAtomPosComponent r2x
            = polymer::PolAtomPosComponent(vector(1.0, 0.0, 0.0), 0, a2);
        polymer::PolAtomPosComponent r3x
            = polymer::PolAtomPosComponent(vector(1.0, 0.0, 0.0), 0, a3);
        polymer::PolAtomPosComponent r1y
            = polymer::PolAtomPosComponent(vector(0.0, 1.0, 0.0), 0, a1);
        polymer::PolAtomPosComponent r2y
            = polymer::PolAtomPosComponent(vector(0.0, 1.0, 0.0), 0, a2);
        polymer::PolAtomPosComponent r3y
            = polymer::PolAtomPosComponent(vector(0.0, 1.0, 0.0), 0, a3);
        polymer::PolAtomPosComponent r1z
            = polymer::PolAtomPosComponent(vector(0.0, 0.0, 1.0), 0, a1);
        polymer::PolAtomPosComponent r2z
            = polymer::PolAtomPosComponent(vector(0.0, 0.0, 1.0), 0, a2);
        polymer::PolAtomPosComponent r3z
            = polymer::PolAtomPosComponent(vector(0.0, 0.0, 1.0), 0, a3);
        std::vector<Observable*> observables_vec
            = { &pe_polymer,
                &w_polymer, &wc_polymer,
                &r1x, &r2x, &r3x, &r1y, &r2y, &r3y, &r1z, &r2z, &r3z };
        ObservableContainer container = ObservableContainer(observables_vec);
        ObservableContainer& obs = container;
        ObservableContainer container2 = ObservableContainer(observables_vec);
        //ObservableContainer& obs2 = container2;
        polymer::RCMComponent rcmx
            = polymer::RCMComponent(vector(1.0, 0.0, 0.0), 0,
                !calc_mean, !calc_err, !print_val);
        polymer::RCMComponent rcmy
            = polymer::RCMComponent(vector(0.0, 1.0, 0.0), 0,
                !calc_mean, !calc_err, !print_val);
        polymer::RCMComponent rcmz
            = polymer::RCMComponent(vector(0.0, 0.0, 1.0), 0,
                !calc_mean, !calc_err, !print_val);
        polymer::BondLength bl1
            = polymer::BondLength(0, 0,
                !calc_mean, !calc_err, !print_val);
        polymer::BondLength bl2
            = polymer::BondLength(0, 1,
                !calc_mean, !calc_err, !print_val);
        std::vector<Observable*> observables_vec3
            = { &pe_polymer,
                &w_polymer, &wc_polymer,
                &rcmx, &rcmy, &rcmz,
                &bl1, &bl2,
                &r1x, &r2x, &r3x, &r1y, &r2y, &r3y, &r1z, &r2z, &r3z};
        ObservableContainer container3 = ObservableContainer(observables_vec3);
        ObservableContainer& obs3 = container3;
        std::vector<Observable*> observables_vec4
            = { &pe_polymer,
                &w_polymer, &wc_polymer,
                &rcmx, &rcmy, &rcmz,
                &bl1, &bl2,
                &r1x, &r2x, &r3x, &r1y, &r2y, &r3y, &r1z, &r2z, &r3z};
        ObservableContainer container4 = ObservableContainer(observables_vec3);
        ObservableContainer& obs4 = container4;
        // setup potentials
        LJPotential pp = LJPotential(settings.epp, settings.sp,
            &pe_polymer, &w_polymer);
        SolventConfigHandler solvent_config = SolventConfigHandler(settings.ss, settings.rho_s, settings.nc);
        AdjustedLJPotential ss = AdjustedLJPotential(settings.ess, settings.ss, settings.rc_ss, solvent_config.box());
        AdjustedLJPotential ps = AdjustedLJPotential(settings.eps, 0.5 * (settings.ss + settings.sp), settings.rc_ps, solvent_config.box());
        // create endpoint records in separate files
        fprintf(stderr, "%s\n", "Creating endpoint records from MD file...");
        geodesic::Manager manager = geodesic::Manager(&pp, &ss, &ps, &container2);
        manager.read_states(settings.cndir, args::sim_name);
        manager.write_geodesic_inputs(settings.geodir_input, args::sim_name);
        std::string ini_path
            = settings.geodir_input + manager.ini_fname(args::sim_name);
        std::string fin_path
            = settings.geodir_output + manager.fin_fname(args::sim_name);
        if (settings.reverse_path){
            std::string temp = ini_path;
            ini_path = fin_path;
            fin_path = temp;
        }
        fprintf(stderr, "%s:\n%s: %s\n%s: %s\n",
            "Endpoint files are created",
            "initial record",
            ini_path.c_str(),
            "final record",
            fin_path.c_str());
        // output MD path file, if necessary
        if (settings.save_md_path) {
            bool overwrite = true;
            std::string mdfile = args::sim_name + "md_path.txt";
            std::string md_file
                = settings.geodir_output + mdfile;
            fprintf(stderr, "%s %s\n",
                "Output MD path file to", md_file.c_str());
            geodesic::Path md_path = manager.MD_path();
            if (settings.reverse_path) {
                md_path.reverse();
            }
            md_path.write(md_file, overwrite);
            manager.output_observables(md_path, settings.geodir_output, args::sim_name + "md");
        }
        // initialize path computers
        geodesic::SLERP slerp_cmp
            = geodesic::SLERP(&pp, &ss, &ps, settings.epsilon);
        obs.add_observable(&slerp_cmp.psi1);
        obs.add_observable(&slerp_cmp.psi2);
        obs.add_observable(&slerp_cmp.delta_psi1);
        obs.add_observable(&slerp_cmp.delta_psi2);
        //geodesic::ShortStep shortstep_cmp
        //    = geodesic::ShortStep(&pp, &ss, &ps, 0.0001);
        //obs2.add_observable(&shortstep_cmp.psi1);
        //obs2.add_observable(&shortstep_cmp.psi2);
        //obs2.add_observable(&shortstep_cmp.delta_psi1);
        //obs2.add_observable(&shortstep_cmp.delta_psi2);
        //obs2.add_observable(&shortstep_cmp.theta1);
        //obs2.add_observable(&shortstep_cmp.theta2);
        //obs2.add_observable(&shortstep_cmp.delta_theta1);
        //obs2.add_observable(&shortstep_cmp.delta_theta2);
        geodesic::SHOVE shove_cmp
            = geodesic::SHOVE(&pp, &ss, &ps, settings.tol, settings.epsilon);
        obs3.add_observable(&shove_cmp.psi1);
        obs3.add_observable(&shove_cmp.psi2);
        geodesic::PLERP plerp_cmp
            = geodesic::PLERP(&pp, &ss, &ps, settings.epsilon);
        // initialize geodesic simulations
        simple::SLERPGeodesicSimulation slerp_sim
            = simple::SLERPGeodesicSimulation(args::sim_name + "SLERP",
                ini_path,
                fin_path,
                settings.geodir_output,     // dir for path file with no forays
                settings.geodir_output,     // dir fpr path file with forays
                settings.geodir_data,       // dir for path observables
                obs,
                &slerp_cmp,
                settings.el,
                settings.should_write_data,
                settings.dtau,    // step in affine parameter
                settings.icalc,   // calculate every icalc steps
                settings.iblock,  // average every iblock calcsteps
                settings.iprint,  // print status every iprint steps
                settings.isave,   // save path + obs every isave steps
                settings.itape,   // save path w/ forays every itape steps
                settings.max_propag_steps,    // max propagation steps
                settings.max_escape_steps);   // max escape steps
        //simple::ShortStepGeodesicSimulation shortstep_sim
        //    = simple::ShortStepGeodesicSimulation(args::sim_name + "ShortStep",
        //        ini_path,
        //        fin_path,
        //        settings.geodir_output,     // dir for path file with no forays
        //        settings.geodir_output,     // dir fpr path file with forays
        //        settings.geodir_data,       // dir for path observables
        //        obs2,
        //        &shortstep_cmp,
        //        settings.el,
        //        settings.should_write_data,
        //        settings.dtau,    // step in affine parameter
        //        settings.icalc,   // calculate every icalc steps
        //        settings.iblock,  // average every iblock calcsteps
        //        settings.iprint,  // print status every iprint steps
        //        settings.isave,   // save path + obs every isave steps
        //        settings.itape,   // save path w/ forays every itape steps
        //        settings.max_propag_steps,    // max propagation steps
        //        settings.max_escape_steps);   // max escape steps
        simple::ShoveGeodesicSimulation shove_sim
            = simple::ShoveGeodesicSimulation(args::sim_name + "SHOVE",
                ini_path,
                fin_path,
                settings.geodir_output,     // dir for path file with no forays
                settings.geodir_output,     // dir fpr path file with forays
                settings.geodir_data,       // dir for path observables
                obs3,
                &shove_cmp,
                settings.el,
                settings.should_write_data,
                settings.dtau,     // max step in cfg space / sigma
                settings.icalc,   // calculate every icalc steps
                settings.iblock,  // average every iblock calcsteps
                settings.iprint,  // print status every iprint steps
                settings.isave,   // save path + obs every isave steps
                settings.itape,   // save path w/ forays every itape steps
                settings.max_propag_steps,    // max propagation steps
                settings.max_escape_steps);   // max escape steps
        simple::PlerpGeodesicSimulation plerp_sim
            = simple::PlerpGeodesicSimulation(args::sim_name + "PLERP",
                ini_path,
                fin_path,
                settings.geodir_output,     // dir for path file with no forays
                settings.geodir_output,     // dir fpr path file with forays
                settings.geodir_data,       // dir for path observables
                obs4,
                &plerp_cmp,
                settings.el,
                settings.should_write_data,
                settings.dtau,    // max step in cfg space / sigma
                settings.icalc,   // calculate every icalc steps
                settings.iblock,  // average every iblock calcsteps
                settings.iprint,  // print status every iprint steps
                settings.isave,   // save path + obs every isave steps
                settings.itape,   // save path w/ forays every itape steps
                settings.max_propag_steps,    // max propagation steps
                settings.max_escape_steps);   // max escape steps
        fprintf(stderr, "%s\n", "Simulation initialized");
        // save simulation settings in a file
        settings.write(settings.geodir_input, args::sim_name);
        fprintf(stderr, "%s\n", "Computing paths...");
        /**********************************************************************/
        slerp_sim.compute_path();
        fprintf(stderr, "%s\n", "SLERP path computed!");
        slerp_sim.write_run_summary();
        // shortstep_sim.compute_path();
        // fprintf(stderr, "%s\n", "ShortStep path computed!");
        // shortstep_sim.write_run_summary();
        /**********************************************************************/
        shove_sim.compute_path();
        fprintf(stderr, "%s\n", "SHOVE path computed!");
        shove_sim.write_run_summary();
        /**********************************************************************/
        plerp_sim.compute_path();
        fprintf(stderr, "%s\n", "PLERP path computed!");
        plerp_sim.write_run_summary();
        /**********************************************************************/
        // get slices of the path
        // geodesic::Path first_half = shove_sim.path().get_slice(2, 1);
        // geodesic::Path second_half = shove_sim.path().get_slice(2, 2);
        // std::string base = settings.geodir_output + args::sim_name;
        // std::string first_ini = base + "_1_ini.cfg";
        // std::string first_fin = base + "_1_fin.cfg";
        // std::string second_ini = base + "_2_ini.cfg";
        // std::string second_fin = base + "_2_fin.cfg";
        // first_half.initial().write(first_ini, true);
        // first_half.final().write(first_fin, true);
        // second_half.initial().write(second_ini, true);
        // second_half.final().write(second_fin, true);
    }
    return res;
}
