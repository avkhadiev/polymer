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
        polymer::PE pe_polymer = polymer::PE(calc_mean, calc_err, print_val);
        polymer::Virial w_polymer = polymer::Virial(calc_mean, calc_err, print_val);
        polymer::ConstraintVirial wc_polymer
            = polymer::ConstraintVirial(calc_mean, calc_err, print_val);
        std::vector<Observable*> observables_vec
            = { &pe_polymer,
                &w_polymer, &wc_polymer };
        ObservableContainer container = ObservableContainer(observables_vec);
        ObservableContainer& obs = container;
        ObservableContainer container2 = ObservableContainer(observables_vec);
        ObservableContainer& obs2 = container2;
        // setup potentials
        LJPotential pp = LJPotential(settings.epp, settings.sp,
            &pe_polymer, &w_polymer);
        SolventConfigHandler solvent_config = SolventConfigHandler(settings.ss, settings.rho_s, settings.nc);
        AdjustedLJPotential ss = AdjustedLJPotential(settings.ess, settings.ss, settings.rc_ss, solvent_config.box());
        AdjustedLJPotential ps = AdjustedLJPotential(settings.eps, 0.5 * (settings.ss + settings.sp), settings.rc_ps, solvent_config.box());
        // create endpoint records in separate files
        fprintf(stderr, "%s\n", "Creating endpoint records from MD file...");
        geodesic::Manager manager = geodesic::Manager(&pp, &ss, &ps);
        manager.read_states(settings.cndir, args::sim_name);
        manager.write_geodesic_inputs(settings.geodir_input, args::sim_name);
        std::string ini_path
            = settings.geodir_input + manager.ini_fname(args::sim_name);
        std::string fin_path
            = settings.geodir_output + manager.fin_fname(args::sim_name);
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
            md_path.write(md_file, overwrite);
            manager.output_observables(md_path, settings.geodir_output, args::sim_name + "md");
        }
        // initialize path computers
        geodesic::SLERP slerp_cmp = geodesic::SLERP(&pp, &ss, &ps, 0.0001);
        obs.add_observable(&slerp_cmp.psi1);
        obs.add_observable(&slerp_cmp.psi2);
        obs.add_observable(&slerp_cmp.delta_psi1);
        obs.add_observable(&slerp_cmp.delta_psi2);
        geodesic::ShortStep shortstep_cmp
            = geodesic::ShortStep(&pp, &ss, &ps, 0.0001);
        obs2.add_observable(&shortstep_cmp.psi1);
        obs2.add_observable(&shortstep_cmp.psi2);
        obs2.add_observable(&shortstep_cmp.delta_psi1);
        obs2.add_observable(&shortstep_cmp.delta_psi2);
        obs2.add_observable(&shortstep_cmp.theta1);
        obs2.add_observable(&shortstep_cmp.theta2);
        obs2.add_observable(&shortstep_cmp.delta_theta1);
        obs2.add_observable(&shortstep_cmp.delta_theta2);
        fprintf(stderr, "%s\n", "Path computers initialized");
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
        simple::ShortStepGeodesicSimulation shortstep_sim
            = simple::ShortStepGeodesicSimulation(args::sim_name + "ShortStep",
                ini_path,
                fin_path,
                settings.geodir_output,     // dir for path file with no forays
                settings.geodir_output,     // dir fpr path file with forays
                settings.geodir_data,       // dir for path observables
                obs2,
                &shortstep_cmp,
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
        fprintf(stderr, "%s\n", "Simulation initialized");
        // save simulation settings in a file
        settings.write(settings.geodir_input, args::sim_name);
        fprintf(stderr, "%s\n", "Computing paths...");
        slerp_sim.compute_path();
        fprintf(stderr, "%s\n", "SLERP path computed!");
        slerp_sim.write_run_summary();
        shortstep_sim.compute_path();
        fprintf(stderr, "%s\n", "ShortStep path computed!");
        shortstep_sim.write_run_summary();
    }
    return res;
}
