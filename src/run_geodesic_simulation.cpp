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
        // normalize tolerance for entire chain - CANCELLED: CONDITION TOO TIGHT GIVEN TOLERANCE
        // settings.epsilon = settings.epsilon / sqrt(settings.nb + 1);
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
        // no observable output except at the end
        settings.isave = 0;
        // setup observables
        std::vector<Observable*> observables_vec = {};
        ObservableContainer container = ObservableContainer(observables_vec);
        ObservableContainer& obs = container;
        ObservableContainer container2 = ObservableContainer(observables_vec);
        std::vector<Observable*> observables_vec3 = {};
        ObservableContainer container3 = ObservableContainer(observables_vec3);
        ObservableContainer& obs3 = container3;
        std::vector<Observable*> observables_vec4 = {};
        ObservableContainer container4 = ObservableContainer(observables_vec3);
        ObservableContainer& obs4 = container4;
        // setup potentials
        LJPotential pp = LJPotential(settings.epp, settings.sp);
        SolventConfigHandler solvent_config = SolventConfigHandler(settings.ss, settings.rho_s, settings.nc);
        AdjustedLJPotential ss = AdjustedLJPotential(settings.ess, settings.ss, settings.rc_ss, solvent_config.box());
        AdjustedLJPotential ps = AdjustedLJPotential(settings.eps, 0.5 * (settings.ss + settings.sp), settings.rc_ps, solvent_config.box());
        // create endpoint records in separate files
        geodesic::Manager manager = geodesic::Manager(&pp, &ss, &ps, &container2);
        manager.read_states(settings.cndir, args::sim_name);
        manager.write_geodesic_inputs(settings.geodir_input, args::sim_name);
        std::string ini_path
            = settings.geodir_input + manager.ini_fname(args::sim_name);
        std::string fin_path
            = settings.geodir_input + manager.fin_fname(args::sim_name);
        if (settings.reverse_path){
            std::string temp = ini_path;
            ini_path = fin_path;
            fin_path = temp;
        }
        //fprintf(stderr, "%s:\n%s: %s\n%s: %s\n",
        //    "Endpoint files are created",
        //    "initial record",
        //    ini_path.c_str(),
        //    "final record",
        //    fin_path.c_str());
        // output MD path file, if necessary
        if (settings.save_md_path) {
            // removed writing path
            //bool overwrite = true;
            std::string mdfile = args::sim_name + "md_path.txt";
            std::string md_file
                = settings.geodir_output + mdfile;
            //fprintf(stderr, "%s %s\n",
            //    "Output MD path file to", md_file.c_str());
            geodesic::Path md_path = manager.MD_path();
            if (settings.reverse_path) {
                md_path.reverse();
            }
            fprintf(stdout, "%5.7f\n", md_path.euc_sep.value);
            fprintf(stdout, "%5.7f\n",
                md_path.euc_sep.value / sqrt((simple::BasePolymer::nb() + 1)));
            fprintf(stdout, "%5.7f\n", md_path.length.value);
            //md_path.write(md_file, overwrite);
            if (settings.should_write_data){
                manager.output_observables(md_path, settings.geodir_output, args::sim_name + "md");
            }
        }
        // initialize path computers
        geodesic::SLERP slerp_cmp
            = geodesic::SLERP(&pp, &ss, &ps, settings.epsilon);
        geodesic::SHOVE shove_cmp
            = geodesic::SHOVE(&pp, &ss, &ps, settings.tol, settings.epsilon);
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
        // save simulation settings in a file
        settings.write(settings.geodir_input, args::sim_name);
        /**********************************************************************/
        slerp_sim.compute_path();
        slerp_sim.write_run_summary();
        /**********************************************************************/
        shove_sim.compute_path();
        shove_sim.write_run_summary();
        /**********************************************************************/
        plerp_sim.compute_path();
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
