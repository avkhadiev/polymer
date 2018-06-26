// 2017 Artur Avkhadiev
/*! \file run_md_simulation.cpp
* Usage: run_md_simulation <name> <config_file> <optional initial_state_file>
*/
#include<string>
#include<vector>
#include "../include/default_macros.h"
#include "../include/settings_parser.h"
// particular potentials and integrators
#include "../include/potential.h"
#include "../include/ljpotential.h"
#include "../include/verlet_integrator.h"
#include "../include/rattle_integrator.h"
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
    std::string initial_state = "";
    std::string config_file;
    std::string sim_name;
    bool is_input_given = false;
}
namespace messages{
    std::string usage = "Usage: run_md_simulation <name> <config_file> <optional initial_state_file>";
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
        // fprintf(stderr, "%s\n", "settings parser set up");
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
        // make observables
        bool calc_mean = true;
        bool calc_err = true;
        bool print_val = true;
        bool remove_linear_momentum = true;
        bool remove_angular_momentum = false;
        observable::Time time = observable::Time();
        polymer::KE ke_polymer = polymer::KE(calc_mean, calc_err, print_val);
        polymer::PE pe_polymer = polymer::PE(calc_mean, calc_err, print_val);
        polymer::Virial w_polymer = polymer::Virial(calc_mean, calc_err, print_val);
        polymer::ConstraintVirial wc_polymer
            = polymer::ConstraintVirial(calc_mean, calc_err, print_val);
        polymer::KineticTemperature temp_kin_polymer
            = polymer::KineticTemperature(
                remove_linear_momentum, remove_angular_momentum,
                calc_mean, calc_err, print_val);
        polymer::LinMomComponent px
            = polymer::LinMomComponent(vector(1.0, 0.0, 0.0),
                !calc_mean, !calc_err, !print_val);
        polymer::LinMomComponent py
            = polymer::LinMomComponent(vector(0.0, 1.0, 0.0),
                !calc_mean, !calc_err, !print_val);
        polymer::LinMomComponent pz
            = polymer::LinMomComponent(vector(0.0, 0.0, 1.0),
                !calc_mean, !calc_err, !print_val);
        polymer::AngMomComponent lx
            = polymer::AngMomComponent(vector(1.0, 0.0, 0.0),
                !calc_mean, !calc_err, !print_val);
        polymer::AngMomComponent ly
            = polymer::AngMomComponent(vector(0.0, 1.0, 0.0),
                !calc_mean, !calc_err, !print_val);
        polymer::AngMomComponent lz
            = polymer::AngMomComponent(vector(0.0, 0.0, 1.0),
                !calc_mean, !calc_err, !print_val);
        std::vector<Observable*> observables_vec
            = {&time,
                &ke_polymer, &pe_polymer,
                &temp_kin_polymer,
                &w_polymer, &wc_polymer,
                &px, &py, &pz,
                &lx, &ly, &lz };
        ObservableContainer container = ObservableContainer(observables_vec);
        ObservableContainer& obs = container;
        bool is_planar;
        if (settings.polymer_planar_conformation){
            is_planar = true;
        }
        else {
            is_planar = false;
        }
        // prepare state configuration
        SolventConfigHandler solvent_config = SolventConfigHandler(settings.ss, settings.rho_s, settings.nc);
        NAtomicConfigHandler polymer_config
            = NAtomicConfigHandler();
        ConfigHandler& cfg = polymer_config;
        // potential setup
        LJPotential pp = LJPotential(settings.epp, settings.sp,
            &pe_polymer, &w_polymer);
        AdjustedLJPotential ss = AdjustedLJPotential(settings.ess, settings.ss, settings.rc_ss, solvent_config.box());
        AdjustedLJPotential ps = AdjustedLJPotential(settings.eps, 0.5 * (settings.ss + settings.sp), settings.rc_ps, solvent_config.box());
        // initialize
        if (args::is_input_given){
            fprintf(stderr, "%s\n",
            "Don't yet know how to work with when initial state configuration is given");
            exit(EXIT_FAILURE);
        }
        else {
            //solvent_config.fcc_positions();
            //solvent_config.ran_velocities(settings.temperature);
        }
        ForceUpdater force_loop = ForceUpdater(&pp, &ss, &ps);
        cfg.set_force_updater(&force_loop);
        if (settings.polymer_planar_conformation) {
            polymer_config.initialize2D(
                settings.temperature
            );
        }
        else {
            polymer_config.initialize3D(
                settings.temperature
            );
        }
        // scale time units if no potential
        if (no_potential){
            size_t ndof =
                std::max(0,
                    (3 * ((simple::BaseState::nm())
                        - (int)remove_linear_momentum
                        - (int)remove_angular_momentum)
                    + 2 * (simple::BaseState::nm()
                        * simple::BasePolymer::nb())
                    )
                );
            size_t natoms = simple::BasePolymer::nb() + 1;
            double omega_avg
                = sqrt(ndof/natoms * settings.temperature)
                    / (settings.d * settings.mp);
            double scale_factor
                = 2 * PI / omega_avg;
            double dt_old = settings.dt;
            double runtime_old = settings.runtime;
            settings.dt = dt_old * scale_factor;
            settings.runtime = runtime_old * scale_factor;
            if (VERBOSE){
                fprintf(stderr, "%s\n", "Zero potential requires change in units:");
                fprintf(stderr, "%s %5.7f %s %5.7f %s %5.7f\n",
                    "dt change from",
                    dt_old,
                    "to",
                    settings.dt,
                    "with scale factor",
                    scale_factor);
                fprintf(stderr, "%s %5.7f %s %5.7f %s %5.7f\n",
                    "runtime change from",
                    runtime_old,
                    "to",
                    settings.runtime,
                    "with scale factor",
                    scale_factor);
                }
        }
        if (DEBUG){
            fprintf(stdout, "%s\n%s",
                "State initialized by the Configuration Handler:",
                cfg.bond_state().to_string(true, true).c_str());
        }
        RattleIntegrator rattle = RattleIntegrator(force_loop,
            settings.tol,
            &ke_polymer, NULL,
            &wc_polymer,
            solvent_config.box(),
            settings.tiny, settings.maxiter);
        VerletIntegrator verlet
            = VerletIntegrator(force_loop, NULL, solvent_config.box());
        Integrator& integrator = rattle;
        // simulation setup
        // fprintf(stderr, "%s\n", "setting up simulation");
        simple::MDSimulation sim
            = simple::MDSimulation(args::sim_name,
                settings.cndir,
                settings.tpdir,
                settings.dtdir,
                args::initial_state,
                cfg,
                integrator,
                obs,
                settings.should_write_data,
                settings.dt,
                settings.icalc,
                settings.iblock,
                settings.iprint,
                settings.isave,
                settings.itape);
        // fprintf(stderr, "%s\n", "simulation set up");
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
        // save simulation settings in a file
        // settings.write(settings.cndir, args::sim_name);
        sim.relax(20.0);
        sim.evolve(settings.runtime);
        sim.write_run_summary();
    }
    return res;
}
