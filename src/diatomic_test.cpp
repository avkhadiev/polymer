// 2017 Artur Avkhadiev
/*! \file diatomic_observables_simulation.cpp
*/
#include<string>
#include<vector>
#include "../include/simple_simulation.h"
#include "../include/diatomic_config_handler.h"
#include "../include/diatomic_observables.h"
#include "../include/dynamic_observables.h"
const double PI = 3.14159265359;
namespace parameters{
    int NM = 1;
    int NB = 1;                          /**> is a diatomic */
    double M = 1.0;                      /**> specifies unit of mass */
    double D = 1.0;                      /**> specifies unit of length */
}
namespace settings{
    std::string sim_name = "diatomic";
    std::string cndir = "/Users/Arthur/stratt/polymer/test/";
    std::string tpdir = "/Users/Arthur/stratt/polymer/test/";
    std::string dtdir = "/Users/Arthur/stratt/polymer/test/";
    std::string infile_random = "random";
    std::string infile = infile_random;
    double runtime;
    double dt = 0.001;                         /**> in units of period T */
    size_t icalc = 10;
    size_t iprint = 1000;
    size_t isave = 100;
    size_t idata = 100;
    size_t itape = 10;
}
namespace rattle{
    double tol = pow(10, -4);                  /**> specified via input */
    double tiny = pow(10, -7);
    int maxiter = pow(10, 7);
}
namespace observables{
    // declare observables
    diatomic::BondPosLength bpl;
    diatomic::BondVelLength bvl = diatomic::BondVelLength();
    diatomic::BondPosProj bpp = diatomic::BondPosProj();
    diatomic::BondVelProj bvp = diatomic::BondVelProj();
    dynamic::IntKE ik = dynamic::IntKE(parameters::M);
    dynamic::LNorm ln = dynamic::LNorm(parameters::M);
    dynamic::LProj lp = dynamic::LProj(parameters::M);
    void setup(DiatomicConfigHandler& cfg){
        double m = parameters::M;
        Vector paxis = cfg.ini_bp_axis();
        Vector vaxis = cfg.ini_bv_axis();
        Vector laxis = lp.L(cfg.atom_state());// outputs angular momentum vector
        laxis = divide(laxis, norm(laxis));
        bpl = diatomic::BondPosLength();
        bvl = diatomic::BondVelLength();
        bpp = diatomic::BondPosProj(paxis);
        bvp = diatomic::BondVelProj(vaxis);
        ik = dynamic::IntKE(m);
        ln = dynamic::LNorm(m);
        lp = dynamic::LProj(m, laxis);
    }
    // define status variables and time of evalation
    container::Unit bpl_unit = {.obs = bpl,
        .status = true, .timelog = true, .average = false,
        .eval_time = container::EvalTime::sim_loop};
    container::Unit bvl_unit = {.obs = bvl,
        .status = true, .timelog = true, .average = false,
        .eval_time = container::EvalTime::sim_loop};
    container::Unit bpp_unit = {.obs = bpp,
        .status = true, .timelog = true, .average = false,
        .eval_time = container::EvalTime::sim_loop};
    container::Unit bvp_unit = {.obs = bvp,
        .status = true, .timelog = true, .average = false,
        .eval_time = container::EvalTime::sim_loop};
    container::Unit ik_unit = {.obs = ik,
        .status = true, .timelog = true, .average = false,
        .eval_time = container::EvalTime::sim_loop};
    container::Unit ln_unit = {.obs = ln,
        .status = true, .timelog = true, .average = false,
        .eval_time = container::EvalTime::sim_loop};
    container::Unit lp_unit = {.obs = lp,
        .status = true, .timelog = true, .average = false,
        .eval_time = container::EvalTime::sim_loop};
    /************************************************************************/
    std::vector<container::Unit> vec
        = {bpl_unit, bvl_unit, bpp_unit, bvp_unit, ik_unit, ln_unit, lp_unit};
}
namespace messages{
    std::string usage
        = "Usage: diatomic_test <name> <infile> <runtime/T> <timestep/T> <tol> <optional icalc iprint isave idata itape>";
    std::string infile
        = "If infile == 'random', diatomic will be initialized in a random orientation; otherwise, the infile should specify a path to the configuration file";
}
int main(int argc, char **argv){
    if (argc < 6){
        fprintf(stderr, "%s\n", messages::usage.c_str());
        fprintf(stderr, "%s\n", messages::infile.c_str());
    }
    else {
        settings::sim_name = argv[1];
        settings::infile = argv[2];                     /*> input config   */
        settings::runtime = atof(argv[3]);              /**>  run time     */
        settings::dt = atof(argv[4]);                   /**>  timestep     */
        rattle::tol = atof(argv[5]);                    /**> tolerance     */
        if (argc == 11) {
            settings::icalc = atoi(argv[6]);     /**> data calc interval   */
            settings::iprint = atoi(argv[7]);    /**> state print interval */
            settings::isave = atoi(argv[8]);     /**> config save interval */
            settings::idata = atoi(argv[9]);     /**> data write interval  */
            settings::itape = atoi(argv[10]);    /**> tape write interval  */
        }
        // state setup
        simple::BaseState::set_nm(parameters::NM);
        simple::BasePolymer::set_nb(parameters::NB);
        simple::BasePolymer::set_m(parameters::M);
        simple::BasePolymer::set_d(parameters::D);
        // initialize in state
        if (settings::infile == settings::infile_random){
            // no input configuration is required; state will be initialized
            // at random
            settings::infile = "";
        }
        DiatomicConfigHandler config = DiatomicConfigHandler();
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
                settings::dt * 2 * PI,
                settings::icalc,
                settings::iprint,
                settings::isave,
                settings::idata,
                settings::itape);
        if (!sim.is_input_given()){
            // The state is now set up.
            // If it was NOT set up from input file,
            //  observables may need the initial-state information. For example,
            //  the orientation of the angular momentum vector initially needs
            //  to be given to the LProj observable.
            // If it was set up from the input file, this information need not
            //  be re-written.
            observables::setup(config);
        }
        // TODO main evolve loop is preceded by relaxation if necessary
        sim.evolve(settings::runtime * 2 * PI);
    }
}
