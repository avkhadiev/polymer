// 2017 Artur Avkhadiev
/*! \file settings_parser.h
* contains all the necessary settings / parameters,
* reads in config files requried to start a simulation
*/
#ifndef POLYMER_SETTINGS_PARSER_H
#define POLYMER_SETTINGS_PARSER_H
#include <cmath>                /**> pow */
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "parsing.h"
#include "default_macros.h"
class SettingsParser{
private:
    std::vector<std::string> get_args(std::string line);
    // read helpers
    void read_potential(std::ifstream& stream);
    void read_state_config(std::ifstream& stream);
    void read_init(std::ifstream& stream);
    void read_integration(std::ifstream& stream);
    void read_geodesic(std::ifstream& stream);
    void read_io(std::ifstream& stream);
    // write helpers
    void write_potential(std::ofstream& stream) const;
    void write_state_config(std::ofstream& stream) const;
    void write_init(std::ofstream& stream) const;
    void write_integration(std::ofstream& stream) const;
    void write_geodesic(std::ofstream& stream) const;
    void write_io(std::ofstream& stream) const;
    std::string potential_header;
    std::string state_config_header;
    std::string init_header;
    std::string integration_header;
    std::string geodesic_header;
    std::string io_header;
public:
    // potential settings
    double epp;                 /**> epsilon polymer-polymer */
    double ess;                 /**> epsilon solvent-solvent */
    double eps;                 /**> epsilon polymer-solvent */
    double sp;                  /**> sigma polymer*/
    double ss;                  /**> sigma solvent*/
    double rc_ss;               /**> cutoff for solvent-solvent*/
    double rc_ps;               /**> cutoff for polymer-solvent*/
    // state configuration settings
    double mp;                  /**> polymer mass */
    double ms;                  /**> solvent mass                           */
    int np;                     /**> number of polymers                     */
    int nb;                     /**> number of bonds                        */
    int nc;                     /**> nsolvents = 4 * nc^3                   */
    double d;                   /**> reduced in units of sp!                */
    // initialization settings
    double rho_s;               /**> density of solvents in reduced units   */
    double temperature;         /**> temperature in reduced units           */
    int polymer_planar_conformation;   /**> 0 = 3d, planar o/w              */
    double polymer_energy;      /**> total initial energy of the polymer    */
    // integration settings
    double runtime;             /**> run time in reduced units              */
    double dt;                  /**> timestep in reduced units              */
    double tol;                 /**> tolerance for rattle                   */
    const double tiny = pow(10, -8.0);
    const double maxiter = pow(10, 4.0);
    // geodesic settings
    std::string geodir_input;   /**> dir for storing geosim inputs          */
    std::string geodir_output;  /**> dir for storing geosim inputs          */
    std::string geodir_data;    /**> dir for storing geosim observables     */
    double el;                  /**> landscape energy mass                  */
    double dtau;                /**> change in progress variable            */
    double epsilon;             /**> how close you need to get              */
    int save_md_path;           /**> 0 = don't save md_path, save o/w       */
    int reverse_path;           /**> 1 = switch initial and final points    */
    const double max_propag_steps = 1000000;
    const double max_escape_steps = 1000;
    // i/o settings
    std::string cndir;          /**> directory for simulation config i/o    */
    std::string dtdir;          /**> directory for observables i/o          */
    std::string tpdir;          /**> directory for state config i/o         */
    bool should_write_data;     /**> should output observables?             */
    size_t icalc;               /**> observable update every ... steps      */
    size_t iprint;              /**> status printout every ... steps        */
    size_t isave;               /**> configuration update every ... steps   */
    size_t iblock;              /**> block avging&writeout every ... steps  */
    size_t itape;               /**> tapefile update every ... steps        */
    // functions
    void read(std::string fin);
    void write(std::string outdir, std::string sim_name) const;
    SettingsParser();
    SettingsParser(std::string fin);
    ~SettingsParser();
};
#endif
