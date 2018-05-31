// 2018 Artur Avkhadiev
/*! \file geodesic_manager.h
*/
#ifndef GEODESIC_MANAGER_H
#define GEODESIC_MANAGER_H
#include <vector>
#include <sys/stat.h>
#include <string>
#include <cassert>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include "geodesic_record.h"
#include "geodesic_path.h"
#include "natomic_config_handler.h"
#include "observable_container.h"
#include "general_observables.h"
namespace geodesic{
    /***************************************************************************
    *                          GEODESIC MANAGER
    ***************************************************************************/
    /**
    * Transforms configuration files into geodesic:: inputs suitable for
    * computing geodesic paths
    * The manager chooses input/output geodesic records from a config file,
    * and outputs those values into files in a given directory
    */
    class Manager {
    public:
        Manager();
        Manager(Potential* polymer_potential,
            Potential* solvent_potential,
            Potential* inter_potential);
        ~Manager();
        /**
        * Read states from cndir + sim_name + "_cn.cfg" and store the sequence
        * if successfully read 1 or more states, will set initial, final records with PE calculated  (if force updater is set)
        */
        void read_states(std::string cndir, std::string sim_name);
        void write_geodesic_inputs(std::string outdir, std::string sim_name) const;
        /**
        * Will turn the sequence of states (from the MD trajectory)
        * into a geodesic::Path with PE calculated along the way
        */
        Path MD_path();
        void output_observables(Path& path,
                                std::string dtdir, std::string name);
        void set_fupd(ForceUpdater& fupd);
        ForceUpdater fupd() const;
        std::string ini_fname(std::string sim_name) const;
        std::string fin_fname(std::string sim_name) const;
        Record initial() const;
        Record final() const;
    protected:
        ForceUpdater _fupd;
        std::vector<simple::BondState> _states;
        bool _states_read;
        Record _initial;
        Record _final;
        double _pe(simple::BondState& state);
    };
} // namespace geodesic
#endif
