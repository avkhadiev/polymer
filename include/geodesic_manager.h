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
        Manager(ForceUpdater* fupd = NULL); /**> need fupd to calculate PE */
        ~Manager();
        /**
        * Read states from cndir + sim_name + "_cn.cfg" and store the sequence
        * if successfully read 1 or more states, will set initial, final records with PE calculated  (if force updater is set)
        */
        void read_states(std::string cndir, std::string sim_name);
        void write_geodesic_inputs(std::string outdir, std::string sim_name) const;
        /**
        * TODO
        * Will turn the sequence of states (from the MD trajectory)
        * into a geodesic::Path with PE calculated along the way
        */
        Path MD_path();
        void set_fupd(ForceUpdater *fupd);
        ForceUpdater fupd() const;
        Record initial() const;
        Record final() const;
    protected:
        NAtomicConfigHandler _config_handler;
        std::vector<simple::BondState> _states;
        bool _states_read;
        Record _initial;
        Record _final;
        void _load_state_to_cfg_manager(simple::BondState& state);
        double _pe(simple::BondState& state);
    };
} // namespace geodesic
#endif
