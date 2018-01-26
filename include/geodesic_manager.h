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
        Manager(ForceUpdater* fupd = NULL);
        ~Manager();
        void read_states(std::string cndir, std::string sim_name,
            ForceUpdater* fupd = NULL);
        void write_geodesic_inputs(std::string outdir) const;
        Record initial() const;
        Record final() const;
    protected:
        std::vector<simple::BondState> _states;
        NAtomicConfigHandler _config_handler;
    };
} // namespace geodesic
#endif
