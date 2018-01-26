// 2018 Artur Avkhadiev
/*! \file geodesic_path_computer.h
*/
#ifndef GEODESIC_PATH_COMPUTER_H
#define GEODESIC_PATH_COMPUTER_H
#include <vector>
#include <sys/stat.h>
#include <string>
#include <cassert>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include "geodesic_record.h"
#include "geodesic_path.h"
#include "geodesic_manager.h"
#include "observable.h"
#include "potential.h"
#include "force_updater.h"
namespace geodesic{
    /***************************************************************************
    *                          GEODESIC PATH COMPUTER
    ***************************************************************************/
    class PathComputer{
    public:
        PathComputer(double landscape_energy, double dt,
            Potential* polymer_potential,
            Potential* solvent_potential,
            Potential* inter_potential);
        ~PathComputer();
        double landscape_energy() const {return _E_l;};
        void set_landscape_energy(double energy) {_E_l = energy;};
        double dt() const {return _dt;};
        void set_dt(double dt) {_dt = dt;};
        /**
        * Path computation methods
        */
        Path naive_SLERP(Record initial, Record final);
        Path perturbative_SLERP(Record initial, Record final);
        Path perturbative_new(Record initial, Record final);
    private:
        ForceUpdater _fupd;
        double _E_l;
        double _dt;
    };
} // namespace geodesic
#endif
