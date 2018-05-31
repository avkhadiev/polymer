// 2017 Artur Avkhadiev
/*! \file simple_polymer.cpp
* Unlike the regular molecule, a simple::Polymer is assumed to be a freely-jointed chain. All atoms have the same mass and all interatomic bond lengths are fixed and equal. A simple::Polymer does not contain any time records. The atoms it contains are of type simple::Atom
* simple::Polymer does not contain records on atomic mass or interatomic distance. Those are recorded in the state.
* a polymer has 2 representations:
*       1. via bond vectors and (RCM, VCM, m, d)
*       2. via atomic positions and velocities in CM frame and (RCM, VCM, m, d)
* These representations are implemented as derived classes from a BasePolymer
* class that contains common characteristics: (RCM, VCM, m, and d)
*/
#ifndef POLYMER_SIMPLE_POLYMER_H
#define POLYMER_SIMPLE_POLYMER_H
#include <string>
#include <sstream>
#include <vector>

#include "vector.h"
#include "simple_atom.h"
#include "simple_bond.h"

#define DEFAULT_NB 0
#define DEFAULT_M 1.0
#define DEFAULT_D 3.0
namespace simple {
    /***************************************************************************
    *                               BASE POLYMER
    ***************************************************************************/
    class BasePolymer {
    private:
        // all polymers in the simulation have the same number of bonds, nb,
        // the same mass, m, and the same bond length, d.
        static int _nb;
        // mass of each atoms
        static double _m;
        // fixed bond length for all bonds
        static double _d;
    public:
        static int nb() {return _nb;};
        static double m() {return _m;};
        static double d() {return _d;};
        static void set_nb(int nb) {_nb = nb;};
        static void set_m(double m) {_m = m;};
        static void set_d(double d) {_d = d;};
        /**
        * Ouptuts an std::string representation of the molecule
        * Depends on the chosen representation (the derived class)
        */
        virtual std::string to_string(bool verbose = true) const = 0;
        BasePolymer();
        virtual ~BasePolymer();
    };
    class AtomPolymer;
    /***************************************************************************
    *                               BOND POLYMER
    ***************************************************************************/
    class BondPolymer :
        public BasePolymer {
        private:
            void subtract_vel_projections();
            friend class AtomPolymer;
            // position of the center of mass of the molecule
            Vector _rcm;
            // velocity of the center of mass of the molecule
            Vector _vcm;
        public:
            std::vector<Bond> bonds;
            const std::vector<Bond> get_bonds() const {return bonds;};
            bool operator==(const BondPolymer &other) const;
            bool operator!=(const BondPolymer &other) const;
            Vector rcm() const {return _rcm;};
            Vector vcm() const {return _vcm;};
            void set_rcm(Vector rcm) {_rcm = rcm;};
            void set_vcm(Vector vcm) {_vcm = vcm;};
            /**
            * Ouptuts an std::string header
            * that constains RCM and VCM vectors of the molecule
            */
            std::string write_cm_str(bool verbose = true) const;
            /**
            * Reads the non-verbose header string of a molecule,
            * saves rcm and vcm
            */
            void read_cm_str(std::string non_verbose_str);
            /**
            * Ouptuts an std::string representation of the molecule
            * in the bond vector representation
            */
            virtual std::string to_string(bool verbose = true) const;
            /**
            * Will remove projections of bond velocities on the bonds.
            */
            BondPolymer(const std::vector<Bond>& new_bonds = {},
                Vector rcm = vector(0.0, 0.0, 0.0),
                Vector vcm = vector(0.0, 0.0, 0.0));
            /**
            * Creates a polymer in the bond representation using the atom representation;
            * Subtracts average velocity from each velocity
            * (total momentum in CM frame is 0)
            * Removes projection of bond velocities on the bonds.
            */
            BondPolymer(const AtomPolymer& polymer);
            ~BondPolymer();
    };
    /***************************************************************************
    *                               ATOM POLYMER
    ***************************************************************************/
    class AtomPolymer :
        public BasePolymer {
        private:
            void subtract_average_vel();
            friend class BondPolymer;
        public:
            std::vector<Atom> atoms;
            const std::vector<Atom>& get_atoms() const {return atoms;};
            bool operator==(const AtomPolymer &other) const;
            bool operator!=(const AtomPolymer &other) const;
            Vector rcm() const;
            Vector vcm() const;
            /**
            * Ouptuts an std::string representation of the molecule
            * in the bond vector representation
            */
            virtual std::string to_string(bool verbose = true) const;
            /**
            * Creates a polymer in atomic representation
            * Subtracts average velocity from each velocity
            * (total momentum in CM frame is 0)
            * Assumes bond velocities projections on the bonds are zero
            */
            AtomPolymer(const std::vector<Atom>& new_atoms = {});
            /**
            * Creates a polymer in the atom representation using the bond representation
            * Removes projection of bond velocities on the bonds.
            */
            AtomPolymer(const BondPolymer& polymer);
            ~AtomPolymer();
    };
    /**
    * Takes a polymer and ouptuts its std::string representation
    */
    ::std::ostream& operator<<(::std::ostream& os,
        const simple::BondPolymer& polymer);
    ::std::ostream& operator<<(::std::ostream& os,
        const simple::AtomPolymer& polymer);
    /**
    * Reads the stream from the current position, advancing the position until after the molecule instance; converts input data to a molecule instance
    */
    BondPolymer string_to_bond_polymer(std::ifstream& input_stream);
    AtomPolymer string_to_atom_polymer(std::ifstream& input_stream);
}   // namespace simple
#endif
