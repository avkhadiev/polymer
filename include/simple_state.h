// 2017 Artur Avkhadiev
/*! \file simple_state.h
*/
#ifndef POLYMER_SIMPLE_STATE_H
#define POLYMER_SIMPLE_STATE_H
#include <vector>
#include <string>
#include <sstream>
#include "simple_polymer.h"
namespace simple {
    class BaseState {
    protected:
        static int _nm;              /*>> number of polymer molecules */
        /**
        * representations of the state will have data members describing
        * positions and velocities of molecules; _time is the time that those
        * records correspond to
        */
        double _time;                /*>> last update time for representation */
        /**
        * Outputs a string (possibly multiline) with the state's molecules in
        * the required representation
        */
        virtual std::string _molecules_str(bool verbose = true) const = 0;
        /**
        * Ouptuts a header information about the state:
        *   number of polymer molecules,
        *   number of bonds in them,
        *   mass of their atoms, and
        *   the interatomic bond length
        */
        std::string _header_str(bool verbose = true) const;
        /**
        * Ouptuts the time that the records of the state correspond to
        */
        std::string _time_str(bool verbose = true) const;
        /**
        * Outputs an std::string representation of the state
        *   Requires the molecule string prepared and passed as argument
        */
        std::string _to_string_helper(std::string molecules_str,
            bool verbose = true, bool output_header = true) const;
        /**
        * Writes an std::string representation of the state
        *   Requires the molecule string prepared and passed as argument
        */
        void _write_to_file_helper(std::string molecules_str,
            std::string outdir,
            std::string fname,
            bool verbose = true,
            bool overwrite = true) const;
    public:
        // getters
        static int nm() {return _nm;};
        static int polymer_nb() {return BasePolymer::nb();};
        static double polymer_mass() {return BasePolymer::m();};
        static double polymer_bondlength() {return BasePolymer::d();};
        double time() const {return _time;};
        // setters
        static void set_nm(int nm) {_nm = nm;};
        /**
        * Base state is an abstract class, so all non-static functions are
        * only called from instances of derived classes
        *     -- representations of the state.
        * If the state time is advanced from a representation, all
        * representation-specific data memebers (like polymer representations)
        * are assumed to be the most recent versions corresponding to that time.
        * Consequently, representation time also has to be updated.
        * On the contrary, if there is more than one representation of the base
        * state, a representation whose data members do not correspond to the
        * most recent time (for example, if the representation is not used to
        * by an integrator to evolve all molecules in time) would have a lag
        * between its representation time and the state time
        */
        void set_time(double time) {_time = time;};
        void advance_time(double timestep) {_time += timestep;};
        /**
        * Reads the non-verbose header string of a state,
        * saves nm, nb, m, and d
        */
        void static read_header(std::string non_verbose_header);
        /**
        * Reads the non-verbose line with time and updates _time data member
        */
        void read_time(std::string non_verbose_time);
        /**
        * Outputs an std::string representation of the state
        *   just call _to_string_helper and pass the output of _molecules_str
        */
        virtual std::string to_string(bool verbose = true,
            bool output_header = true) const = 0;
        /**
        * If the file outdir/fname.cfg does not exist, creates it.
        * Otherwise, if overwrite = false, appends the state to the end of the file; if overwrite=true, overwrites the file; adds header information.
        *   just call _write_to_file_helper and pass the output of
        * _molecules_str
        */
        virtual void write_to_file(std::string outdir,
            std::string fname,
            bool verbose = false,
            bool overwrite = false) = 0;
        BaseState(double time = 0.0);
        ~BaseState();
    };
    class AtomState;
    class BondState :
        public BaseState {
    protected:
        virtual std::string _molecules_str(bool verbose = true) const;
    public:
        std::vector<BondPolymer> polymers;
        bool operator==(const BondState &other) const;
        bool operator!=(const BondState &other) const;
        virtual std::string to_string(bool verbose = true,
            bool output_header = true) const;
        virtual void write_to_file(std::string outdir,
            std::string fname,
            bool verbose = false,
            bool overwrite = false);
        /**
        * updates polymers of the state according to polymers in the atom
        * representation
        */
        void update(const AtomState &s);
        BondState(const std::vector<BondPolymer>& new_polymers = {},
            double time = 0.0);
        BondState(const AtomState &s);
        ~BondState();
    };
    class AtomState :
        public BaseState {
    protected:
        virtual std::string _molecules_str(bool verbose = true) const;
    public:
        std::vector<AtomPolymer> polymers;
        bool operator==(const AtomState &other) const;
        bool operator!=(const AtomState &other) const;
        virtual std::string to_string(bool verbose = true,
            bool output_header = true) const;
        virtual void write_to_file(std::string outdir,
            std::string fname,
            bool verbose = false,
            bool overwrite = false);
        /**
        * updates polymers of the state according to polymers in the bond
        * representation
        */
        void update(const BondState &s);
        AtomState(const std::vector<AtomPolymer>& new_polymers = {},
            double time = 0.0);
        AtomState(const BondState &s);
        ~AtomState();
    };
    ::std::ostream& operator<<(::std::ostream& os, const BondState& s);
    ::std::ostream& operator<<(::std::ostream& os, const AtomState& s);
    /**
    * Reads the stream from the current position, advancing the position until
    * after the state struct; converts input data to a state struct.
    * assumes no header!
    */
    BondState string_to_bond_state(std::ifstream& input_stream);
    AtomState string_to_atom_state(std::ifstream& input_stream);
    /**
    * Opens a file outdir/fname.cfg for reading and reads states from the file
    * into a vector of states. The file should contain non-verbose string
    * representations of states.
    */
    void read_states_from_file(std::string indir,
        std::string fname, std::vector<BondState>& states);
    void read_states_from_file(std::string indir,
        std::string fname, std::vector<AtomState>& states);
}   // namespace simple
#endif
