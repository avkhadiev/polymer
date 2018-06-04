// 2018 Artur Avkhadiev
/*! \file geodesic_record.h
*/
#ifndef GEODESIC_RECORD_H
#define GEODESIC_RECORD_H
#include <string>
#include <fstream>
#include "parsing.h"
#include "simple_state.h"
#include "config_handler.h"
namespace geodesic{
    /***************************************************************************
    *                               GEODESIC RECORD
    ***************************************************************************/
    class Record {
    public:
        Record();
        Record(simple::BondState state, double pe);
        Record(simple::AtomState state, double pe);
        Record(std::ifstream& input_stream,
            bool read_state_header = true,
            bool atom_state = false);        /**> read bond state by default */
        /**> state header always read */
        Record(std::string file, bool atom_state = false);
        ~Record();
        // can't really compare records, because records have states while completely ignoring time- and velocity- related information in them
        // these operators are only needed for testing with the Google suite
        bool operator==(const Record &other) const;
        bool operator!=(const Record &other) const;
        ConfigHandler cfg_handler;
        simple::BondState& bond_state();
        simple::AtomState& atom_state();
        // simple::BondState state() const;
        // simple::BondState& modifiable_state();
        double pe() const {return _pe;};
        void set_pe(double pe) {_pe = pe;};
        /* Record output. By default, outputs bond_state representation */
        std::string to_string(bool output_header,
            bool verbose = false,
            bool atom_state = false) const;
        /* writes to file stream can output state header
        By default, outputs bond_state representation */
        void write(std::ofstream& output_stream,
            bool output_header,
            bool verbose = false,
            bool atom_state = false) const;
        /* writes to file, truncates and outputs header if overwrite = true
        By default, outputs bond_state representation */
        void write(std::string file, bool overwrite,
            bool verbose = false,
            bool atom_state = false) const;
    private:
        double _pe;
    };
    ::std::ostream& operator<<(::std::ostream& os, Record& record);
} // namespace geodesic
#endif
