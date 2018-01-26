// 2018 Artur Avkhadiev
/*! \file geodesic_record.h
*/
#ifndef GEODESIC_RECORD_H
#define GEODESIC_RECORD_H
#include <string>
#include <fstream>
#include "parsing.h"
#include "simple_state.h"
namespace geodesic{
    /***************************************************************************
    *                               GEODESIC RECORD
    ***************************************************************************/
    class Record {
    public:
        Record();
        Record(simple::BondState state, double pe);
        Record(std::ifstream& input_stream, bool read_state_header = true);
        Record(std::string file);       /**> state header always read */
        ~Record();
        bool operator==(const Record &other) const;
        bool operator!=(const Record &other) const;
        simple::BondState state() const {return _state;};
        double pe() const {return _pe;};
        /* Record output */
        std::string to_string(bool output_header) const;
        /* writes to file stream can output state header */
        void write(std::ofstream& input_stream, bool output_header) const;
        /* writes to file, truncates and outputs header if overwrite = true */
        void write(std::string file, bool overwrite) const;
    private:
        simple::BondState _state;
        double _pe;
    };
    ::std::ostream& operator<<(::std::ostream& os, const Record& record);
} // namespace geodesic
#endif
