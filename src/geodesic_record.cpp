// 2018 Artur Avkhadiev
/*! \file geodesic_record.cpp
*/
#include "../include/geodesic_record.h"
namespace geodesic{
    /***************************************************************************
    *                               GEODESIC RECORD
    ***************************************************************************/
    Record::Record():
        _state(),
        _pe(){}
    Record::Record(simple::BondState state, double pe):
        _state(state),
        _pe(pe){}
    Record::Record(std::ifstream& input_stream, bool read_state_header):
        _state(),
        _pe()
    {
        std::string line;
        // read header information first, if it is provided
        if (read_state_header){
            std::getline(input_stream, line);
            simple::BaseState::read_header(line);
        }
        _state = simple::string_to_bond_state(input_stream);
        // read potential energy
        std::getline(input_stream, line);
        _pe = std::stof(line);
    }
    Record::Record(std::string file) :
        _state(),
        _pe(0.0)
    {
        std::ifstream readout;
        readout.open(file, std::ifstream::in);
        if (!readout.is_open()) {
            std::string err_msg = "Record: unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), file.c_str());
            perror("open");
        }
        else{
            bool read_state_header = true;
            *this = Record(readout, read_state_header);
        }
    }
    Record::~Record(){
    }
    bool Record::operator==(const Record &other) const{
        //bool states_equal = (state() == other.state());
        //bool energies_equal = (pe() == other.pe());
        //return states_equal && energies_equal;
        // FIXME cannot compare "states" because they include both time and
        // velocities, whereas
        // configuration space does not retain this information
        return true;
    }
    bool Record::operator!=(const Record &other) const{
        return !(*this == other);
    }
    std::string Record::to_string(bool output_header, bool verbose) const {
        std::string record_string
            = state().to_string(verbose, output_header) + std::to_string(pe());
        return record_string;
    }
    ::std::ostream& operator<<(::std::ostream& os, const Record& record){
        bool output_header = true;
        // this function is used for debugging and verbosity helps
        bool verbose = true;
        return os << record.to_string(output_header, verbose).c_str();
    }
    void Record::write(std::ofstream& writeout, bool output_header) const {
        writeout << to_string(output_header) << std::endl;
    }
    void Record::write(std::string fout, bool overwrite) const{
        std::ofstream writeout;
        bool overwritten = false;
        // open writeout for output operations and
        // set the stream's position indicator to the end of the stream before each output operation.
        if (overwrite) {
            // if overwrite is allowed, try to open in truncate mode
            writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
            overwritten = true;
        }
        if (!overwrite) {
            // if overwrite is forbidden, try to open the file in append mode
            writeout.open(fout, std::ofstream::out | std::ofstream::app);
            if (!writeout.is_open()) {
            // if file does not exist, open in truncate mode
                writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
                overwritten = true;
            }
        }
        // now file has to be opened
        if (writeout.is_open()) {
            // if file could be opened...
            // if overwrote the file, output header
            bool output_header;
            if (overwritten) {
                output_header = true;
            }
            else{
                output_header = false;
            }
            write(writeout, output_header);
        }
        else {
            // if file still could not be opened
            std::string err_msg = "write_state_to_file: unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fout.c_str());
            perror("open");
        }
        writeout.close();
    }
} // namespace geodesic
