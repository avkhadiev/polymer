// 2018 Artur Avkhadiev
/*! \file geodesic_record.cpp
*/
#include "../include/geodesic_record.h"
namespace geodesic{
    /***************************************************************************
    *                               GEODESIC RECORD
    ***************************************************************************/
    Record::Record():
        cfg_handler(),
        _pe(){}
    Record::Record(simple::BondState state, double pe):
        cfg_handler(state),
        _pe(pe){}
    Record::Record(simple::AtomState state, double pe):
        cfg_handler(state),
        _pe(pe){}
    Record::Record(std::ifstream& input_stream,
        bool read_state_header,
        bool atom_state):
        Record::Record()
    {
        std::string line;
        // read header information first, if it is provided
        if (read_state_header){
            std::getline(input_stream, line);
            simple::BaseState::read_header(line);
        }
        // if reading in atomic state representation
        if (atom_state){
            simple::AtomState state
                = simple::string_to_atom_state(input_stream);
            cfg_handler.set_state(state);
        }
        // else reading in bond state representation
        else {
            simple::BondState state
                = simple::string_to_bond_state(input_stream);
            cfg_handler.set_state(state);
        }
        // read potential energy
        std::getline(input_stream, line);
        _pe = std::stof(line);
    }
    Record::Record(std::string file, bool atom_state) :
        Record::Record()
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
            *this = Record(readout, read_state_header, atom_state);
        }
    }
    Record::~Record(){
    }
    bool Record::operator==(const Record &other) const{
        simple::AtomState this_state = cfg_handler.atom_state();
        simple::AtomState other_state = other.cfg_handler.atom_state();
        // check polymer equality
        bool polymers_equal = true;
        simple::AtomPolymer this_polymer;
        simple::AtomPolymer other_polymer;
        for(int i = 0; i < simple::BaseState::nm(); ++i){
            this_polymer = this_state.polymers.at(i);
            other_polymer = this_state.polymers.at(i);
            for(int j = 0; j < simple::BasePolymer::nb() + 1; ++j){
                polymers_equal = polymers_equal &&
                (this_polymer.atoms.at(j).position
                    == other_polymer.atoms.at(j).position);
            }
        }
        bool solvents_equal = true;
        for(int i = 0; i < simple::BaseState::nsolvents(); ++i){
            solvents_equal = solvents_equal &&
                (this_state.solvents.at(i).atom.position
                    == other_state.solvents.at(i).atom.position);
        }
        bool pe_equal = (pe() == other.pe());
        return pe_equal && polymers_equal && solvents_equal;
    }
    bool Record::operator!=(const Record &other) const{
        return !(*this == other);
    }
    simple::BondState& Record::bond_state(){
        return cfg_handler.bond_state();
    }
    simple::AtomState& Record::atom_state(){
        return cfg_handler.atom_state();
    }
    std::string Record::to_string(bool output_header,
        bool verbose,
        bool atom_state) const{
        std::string state_string;
        // if using atom_state representation
        if (atom_state) {
            state_string
                = cfg_handler.atom_state().to_string(verbose, output_header);
        }
        else {
            state_string
                = cfg_handler.bond_state().to_string(verbose, output_header);
        }
        // else using bond_state representation
        std::string record_string = state_string + std::to_string( pe() );
        return record_string;
    }
    ::std::ostream& operator<<(::std::ostream& os, Record& record){
        bool output_header = true;
        // this function is used for debugging and verbosity helps
        bool verbose = true;
        return os << record.to_string(output_header, verbose).c_str();
    }
    void Record::write(std::ofstream& writeout,
        bool output_header, bool verbose, bool atom_state) const{
        writeout << to_string(output_header, verbose, atom_state) << std::endl;
    }
    void Record::write(std::string fout,
        bool overwrite, bool verbose, bool atom_state) const{
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
            write(writeout, output_header, verbose, atom_state);
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
