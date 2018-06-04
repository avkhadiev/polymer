// 2018 Artur Avkhadiev
/*! \file geodesic_path.cpp
*/
#include "../include/geodesic_path.h"
namespace geodesic{
    /***************************************************************************
    *                               GEODESIC PATH
    ***************************************************************************/
    /***************************************************************************
    *                        CONSTRUCTORS & DESTRUCTORS
    ***************************************************************************/
    Path::Path() :
        _initial(),
        _final(),
        _path(),
        // print_inst_value, e_format
        length(true, false){}
    Path::Path(Record initial, Record final) :
        _initial(initial),
        _final(final),
        _path(),
        // print_inst_value, e_format
        length(true, false)
    {
        _path.clear();
        //this is a ``bare'', empty path with just two endpoints ---
        // the length is to be computed; so far, it is zero
        //recompute_set_theta();
        length.value = 0.0;
    }
    Path::Path(std::list<Record>& records) :
        _initial(records.front()),
        _final(records.back()),
        _path(records),
        // print_inst_value, e_format
        length(true, false){
            // omit boundary values in path record
            _path.pop_front();
            _path.pop_back();
            recompute_length();
        }
    Path::~Path() {
        _path.clear();
    }
    Path::Path(std::string fin) :
        _initial(),
        _final(),
        _path(),
        // print_inst_value, e_format
        length(true, false)
    {
        std::ifstream readout;
        readout.open(fin, std::ifstream::in);
        if (!readout.is_open()) {
            std::string err_msg = "Path: unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fin.c_str());
            perror("open");
        }
        else{
            // read header information first
            _read_header(readout);
            // read all states
            int c;                        /*>> character for peeking */
            bool read_header = false;
            _path.clear();
            size_t counter = 0;
            while (!readout.eof() && !readout.fail()) {
                // read states into a vector of states until EOF is reached
                // or reading fails for some reason
                _path.push_back(Record(readout, read_header));
                ++counter;
                // peek next character; if EOF is reached, eofbit will be set
                // and while loop will be terminated
                c = readout.peek();
                // if failbit was set, something is wrong!
                if (readout.fail()) {
                    std::string err_msg = "Path failbit was set when reading file at";
                    fprintf(stderr, "%s %s\n", err_msg.c_str(), fin.c_str());
                    perror("readout ifstream:");
                }
            }
            readout.close();
            if (_path.size() > 1){
                // initialize other members
                _initial = _path.front();
                _final = _path.back();
                // omit boundary values from path record;
                // they are saved separately
                _path.pop_front();
                _path.pop_back();
                recompute_length();
            }
            else {
                std::string err_msg = "Path: less than two records provided";
                throw std::invalid_argument(err_msg);
            }
        }
    }
    /***************************************************************************
    *                           BINARY OPERATORS
    ***************************************************************************/
    // A simplistic implementation of operator= (see better implementation below)
    Path& Path::operator= (const Path &other)
    {
        _initial = other.initial();
        _final = other.final();
        // TODO compare length
        _path = other._path;
        return *this;
    }
    /***************************************************************************
    *                                  GETTERS
    ***************************************************************************/
    size_t Path::nrecords() const {
        return _path.size() + 2;           /**> include boundary values */
    }
    Record Path::current_tail() const {
        // current tail is the record in the path that preceeds the final
        // configuration
        Record current_tail;
        if (_path.size() > 0){
            // usually it is the back of the linked list of intermediate values
            current_tail = _path.back();
        }
        else {
            // but, if there are no intermediate values,
            // then the initial value is the current tail
            current_tail = initial();
        }
        return current_tail;
    }
    /***************************************************************************
    *         FUNCTIONS AND A NESTED CLASS DEALING WITH PATH LENGTH
    ***************************************************************************/
    void Path::_increment_length( const Record& last_record,
                                  const Record& next_record){
        length.update(last_record, next_record);
    }
    void Path::recompute_length(){
        //length._zero();
        if (_path.size() == 0){
            // no intermediate values
            _increment_length(initial(), final());
        }
        else {
            // 1 intermeditate value
            _increment_length(initial(), _path.front());
            _increment_length(_path.back(), final());
            if (_path.size() > 1){
                // several intermediate valuues
                std::list<Record>::const_iterator rec;
                for (rec = _path.begin();
                     rec != std::prev(_path.end());
                     std::advance(rec, 1)){
                    _increment_length(*rec, *(std::next(rec)));
                }
            }
        }
    }
    void Path::append(const Record& new_record){
        _increment_length(current_tail(), new_record);
        _path.push_back(new_record);
    }
    void Path::reverse(){
        // reverse the path
        _path.reverse();
        // swap initial and final records
        Record temp;
        _initial = temp;
        _initial = final();
        _final = temp;
    }
    void Path::merge(const Path& new_path){
        // TODO make sure time flows in the right direction
        fprintf(stderr, "%s\n", "geodesic::Path::merge not implemented yet");
    }
    void Path::write(std::string fout, bool overwrite){
        std::ofstream writeout;
        // open writeout for output operations and
        // set the stream's position indicator to the end of the stream before each output operation.
        if (overwrite) {
            // if overwrite is allowed, try to open in truncate mode
            writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
        }
        if (!overwrite) {
            // if overwrite is forbidden, try to open the file in append mode
            writeout.open(fout, std::ofstream::out | std::ofstream::app);
            if (!writeout.is_open()) {
                // if file does not exist, open in truncate mode
                writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
                overwrite = true;
            }
        }
        // now file has to be opened
        if (writeout.is_open()) {
            // if file could be opened...
            if (overwrite) {
                writeout << _header_str() << std::endl;
                overwrite = false;
            }
            std::list<Record>::const_iterator rec;
            initial().write(writeout, overwrite);
            for (rec = _path.begin(); rec != _path.end(); std::advance(rec, 1)){
                (*rec).write(writeout, overwrite);
            }
            final().write(writeout, overwrite);
        }
        else {
            // if file still could not be opened
            std::string err_msg = "geodesic::Path::write(): unable to open file at";
            fprintf(stderr, "%s %s\n", err_msg.c_str(), fout.c_str());
            perror("open");
        }
        writeout.close();
    }
    std::string Path::_header_str() const{
        bool verbose = false;
        std::string header
            //= length.to_string() + "\n"
            = initial().atom_state().header_str(verbose);
        return header;
    }
    void Path::_read_header(std::ifstream& readout){
        std::string line;
        /*state header string                                     */
        std::getline(readout, line);
        simple::BaseState::read_header(line);
    }
    std::list<Record> Path::get_path() const{
        std::list<Record> full_path(_path);
        full_path.push_back(final());
        full_path.push_front(initial());
        return full_path;
    }
} // namespace geodesic
