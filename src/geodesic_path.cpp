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
        length(true, false),
        euc_sep(false, false){}
    Path::Path(Record initial, Record final) :
        _initial(initial),
        _final(final),
        _path(),
        // print_inst_value, e_format
        length(true, false),
        euc_sep(false, false)
    {
        _path.clear();
        //this is a ``bare'', empty path with just two endpoints ---
        // the length is to be computed; so far, it is zero
        //recompute_set_theta();
        length.value = 0.0;
        euc_sep.value = 0.0;
        euc_sep.update(_initial, _final);
    }
    Path::Path(std::list<Record>& records) :
        _initial(records.front()),
        _final(records.back()),
        _path(records),
        // print_inst_value, e_format
        length(true, false),
        euc_sep(false, false){
            // omit boundary values in path record
            _path.pop_front();
            _path.pop_back();
            recompute_length();
            //fprintf(stderr, "%s\n", _header_str().c_str());
        }
    Path::~Path() {
        _path.clear();
    }
    Path::Path(std::string fin) :
        _initial(),
        _final(),
        _path(),
        // print_inst_value, e_format
        length(true, false),
        euc_sep(false, false)
    {
        std::ifstream readout;
        readout.open(fin.c_str(), std::ifstream::in);
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
    Path& Path::operator= (const Path &other)
    {
        _initial = other.initial();
        _final = other.final();
        _path = other._path;
        euc_sep.value = other.euc_sep.value;
        length.value = other.length.value;
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
        length.value = 0.0;
        euc_sep.value = 0.0;
        euc_sep.update(initial(), final());
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
            writeout.open(fout.c_str(), std::ofstream::out | std::ofstream::trunc);
        }
        if (!overwrite) {
            // if overwrite is forbidden, try to open the file in append mode
            writeout.open(fout, std::ofstream::out | std::ofstream::app);
            if (!writeout.is_open()) {
                // if file does not exist, open in truncate mode
                writeout.open(fout.c_str(), std::ofstream::out | std::ofstream::trunc);
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
    Path Path::get_slice(size_t n_slices, size_t mth_slice){
        if (n_slices > nrecords() - 1){
            fprintf(stderr, "%s\n", "Path::get_slice number of slices should be no greater than one less the number of records. Returning an empty path");
            return Path();
        }
        if ((mth_slice > n_slices) || (mth_slice < 1)) {
            fprintf(stderr, "%s\n", "Path::get_slice number of the slice should be a positive number no greater than the number of slices. Returning an empty path");
            return Path();
        }
        size_t overlap = n_slices - 1;
        double check = (double)(nrecords() + overlap);
        size_t nrecords_per_slice
            = (size_t)(std::ceil(check / ((double)n_slices)));
        fprintf(stderr, "%s %zu %s %zu %s %zu\n",
            "nrecords", nrecords(), "overlap", overlap,
            "nrecords_per_slice", nrecords_per_slice);
        // these indices index the path excluding the endpoints
        // if the entire path, including endpoints, is 0-index, then
        // the range has to be from 1 to nrecords - 2;
        size_t start_index = (mth_slice - 1) * (nrecords_per_slice - 1);
        size_t finish_index = start_index + nrecords_per_slice - 1;
        // if this is the first or the last slice,
        // will have to start later or end earlier
        if (start_index == 0) ++start_index;
        if (finish_index >= nrecords() - 1){
            finish_index = nrecords() - 2;
        }
        // find initial and final iterators for the list of records
        fprintf(stderr, "start %zu finish %zu\n", start_index, finish_index);
        std::list<Record>::const_iterator begin;
        std::list<Record>::const_iterator end;
        size_t running_index = 1;
        for (std::list<Record>::const_iterator rec = _path.begin();
             rec != _path.end();
             std::advance(rec, 1)){
            if (running_index == start_index){
                begin = rec;
            }
            else if (running_index == finish_index){
                end = rec;
                fprintf(stderr, "%s\n", "got final index");
                break;
            }
            ++running_index;
        }
        fprintf(stderr, "%s %zu\n", "running index", running_index);
        // use range constructor to make a new_path,
        // add the endpoints if necessary
        fprintf(stderr, "%s\n", "constructing new record");
        std::list<Record> new_path(begin, end);
        if (start_index == 1) new_path.push_front(initial());
        if (finish_index == nrecords() - 2) new_path.push_back(final());
        fprintf(stderr, "%s\n", "making a path");
        Path sliced_path = Path(new_path);
        fprintf(stderr, "%s with %zu records\n", "Returning sliced path with",
            sliced_path.nrecords());
        return sliced_path;
    }
    std::string Path::_header_str() const{
        bool verbose = false;
        std::string header
            = euc_sep.print_value() + "\n"
            + initial().atom_state().header_str(verbose);
        return header;
    }
    void Path::_read_header(std::ifstream& readout){
        std::string line;
        std::getline(readout, line);
        std::istringstream ss(line.c_str());
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> args(begin, end);
        euc_sep.value = atof(args.at(0).c_str());
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
