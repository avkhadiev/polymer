// 2018 Artur Avkhadiev
/*! \file geodesic_path.cpp
*/
#include "../include/geodesic_path.h"
namespace geodesic{
    /***************************************************************************
    *                               GEODESIC PATH
    ***************************************************************************/
    Path::Path() :
        _initial(),
        _final(),
        _ell(0.0),
        _path(){}
    Path::Path(Record initial, Record final) :
        _initial(initial),
        _final(final),
        _ell(0.0),
        _path()
    {
        _path.clear();
    }
    Path::Path(std::list<Record>& records) :
        _initial(records.front()),
        _final(records.back()),
        _ell(0.0),
        _path(records){
            // omit boundary values in path record
            _path.pop_front();
            _path.pop_back();
            // update path length
            if (_check_time_flow()){
                recompute_length();
            }
            else {
                std::string err_msg = "Path: time direction not unique";
                throw std::invalid_argument(err_msg);
            }
        }
    Path::~Path() {
        _path.clear();
    }
    Path::Path(std::string fin) :
        _initial(),
        _final(),
        _ell(0.0),
        _path()
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
                // update path length
                if (_check_time_flow()){
                    recompute_length();
                }
                else {
                    std::string err_msg = "Path: time direction not unique";
                    throw std::invalid_argument(err_msg);
                }
            }
            else {
                std::string err_msg = "Path: less than two records provided";
                throw std::invalid_argument(err_msg);
            }
        }
    }
    size_t Path::nrecords() const {
        return _path.size() + 2;           /**> include boundary values */
    }
    // A simplistic implementation of operator= (see better implementation below)
    Path& Path::operator= (const Path &other)
    {
        _initial = other.initial();
        _final = other.final();
        _ell = other.kinematic_length();
        _path = other._path;
        return *this;
    }
    const Record Path::current_tail() const {
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
    double Path::time() const {
        // total travel time (non-negative quantity)
        return abs(final().state().time() - initial().state().time());
    }
    bool Path::_check_time_flow(const Record& preceding_record, const Record& next_record) const{
        bool is_time_flow_ok = true;
        double preceding_record_t = preceding_record.state().time();
        double next_record_t = next_record.state().time();
        double initial_t = initial().state().time();
        double final_t = final().state().time();
        double time_diff = final_t - initial_t;
        if (time_diff > 0){
            /* if t_f - t_i > 0, then current_tail.t < new_record.t < final.t */
            time_diff
                = (preceding_record_t < next_record_t)
                && (next_record_t < final_t);
        }
        else if (time_diff < 0){
            /* if t_f - t_i < 0, then current_tail.t > new_record.t > final.t */
            time_diff
                = (preceding_record_t > next_record_t)
                && (next_record_t > final_t);
        }
        else {
            /* if t_f - t_i ==0, there can be no intermediate values */
            is_time_flow_ok = false;
        }
        return is_time_flow_ok;
    }
    bool Path::_check_time_flow() const{
        bool is_time_flow_ok = true;
        if (!_path.empty()){
            Record preceding_record = initial();
            std::list<Record>::const_iterator rec;
            bool update;
            for (rec = _path.begin(); rec != _path.end(); std::advance(rec, 1)){
                update =  _check_time_flow(preceding_record, *rec);
                is_time_flow_ok = is_time_flow_ok && update;
                if (is_time_flow_ok){
                    continue;
                }
                else {
                    break;
                }
            }
        }
        return is_time_flow_ok;
    }
    void Path::_increment_length( const Record& preceding_record, const Record& next_record){
        // FIXME
        fprintf(stderr, "%s\n", "geodesic::Path::_increment_length:: is to be implemented in daugther classes; this function is only enabled for testing");
    }
    void Path::recompute_length(){
        _ell = 0.0;
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
        if (_check_time_flow(current_tail(), new_record)){
            _increment_length(current_tail(), new_record);
            _path.push_back(new_record);
        }
        else {
            std::string err_msg = "Path::append: wrong time direction";
            throw std::invalid_argument(err_msg);
        }
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
            writeout << std::endl;
            for (rec = _path.begin(); rec != _path.end(); std::advance(rec, 1)){
                (*rec).write(writeout, overwrite);
                writeout << std::endl;
            }
            final().write(writeout, overwrite);
            writeout << std::endl;
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
        std::string header = std::to_string(kinematic_length()) + " "
            + std::to_string(nrecords()) + " "
            + std::to_string(time()) + "\n"
            + initial().state().header_str(verbose);
        return header;
    }
    void Path::_read_header(std::ifstream& readout){
        std::string line;
        /* (1st line) kinematic_length nrecords time */
        std::getline(readout, line);
        std::istringstream ss(line.c_str());
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> words(begin, end);
        _ell = atof(words.at(0).c_str());
        /* (2nd line) state header string            */
        std::getline(readout, line);
        simple::BaseState::read_header(line);
    }
} // namespace geodesic
