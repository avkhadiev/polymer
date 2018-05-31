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
    //Path::Length::Length(){}
    //    _calculate_diag_length(true),
    //    _I0_over_M(pow(simple::BasePolymer::d(), 2.0)) {}
    //Path::Length::Length(bool calculate_diag_length) :
    //    Path::Length::Length()
    //{
    //    _calculate_diag_length = calculate_diag_length;
    //}
    //Path::Length::Length(double g_pol_trans,
    //                    double g_pol_rot,
    //                    double g_pol_rot_diag,
    //                    double g_sol,
    //                    bool calculate_diag_length) :
    //    Path::Length::Length(calculate_diag_length)
    //{
    //    _calculate_diag_length = calculate_diag_length;
    //    _g_pol_trans = g_pol_trans;
    //    _g_pol_rot_full = g_pol_rot;
    //    _g_pol_rot_diag = g_pol_rot_diag;
    //    _g_sol = g_sol;
    //}
    //Path::Length::~Length(){};
    //void Path::Length::_zero_sq_increments(){
    //    _g_pol_trans_sq    = 0.0;
    //    _g_pol_rot_full_sq = 0.0;
    //    _g_pol_rot_diag_sq = 0.0;
    //    _g_sol_sq          = 0.0;
    //}
    //void Path::Length::_zero(){
    //    _ell_full       = 0.0;
    //    _ell_diag       = 0.0;
    //    _g_pol_diag     = 0.0;
    //    _g_pol_full     = 0.0;
    //    _g_pol_trans    = 0.0;
    //    _g_pol_rot_full = 0.0;
    //    _g_pol_rot_diag = 0.0;
    //    _g_sol          = 0.0;
    //    _zero_sq_increments();
    //}
    //std::string Path::Length::to_string() const {
    //    std::string str
    //        = std::to_string(ell())
    //        + " "
    //        + std::to_string(ell_diag())
    //        + " "
    //        + std::to_string(g_polymer_trans())
    //        + " "
    //        + std::to_string(g_polymer_rot())
    //        + " "
    //        + std::to_string(g_polymer_rot_diag())
    //        + " "
    //        + std::to_string(g_solvent());
    //    return str;
    //}
    //void Path::Length::_update_from_string(std::string line){
    //    std::istringstream ss(line.c_str());
    //    std::istream_iterator<std::string> begin(ss);
    //    std::istream_iterator<std::string> end;
    //    std::vector<std::string> words(begin, end);
    //    _g_pol_trans = atof(words.at(0).c_str());
    //    _g_pol_rot_full = atof(words.at(1).c_str());
    //    _g_pol_rot_diag = atof(words.at(2).c_str());
    //    _g_sol = atof(words.at(3).c_str());
    //}
    //void Path::Length::_increment(  simple::BondState& last_state,
    //                                simple::BondState& next_state){
    //    _zero_sq_increments();
    //    _calculate_square_increments(last_state, next_state);
    //    // update components
    //    _g_pol_trans += sqrt(_g_pol_trans_sq);
    //    // fprintf(stderr, "g_pol_rot_full_sq = %f\n", _g_pol_rot_full_sq);
    //    // "increment" for the rotational component can be negative
    //    // due to the cross-terms
    //    if (_g_pol_rot_full_sq >= 0){
    //        _g_pol_rot_full += sqrt(_g_pol_rot_full_sq);
    //        _g_pol_full += sqrt(_g_pol_trans_sq + _g_pol_rot_full_sq);
    //    }
    //    else{
    //        fprintf(stderr, "%s\n", "Ouch component increment is negative");
    //        _g_pol_rot_full -= sqrt(-_g_pol_rot_full_sq);
    //        if (_g_pol_trans_sq + _g_pol_rot_full_sq >= 0)
    //            _g_pol_full += sqrt(_g_pol_trans_sq + _g_pol_rot_full_sq);
    //        else {
    //            _g_pol_full -= sqrt(- _g_pol_trans_sq - _g_pol_rot_full_sq);
    //        }
    //    }
    //    _g_sol += sqrt(_g_sol_sq);
    //    double ell_full = _g_pol_trans_sq + _g_pol_rot_full_sq + _g_sol_sq;
    //    if (ell_full >= 0){
    //        _ell_full += sqrt(ell_full);
    //    }
    //    else {
    //        _ell_full -= sqrt(- ell_full);
    //    }
    //    if (_calculate_diag_length) {
    //        _g_pol_rot_diag += sqrt(_g_pol_rot_diag_sq);
    //        _g_pol_diag += sqrt(_g_pol_trans_sq + _g_pol_rot_diag_sq);
    //        _ell_diag += sqrt(_g_pol_trans_sq + _g_pol_rot_diag_sq + _g_sol_sq);
    //    }
    //}
    //void Path::Length::_calculate_square_increments(
    //    simple::BondState& last_state,
    //    simple::BondState& next_state){
    //    _calculate_g_pol_trans_sq_increment(last_state, next_state);
    //    _calculate_g_pol_rot_sq_increment(last_state, next_state);
    //    _calculate_g_sol_sq_increment(last_state, next_state);
    //}
    //void Path::Length::_calculate_g_pol_trans_sq_increment(
    //    simple::BondState& last_state,
    //    simple::BondState& next_state){
    //    int nm = simple::BaseState::nm();
    //    Vector r_last, r_next;
    //    _g_pol_trans_sq = 0.0;
    //    for(int ip = 0; ip < nm; ++ip){
    //        r_last = last_state.polymers.at(ip).rcm();
    //        r_next = next_state.polymers.at(ip).rcm();
    //        _g_pol_trans_sq += _trans_sq_increment(r_last, r_next);
    //    }
    //}
    //void Path::Length::_calculate_g_pol_rot_sq_increment(
    //    simple::BondState& last_state,
    //    simple::BondState& next_state){
    //    int nm = simple::BaseState::nm();
    //    int nb = simple::BondPolymer::nb();
    //    simple::BondPolymer last_polymer, next_polymer;
    //    Vector omega_last_i, omega_next_i;
    //    Vector omega_last_j, omega_next_j;
    //    _g_pol_rot_full_sq = 0.0;
    //    if (_calculate_diag_length) _g_pol_rot_diag_sq = 0.0;
    //    double square_increment;
    //    double Dij;
    //    for(int ip = 0; ip < nm; ++ip){             // loop over polymers
    //        last_polymer = last_state.polymers.at(ip);
    //        next_polymer = next_state.polymers.at(ip);
    //        for(int ib = 0; ib < nb; ++ib){         // double loop over bonds
    //            omega_last_i = last_polymer.bonds.at(ib).position;
    //            omega_next_i = next_polymer.bonds.at(ib).position;
    //            for(int jb = ib; jb < nb; ++jb){
    //                omega_lafadsfst_j = last_polymer.bonds.at(jb).position;
    //                omega_next_j = next_polymer.bonds.at(jb).position;
    //                Dij = simple::BaseState::Dij(ib + 1, jb + 1);
    //                square_increment = Dij * _rot_sq_increment(
    //                        omega_last_i, omega_next_i,
    //                        omega_last_j, omega_next_j);
    //                if (ib != jb) {
    //                    _g_pol_rot_full_sq += 2 * square_increment;
    //                }
    //                else{
    //                    _g_pol_rot_full_sq += square_increment;
    //                    if (_calculate_diag_length){
    //                        _g_pol_rot_diag_sq += square_increment;
    //                    }
    //                }
    //            }
    //        }
    //    }
    //    _g_pol_rot_full_sq = _g_pol_rot_full_sq * _I0_over_M;
    //    if (_calculate_diag_length){
    //        _g_pol_rot_diag_sq = _g_pol_rot_diag_sq * _I0_over_M;
    //    }
    //}
    //void Path::Length::_calculate_g_sol_sq_increment(
    //    simple::BondState& last_state,
    //    simple::BondState& next_state){
    //    int ns = simple::BaseState::nsolvents();
    //    Vector r_last, r_next;
    //    _g_sol_sq = 0.0;
    //    for(int is = 0; is < ns; ++is){
    //        r_last = last_state.solvents.at(is).r();
    //        r_next = next_state.solvents.at(is).r();
    //        _g_sol_sq += _trans_sq_increment(r_last, r_next);
    //    }
    //}
    ///** returns (r_next - r_last)^2 */
    //double Path::Length::_trans_sq_increment(
    //    const Vector &r_last, const Vector &r_next){
    //    return normsq(subtract(r_next, r_last));
    //}
    ///**
    //* returns D_ij del Psi_i del Psi_j p_ij, where
    //*   del Psi = arccos(omega_next * omega_last)
    //*   p_ij = del_Omega_i * del_Omega_j, where
    //*       del Omega = omega_next - omega_last
    //*/
    //double Path::Length::_rot_sq_increment(
    //    const Vector& omega_last_i, const Vector& omega_next_i,
    //    const Vector& omega_last_j, const Vector& omega_next_j){
    //    double del_Psi_i = acos(dot(omega_last_i, omega_next_i));
    //    double del_Psi_j = acos(dot(omega_last_j, omega_next_j));
    //    Vector del_omega_i = subtract(omega_next_i, omega_last_i);
    //    Vector del_omega_j = subtract(omega_next_j, omega_last_j);
    //    double p_ij
    //        = dot(del_omega_i, del_omega_j)
    //        / (norm(del_omega_i) * norm(del_omega_j));
    //    return del_Psi_i * del_Psi_j * p_ij;
    //}
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
            = initial().state().header_str(verbose);
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
