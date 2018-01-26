// 2017 Artur Avkhadiev
/*! \file observable.cpp
*/
#include <vector>
#include <fstream>
#include <sys/stat.h>                       /**> checks for file existence */
#include <limits>                           /**> controls output precision */
#include <string>
#include <../include/parsing.h>
#include "../include/vector.h"
#include "../include/observable.h"
typedef std::numeric_limits< double > dbl;
namespace observable{
    std::string to_string(const name_t &name, bool verbose){
        std::string s;
        if (verbose) {
            s = "NAME: " + name.full + " (" + name.abridged + ")\n";
            s += "AXIS NAME: " + name.latex + " (" + name.units + ")";
        }
        else {
            s = name.full + "\n"
                + name.abridged + "\n"
                + name.latex + "\n"
                + name.units;
        }
        return s;
    }
    ::std::ostream& operator<<(::std::ostream& os, const observable::name_t& name){
        return os << to_string(name).c_str();
    }
    std::string to_string(const update_time_t &update_time, bool verbose){
        std::string s = "";
        if (verbose){
            s += "Update Time Rule: ";
            switch (update_time) {
            case observable::FORCE_LOOP:
                s += "force loop";
                break;
            case observable::INTEGRATION_STEP:
                s += "integration step";
                break;
            case observable::MAIN_LOOP:
                s += "mainloop";
                break;
            default:
                s += "UNCRECOGNIZED";
            }
        }
        else {
            s += std::to_string((int)update_time);
        }
        return s;
    }
    ::std::ostream& operator<<(::std::ostream& os, const observable::update_time_t& update_time){
        return os << to_string(update_time).c_str();
    }
    std::string to_string(const observable::calculate_avg_t &estimate, bool verbose){
        std::string s = "";
        if (verbose) {
            std::string estimate_mean = "ESIMATE MEAN: ";
            std::string estimate_meansq = "ESTIMATE MEANSQ: ";
            std::string yes = "YES";
            std::string no = "NO";
            estimate.mean?
                estimate_mean += yes :
                estimate_mean += no;
            estimate.meansq?
                estimate_meansq += yes :
                estimate_meansq += no;
            s += estimate_mean + "\n" + estimate_meansq;
        }
        else {
            s += std::to_string((int)estimate.mean)
                + " " + std::to_string((int)estimate.meansq);
        }
        return s;
    }
    ::std::ostream& operator<<(::std::ostream& os, const observable::calculate_avg_t& estimate){
        return os << to_string(estimate);
    };
} // namespace observable
Observable::Observable(observable::name_t name,
    observable::update_time_t update_time,
    observable::calculate_avg_t calc_instructions,
    bool print_inst_val,   // print out instantaneous value in log & state?
    bool e_format          // use precision notation when writing out?
) :
    _name(name),
    _update_time(update_time),
    _calculate(calc_instructions),
    _print_inst_val(print_inst_val),
    _e_format(e_format),
    _msd_additive_const(0.0),
    _block_number(0),
    _is_block_averaged(false),
    value(0.0)
{
    // parse observable name
    if (_name.full == ""){
        _name.full = "Unknown Observable";
    }
    if (_name.abridged == "") {
        _name.abridged = _name.full;
    }
    _name.abridged = parse_string(_name.abridged);
    // zero run and block averages
    _zero_run();
    _zero_block();
}
std::string Observable::print_value() const {
    if (e_format()){
        std::ostringstream out;
        out.precision(dbl::max_digits10);
        out << value;
        return out.str();
    }
    else {
        return std::to_string(value);
    }
};
std::string Observable::status_string() const{
    std::string s = _name.abridged + " " + print_value();
    return s;
}
void Observable::read_status(std::string line) {
    std::istringstream ss(line.c_str());
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> words(begin, end);
    value = atof(words.at(1).c_str());
}
std::string Observable::to_string(bool verbose) const{
    std::string s = observable::to_string(_name, verbose) + "\n";
    s += observable::to_string(_update_time, verbose) + "\n";
    s += observable::to_string(_calculate, verbose) + "\n";
    if (verbose) s += "instantaneous value: ";
    s += print_value() + "\n";
    if (_calculate.mean || _calculate.meansq){
        if (_calculate.mean) {
            if (verbose) s += "run mean, mean_acc, norm: ";
            s += std::to_string(_run.mean) + " "
                + std::to_string(_run.mean_acc) + " "
                + std::to_string(_run.norm) + "\n";
            if (verbose) {
                s += "block (" + std::to_string(block_number()) + ") "
                    + "mean, mean_acc, norm: ";
            }
            else {
                s += std::to_string(block_number()) + " ";
            }
            s += std::to_string(_block.mean) + " "
                + std::to_string(_block.mean_acc) + " "
                + std::to_string(_block.norm) + "\n";
        }
        if (_calculate.meansq) {
            if (verbose) s += "run meansq, meansq_acc, norm: ";
            s += std::to_string (_run.meansq) + " "
                + std::to_string(_run.meansq_acc) + " "
                + std::to_string(_run.norm) + "\n";
            if (verbose) {
                s += "block (" + std::to_string(block_number()) + ") "
                    + "meansq, meansq_acc, norm: ";
            }
            else {
                s += std::to_string(block_number()) + " ";
            }
            s += std::to_string(_block.meansq) + " "
                + std::to_string(_block.meansq_acc) + " "
                + std::to_string(_block.norm) + "\n";
        }
        if (_calculate.mean || _calculate.meansq) {
            if (verbose) s += "additive const: ";
            s += std::to_string(_msd_additive_const) + "\n";
            if (verbose) s += "run, block fluctuation: ";
            s += std::to_string(_run.fluctuation)
                + " " + std::to_string(_block.fluctuation) + "\n";
        }
    }
    if (verbose) s += "print instantaneous value to log: ";
    s += std::to_string((int)print_inst_val()) + "\n";
    if (verbose) s += "use precision notation when printing: ";
    s += std::to_string((int)e_format());
    return s;
}
::std::ostream& operator<<(::std::ostream& os, const Observable& obs){
    return os << obs.to_string();
}
void Observable::_zero_block() {
    _block.mean = 0.0;
    _block.meansq = 0.0;
    _block.fluctuation = 0.0;
    _block.mean_acc = 0.0;
    _block.meansq_acc = 0.0;
    _block.norm = 0;
}
void Observable::start_new_block(){
    _zero_block();
    _is_block_averaged = false;
    ++_block_number;
}
void Observable::update_block(){
    if (_calculate.mean || _calculate.meansq){
        ++_block.norm;
        if (_calculate.mean) _block.mean_acc += value;
        if (_calculate.meansq) _block.meansq_acc += pow(value, 2.0);
    }
}
void Observable::_average_block(){
    _is_block_averaged = true;
    if (_calculate.mean) {
        _block.mean = _block.mean_acc / (double(_block.norm));
    }
    if (_calculate.meansq) {
        _block.meansq = _block.meansq_acc / (double(_block.norm));
    }
    if (_calculate.mean && _calculate.meansq){
        // s^2(A) = 1/(norm-1) \sum (A - <A>block)^2
        _block.fluctuation
            = _block.meansq - pow(_block.mean, 2.0);
        // to calculate norm/(norm - 1), make sure norm > 1
        if (_block.norm < 2){
            std::string err_msg
                = "Observable: number of points in block is less than 2; can't get an unbiased estimate of variance (_block.norm = "
                + std::to_string(_block.norm)
                + ")\n";
            throw std::invalid_argument(err_msg);
        }
        double normalize = ((double)_block.norm)/((double)(_block.norm - 1.0));
        _block.fluctuation
            = _msd_additive_const + normalize * (_block.fluctuation);
        }
}
void Observable::_zero_run(){
    _run.mean = 0.0;
    _run.meansq = 0.0;
    _run.fluctuation = 0.0;
    _run.mean_acc = 0.0;
    _run.meansq_acc = 0.0;
    _run.norm = 0;
    // ultimate error in the mean (square root of (variance/Nblocks))
    _run_err = 0.0;
}
void Observable::start_new_run(){
    _zero_run();
    _zero_block();
    _block_number = 0;
}
void Observable::_update_run(){
    if (_calculate.mean || _calculate.meansq){
        ++_run.norm;
        if (_calculate.mean) _run.mean_acc += _block.mean;
        // <A^2>_run = 1/nblock \sum <A>_b^2
        if (_calculate.meansq) _run.meansq_acc += pow(_block.mean, 2.0);
    }
}
void Observable::average_block_update_run(){
    // check block for accumulation
    if (_block.norm < 1){
        std::string err_msg
            = "Observable: can't average with block norm = "
            + std::to_string(_block.norm)
            + "\n";
        throw std::invalid_argument(err_msg);
    }
    if (!is_block_averaged()){
        _average_block();               // calculate averages
        _update_run();                  // use averages to update run acczero_
    }
    else {
        std::string err_msg
            = "Observable: attempt to average same block twice\n";
        throw std::invalid_argument(err_msg);
    }
}
void Observable::average_run(){
    // check block for accumulation
    if (_run.norm < 1){
        std::string err_msg
            = "Observable: can't average with block norm = "
            + std::to_string(_block.norm)
            + "\n";
        throw std::invalid_argument(err_msg);
    }
    if (_calculate.mean) {
        _run.mean = _run.mean_acc / (double(_run.norm));
    }
    if (_calculate.meansq) {
        _run.meansq = _run.meansq_acc / (double(_run.norm));
    }
    if (_calculate.mean && _calculate.meansq){
        // s^2(<A>block) = 1/(NBLOCKS-1) \sum (<A>block - <A>run)^2
        _run.fluctuation
            = _run.meansq - pow(_run.mean, 2.0);
        // to calculate norm/(norm - 1), make sure norm > 1
        if (_run.norm < 2){
            std::string err_msg
                = "Observable: number of blocks used to calculate run average is less than 2; can't get an unbiased estimate of variance (_run.norm = "
                + std::to_string(_run.norm)
                + ")\n";
            throw std::invalid_argument(err_msg);
        }
        double normalize = ((double)_run.norm)/((double)(_run.norm - 1.0));
        _run.fluctuation
            = _msd_additive_const + normalize * (_run.fluctuation);
        if (_run.fluctuation >= -std::numeric_limits<double>::epsilon()) {
            // normalize and get estimated errors
            // err = sqrt(s^2(<A>block)/Nblocks)
            _run_err = sqrt(_run.fluctuation / ((double)_run.norm));
        }
        else {
            fprintf(stderr, "%s: %s (%f), %s\n",
                _name.full.c_str(),
                "Estimate of run variance is < 0",
                _run.fluctuation,
                "can't calculate run error");
        }
    }
}
// getters
size_t Observable::block_number() const {
    return _block_number;
}
double Observable::block_mean() const{
    if (!_calculate.mean) {
        fprintf(stderr, "%s %s: %s\n", "Observable",
            short_name().c_str(),
            "block mean requested but mean is not tracked");
    }
    return _block.mean;
}
double Observable::block_meansq() const{
    if (!_calculate.meansq) {
        fprintf(stderr, "%s %s: %s\n", "Observable",
            short_name().c_str(),
            "block meansq requested but meansq is not tracked");
    }
    return _block.meansq;
}
double Observable::block_fluctuation() const {
    if (!(_calculate.mean && _calculate.meansq)) {
        fprintf(stderr, "%s %s: %s\n", "Observable",
            short_name().c_str(),
            "block fluctuation requested but either mean or meansq are not tracked");
    }
    return _block.fluctuation;
}
double Observable::run_mean() const {
    if (!_calculate.mean) {
        fprintf(stderr, "%s %s: %s\n", "Observable",
            short_name().c_str(),
            "run mean requested but mean is not tracked");
    }
    return _run.mean;
}
double Observable::run_meansq() const{
    if (!_calculate.meansq) {
        fprintf(stderr, "%s %s: %s\n", "Observable",
            short_name().c_str(),
            "run meansq requested but meansq is not tracked");
    }
    return _run.meansq;
}
double Observable::run_fluctuation() const {
    if (!(_calculate.mean && _calculate.meansq)) {
        fprintf(stderr, "%s %s: %s\n", "Observable",
            short_name().c_str(),
            "run fluctuation requested but either mean or meansq are not tracked");
    }
    return _run.fluctuation;
}
double Observable::run_err() const {
    if (!(_calculate.mean && _calculate.meansq)) {
        fprintf(stderr, "%s %s: %s\n", "Observable",
            short_name().c_str(),
            "run error requested but either mean or meansq are not tracked");
    }
    return _run_err;
}
//void Observable::prepare_ofstream(std::string datadir,
//    std::string sim_name, std::ofstream& writeout, bool overwrite) const{
//    std::string fout;
//    fout = datadir
//        + parse_string(sim_name)
//        + "_"
//        + parse_string(fname())
//        + ".csv";
//    // if overwrite is allowed, add header
//    if (overwrite) {
//        // if overwrite is allowed, open in truncate mode
//        writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
//    }
//    if (!overwrite) {
//        // if overwrite is forbidden, to open the file in append mode
//        writeout.open(fout, std::ofstream::out | std::ofstream::app);
//    }
//    if (!writeout.is_open()){
//        // if file still could not be opened
//        std::string err_msg = "prepare_ofstream: unable to open file at";
//        fprintf(stderr, "%s %s. %s \n",
//            err_msg.c_str(), fout.c_str(), "file stream will be closed.");
//        perror("open");
//        writeout.close();
//    }
//    return;
//}
//void Observable::writeout(std::ofstream& fs, bool add_header){
//    std::string delim = ",";                /**> delimeter for CSV data */
//    if (fs.is_open()) {
//        if (add_header) {
//            fs << "value" << std::endl;
//        }
//        // then output the data
//        fs.precision(dbl::max_digits10);
//        fs << value() << std::endl;
//    }
//    else {
//        // if file still could not be opened
//        std::string err_msg = "writeout: unable to open file";
//        fprintf(stderr, "%s\n", err_msg.c_str());
//        perror("open");
//    }
//    fs.close();
//}
//void Observable::update(std::string obs_string){
//    std::istringstream ss(obs_string.c_str());
//    std::istream_iterator<std::string> begin(ss);
//    std::istream_iterator<std::string> end;
//    std::vector<std::string> words(begin, end);
//    // last space-separated substring element is the value of the observable
//    // (in case observable name contains more than one word)
//    update(atof(words.back().c_str()));
//}
