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
/**
* Base Class Observable
*/
Observable::Observable(std::string name,
    std::string fname,
    std::string units,
    std::string axis_name,
    double value) :
    _name(name),
    _fname(fname),
    _units(units),
    _axis_name(axis_name),
    _value(value){
    if (_axis_name == "") _axis_name = _name;
}
void Observable::prepare_ofstream(std::string datadir,
    std::string sim_name, std::ofstream& writeout, bool overwrite) const{
    std::string fout;
    fout = datadir
        + parse_string(sim_name)
        + "_"
        + parse_string(fname())
        + ".csv";
    // if overwrite is allowed, add header
    if (overwrite) {
        // if overwrite is allowed, open in truncate mode
        writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
    }
    if (!overwrite) {
        // if overwrite is forbidden, to open the file in append mode
        writeout.open(fout, std::ofstream::out | std::ofstream::app);
    }
    if (!writeout.is_open()){
        // if file still could not be opened
        std::string err_msg = "prepare_ofstream: unable to open file at";
        fprintf(stderr, "%s %s. %s \n",
            err_msg.c_str(), fout.c_str(), "file stream will be closed.");
        perror("open");
        writeout.close();
    }
    return;
}
std::string Observable::print_value() const {
    std::ostringstream out;
    out.precision(dbl::max_digits10);
    out << value();
    return out.str();
};
void Observable::writeout(std::ofstream& fs, bool add_header){
    std::string delim = ",";                /**> delimeter for CSV data */
    if (fs.is_open()) {
        if (add_header) {
            fs << "value" << std::endl;
        }
        // then output the data
        fs.precision(dbl::max_digits10);
        fs << value() << std::endl;
    }
    else {
        // if file still could not be opened
        std::string err_msg = "writeout: unable to open file";
        fprintf(stderr, "%s\n", err_msg.c_str());
        perror("open");
    }
    fs.close();
}
void Observable::update(std::string obs_string){
    std::istringstream ss(obs_string.c_str());
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> words(begin, end);
    // last space-separated substring element is the value of the observable
    // (in case observable name contains more than one word)
    update(atof(words.back().c_str()));
}
::std::ostream& operator<<(::std::ostream& os, const Observable& obs) {
    return os << obs.to_string().c_str();
};
/**
* Class AvgObservable with Virtual Inheritance of Observable
*/
AvgObservable::AvgObservable(Observable& inst_obs, double acc,
    std::string name,
    std::string fname,
    std::string units,
    std::string axis_name) :
    /**
    * this base class constructor is irrelevant and will not be called by
    * a concrete observable due to virtual inheritance
    */
    Observable(name, fname, units, axis_name),
    _inst_obs(inst_obs),
    _acc(acc){}
void AvgObservable::update(std::string obs_string){
    std::istringstream ss(obs_string.c_str());
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> words(begin, end);
    _acc = atof(words.at( words.size() - 2 ).c_str());
    Observable::update(atof(words.back().c_str()));
}
std::string AvgObservable::print_value() const {
    std::ostringstream out;
    out.precision(dbl::max_digits10);
    out << _acc;
    out << " ";
    out << value();
    return out.str();
};
/**
* Class TimeLogObservable with Virtual Inheritance of Observable
*/
TimeLogObservable::TimeLogObservable(std::string name, std::string fname,
    std::string units, std::string axis_name) :
    /**
    * this base class constructor is irrelevant and will not be called by
    * a concrete observable due to virtual inheritance
    */
    Observable (name, fname, units, axis_name) {}
void TimeLogObservable::writeout(std::ofstream& fs, bool add_header){
    std::string delim = ",";                /**> delimeter for CSV data */
    if (fs.is_open()) {
        if (add_header) {
            fs << "time" << delim;
            fs << "value" << std::endl;
        }
        // then output the data
        for (const Record& rec: _records) {
            fs.precision(dbl::max_digits10);
            fs << rec.time << delim;
            fs << rec.value << std::endl;
        }
    }
    else {
        // if file still could not be opened
        std::string err_msg = "writeout: unable to open file";
        fprintf(stderr, "%s\n", err_msg.c_str());
        perror("open");
    }
    fs.close();
    // now that everything is written, flush records to free memory and avoid
    // overwriting old data twice next time
    _records.clear();
}
