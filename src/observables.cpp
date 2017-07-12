// 2017 Artur Avkhadiev
/*! \file observables.cpp
*/
#include <map>
#include <utility>      /* std::pair, std::make_pair */
#include <vector>
#include <fstream>
#include <iomanip>      /* std::setw */
#include <string>
#include "../include/vector.h"
#include "../include/observables.h"
std::string scalar_observable_to_string(ScalarObservable so) {
    std::string so_str;
    std::string so_name = so.name;
    std::string so_units = so.units;
    so_str = so.name + ", " + so.units;
    return so_str;
}
std::string vector_observable_to_string(VectorObservable vo) {
    std::string vo_str;
    std::string vo_name = vo.name;
    std::string vo_units = vo.units;
    vo_str = vo.name + ", " + vo.units;
    return vo_str;
}
::std::ostream& operator<<(::std::ostream& os, const ScalarObservable& so) {
    return os << scalar_observable_to_string(so).c_str();
};
::std::ostream& operator<<(::std::ostream& os, const VectorObservable& vo) {
    return os << vector_observable_to_string(vo).c_str();
};
/**
* Given the name and the units of an observable, initialize the necessary
* struct; the value_time vector of value, time pairs will be empty;
*/
ScalarObservable declare_scalar_observable(std::string name, std::string units) {
    std::vector<std::pair<double, double> > empty_vector;
    ScalarObservable so = {.value_time = empty_vector,
        .accumulator = 0.0,
        .name = name,
        .units = units};
    return so;
}
VectorObservable declare_vector_observable(std::string name, std::string units) {
    std::vector<std::pair<Vector, double> > empty_vector;
    VectorObservable vo = {.value_time = empty_vector,
        .accumulator = vector(0.0, 0.0, 0.0),
        .name = name,
        .units = units};
    return vo;
}
void write_vector_observable_to_file(VectorObservable& vo,
    std::string outdir, std::string sim_name, bool overwrite){
    // add header files?
    bool add_header;
    if (overwrite) {
        add_header = true;
    }
    else {
        add_header = false;
    }
    // integer for field width
    int indent = 15;
    // stores path to output file
    std::string fout;
    fout = outdir + sim_name + "_" + vo.name + ".dat";
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
        // if file does not exist, open in truncate mode; allow to add header
            add_header = true;
            writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
        }
    }
    // now file has to be opened
    if (writeout.is_open()) {
        // if file could be opened...
        // first write the header file if necessary
        if (add_header) {
            // FIRST LINE
            writeout << std::left << std::setw(indent) << "# Time";
            writeout << std::right << std::setw(indent) << vector_observable_to_string(vo).c_str() << std::endl;
            // SECOND LINE
            writeout << std::left << std::setw(indent) << "#";
            writeout << std::left << std::setw(indent) << "x";
            writeout << std::left << std::setw(indent) << "y";
            writeout << std::left << std::setw(indent) << "z" << std::endl;
        }
        // then output the data
        std::pair<Vector, double> value_time;
        for (int i = 0; i < vo.value_time.size(); ++i) {
            value_time = vo.value_time.at(i);
            writeout << std::left << std::setw(indent) << value_time.second;
            writeout << std::left << std::setw(indent) << value_time.first.x;
            writeout << std::left << std::setw(indent) << value_time.first.y;
            writeout << std::left << std::setw(indent) << value_time.first.z << std::endl;
        }
    }
    else {
        // if file still could not be opened
        std::string err_msg = "write_scalar_observable_to_file: unable to open file at";
        fprintf(stderr, "%s %s\n", err_msg.c_str(), fout.c_str());
        perror("open");
    }
    writeout.close();
}
void write_scalar_observable_to_file(ScalarObservable& so,
    std::string outdir, std::string sim_name, bool overwrite){
    // add header files?
    bool add_header;
    if (overwrite) {
        add_header = true;
    }
    else {
        add_header = false;
    }
    // integer for field width
    int indent = 15;
    // stores path to output file
    std::string fout;
    fout = outdir + sim_name + "_" + so.name + ".dat";
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
        // if file does not exist, open in truncate mode; allow to add header
            add_header = true;
            writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
        }
    }
    // now file has to be opened
    if (writeout.is_open()) {
        // if file could be opened...
        // first write the header file if necessary
        if (add_header) {
            // FIRST LINE
            writeout << std::left << std::setw(indent) << "# Time";
            writeout << std::left << std::setw(indent) << scalar_observable_to_string(so).c_str() << std::endl;
        }
        // then output the data
        std::pair<double, double> value_time;
        for (int i = 0; i < so.value_time.size(); ++i) {
            value_time = so.value_time.at(i);
            writeout << std::left << std::setw(indent) << value_time.second;
            writeout << std::left << std::setw(indent) << value_time.first << std::endl;
        }
    }
    else {
        // if file still could not be opened
        std::string err_msg = "write_scalar_observable_to_file: unable to open file at";
        fprintf(stderr, "%s %s\n", err_msg.c_str(), fout.c_str());
        perror("open");
    }
    writeout.close();
}
