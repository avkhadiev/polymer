// 2017 Artur Avkhadiev
/*! \file observables.cpp
*/
#include <map>
#include <utility>      /* std::pair, std::make_pair */
#include <vector>
#include <fstream>
#include <iomanip>      /* std::setw */
#include <string>
#include <../include/parsing.h>
#include "../include/vector.h"
#include "../include/observables.h"
std::string scalar_observable_to_string(ScalarObservable so) {
    std::string so_str;
    std::string so_axis_name = so.axis_name;
    std::string so_units = so.units;
    so_str = so_axis_name + ", " + so.units;
    return so_str;
}
std::string vector_observable_to_string(VectorObservable vo) {
    std::string vo_str;
    std::string vo_axis_name = vo.axis_name;
    std::string vo_units = vo.units;
    vo_str = vo_axis_name + ", " + vo.units;
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
ScalarObservable declare_scalar_observable(std::string name,
    std::string units,
    std::string axis_name) {
    std::vector<std::pair<double, double> > empty_vector;
    if (axis_name == "") {
        axis_name = name;
    }
    ScalarObservable so = {.value_time = empty_vector,
        .accumulator = 0.0,
        .name = name,
        .axis_name = axis_name,
        .units = units};
    return so;
}
VectorObservable declare_vector_observable(std::string name,
    std::string units,
    std::string axis_name) {
    std::vector<std::pair<Vector, double> > empty_vector;
    if (axis_name == "") {
        axis_name = name;
    }
    VectorObservable vo = {.value_time = empty_vector,
        .accumulator = vector(0.0, 0.0, 0.0),
        .name = name,
        .axis_name = axis_name,
        .units = units};
    return vo;
}
void clear_observable_records(ScalarObservable& so){
    so.value_time.clear();
};
void clear_observable_records(VectorObservable& vo){
    vo.value_time.clear();
};
void write_vector_observable_to_file(VectorObservable& vo,
    std::string outdir, std::string sim_name, bool overwrite){
    // add header files?
    bool add_header;
    // CSV file delimeter
    std::string delim = ",";
    if (overwrite) {
        add_header = true;
    }
    else {
        add_header = false;
    }
    // stores path to output file
    std::string fout;
    fout = outdir
        + parse_string(sim_name)
        + "_"
        + parse_string(vo.name)
        + ".csv";
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
            std::string s = parse_string(vo.name);
            writeout << "time" << delim;
            writeout << s.c_str() << "_x" << delim;
            writeout << s.c_str() << "_y" << delim;
            writeout << s.c_str() << "_z" << delim;
            writeout << "axis_name" << delim;
            writeout << "units" << std::endl;
        }
        // then output the data
        std::pair<Vector, double> value_time;
        for (int i = 0; i < vo.value_time.size(); ++i) {
            value_time = vo.value_time.at(i);
            writeout << value_time.second << delim;
            writeout << value_time.first.x << delim;
            writeout << value_time.first.y << delim;
            writeout << value_time.first.z << delim;
            writeout << vo.axis_name << delim;
            writeout << vo.units << std::endl;
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
    // CSV file delimeter
    std::string delim = ",";
    if (overwrite) {
        add_header = true;
    }
    else {
        add_header = false;
    }
    // integer for field width
    // int indent = 15;
    // stores path to output file
    std::string fout;
    fout = outdir
        + parse_string(sim_name)
        + "_"
        + parse_string(so.name)
        + ".csv";
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
            writeout << "time" << delim;
            writeout << parse_string(so.name) << delim;
            writeout << "axis_name" << delim;
            writeout << "units" << std::endl;
        }
        // then output the data
        std::pair<double, double> value_time;
        for (int i = 0; i < so.value_time.size(); ++i) {
            value_time = so.value_time.at(i);
            writeout << value_time.second << delim;
            writeout << value_time.first << delim;
            writeout << so.axis_name << delim;
            writeout << so.units << std::endl;
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
