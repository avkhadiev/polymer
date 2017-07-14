// 2017 Artur Avkhadiev
/*! \file observables.h
*/
#ifndef POLYMER_OBSERVABLES_H
#define POLYMER_OBSERVABLES_H
#include <string>
#include <utility>      /* std::pair, std::make_pair */
#include <vector>
#include "vector.h"
struct scalar_observable_t {
    std::vector<std::pair<double, double> > value_time;
    double accumulator;
    std::string name;
    std::string axis_name;
    std::string units;
    scalar_observable_t& operator=(const scalar_observable_t& so)
    {
        value_time = so.value_time;
        accumulator = so.accumulator;
        name = so.name;
        axis_name = so.axis_name;
        units = so.units;
        return *this;
    }
};
typedef struct scalar_observable_t ScalarObservable;
struct vector_observable_t {
    std::vector<std::pair<Vector, double> > value_time;
    Vector accumulator;
    std::string name;
    std::string axis_name;
    std::string units;
    vector_observable_t& operator=(const vector_observable_t& vo)
    {
        value_time = vo.value_time;
        accumulator = vo.accumulator;
        name = vo.name;
        axis_name = vo.axis_name;
        units = vo.units;
        return *this;
    }
};
typedef struct vector_observable_t VectorObservable;
/**
* Given the observable, return a string with "Axis name, Units"
*/
std::string scalar_observable_to_string(ScalarObservable so);
std::string vector_observable_to_string(VectorObservable vo);
::std::ostream& operator<<(::std::ostream& os, const ScalarObservable& so);
::std::ostream& operator<<(::std::ostream& os, const VectorObservable& vo);
/**
* Given the name and the units of an observable, initialize the necessary
* struct
*/
ScalarObservable declare_scalar_observable(std::string name,
    std::string units = "",
    std::string axis_name = "");
VectorObservable declare_vector_observable(std::string name,
    std::string units = "",
    std::string axis_name = "");
/**
* Clear the value_time vector in scalar observable
*/
void clear_observable_records(ScalarObservable *so);
/**
* Clear the value_time vector in vector observable
*/
void clear_observable_records(VectorObservable *vo);
/**
* Given the observable, the output directory and the name of the simulation,
* output the observable into outdir/sim_name_observable_name.dat.
* if overwrite = false:
*   if the file already exists, appends data to already written data;
*   if the file does not exist, creates it.
* if overwrite = true:
*   if the file already exists, overwrites it;
*   if the file does not exist, creates it.
*/
void write_scalar_observable_to_file(ScalarObservable& so,
    std::string outdir, std::string sim_name, bool overwrite = false);
void write_vector_observable_to_file(VectorObservable& vo,
    std::string outdir, std::string sim_name, bool overwrite = false);
#endif
