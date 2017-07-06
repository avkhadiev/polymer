// 2017 Artur Avkhadiev
/*! \file observables.h
*/
#ifndef POLYMER_OBSERVABLES_H
#define POLYMER_OBSERVABLES_H
#include <string>
#include <utility>      /* std::pair, std::make_pair */
#include <vector>
#include "vector.h"
    typedef struct scalar_observable_t {
        std::vector<std::pair<double, double> > value_time;
        std::string name;
        std::string units;
    } ScalarObservable;
    typedef struct vector_observable_t {
        std::vector<std::pair<Vector, double> > value_time;
        std::string name;
        std::string units;
    } VectorObservable;
    /**
    * Given the observable, return a string with "Name, Units"
    */
    std::string scalar_observable_to_string(ScalarObservable so);
    std::string vector_observable_to_string(VectorObservable vo);
    ::std::ostream& operator<<(::std::ostream& os, const ScalarObservable& so);
    ::std::ostream& operator<<(::std::ostream& os, const VectorObservable& vo);
    /**
    * Given the name and the units of an observable, initialize the necessary
    * struct
    */
    ScalarObservable declare_scalar_observable(std::string name, std::string units);
    VectorObservable declare_vector_observable(std::string name, std::string units);
#endif
