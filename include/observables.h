// 2017 Artur Avkhadiev
/*! \file observables.h
*/
#ifndef POLYMER_OBSERVABLES_H
#define POLYMER_OBSERVABLES_H
#include <string>
#include <utility>      // std::pair, std::make_pair
#include <vector>
#include "vector.h"
    typedef struct scalar_observable_t {
        std::vector<std::pair<double, double>> value_time;
        std::string name;
        std::string units;
    } ScalarObservable;
    typedef struct vector_observable_t {
        std::vector<std::pair<Vector, double>> value_time;
        std::string name;
        std::string units;
    } VectorObservable;
    std::string get_scalar_observable_str(ScalarObservable so);
    std::string get_vector_observable_str(VectorObservable vo);
#endif
