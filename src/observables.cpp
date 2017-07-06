// 2017 Artur Avkhadiev
/*! \file observables.cpp
*/
#include <map>
#include <utility>      /* std::pair, std::make_pair */
#include <vector>
#include <string>
#include "../include/vector.h"
#include "../include/observables.h"
std::string scalar_observable_to_string(ScalarObservable vo) {
    std::string vo_str;
    std::string vo_name = vo.name;
    std::string vo_units = vo.units;
    vo_str = vo.name + ", " + vo.units;
    return vo_str;
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
    ScalarObservable so = {.value_time = empty_vector, .name = name, .units = units};
    return so;
}
VectorObservable declare_vector_observable(std::string name, std::string units) {
    std::vector<std::pair<Vector, double> > empty_vector;
    VectorObservable vo = {.value_time = empty_vector, .name = name, .units = units};
    return vo;
}
