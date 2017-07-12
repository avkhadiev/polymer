// 2017 Artur Avkhadiev
/*! \file ljpotential_observable_container.cpp
*/
#include <map>
#include <utility>      /* std::pair, std::make_pair */
#include <string>
#include "../include/observables.h"
#include "../include/ljpotential_observable_container.h"
std::string pair_potential_name = "v_{ij}";
std::string pair_virial_name = "w_{ij}";
std::string pair_fstrength_name = "f_{ij}";
std::string pair_potential_units = "\\epsilon^{-1}";
std::string pair_virial_units = "\\epsilon^{-1}";
std::string pair_fstrength_units = "\\sigma \\epsilon^{-1}";
LJPotentialObservableContainer::LJPotentialObservableContainer()
{
    _pair_potential = declare_scalar_observable(pair_potential_name,
        pair_potential_units);
    _pair_virial = declare_scalar_observable(pair_virial_name,
        pair_virial_units);
    _pair_fstrength = declare_scalar_observable(pair_fstrength_name,
        pair_fstrength_units);
    add_scalar_observable(&_pair_potential);
    add_scalar_observable(&_pair_virial);
    add_scalar_observable(&_pair_fstrength);
};
LJPotentialObservableContainer::~LJPotentialObservableContainer(){
};
ScalarObservable const *LJPotentialObservableContainer::get_pair_potential() const {
    return &_pair_potential;
}
ScalarObservable const *LJPotentialObservableContainer::get_pair_virial() const {
    return &_pair_virial;
}
ScalarObservable const *LJPotentialObservableContainer::get_pair_fstrength() const {
    return &_pair_fstrength;
}
