// 2017 Artur Avkhadiev
/*! \file ljpotential_observable_container.cpp
*/
#include <map>
#include <utility>      /* std::pair, std::make_pair */
#include <string>
#include "../include/observables.h"
#include "../include/ljpotential_observable_container.h"
LJPotentialObservableContainer::LJPotentialObservableContainer()
{
    _pair_potential = declare_scalar_observable("Pair Potential",
        "\\epsilon",
        "v_{ij}");
    _pair_virial = declare_scalar_observable("Pair Virial",
        "\\epsilon",
        "w_{ij}");
    _pair_fstrength = declare_scalar_observable("Force Strength",
        "\\epsilon \\sigma^{-1}",
        "f_{ij}");
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
