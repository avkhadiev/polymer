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
    _pair_potential = declare_scalar_observable("Potential",
        "\\epsilon",
        "v_{ij}");
    _pair_virial = declare_scalar_observable("Virial",
        "\\epsilon",
        "w_{ij}");
    _pair_fstrength = declare_scalar_observable("Force",
        "\\epsilon \\sigma^{-1}",
        "f_{ij}");
    _distance = declare_scalar_observable("Distance",
        "\\sigma",
        "r_{ij}");
    add_scalar_observable(&_pair_potential);
    add_scalar_observable(&_pair_virial);
    add_scalar_observable(&_pair_fstrength);
    add_scalar_observable(&_distance);
};
LJPotentialObservableContainer::~LJPotentialObservableContainer(){
};
