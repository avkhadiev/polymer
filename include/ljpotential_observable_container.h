// 2017 Artur Avkhadiev
/*! \file ljpotential_observable_container.h
*/
#ifndef LJPOTENTIAL_OBSERVABLE_CONTAINER_H
#define LJPOTENTIAL_OBSERVABLE_CONTAINER_H
#include <map>
#include "observables.h"
#include "observable_container.h"
    class LJPotentialObservableContainer :
        public ObservableContainer {
    private:
        // scalar observables
        ScalarObservable _pair_potential;
        ScalarObservable _pair_virial;
        ScalarObservable _pair_fstrength;
        ScalarObservable _distance;
    public:
        // a constructor and a destructor
        LJPotentialObservableContainer();
        ~LJPotentialObservableContainer();
    };
#endif
