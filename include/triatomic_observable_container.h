// 2017 Artur Avkhadiev
/*! \file triatomic_observable_container.h
*/
#ifndef POLYMER_TRIATOMIC_OBSERVABLE_CONTAINER_H
#define POLYMER_TRIATOMIC_OBSERVABLE_CONTAINER_H
#include <map>
#include "observables.h"
#include "observable_container.h"
    class TriatomicObservableContainer :
        public ObservableContainer {
    private:
        // scalar observables
        ScalarObservable _bonds[3];
        ScalarObservable _kinetic_energy;
        ScalarObservable _potential_energy;
        ScalarObservable _bond_angle;
        // vector observables
        VectorObservable _linear_momentum;
        VectorObservable _angular_momentum;
    public:
        // getters
        // default argument measure_time = -1 returns the most recent observation
        const ScalarObservable **get_bonds() const;
        const ScalarObservable *get_kinetic_energy() const;
        const ScalarObservable *get_potential_energy() const;
        const ScalarObservable *get_bond_angle() const;
        const VectorObservable *get_linear_momentum() const;
        const VectorObservable *get_angular_momentum() const;
        // a constructor and a destructor
        TriatomicObservableContainer();
        ~TriatomicObservableContainer();
    };
#endif
