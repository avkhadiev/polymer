// 2017 Artur Avkhadiev
/*! \file triatomic_observable_container
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
        std::map<std::string, ScalarObservable *> _scalar_observables;
        // vector observables
        VectorObservable _linear_momentum;
        VectorObservable _angular_momentum;
        std::map<std::string, VectorObservable *> _vector_observables;
    public:
        // getters
        // default argument measure_time = -1 returns the most recent observation
        const ScalarObservable **get_bonds(double measure_time = -1);
        const ScalarObservable *get_kinetic_energy(double measure_time = -1);
        const ScalarObservable *get_potential_energy(double measure_time = -1);
        const VectorObservable *get_linear_momentum(double measure_time = -1);
        const VectorObservable *get_angular_momentum(double measure_time = -1);
        virtual std::map<std::string, ScalarObservable *> get_scalar_observables(std::vector<std::string> names = {});
        virtual std::map<std::string, VectorObservable *> get_vector_observables(std::vector<std::string> names = {});
        virtual void writeout_observables(std::vector<std::string> names = {});
        // setters
        virtual void update_scalar_observables(std::map<std::string, std::pair<double, double>> time_value);
        virtual void update_vector_observables(std::map<std::string, std::pair<Vector, double>> time_value);
        // a constructor and a destructor
        TriatomicObservableContainer();
        ~TriatomicObservableContainer();
    };
#endif
