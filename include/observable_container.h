// 2017 Artur Avkhadiev
/*! \file observable_container.h
*/
#ifndef POLYMER_OBSERVABLE_CONTAINER_H
#define POLYMER_OBSERVABLE_CONTAINER_H
#include <observables.h>
#include <map>
    class ObservableContainer {
    public:
        virtual std::map<std::string, ScalarObservable *> get_scalar_observables(std::vector<std::string> names = {}) = 0;
        virtual std::map<std::string, VectorObservable *> get_vector_observables(std::vector<std::string> names = {}) = 0;
        virtual void update_scalar_observables(std::map<std::string, std::pair<double, double>> time_value) = 0;
        virtual void update_vector_observables(std::map<std::string, std::pair<Vector, double>> time_value) = 0;
        virtual void writeout_observables(std::vector<std::string> names = {}) = 0;
        virtual ~ObservableContainer();
    };
#endif
