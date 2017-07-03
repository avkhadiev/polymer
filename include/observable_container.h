// 2017 Artur Avkhadiev
/*! \file observable_container.h
*/
#ifndef POLYMER_OBSERVABLE_CONTAINER_H
#define POLYMER_OBSERVABLE_CONTAINER_H
#include <observables.h>
#include <map>
    class ObservableContainer {
    public:
        virtual void get_scalar_observables() = 0;
        virtual void get_vector_observables() = 0;
        virtual void update_scalar_observables(std::map<std::string, double> updates) = 0;
        virtual void update_vector_observables(std::map<std::string, Vector> updates) = 0;
        virtual void writeout_observables() = 0;
        virtual ~ObservableContainer();
    };
#endif
