// 2017 Artur Avkhadiev
/*! \file observable_container.h
*/
#ifndef POLYMER_OBSERVABLE_CONTAINER_H
#define POLYMER_OBSERVABLE_CONTAINER_H
#include "observables.h"
#include <map>
#include <set>
    class ObservableContainer {
    protected:
        std::map<std::string, ScalarObservable> scalar_observables;
        std::map<std::string, VectorObservable> vector_observables;
        std::vector<std::string> observable_names;
        // update records
        void update_scalar(std::string name,
            std::pair<double, double> value_time);
        void update_vector(std::string name,
            std::pair<Vector, double> value_time);
        // clear records
        void clear_scalar(std::string name);
        void clear_vector(std::string name);
    public:
        const std::vector<std::string>& get_observable_names() const;
        void add_scalar(ScalarObservable& so);
        void add_vector(VectorObservable& vo);
        ScalarObservable& get_scalar(std::string name);
        VectorObservable& get_vector(std::string name);
        // zeros accumulators of all observables
        void zero_accumulators();
        // clears all records
        void clear();
        // update observables via accumulators
        void update(double time);                       // updates all
        void update(std::string name, double time);     // searches for name
        void writeout(std::string outdir,
            std::string sim_name,
            bool overwrite = false);
        ObservableContainer();
        ~ObservableContainer();
    };
#endif
