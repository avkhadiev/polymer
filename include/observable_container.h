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
        std::map<std::string, ScalarObservable *> _scalar_observables;
        std::map<std::string, VectorObservable *> _vector_observables;
        std::vector<std::string> _observable_names;
    public:
        virtual void add_scalar_observable(ScalarObservable *so);
        virtual void add_scalar_observables(std::vector<ScalarObservable *> so_vec);
        virtual void add_vector_observable(VectorObservable *vo);
        virtual void add_vector_observables(std::vector<VectorObservable *> vo_vec);
        virtual ScalarObservable *get_scalar_observable(std::string name);
        virtual VectorObservable *get_vector_observable(std::string name);
        virtual void update_scalar_observable(std::string name, std::pair<double, double> value_time);
        virtual void update_vector_observable(std::string name, std::pair<Vector, double> value_time);
        virtual void writeout_observables_to_file(std::vector<std::string> names,
            std::string outdir,
            std::string sim_name,
            bool overwrite = false);
        virtual ~ObservableContainer();
    };
#endif
