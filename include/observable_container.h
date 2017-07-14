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
        std::vector<std::string> get_observable_names() const;
        // adds scalar observable pointer to the dictionary keyed on
        // scalar observables' names; appends the name to the vector of names
        // of all observables (vector and scalar)
        virtual void add_scalar_observable(ScalarObservable *so);
        // applies add_scalar_observable to all observables in the vector
        virtual void add_scalar_observables(std::vector<ScalarObservable *> so_vec);
        // adds vector observable pointer to the dictionary keyed on
        // vector observables' names; appends the name to the vector of names
        // of all observables (vector and scalar)
        virtual void add_vector_observable(VectorObservable *vo);
        // applies add_vector_observable to all observables in the vector
        virtual void add_vector_observables(std::vector<VectorObservable *> vo_vec);
        // returns a pointer to scalar observable given its name
        virtual ScalarObservable *get_scalar_observable(std::string name);
        // returns a pointer to vector observable given its name
        virtual VectorObservable *get_vector_observable(std::string name);

        // returns a pointer to scalar observable accumulator given its name
        virtual double *get_scalar_observable_accumulator(std::string name);
        // returns a pointer to vector observable accumulator given its name
        virtual Vector *get_vector_observable_accumulator(std::string name);
        // zeroes scalar observable accumulator given its name
        virtual void zero_scalar_observable_accumulator(std::string name);
        // zeroes vector observable accumulator given its name
        virtual void zero_vector_observable_accumulator(std::string name);
        // zeros accumulators of all observables
        virtual void zero_accumulators(std::vector<std::string> names = {});
        // clears all records of the scalar observable with the given name
        virtual void clear_scalar_observable_records(std::string name);
        // clears all records of the vector observable with the given name
        virtual void clear_vector_observable_records(std::string name);
        // clears all recods of observables with names in the given vector of
        // names
        virtual void clear_observables_records(std::vector<std::string> names={});
        // pushes the given pair (value, time) to a vector of observed values
        // for a scalar observable of the given name
        virtual void update_scalar_observable(std::string name, std::pair<double, double> value_time);
        // pushes the given pair (value, time) to a vector of observed values
        // for a vector observable of the given name
        virtual void update_vector_observable(std::string name, std::pair<Vector, double> value_time);
        // updates an observable by name by pushing its accumulator paired with
        // the given time on the back of the vector of observed values.
        virtual void update_observable_through_accumulator(std::string name, double time);
        // applies update_observable_through_accumulator to all observables
        // corresponding to a vector of given names; if the vector is empty,
        // applies to all observables
        virtual void update_observables_through_accumulators(std::vector<std::string> names, double time);
        // writes out observables to oudir/sim_name_<observable_name>.dat
        // if the vector of names is empty, outputs all observables.
        virtual void writeout_observables_to_file(std::vector<std::string> names,
            std::string outdir,
            std::string sim_name,
            bool overwrite = false);
        ObservableContainer();
        virtual ~ObservableContainer();
    };
#endif
