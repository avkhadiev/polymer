// 2017 Artur Avkhadiev
/*! \file observables.h
*/
#ifndef POLYMER_OBSERVABLE_H
#define POLYMER_OBSERVABLE_H
#include <string>
#include <fstream>
#include <utility>      /* std::pair, std::make_pair */
#include <cmath>        /* pow */
#include <stdexcept>
#include "vector.h"
#include "parsing.h"
#include "simple_state.h"
namespace observable {
    typedef struct {            // names for the observable
        std::string full;       // observable full name, descriptive and telling
        std::string abridged;   // observable short name, a CSV file column name
        std::string latex;      // observable LaTeX name, for plotting
        std::string units;      // units of measurement, for plotting
    } name_t;
    typedef enum {              // when to update an observable?
        FORCE_LOOP,             // during force calculations
        INTEGRATION_STEP,       // during the integration step
        MAIN_LOOP               // after the integration step
    } update_time_t;
    typedef struct {            // estimate statistical properties of observable?
        bool mean;              // estimate the first moment?
        bool meansq;            // estimate the second moment?
    } calculate_avg_t;
    std::string to_string(const name_t& s, bool verbose = true);
    std::string to_string(const update_time_t& s, bool verbose = true);
    std::string to_string(const calculate_avg_t& s, bool verbose = true);
    ::std::ostream& operator<<(::std::ostream& os, const name_t& s);
    ::std::ostream& operator<<(::std::ostream& os, const update_time_t& s);
    ::std::ostream& operator<<(::std::ostream& os, const calculate_avg_t& s);
} // namespace observable
class Observable {
protected:
    observable::name_t _name;           // potentially parsed input argument
    observable::update_time_t _update_time;
    observable::calculate_avg_t _calculate;
    bool _print_inst_val;               // print inst val in log&state updates?
    bool _e_format;                     // use precision notation for writeout?
    // structures for proper averaging
    double _msd_additive_const;         // for mean-square-deviation, 0.0 dflt
    typedef struct {
        double mean;                    // mean value
        double meansq;                  // mean of the square value
        double fluctuation;             // fluctuation value
        double mean_acc;                // mean value accumulator
        double meansq_acc;              // mean of the square accumulator
        size_t norm;                    // normalization counter
    } averages_t;
    averages_t _run;                    // run averages
    averages_t _block;                  // block averages
    double _run_err;
    size_t _block_number;
    void _average_block();              // finish block averages
    void _update_run();                 // update run averages
    bool _is_block_averaged;            // has block been averaged?
    void _zero_block();                 // reset all block averages
    void _zero_run();                   // reset all block averages
    friend class ObservableContainer;   // needs access to average information
public:
    double value;                       // instantaneous value of observable
    std::string print_value() const;    // uses precision notation if needed
    void zero_value() {value = 0.0;};
    // getters --- general infromation
    std::string long_name() const {return _name.full;};
    std::string short_name() const {return _name.abridged;};
    std::string tex_name() const {return _name.latex;};
    std::string units() const {return _name.units;};
    observable::update_time_t update_time() const {return _update_time;};
    bool should_average() const {
        return (_calculate.mean || _calculate.meansq);
    };
    bool calculate_mean() const {return _calculate.mean;};
    bool calculate_meansq() const {return _calculate.meansq;};
    bool calculate_fluctuations() const {
        return (_calculate.mean && _calculate.meansq);
    };
    virtual bool update_method_specified() const {return false;};
    bool print_inst_val() const {return _print_inst_val;};
    double msd_additive_const() const {return _msd_additive_const;};
    bool e_format() const {return _e_format;};
    bool is_block_averaged() const {return _is_block_averaged;};
    // getters --- values
    // this getters may raise an exception if requested information does not m
    // make sense for the given observable --- e.g., if fluctuations are not
    // computed but user requests fluctuations
    size_t block_number() const;
    double block_mean() const;
    double block_meansq() const;
    double block_fluctuation() const;
    double run_mean() const;
    double run_meansq() const;
    double run_fluctuation() const;
    double run_err() const;
    // setters
    void set_print_inst_val(bool print) {_print_inst_val = print;};
    void set_e_format(bool e_format) {_e_format = e_format;};
    void set_msd_add_const(bool constant) {_msd_additive_const = constant;};
    // averaging methods
    // zero block, set average = false, increment block number
    void start_new_block();
    void update_block();                // update all block average counters
    // AVERAGE BLOCK:
    //  must have # instvalues > 1
    //  every block is averaged only once. is_block_averaged keeps track
    // UPDATE RUN: update run average counters
    void average_block_update_run();
    void start_new_run();               // zero_run counter, set block # = 0
    void zero_run();                    // zero run averages
    // perform run averages (can be done many times)
    // must have # blocks > 1
    void average_run();
    // I/O
    std::string to_string(bool verbose = true) const;
    // outputs instantaneous value
    std::string status_string() const;
    void read_status(std::string line);
    // update method: does not do anything; can be specified in derived classes
    virtual void update(const simple::AtomState& state){};
    Observable(observable::name_t name,
        observable::update_time_t update_time,
        observable::calculate_avg_t calc_instructions,
        bool print_inst_val,   // print out instantaneous value in log & state?
        bool e_format = false  // use precision notation when writing out?
    );
    virtual ~Observable(){};
};
::std::ostream& operator<<(::std::ostream& os, const Observable& obs);
#endif
