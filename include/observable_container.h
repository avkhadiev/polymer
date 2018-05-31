// 2017 Artur Avkhadiev
/*! \file observable_container.h
*/
#ifndef POLYMER_OBSERVABLE_CONTAINER_H
#define POLYMER_OBSERVABLE_CONTAINER_H
#include <fstream>
#include <string>
#include <vector>
#include <functional>         /** std::ref */
#include <queue>              /** std::queue */
#include <limits>             /**> controls output precision */
#include <stdexcept>
#include "observable.h"
#include "simple_state.h"
namespace observable{
    class BlockRecord {
    private:
        Observable *_obs;
        std::queue<const double> _inst_values;
    public:
        Observable *obs() {return _obs;};
        friend class ObservableContainer;
        // getters
        size_t number() const {return _obs->block_number();};
        double mean() const {return _obs->block_mean();};
        double meansq() const {return _obs->block_meansq();};
        double fluctuation() const {return _obs->block_fluctuation();};
        bool e_format() const {return _obs->e_format();};
        bool is_averaged() const {return _obs->is_block_averaged();};
        std::string mean_str() const;
        std::string meansq_str() const;
        std::string fluctuation_str() const;
        // handling records
        bool records_empty() const {return _inst_values.empty();};
        size_t nrecords() const {return _inst_values.size();};
        void add_newest_record() {_inst_values.push(_obs->value);};
        // IO
        void writeout_pop_oldest_record(std::ofstream& output);
        void clear_records();
        // constructor-destructor
        BlockRecord(Observable* obs);
        ~BlockRecord() {clear_records();};
    };
} // namespace observable
class ObservableContainer {
private:
    size_t _blockstep;
    bool _average_data;
    std::vector<observable::BlockRecord> _observable_logs;
    // helper I/O
    const std::string inst_datafile = "inst_data";
    const std::string block_datafile = "block_data";
    const std::string meta_datafile = "meta_data";
    const std::string run_summary = "run_summary";
    void _prepare_inst_datafile(std::string datadir,
        std::string sim_name);
    void _prepare_block_datafile(std::string datadir,
        std::string sim_name);
    void _prepare_run_summary_file(std::string datadir,
        std::string sim_name);
    // writeout a meta-data file with observable names, axis_names, etc. in CSV
    void _record_meta_data(std::string datadir,
        std::string sim_name);
    // given a block record, writeout its contents and erase data
    void _writeout_block(std::ofstream& inst, std::ofstream& block);
    void _write_inst_data(std::ofstream& inst);
    void _write_block_data(std::ofstream& block);
    void _zero_instantaneous_values();
    size_t _observables_to_average;
public:
    void set_average_data(bool should) {_average_data = should;};
    bool should_average() const {return _average_data;};
    // updating observables
    // will calculate any observables to be calculated in the MAIN_LOOP
    virtual void calculate_observables(const simple::AtomState& state);
    // will record all observables' instantaneous values and then zero them
    void add_observable(Observable* observable);
    void record_observables();
    void zero_accumulators();
    // will erase all recorded data
    void clear_records();
    // if should_write_data: setup datafile streams --- instantaneous (& average, if neccessary)
    // set blockstep == 0
    // call start_new_run for all observable
    void run_begin(std::string dtdir, std::string sim_name, bool should_write_data);
    // perform final run averaging for all observables
    //  ensure _average_data == true
    //  ensure blockstep > 1 (can't average otherwise)
    void run_end();
    // if averaging is enabled,
    //  increment blockstep
    //  begin blocks for all observables,
    //      regardless of their averaging instructions: this ensures all records
    //      made at the same time have the same block number
    //  set all block record averages to 0.0 (done implicitly);
    void block_begin();
    // if averaging is enabled,
    //  update blocks for all observables that require averaging
    void update_block();
    // if averaging is enabled,
    //  average blocks and update run averages
    //      for all observables that require averaging
    //  update all block records of averages, set block record as averaged
    //  if blockstep > 1, also perform run averaging (to reflect newest changes)
    void block_end();
    // counters
    size_t nobservables() const {return _observable_logs.size();};
    // I/O
    // for observables that request instantaneous value in the status
    // output the instantaneous value
    std::string status_string(bool verbose = false);
    void print_status(std::ofstream& output, bool verbose = false);
    void read_status(std::ifstream& intput);
    // prepare instantaneous datafile by writing out the column headers for CSV
    // if averaging is enabled, prepare average datafile similarly
    // record a meta-data file with observable names, axis_name, and other data
    void prepare_datafiles(std::string datadir,
        std::string sim_name);
    // open instantaneous data file.
    // if necessary, open block data file
    // loop over all block recorsd
    // writeout data with write_block
    void write_data(std::string datadir, std::string sim_name);
    void write_run_summary(std::string datadir, std::string sim_name);
    // constructors and destructors
    // make a Log for each observable,
    // build a Log vector
    // disable averaging (by default, simulation will enable if necessary)
    // set blockstep = 0;
    ObservableContainer(std::vector<Observable*> observables);
    // delete all block records for each log
    ~ObservableContainer();
};
#endif
