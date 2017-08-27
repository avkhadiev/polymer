// 2017 Artur Avkhadiev
/*! \file observables.h
*/
#ifndef POLYMER_OBSERVABLE_H
#define POLYMER_OBSERVABLE_H
#include <string>
#include <fstream>
#include <utility>      /* std::pair, std::make_pair */
#include "vector.h"
#include "simple_state.h"
class Observable {
protected:
    std::string _name;
    std::string _fname;
    std::string _units;
    std::string _axis_name;
    double _value;
public:
    std::string name() const {return _name;};
    std::string fname() const {return _fname;};
    std::string units() const {return _units;};
    std::string axis_name() const {return _axis_name;};
    double value() const {return _value;};
    virtual std::string print_value() const;
    virtual std::string to_string() const
        {return fname() + " " + print_value();};
    void prepare_ofstream(std::string datadir,
        std::string sim_name,
        std::ofstream& writeout,
        bool overwrite = false) const;
    virtual void writeout(std::ofstream& fs, bool add_header);
    // derived classes have designated values, times, and accumulators;
    // these functions update these quantities; they may also change any
    // additional data members introduced in concrete derived classes
    virtual void zero() {_value = 0.0;};
    virtual void add_record(double time){};
    void update(double value) {_value = value;};
    virtual void update(std::string obs_string);  //*> update from to_string */
    virtual void calc_avg(size_t calcsteps){}     //*> update average        */
    virtual void update(const simple::Atom& atom){};
    virtual void update(const simple::Bond& bond){};
    virtual void update(const simple::AtomPolymer& polymer){};
    virtual void update(const simple::BondPolymer& polymer){};
    virtual void update(const simple::AtomState& state){};
    virtual void update(const simple::BondState& state){};
    Observable(std::string name,
        std::string fname,
        std::string units = "",
        std::string axis_name = "",
        double value = 0.0);
    virtual ~Observable(){};
};
class AvgObservable :
    public virtual Observable{
protected:
    Observable& _inst_obs;       /**> corresponding instantaneous observable */
    double _acc;
    virtual void zero() {Observable::zero(); _acc = 0.0;};
public:
    double acc() const {return _acc;};
    virtual std::string print_value() const;
    virtual void update(std::string obs_string);  //*> update from to_string */
    using Observable::update;
    void calc_avg(size_t ncalcsteps){
        _acc += _inst_obs.value();
        update(acc()/((double)ncalcsteps));
    }
    AvgObservable(Observable& inst_obs, double acc = 0.0,
        std::string name = "virtual name",
        std::string fname = "vfname",
        std::string units = "virtual units",
        std::string axis_name = "virtual axis_name");
    virtual ~AvgObservable(){};
};
class TimeLogObservable:
    public virtual Observable{
private:
    typedef struct record_t {
        double time;
        double value;
    } Record;
    std::vector<Record> _records; /*>> time-value pairs */
public:
    virtual void add_record(double time){
        Record rec = {.time = time, .value = value()};
        _records.push_back(rec);
    }
    // writes out records and deletes them from memory
    virtual void writeout(std::ofstream& fs, bool add_header);
    TimeLogObservable(std::string name = "virtual name",
        std::string fname = "vfname",
        std::string units = "virtual units",
        std::string axis_name = "virtual axis_name");
    virtual ~TimeLogObservable(){};
};
::std::ostream& operator<<(::std::ostream& os, const Observable& obs);
#endif
