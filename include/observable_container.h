// 2017 Artur Avkhadiev
/*! \file observable_container.h
*/
#ifndef POLYMER_OBSERVABLE_CONTAINER_H
#define POLYMER_OBSERVABLE_CONTAINER_H
#include <fstream>
#include <string>
#include <vector>
#include "observable.h"
#include "simple_state.h"
namespace container{
    typedef enum eval_time_t {
        force_loop,             /**> observable evaluated by the force loop */
        integrator,             /**> observable evaluated by the integrator */
        sim_loop                /**> observable evaluated after integration */
    } EvalTime;
    typedef struct unit_t {
        Observable& obs;        /**> observables reside outside the container */
        bool status;            /**> display observable value in status? */
        bool timelog;           /**> does the observable have a time log?*/
        bool average;           /**> is observable an average quantity?  */
        EvalTime eval_time;     /**> when to evaluate the observable?    */
    } Unit;
} // namespace container
class ObservableContainer {
private:
    std::vector<container::Unit> _observables;
    std::vector<container::Unit> _status_variables;
public:
    // I/O
    std::string config_string();           
    std::string status_string();
    void print_status(std::ofstream& output);
    void read_status(std::ifstream& input);
    virtual void write_data(std::string datadir,
            std::string sim_name,
            bool overwrite = false);
    virtual void update(const simple::AtomState& state, size_t calcstep);
    virtual void update(const simple::BondState& state, size_t calcstep);
    ObservableContainer(std::vector<container::Unit>& observables);
    ~ObservableContainer(){};
};
#endif
