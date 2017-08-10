// 2017 Artur Avkhadiev
/*! \file observable_container.cpp
*/
#include <map>
#include <vector>
#include <stdexcept>               /* std::out_of_range, std::invalid_argument */
#include <utility>                 /* std::pair, std::make_pair */
#include "../include/observables.h"
#include "../include/parsing.h"
#include "../include/observable_container.h"
ObservableContainer::ObservableContainer(){
}
ObservableContainer::~ObservableContainer(){
}
const std::vector<std::string>& ObservableContainer::get_observable_names() const {
    return observable_names;
}
void ObservableContainer::add_scalar(ScalarObservable& so) {
    scalar_observables.emplace(so.name, so);
    observable_names.push_back(so.name);
}
void ObservableContainer::add_vector(VectorObservable& vo) {
    vector_observables.emplace(vo.name, vo);
    observable_names.push_back(vo.name);
}
ScalarObservable& ObservableContainer::get_scalar(std::string name) {
    try
    {
        return scalar_observables.at(name);
    }
    catch (std::out_of_range& e)
    {
        std::string err_msg = "get_scalar_observable: \""
            + name
            + "\" observable does not exist";
        throw std::invalid_argument(err_msg);
    }
}
VectorObservable& ObservableContainer::get_vector(std::string name) {
    try
    {
        return vector_observables.at(name);
    }
    catch (std::out_of_range& e)
    {
        std::string err_msg = "get_vector_observable: \""
            + name
            + "\" observable does not exist";
        throw std::invalid_argument(err_msg);
    }
}
void ObservableContainer::clear_scalar(std::string name){
    ScalarObservable& so = get_scalar(name);
    clear_observable_records(so);
}
void ObservableContainer::clear_vector(std::string name){
    VectorObservable& vo = get_vector(name);
    clear_observable_records(vo);
}
void ObservableContainer::clear(){
    for (std::string& name : observable_names){
        if (scalar_observables.find(name) != scalar_observables.end()){
            clear_scalar(name);
        }
        // update a vector observable if such a name exists
        if (vector_observables.find(name) != vector_observables.end()){
            clear_vector(name);
        }
    }
}
void ObservableContainer::update_scalar(std::string name,
    std::pair<double, double> value_time) {
    ScalarObservable &so = get_scalar(name);
    so.value_time.push_back(value_time);
}
void ObservableContainer::update_vector(std::string name,
    std::pair<Vector, double> value_time) {
    VectorObservable& vo = get_vector(name);
    vo.value_time.push_back(value_time);
}
void ObservableContainer::update(std::string name, double time){
    if (scalar_observables.find(name) != scalar_observables.end()){
        update_scalar(name,
            std::make_pair(scalar_observables.at(name).accumulator, time));
    }
    // update a vector observable if such a name exists
    if (vector_observables.find(name) != vector_observables.end()){
        update_vector(name,
            std::make_pair(vector_observables.at(name).accumulator, time));
    }
}
void ObservableContainer::update(double time){
    for(std::string& name : observable_names){
        update(name, time);
    }
}
void ObservableContainer::writeout(std::string outdir,
    std::string sim_name,
    bool overwrite)
{
    for (std::string& name : observable_names){
        // print scalar observable if key exists
        if (scalar_observables.find(name) != scalar_observables.end()){
            write_scalar_observable_to_file(scalar_observables.at(name),
                outdir,
                sim_name,
                overwrite);
        }
        // print vector observable if key exists
        if (vector_observables.find(name) != vector_observables.end()){
            write_vector_observable_to_file(vector_observables.at(name),
                outdir,
                sim_name,
                overwrite);
        }
    }
}
