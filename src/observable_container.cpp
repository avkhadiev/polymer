// 2017 Artur Avkhadiev
/*! \file observable_container.cpp
*/
#include "../include/parsing.h"
#include "../include/observable_container.h"
typedef std::numeric_limits< double > dbl;
namespace observable{
    BlockRecord::BlockRecord(Observable* obs) :
        _obs(obs),
        _inst_values(){}
    void BlockRecord::writeout_pop_oldest_record(std::ofstream& output){
        std::string s;
        if (!records_empty()){
            if (e_format()){
                std::ostringstream out;
                out.precision(dbl::max_digits10);
                out << _inst_values.front();
                s = out.str();
            }
            else {
                s = std::to_string(_inst_values.front());
            }
            _inst_values.pop();
            output << s.c_str();
        }
        else {
            fprintf(stderr, "%s\n",
                "Block Record: can't writeout, record empty");
        }
    }
    std::string BlockRecord::mean_str() const {
        std::string s;
        if (e_format()){
            std::ostringstream out;
            out.precision(dbl::max_digits10);
            out << mean();
            s = out.str();
        }
        else {
            s = std::to_string(mean());
        }
        return s;
    }
    std::string BlockRecord::meansq_str() const{
        std::string s;
        if (e_format()){
            std::ostringstream out;
            out.precision(dbl::max_digits10);
            out << meansq();
            s = out.str();
        }
        else {
            s = std::to_string(meansq());
        }
        return s;
    }
    std::string BlockRecord::fluctuation_str() const{
        std::string s;
        if (e_format()){
            std::ostringstream out;
            out.precision(dbl::max_digits10);
            out << fluctuation();
            s = out.str();
        }
        else {
            s = std::to_string(fluctuation());
        }
        return s;
    }
    void BlockRecord::clear_records(){
        while(!_inst_values.empty()){
            _inst_values.pop();
        }
    }
} // namespace observable
ObservableContainer::ObservableContainer(std::vector<Observable*> observables) :
    _blockstep(0),
    _average_data(false),
    _observable_logs(),
    _observables_to_average(0)
{
    _observable_logs.clear();
    for (Observable *obs : observables){
        if (obs->should_average()) {
            _average_data = true;
            ++_observables_to_average;
        }
        _observable_logs.push_back(observable::BlockRecord(obs));
    }
}
ObservableContainer::~ObservableContainer(){
    _observable_logs.clear();
}
void ObservableContainer::calculate_observables(const simple::AtomState& state){
    for (observable::BlockRecord& record : _observable_logs){
        if (record.obs()->update_time() == observable::MAIN_LOOP){
            record.obs()->update(state);
        }
    }
}
void ObservableContainer::record_observables(){
    for (observable::BlockRecord& record : _observable_logs){
        record.add_newest_record();
        record.obs()->zero_value();
    }
}
void ObservableContainer::run_begin(std::string dtdir, std::string sim_name, bool should_write_data){
    if (should_write_data) prepare_datafiles(dtdir, sim_name);
    _blockstep = 0;
    for (observable::BlockRecord& record : _observable_logs){
        record.obs()->start_new_run();
    }
}
void ObservableContainer::run_end(){
    if (_average_data) {
        if (_blockstep > 1){
            for (observable::BlockRecord& record : _observable_logs){
                if (record.obs()->should_average()){
                    record.obs()->average_run();
                }
            }
        }
        else {
            fprintf(stderr, "%s %zu, %s\n",
                "blockstep = ", _blockstep, "can't average run");
        }
    }
}
void ObservableContainer::block_begin(){
    if (_average_data) {
        ++_blockstep;
        for (observable::BlockRecord& record : _observable_logs){
            record.obs()->start_new_block();
        }
    }
}
void ObservableContainer::update_block(){
    if (_average_data) {
        for (observable::BlockRecord& record : _observable_logs){
            if (record.obs()->should_average()){
                record.obs()->update_block();   // updates *observable*
            }
        }
    }
}
void ObservableContainer::block_end(){
    if (_average_data) {
        for (observable::BlockRecord& record : _observable_logs){
            if (record.obs()->should_average()
                && !record.obs()->is_block_averaged()){
                record.obs()->average_block_update_run();
                if (_blockstep > 1){
                    // saves new run averages
                    // with the new block taken into account
                    record.obs()->average_run();
                }
            }
        }
    }
}
void ObservableContainer::prepare_datafiles(std::string datadir,
    std::string sim_name){
    _prepare_inst_datafile(datadir, sim_name);
    if (_average_data) {
        _prepare_block_datafile(datadir, sim_name);
        _prepare_run_summary_file(datadir, sim_name);
    }
    _record_meta_data(datadir, sim_name);
}
void ObservableContainer::_prepare_inst_datafile(std::string datadir,
    std::string sim_name){
    std::string fout;
    fout = datadir
          + parse_string(sim_name)
          + "_"
          + inst_datafile
          + ".csv";
    std::ofstream writeout;
    writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
    if (writeout.is_open()){
        if (nobservables() > 0) {
            writeout << "block" << ",";
            for (size_t i = 0; i < nobservables(); ++i){
                writeout << _observable_logs.at(i).obs()->short_name();
                if (i < nobservables() - 1){
                    writeout << ",";
                }
                else {
                    writeout << std::endl;
                }
            }
        }
    }
    else {
        // if file still could not be opened
        std::string err_msg = "prepare_inst_datafile: unable to open file at";
        fprintf(stderr, "%s %s. %s \n",
              err_msg.c_str(), fout.c_str(), "file stream will be closed.");
        perror("open");
        writeout.close();
    }
}
void ObservableContainer::_prepare_block_datafile(std::string datadir,
    std::string sim_name){
    std::string fout;
    fout = datadir
          + parse_string(sim_name)
          + "_"
          + block_datafile
          + ".csv";
    std::ofstream writeout;
    writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
    if (writeout.is_open()){
        if (_observables_to_average > 0) {
            writeout << "block" << ",";
            bool something_written;
            size_t nwritten = 0;
            for (size_t i = 0; i < nobservables(); ++i){
                something_written = false;
                if (_observable_logs.at(i).obs()->should_average()){
                    ++nwritten;
                    if (_observable_logs.at(i).obs()->calculate_mean()){
                        writeout
                            << _observable_logs.at(i).obs()->short_name()
                            << "_mean";
                        something_written = true;
                    }
                    if ( _observable_logs.at(i).obs()->calculate_meansq()){
                        if (something_written){
                            writeout << ",";
                        }
                        writeout
                            << _observable_logs.at(i).obs()->short_name()
                            << "_meansq";
                        something_written = true;
                    }
                    if ( _observable_logs.at(i).obs()->calculate_fluctuations()){
                        if (something_written){
                            writeout << ",";
                        }
                        writeout
                            << _observable_logs.at(i).obs()->short_name()
                            << "_fluc";
                        something_written = true;
                    }
                    if ((nwritten < _observables_to_average) && something_written){
                        writeout << ",";
                    }
                }
            }
            writeout << std::endl;
        }
    }
    else {
        // if file still could not be opened
        std::string err_msg = "prepare_block_datafile: unable to open file at";
        fprintf(stderr, "%s %s. %s \n",
              err_msg.c_str(), fout.c_str(), "file stream will be closed.");
        perror("open");
        writeout.close();
    }
}
void ObservableContainer::_prepare_run_summary_file(std::string datadir,
    std::string sim_name){
    std::string fout;
    fout = datadir
          + parse_string(sim_name)
          + "_"
          + run_summary
          + ".csv";
  std::ofstream writeout;
  writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
  if (writeout.is_open()){
      if (_observables_to_average > 0) {
          writeout << "nblocks" << ",";
          bool something_written;
          size_t nwritten = 0;
          for (size_t i = 0; i < nobservables(); ++i){
              something_written = false;
              if (_observable_logs.at(i).obs()->should_average()){
                  ++nwritten;
                  if (_observable_logs.at(i).obs()->calculate_mean()){
                      writeout
                          << _observable_logs.at(i).obs()->short_name()
                          << "_mean";
                      something_written = true;
                  }
                  if ( _observable_logs.at(i).obs()->calculate_meansq()){
                      if (something_written){
                          writeout << ",";
                      }
                      writeout
                          << _observable_logs.at(i).obs()->short_name()
                          << "_meansq";
                      something_written = true;
                  }
                  if ( _observable_logs.at(i).obs()->calculate_fluctuations()){
                      if (something_written){
                          writeout << ",";
                      }
                      writeout
                          << _observable_logs.at(i).obs()->short_name()
                          << "_fluc"
                          << ","
                          << _observable_logs.at(i).obs()->short_name()
                          << "_err";
                      something_written = true;
                  }
                  if ((nwritten < _observables_to_average) && something_written){
                      writeout << ",";
                  }
              }
          }
          writeout << std::endl;
      }
  }
  else {
      // if file still could not be opened
      std::string err_msg = "prepare_run_summary_file: unable to open file at";
      fprintf(stderr, "%s %s. %s \n",
            err_msg.c_str(), fout.c_str(), "file stream will be closed.");
      perror("open");
      writeout.close();
  }
}
void ObservableContainer::_record_meta_data(std::string datadir,
    std::string sim_name){
    std::string fout;
    fout = datadir
          + parse_string(sim_name)
          + "_"
          + meta_datafile
          + ".csv";
    std::ofstream writeout;
    writeout.open(fout, std::ofstream::out | std::ofstream::trunc);
    if (writeout.is_open()){
        if (nobservables() > 0) {
            writeout << "long_name" << ",";
            writeout << "short_name" << ",";
            writeout << "tex_name" << ",";
            writeout << "units" << std::endl;
            for (size_t i = 0; i < nobservables(); ++i){
                writeout
                    << _observable_logs.at(i).obs()->long_name()
                    << ","
                    << _observable_logs.at(i).obs()->short_name()
                    << ","
                    << _observable_logs.at(i).obs()->tex_name()
                    << ","
                    << _observable_logs.at(i).obs()->units()
                    << std::endl;
            }
        }
    }
    else {
        // if file still could not be opened
        std::string err_msg = "prepare_meta_datafile: unable to open file at";
        fprintf(stderr, "%s %s. %s \n",
              err_msg.c_str(), fout.c_str(), "file stream will be closed.");
        perror("open");
        writeout.close();
    }
}
void ObservableContainer::write_data(std::string datadir, std::string sim_name){
    std::ofstream writeout;
    std::string fout;
    // write instantaneous data
    fout = datadir
          + parse_string(sim_name)
          + "_"
          + inst_datafile
          + ".csv";
    writeout.open(fout, std::ofstream::out | std::ofstream::app);
    if (writeout.is_open()){
        _write_inst_data(writeout);
    }
    else {
        std::string err_msg = "write_data: unable to open file at";
        fprintf(stderr, "%s %s. %s \n",
              err_msg.c_str(), fout.c_str(), "file stream will be closed.");
        perror("open");
    }
    writeout.close();
    // write average data, if possible
    if (_average_data){
        fout = datadir
              + parse_string(sim_name)
              + "_"
              + block_datafile
              + ".csv";
        writeout.open(fout, std::ofstream::out | std::ofstream::app);
        if (writeout.is_open()){
            _write_block_data(writeout);
        }
        else {
            std::string err_msg = "write_data: unable to open file at";
            fprintf(stderr, "%s %s. %s \n",
                  err_msg.c_str(), fout.c_str(), "file stream will be closed.");
            perror("open");
        }
        writeout.close();
    }
}
void ObservableContainer::clear_records(){
    for (observable::BlockRecord block : _observable_logs){
        block.clear_records();
    }
}
void ObservableContainer::zero_accumulators(){
    _blockstep = 0;
    for (observable::BlockRecord block : _observable_logs){
        block.obs()->_zero_block();
        block.obs()->_zero_run();
    }
}
void ObservableContainer::_write_inst_data(std::ofstream& writeout){
    if (writeout.is_open()) {
        if (nobservables() > 0){
            while (!_observable_logs.at(0).records_empty()){
                writeout << _blockstep << ",";
                for (size_t i = 0; i < nobservables(); ++i){
                    _observable_logs.at(i).writeout_pop_oldest_record(writeout);
                    if (i < nobservables() - 1){
                        writeout << ",";
                    }
                    else {
                        writeout << std::endl;
                    }
                }
            }
        }
    }
    else {
       // if file still could not be opened
       std::string err_msg = "_write_inst_data: unable to open file";
       fprintf(stderr, "%s\n", err_msg.c_str());
       perror("open");
    }
}
void ObservableContainer::_write_block_data(std::ofstream& writeout){
    if (writeout.is_open()) {
        if (_observables_to_average > 0) {
            writeout << _blockstep << ",";
            bool something_written;
            size_t nwritten = 0;
            for (size_t i = 0; i < nobservables(); ++i){
                something_written = false;
                if (_observable_logs.at(i).obs()->should_average()){
                    ++nwritten;
                    if (_observable_logs.at(i).obs()->calculate_mean()){
                        writeout
                            << _observable_logs.at(i).mean_str();
                        something_written = true;
                    }
                    if ( _observable_logs.at(i).obs()->calculate_meansq()){
                        if (something_written){
                            writeout << ",";
                        }
                        writeout
                            << _observable_logs.at(i).meansq_str();
                        something_written = true;
                    }
                    if ( _observable_logs.at(i).obs()->calculate_fluctuations()){
                        if (something_written){
                            writeout << ",";
                        }
                        writeout
                            << _observable_logs.at(i).fluctuation_str();
                        something_written = true;
                    }
                    if ((nwritten < _observables_to_average) && something_written){
                        writeout << ",";
                    }
                }
            }
            writeout << std::endl;
        }
    }
    else {
       // if file still could not be opened
       std::string err_msg = "_write_block_data: unable to open file";
       fprintf(stderr, "%s\n", err_msg.c_str());
       perror("open");
    }
}
void ObservableContainer::write_run_summary(std::string datadir, std::string sim_name){
    std::ofstream writeout;
    std::string fout;
    // write average data, if possible
    if (_average_data){
        fout = datadir
              + parse_string(sim_name)
              + "_"
              + run_summary
              + ".csv";
        writeout.open(fout, std::ofstream::out | std::ofstream::app);
        if (writeout.is_open()){
            if (_observables_to_average > 0) {
                writeout << _blockstep << ",";
                bool something_written;
                size_t nwritten = 0;
                for (size_t i = 0; i < nobservables(); ++i){
                    something_written = false;
                    if (_observable_logs.at(i).obs()->should_average()){
                        ++nwritten;
                        if (_observable_logs.at(i).obs()->calculate_mean()){
                            writeout
                                << _observable_logs.at(i).obs()->run_mean();
                            something_written = true;
                        }
                        if ( _observable_logs.at(i).obs()->calculate_meansq()){
                            if (something_written){
                                writeout << ",";
                            }
                            writeout
                                << _observable_logs.at(i).obs()->run_meansq();
                            something_written = true;
                        }
                        if ( _observable_logs.at(i).obs()->calculate_fluctuations()){
                            if (something_written){
                                writeout << ",";
                            }
                            writeout
                                << _observable_logs.at(i).obs()->run_fluctuation()
                                << ","
                                << _observable_logs.at(i).obs()->run_err();
                            something_written = true;
                        }
                        if ((nwritten < _observables_to_average) && something_written){
                            writeout << ",";
                        }
                    }
                }
                writeout << std::endl;
            }
        }
        else {
            std::string err_msg = "write_run_summary: unable to open file at";
            fprintf(stderr, "%s %s. %s \n",
                  err_msg.c_str(), fout.c_str(), "file stream will be closed.");
            perror("open");
        }
        writeout.close();
    }
}
// ObservableContainer::ObservableContainer
//    (std::vector<container::Unit>& observables) :
//    _observables(observables){
//    _status_variables.clear();
//    for (const container::Unit& unit: _observables){
//        if (unit.status == true) _status_variables.push_back(unit);
//    }
//}
//// I/O
//std::string ObservableContainer::config_string(){
//    std::string config = "";
//    for (const container::Unit& unit: _observables){
//        config += unit.obs.to_string() + "\n";
//    }
//    return config;
//}
std::string ObservableContainer::status_string(bool verbose) {
    std::string status = "";
    for (observable::BlockRecord& record : _observable_logs){
        if (record.obs()->print_inst_val()) {
            status += record.obs()->status_string() + "\n";
        }
    }
    return status;
}
void ObservableContainer::print_status(std::ofstream& output, bool verbose) {
    output << status_string(verbose).c_str();
}
void ObservableContainer::read_status(std::ifstream& input){
    std::string line;
    for (observable::BlockRecord& record : _observable_logs){
        if (record.obs()->print_inst_val()) {
            std::getline(input, line);
            record.obs()->read_status(line);
        }
    }
}
//void ObservableContainer::write_data(std::string datadir,
//    std::string sim_name, bool overwrite){
//        std::ofstream fs;
//        for (const container::Unit& unit: _observables){
//            unit.obs.prepare_ofstream(datadir, sim_name, fs, overwrite);
//            unit.obs.writeout(fs, overwrite);
//        }
//}
//void ObservableContainer::update(const simple::AtomState& state, size_t calcstep){
//    for (const container::Unit& unit: _observables){
//        if (unit.eval_time == container::EvalTime::sim_loop){
//            if (unit.average) {
//                // average observables are not zeroed -- accumulators should not
//                // be reset.
//                unit.obs.calc_avg(calcstep);
//            }
//            else {
//                unit.obs.zero();
//                unit.obs.update(state);
//            }
//        }
//        if (unit.timelog) unit.obs.add_record(state.time());
//    }
//}
