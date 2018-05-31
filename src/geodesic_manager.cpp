// 2018 Artur Avkhadiev
/*! \file geodesic_manager.cpp
*/
#include "../include/geodesic_manager.h"
namespace geodesic{
    /***************************************************************************
    *                               GEODESIC MANAGER
    ***************************************************************************/
    Manager::Manager() :
        _fupd(),
        _states(),
        _states_read(false),
        _initial(),
        _final(){}
    Manager::Manager(Potential* polymer_potential,
                     Potential* solvent_potential,
                     Potential* inter_potential) :
        _fupd(ForceUpdater( polymer_potential,
                            solvent_potential,
                            inter_potential)),
        _states(),
        _states_read(false),
        _initial(),
        _final(){}
    Manager::~Manager(){}
    ForceUpdater Manager::fupd() const {
        return _fupd;
    }
    void Manager::set_fupd(ForceUpdater &fupd){
        _fupd = fupd;
    }
    double Manager::_pe(simple::BondState& state){
        simple::AtomState atom_state = simple::AtomState(state);
        double pe = _fupd.calc_pot_energy(atom_state);
        return pe;
    }
    void Manager::read_states(std::string cndir, std::string sim_name){
        std::string fname = sim_name + "_cn.txt";
        _states.clear();
        simple::read_states_from_file(cndir, fname, _states);
        // read first and last state
        if (_states.size() > 0){
            _states_read = true;
            _initial = Record(_states.front(), _pe(_states.front()));
            _final = Record(_states.back(), _pe(_states.back()));
        }
        else {
            _states_read = false;
            fprintf(stderr, "%s\n", "geodesic::Manager::read_states: no states were read!");
        }
    }
    Record Manager::initial() const{
        if (!_states_read){
            fprintf(stderr, "%s\n", "geodesic::Manager: no states have been read yet, but boundary value requested");
        }
        return _initial;
    }
    Record Manager::final() const{
        if (!_states_read){
            fprintf(stderr, "%s\n", "geodesic::Manager: no states have been read yet, but boundary value requested");
        }
        return _final;
    }
    std::string Manager::ini_fname(std::string sim_name) const {
        return sim_name + "_ini.cfg";
    }
    std::string Manager::fin_fname(std::string sim_name) const {
        return sim_name + "_fin.cfg";
    }
    void Manager::write_geodesic_inputs(std::string outdir, std::string sim_name) const{
        bool overwrite = true;
        std::string ini_path = outdir + ini_fname(sim_name);
        std::string fin_path = outdir + fin_fname(sim_name);;
        initial().write(ini_path, overwrite);
        final().write(fin_path, overwrite);
    }
    Path Manager::MD_path(){
        Path MD_path = Path();
        if (_states_read){
            std::list<Record> records;
            // build a list of records from a vector of states
            for(simple::BondState& state : _states){
                records.push_back(Record(state, _pe(state)));
            }
            MD_path = Path(records);
        }
        else {
            fprintf(stderr, "%s\n", "geodesic::Manager: no states have been read yet, but MD path requested");
        }
        return MD_path;
    }
    // only works with the first 2 links in the polymer
    void Manager::output_observables(Path& path,
        std::string dtdir, std::string name){
        // make observables
        bool compute_mean = false;
        bool compute_err = false;
        bool print_val = true;
        bool e_format = false;
        // Length
        Length length = Length(print_val, e_format);
        // omega projections
        std::vector<simple::Bond> ini_bonds, fin_bonds;
        ini_bonds = path.initial().state().get_polymers().at(0).get_bonds();
        fin_bonds = path.final().state().get_polymers().at(0).get_bonds();
        Vector ini_omega1, ini_omega2, fin_omega1, fin_omega2;
        ini_omega1 = ini_bonds.at(0).position;
        ini_omega2 = ini_bonds.at(1).position;
        fin_omega1 = fin_bonds.at(0).position;
        fin_omega2 = fin_bonds.at(1).position;
        OmegaProj omega_proj1
            = OmegaProj(ini_omega1, fin_omega1, 1,
                compute_mean, compute_err, print_val, e_format);
        OmegaProj omega_proj2
            = OmegaProj(ini_omega2, fin_omega2, 2,
                compute_mean, compute_err, print_val, e_format);
        // angles
        Psi psi1 = Psi(1, fin_omega1, print_val, e_format);
        Psi psi2 = Psi(2, fin_omega2, print_val, e_format);
        DeltaTheta delta_theta1
            = DeltaTheta(1, fin_omega1,
                compute_mean, compute_err, print_val, e_format);
        DeltaTheta delta_theta2
            = DeltaTheta(2, fin_omega2,
                compute_mean, compute_err, print_val, e_format);
        // observable container
        std::vector<Observable*> observables_vec
            = { &length,
                &omega_proj1, &omega_proj2,
                &psi1, &psi2, &delta_theta1, &delta_theta2};
        ObservableContainer obs = ObservableContainer(observables_vec);
        // traverse path, compute observables
        std::list<Record> full_path = path.get_path();
        std::list<Record>::const_iterator rec;
        bool should_write_data = true;
        bool verbose = true;
        size_t idata = 1000;
        size_t iprint = 1000;
        size_t step = 0;
        obs.run_begin(dtdir, name, should_write_data);
        for (rec = full_path.begin();
             rec != std::prev(full_path.end());
             std::advance(rec, 1))
        {
            ++step;
            length.update(*rec, *(std::next(rec)));
            omega_proj1.update(*(std::next(rec)));
            omega_proj2.update(*(std::next(rec)));
            psi1.update(*(std::next(rec)));
            psi2.update(*(std::next(rec)));
            delta_theta1.update(*rec, *(std::next(rec)));
            delta_theta2.update(*rec, *(std::next(rec)));
            obs.record_observables();
            if (step % iprint == 0) {
                fprintf(stdout, "%s", obs.status_string(verbose).c_str());
            }
            if (step % idata == 0) obs.write_data(dtdir, name);
        }
        if (step % idata != 0) obs.write_data(dtdir, name);
    }
} // namespace geodesic
