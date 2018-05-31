// 2018 Artur Avkhadiev
/*! \file geodesic_path_computer.cpp
*/
#include "../include/geodesic_path_computer.h"
namespace geodesic{
    /***************************************************************************
    *                 GEODESIC PATH COMPUTER (BASE CLASS)
    ***************************************************************************/
    PathComputer::PathComputer(Potential* polymer_potential,
        Potential* solvent_potential,
        Potential* inter_potential,
        double epsilon) :
    _fupd(ForceUpdater(polymer_potential, solvent_potential, inter_potential)),
    _epsilon(epsilon){
    }
    PathComputer::~PathComputer(){};
    void PathComputer::update_PE(Record& rec){
        simple::AtomState atom_state = simple::AtomState(rec.state());
        rec.set_pe(_fupd.calc_pot_energy(atom_state));
    }
    /***************************************************************************
    *                       GEODESIC SLERP PATH COMPUTER
    ***************************************************************************/
    SLERP::SLERP(Potential* polymer_potential,
        Potential* solvent_potential,
        Potential* inter_potential,
        double epsilon) :
    PathComputer(polymer_potential,solvent_potential,inter_potential,epsilon),
    psi1(1, vector(0.0, 0.0, 0.0), true, true),
    psi2(2, vector(0.0, 0.0, 0.0), true, true){}
    SLERP::~SLERP(){};
    SLERP::Omega::Omega() :
        _should_move(true),
        _cur(vector(0.0, 0.0, 0.0)),
        _fin(vector(0.0, 0.0, 0.0)),
        _cosPsi(0.0),
        _sinPsi(0.0),
        _psi(0.0){}
    SLERP::Omega::Omega(const simple::Bond& cur, const simple::Bond& fin) :
        _should_move(true),
        _cur(cur.position),
        _fin(fin.position),
        _cosPsi(dot(_cur, _fin) / (norm(_cur) * norm(_fin))),
        _sinPsi(norm(cross(_cur, _fin)) / (norm(_cur) * norm(_fin))),
        _psi(acos(_cosPsi)){}
    void SLERP::Omega::_recompute_angles() {
        _cosPsi = dot(_cur, _fin) / (norm(_cur) * norm(_fin));
        _sinPsi = norm(cross(_cur, _fin)) / (norm(_cur) * norm(_fin));
        _psi =  acos(_cosPsi);
    }
    SLERP::Omega::~Omega(){};
    double SLERP::Omega::_a(double dtau) const{
        return sin(_psi * (1 - dtau)) / _sinPsi;
    }
    double SLERP::Omega::_b(double dtau) const{
        return sin(_psi * dtau) / _sinPsi;
    }
    Vector SLERP::Omega::at(double dtau) const {
        return add(multiply(_cur, _a(dtau)), multiply(_fin, _b(dtau)));
    }
    void SLERP::Omega::update(Vector cur) {
        _cur = cur;
        _recompute_angles();
    }
    void SLERP::Omega::update(const simple::Bond& bond) {
        update(bond.position);
    }
    void SLERP::_update_bond_from_link(simple::Bond& bond, const Omega& omega){
        bond.position = omega.cur();
    }
    void SLERP::_update_record_from_links(Record& record,
        const std::vector<Omega>& links){
        simple::BondPolymer& polymer = record.modifiable_state().polymers.at(0);
        for(size_t i = 0; i < polymer.nb(); ++i){
            _update_bond_from_link(polymer.bonds.at(i), links.at(i));
        }
    }
    void SLERP::_update_links_from_record(std::vector<Omega>& links,
        const Record& record){
        simple::BondPolymer& polymer = record.state().polymers.at(0);
        for(size_t i = 0; i < links.size(); ++i){
            links.at(i).update(polymer.bonds.at(i));
        }
    }
    bool SLERP::Omega::is_close(double epsilon) const {
        return abs(1.0 - _cosPsi) < epsilon;
    }
    bool SLERP::_is_path_complete(std::vector<Omega>& links, double epsilon){
        bool is_path_complete = true;
        bool is_link_close;
        for(Omega& link : links){
            if (link.is_close(epsilon)){
                is_link_close = true;
                link.set_should_move(false);
            }
            else{
                is_link_close = false;
                link.set_should_move(true);
            }
            is_path_complete = is_path_complete && is_link_close;
        }
        return is_path_complete;
    }
    bool SLERP::move(Record &cur, const Record &fin, double dtau){
        bool was_moved = false;
        size_t nb = simple::BasePolymer::nb();
        std::vector<Omega> links;
        links.resize(nb);
        std::vector<simple::Bond>& ini_bonds
            = cur.modifiable_state().polymers.at(0).bonds;
        const std::vector<simple::Bond>& fin_bonds
            = fin.state().get_polymers().at(0).get_bonds();
        /* initialize a vector of links to propagate path */
        for (size_t i = 0; i < nb; ++i){
            links.at(i) = Omega(ini_bonds.at(i), fin_bonds.at(i));
        }
        /* while checking, will mark which links should be moved */
        if (!_is_path_complete(links, _epsilon)){
            for (size_t i = 0; i < nb; ++i){
                if (links.at(i).should_move()) {
                    links.at(i).update(links.at(i).at(dtau));
                }
            }
            was_moved = true;
            _update_record_from_links(cur, links);
            psi1.value = links.at(0).psi();
            psi2.value = links.at(1).psi();
        }
        else {
            was_moved = false;
        }
        return was_moved;
    }
    void SLERP::escape(Record &cur, const Record &fin, double param){
        fprintf(stderr, "%s\n", "geodesic::SLERP::escape(): not yet implemented");
    }
    /***************************************************************************
    *                 GEODESIC SHORT-STEP PATH COMPUTER
    ***************************************************************************/
    ShortStep::ShortStep(
        Potential* polymer_potential,
        Potential* solvent_potential,
        Potential* inter_potential,
        double epsilon) :
    SLERP(polymer_potential, solvent_potential, inter_potential, epsilon),
    theta1(1, true, true),
    theta2(2, true, true),
    delta_theta1(1, vector(0.0, 0.0, 0.0), false, false, true, true),
    delta_theta2(2, vector(0.0, 0.0, 0.0), false, false, true, true){}
    ShortStep::~ShortStep(){};
    ShortStep::Omega::Omega() :
        SLERP::Omega::Omega(),
        _theta_computed(false),
        /**> nhat = cur x fin / |cur x fin|*/
        _nhat(vector(0.0, 0.0, 0.0)),
        /**> uhat = 1/sinPsi (fin - cur * cosPsi)*/
        _uhat(vector(0.0, 0.0, 0.0)){}
    ShortStep::Omega::Omega(const simple::Bond& cur, const simple::Bond& fin) :
        SLERP::Omega::Omega(cur, fin),
        _theta_computed(false),
        /**> nhat = cur x fin / |cur x fin|*/
        _nhat(divide(cross(_cur, _fin), _sinPsi)),
        /**> uhat = 1/sinPsi (fin - cur * cosPsi)*/
        _uhat(divide(subtract(_fin, multiply(_cur, _cosPsi)), _sinPsi)){}
    ShortStep::Omega::~Omega(){};
    Vector ShortStep::Omega::at(double dtau) const {
        double dtheta = _dtheta(dtau);
        Vector in_plane = multiply(SLERP::Omega::at(dtau), cos(dtheta));
        Vector out_plane = multiply(_nhat, sin(dtheta));
        return add(in_plane, out_plane);
    }
    void ShortStep::Omega::_recompute_nhat(){
        _nhat = divide(cross(_cur, _fin), _sinPsi);
    }
    void ShortStep::Omega::_recompute_uhat(){
        _uhat = divide(subtract(_fin, multiply(_cur, _cosPsi)), _sinPsi);
    }
    void ShortStep::Omega::update(Vector cur) {
        SLERP::Omega::update(cur);
        _recompute_nhat();
        _recompute_uhat();
    }
    void ShortStep::Omega::update(const simple::Bond& bond) {
        update(bond.position);
    }
    void ShortStep::_update_bond_from_link(simple::Bond& bond,
        const Omega& omega){
        bond.position = omega.cur();
    }
    void ShortStep::_update_record_from_links(Record& record,
        const std::vector<Omega>& links){
        simple::BondPolymer& polymer = record.modifiable_state().polymers.at(0);
        for(size_t i = 0; i < polymer.nb(); ++i){
            _update_bond_from_link(polymer.bonds.at(i), links.at(i));
        }
    }
    void ShortStep::_update_links_from_record(std::vector< Omega>& links,
        const Record& record){
        simple::BondPolymer& polymer = record.state().polymers.at(0);
        for(size_t i = 0; i < links.size(); ++i){
            links.at(i).update(polymer.bonds.at(i));
        }
    }
    bool ShortStep::_is_path_complete(std::vector< Omega>& links,
        double epsilon){
        bool is_path_complete = true;
        bool is_link_close;
        for(Omega& link : links){
            if (link.is_close(epsilon)){
                is_link_close = true;
                link.set_should_move(false);
            }
            else{
                is_link_close = false;
                link.set_should_move(true);
            }
            is_path_complete = is_path_complete && is_link_close;
        }
        return is_path_complete;
    }
    void ShortStep::_compute_thetas(std::vector<Omega>& links){
        if (links.size() == 2){
            Omega& omega1 = links.at(0);
            Omega& omega2 = links.at(1);
            double num1, num2, den1, den2;
            num1 = dot(omega1.uhat(), omega2.cur());
            num2 = dot(omega2.uhat(), omega1.cur());
            den1 = dot(omega1.nhat(), omega2.cur());
            den2 = dot(omega2.nhat(), omega1.cur());
            omega1.set_theta(- num1 / den1 * omega1.psi());
            omega2.set_theta(- num2 / den2 * omega2.psi());
            if ((abs(den1) < 0.001) || (abs(den2) < 0.001)){
                fprintf(stderr, "%s %f\n", "Numerator 1", num1);
                fprintf(stderr, "%s %f\n", "Denominator 1", den1);
                fprintf(stderr, "%s %f\n", "Theta 1", omega1.theta());
                fprintf(stderr, "%s %f\n", "Numerator 2", num2);
                fprintf(stderr, "%s %f\n", "Denominator 2", den2);
                fprintf(stderr, "%s %f\n", "Theta 2", omega2.theta());
                //fprintf(stderr, "%s: %s\n",
                //    "omega1 cur", vector_to_string(omega1.cur()).c_str());
                //fprintf(stderr, "%s %s\n",
                //    "omega1 fin", vector_to_string(omega1.fin()).c_str());
                //Vector nhat1 = cross(omega1.cur(), omega1.fin());
                //nhat1 = divide(nhat1, norm(nhat1));
                //fprintf(stderr, "%s: %s\n",
                //    "nhat1", vector_to_string(nhat1).c_str());
                //fprintf(stderr, "%s: %s\n",
                //    "omega2 cur", vector_to_string(omega2.cur()).c_str());
                //fprintf(stderr, "%s %s\n",
                //    "omega2 fin", vector_to_string(omega2.fin()).c_str());
                //Vector nhat2 = cross(omega2.cur(), omega2.fin());
                //nhat2 = divide(nhat2, norm(nhat2));
                //fprintf(stderr, "%s: %s\n",
                //    "nhat2", vector_to_string(nhat2).c_str()
                //);
                //fprintf(stderr, "%s %f\n",
                //    "nhat1 dot omega2", dot(nhat1, omega2.cur()));
                //fprintf(stderr, "%s %f\n",
                //    "nhat2 dot omega1", dot(nhat2, omega1.cur()));
            }
        }
        else{
            fprintf(stderr, "%s\n", "ShortStep::_compute_thetas: don't yet know how to deal with nb \\neq 2!");
        }
    }
    bool ShortStep::move(Record& cur, const Record& fin, double dtau){
        bool was_moved = false;
        size_t nb = simple::BasePolymer::nb();
        std::vector<Omega> links;
        links.resize(nb);
        std::vector<simple::Bond>& ini_bonds
            = cur.modifiable_state().polymers.at(0).bonds;
        const std::vector<simple::Bond>& fin_bonds
            = fin.state().get_polymers().at(0).get_bonds();
        /* initialize a vector of links to propagate path */
        for (size_t i = 0; i < nb; ++i){
            links.at(i) = Omega(ini_bonds.at(i), fin_bonds.at(i));
        }
        /* while checking, will mark which links should be moved */
        if (!_is_path_complete(links, _epsilon)){
            _compute_thetas(links);
            for (size_t i = 0; i < nb; ++i){
                if (links.at(i).should_move() && links.at(i).theta_computed()){
                    links.at(i).update(links.at(i).at(dtau));
                }
                else if (!links.at(i).theta_computed()){
                    fprintf(stderr, "%s\n", "ShortStep::_move: theta was not computed");
                }
            }
            was_moved = true;
            _update_record_from_links(cur, links);
            psi1.value = links.at(0).psi();
            delta_psi1.value = links.at(0).psi() * dtau;
            psi2.value = links.at(1).psi();
            delta_psi2.value = links.at(1).psi() * dtau;
            double dtheta;
            if (links.at(0).theta_computed()){
                theta1.value = links.at(0).theta();
                dtheta = links.at(0).theta() * dtau;
                delta_theta1.value = dtheta;
            }
            if (links.at(1).theta_computed()){
                theta2.value = links.at(1).theta();
                dtheta = links.at(1).theta() * dtau;
                delta_theta2.value = dtheta;
            }
        }
        else {
            was_moved = false;
        }
        return was_moved;
    }
    void ShortStep::escape(Record &cur, const Record &fin, double param){
        fprintf(stderr, "%s\n", "geodesic::ShortStep::escape(): not yet implemented");
    }
} // namespace geodesic
