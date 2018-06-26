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
        simple::AtomState atom_state = rec.atom_state();
        rec.set_pe(_fupd.calc_pot_energy(atom_state));
    }
    bool PathComputer::_is_path_complete(Record &cur, const Record &fin,
        double epsilon){
        simple::AtomState& cur_state = cur.atom_state();
        const simple::AtomState& fin_state = fin.atom_state();
        bool is_path_complete = true;
	double diffsq = 0.0;
        // check (single polymer) positions
        simple::AtomPolymer& cur_polymer = cur_state.polymers.at(0);
        const simple::AtomPolymer fin_polymer = fin_state.polymers.at(0);
        Vector cur_pos;
        Vector fin_pos;
        for (size_t i = 0; i < simple::BasePolymer::nb() + 1; ++i){
            cur_pos = cur_polymer.atoms.at(i).position;
            fin_pos = fin_polymer.atoms.at(i).position;
            diffsq += normsq(subtract(cur_pos, fin_pos)); 
        }
	double diff = sqrt(diffsq);
	if (diff > epsilon){
	    is_path_complete = false;
	}
        //if(is_path_complete){
        //    fprintf(stderr, "%s:\n%s\n",
        //        "TERMINATED PATH",
        //        cur.atom_state().to_string(true, false).c_str());
        //}
        return is_path_complete;
    }
    /***************************************************************************
    *                       GEODESIC SLERP PATH COMPUTER
    ***************************************************************************/
    SLERP::SLERP(Potential* polymer_potential,
        Potential* solvent_potential,
        Potential* inter_potential,
        double epsilon) :
    PathComputer(polymer_potential,solvent_potential,inter_potential,epsilon){}
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
        simple::BondPolymer& polymer = record.bond_state().polymers.at(0);
        for(size_t i = 0; i < polymer.nb(); ++i){
            _update_bond_from_link(polymer.bonds.at(i), links.at(i));
        }
    }
    void SLERP::_update_links_from_record(std::vector<Omega>& links,
        const Record& record){
        const simple::BondPolymer& polymer = record.bond_state().polymers.at(0);
        for(size_t i = 0; i < links.size(); ++i){
            links.at(i).update(polymer.bonds.at(i));
        }
    }
    bool SLERP::Omega::is_close(double epsilon) const {
        return abs(1.0 - _cosPsi) < pow(epsilon, 2.0);
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
    bool SLERP::move(Record &cur, const Record &fin, double dr){
        bool keep_going = false;
        cur.atom_state().update(cur.bond_state());
        // CALCULATE DTAU
        simple::AtomState& cur_state = cur.atom_state();
        const simple::AtomState& fin_state = fin.atom_state();
        // working with a single polymer
        simple::AtomPolymer& cur_polymer = cur_state.polymers.at(0);
        const simple::AtomPolymer fin_polymer = fin_state.polymers.at(0);
        Vector cur_pos, fin_pos, cur_direction;
        double r_normsq = 0.0;
        // assign unit velocities to each atoms in the direction towards final
        // position, and find the atom farthest from its destination
        // zero all forces
        for (size_t i = 0; i < simple::BasePolymer::nb() + 1; ++i){
            cur_pos = cur_polymer.atoms.at(i).position;
            fin_pos = fin_polymer.atoms.at(i).position;
            cur_direction = subtract(fin_pos, cur_pos);
            r_normsq += normsq(cur_direction);
        }
        double r_norm = sqrt(r_normsq);
        if (dr > r_norm) {
            dr = r_norm;
        }
        double dtau = dr / r_norm;
        // MOVE LINKS
        size_t nb = simple::BasePolymer::nb();
        std::vector<Omega> links;
        links.resize(nb);
        std::vector<simple::Bond>& cur_bonds
            = cur.bond_state().polymers.at(0).bonds;
        const std::vector<simple::Bond>& fin_bonds
            = fin.bond_state().get_polymers().at(0).get_bonds();
        /* initialize a vector of links to propagate path */
        for (size_t i = 0; i < nb; ++i){
            links.at(i) = Omega(cur_bonds.at(i), fin_bonds.at(i));
            links.at(i).update(links.at(i).at(dtau));
        }
        _update_record_from_links(cur, links);
        cur.atom_state().update(cur.bond_state());
        if (dtau < 1.0){
            keep_going = true;
        }
        else {
            keep_going = false;
            //fprintf(stderr, "%s:\n%s\n",
            //    "TERMINATED PATH",
            //    cur.atom_state().to_string(true, false).c_str());
        }
        return keep_going;
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
    SLERP(polymer_potential, solvent_potential, inter_potential, epsilon){}
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
        simple::BondPolymer& polymer = record.bond_state().polymers.at(0);
        for(size_t i = 0; i < polymer.nb(); ++i){
            _update_bond_from_link(polymer.bonds.at(i), links.at(i));
        }
    }
    void ShortStep::_update_links_from_record(std::vector< Omega>& links,
        const Record& record){
        const simple::BondPolymer& polymer = record.bond_state().polymers.at(0);
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
                if (DEBUG){
                    fprintf(stderr, "%s %f\n", "Numerator 1", num1);
                    fprintf(stderr, "%s %f\n", "Denominator 1", den1);
                    fprintf(stderr, "%s %f\n", "Theta 1", omega1.theta());
                    fprintf(stderr, "%s %f\n", "Numerator 2", num2);
                    fprintf(stderr, "%s %f\n", "Denominator 2", den2);
                    fprintf(stderr, "%s %f\n", "Theta 2", omega2.theta());
                }
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
            = cur.bond_state().polymers.at(0).bonds;
        const std::vector<simple::Bond>& fin_bonds
            = fin.bond_state().get_polymers().at(0).get_bonds();
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
            double dtheta;
            if (links.at(0).theta_computed()){
                dtheta = links.at(0).theta() * dtau;
            }
            if (links.at(1).theta_computed()){
                dtheta = links.at(1).theta() * dtau;
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
    SHOVE::SHOVE(Potential* polymer_potential,
        Potential* solvent_potential,
        Potential* inter_potential,
        double tol,
        double epsilon) :
        PathComputer(   polymer_potential,
                        solvent_potential,
                        inter_potential,
                        epsilon ),
        _integrator(tol){
    }
    SHOVE::~SHOVE(){};
    void SHOVE::_assign_velocities(
        Record &cur,
        const Record &fin,
        double dR) {
        size_t natoms = simple::BasePolymer::nb() + 1;
        std::vector<simple::Atom>& cur_atoms
            = cur.atom_state().polymers.at(0).atoms;
        const std::vector<simple::Atom> fin_atoms
            = fin.atom_state().get_polymers().at(0).get_atoms();
        Vector cur_dir;
        Vector cur_pos, fin_pos;
        double R_norm_sq = 0.0;
        for (size_t j = 0; j < natoms; ++j){
            cur_pos = cur_atoms.at(j).position;
            fin_pos = fin_atoms.at(j).position;
            cur_dir = subtract(fin_pos, cur_pos);
            cur_atoms.at(j).velocity = cur_dir;
            R_norm_sq += normsq(cur_dir);
        }
        double R_norm = sqrt(R_norm_sq);
        if (dR > R_norm) {
            dR = R_norm;
        }
        double weight = dR / R_norm;
        for (size_t j = 0; j < natoms; ++j){
            cur_atoms.at(j).velocity
                = multiply(cur_atoms.at(j).velocity, weight);
        }
        Vector rb;      // vector whose outerpoduct w itself gives the matrix
        Vector dr;      // difference of deltas of j+1 atom and jth atom
        double x, y, z;
        for (size_t j = 0; j < natoms - 1; ++j){
            // get current dr = delta r_{j+1} - delta r_{j}
            // will be continuously updated
            dr = subtract(cur_atoms.at(j+1).velocity, cur_atoms.at(j).velocity);
            // this unit bond vector is fixed in time
            rb = subtract(cur_atoms.at(j+1).position, cur_atoms.at(j).position);
            rb = divide(rb, norm(rb));
            x = rb.x * rb.x * dr.x + rb.x * rb.y * dr.y + rb.x * rb.z * dr.z;
            y = rb.y * rb.x * dr.x + rb.y * rb.y * dr.y + rb.y * rb.z * dr.z;
            z = rb.z * rb.x * dr.x + rb.z * rb.y * dr.y + rb.z * rb.z * dr.z;
            dr = vector(x, y, z);
            cur_atoms.at(j).velocity   += multiply(dr, 0.5);
            cur_atoms.at(j+1).velocity -= multiply(dr, 0.5);
            //fprintf(stderr, "%s:%s\n",
            //    "Delta R",
            //    vector_to_string(deltaR.at(j)).c_str());
        }
    }
    bool SHOVE::move(Record &cur, const Record &fin, double dr){
        // this function will record which atoms need to be moved at all
        // in the vector _move_atom
        bool keep_going = false;
        cur.bond_state().update(cur.atom_state());
        // compute new values of psi
        std::vector<simple::Bond> ini_bonds
            = cur.bond_state().get_polymers().at(0).get_bonds();
        std::vector<simple::Bond> fin_bonds
            = fin.bond_state().get_polymers().at(0).get_bonds();
        Vector cur_omega1, cur_omega2, fin_omega1, fin_omega2;
        cur_omega1 = ini_bonds.at(0).position;
        fin_omega1 = fin_bonds.at(0).position;
        cur_omega2 = ini_bonds.at(1).position;
        fin_omega2 = fin_bonds.at(1).position;
        if (!_is_path_complete(cur, fin, _epsilon)){
            _assign_velocities(cur, fin, dr);
            _integrator.move(1.0, cur.atom_state());
            keep_going = true;
        }
        return keep_going;
    }
    void SHOVE::escape(Record &cur, const Record &fin, double param){
        fprintf(stderr, "%s\n", "geodesic::SLERP::escape(): not yet implemented");
    }
    PLERP::PLERP(Potential* polymer_potential,
        Potential* solvent_potential,
        Potential* inter_potential,
        double epsilon) :
        PathComputer(   polymer_potential,
                        solvent_potential,
                        inter_potential,
                        epsilon ),
        deltaR() {
        deltaR.resize(simple::BasePolymer::nb() + 1);
    }
    PLERP::~PLERP(){}
    bool PLERP::move(Record &cur, const Record &fin, double dR){
        bool keep_going = false;
        if (!_is_path_complete(cur, fin, _epsilon)){
            move_at_once(cur, fin, dR);
            keep_going = true;
        }
        cur.bond_state().update(cur.atom_state());
        return keep_going;
    }
    void PLERP::move_at_once(Record &cur, const Record &fin, double dR){
        _move0(cur, fin, dR);
        _move1(cur);
        _move2(cur);
    }
    void PLERP::move_by_bond(Record &cur, const Record &fin, double dR){
        // gets naive displacement for dR
        _move0(cur, fin, dR);
        size_t natoms = simple::BasePolymer::nb() + 1;
        std::vector<simple::Atom>& cur_atoms
            = cur.atom_state().polymers.at(0).atoms;
        Vector r1, r2;
        for (size_t j = 0; j < natoms - 1; ++j){
            r1 = cur_atoms.at(j).position;
            r2 = cur_atoms.at(j+1).position;
            // apply jth projection to the naive displacement
            // affects atoms j, j+1, modifies del r_j, del r_{j+1}
            _project_dR(j, _rb(r1, r2));
            // calculate the proposed positions j, j+1 using projected naive
            // displacements
            r1 += deltaR.at(j);
            r2 += deltaR.at(j+1);
            // given the proposed positions, correct the displacements yet again
            // to preserve the constratins of the link
            // and its closest neighbors.
            // modifies del r_j, del r_{j+1}
            _correct(j, r1, r2, cur);
            // use the modified displacements to apply the changes
            cur_atoms.at(j).position += deltaR.at(j);
            cur_atoms.at(j+1).position += deltaR.at(j+1);
        }
    }
    void PLERP::_move0(const Record &cur, const Record &fin, double dR){
        size_t natoms = simple::BasePolymer::nb() + 1;
        const std::vector<simple::Atom> cur_atoms
            = cur.atom_state().get_polymers().at(0).get_atoms();
        const std::vector<simple::Atom> fin_atoms
            = fin.atom_state().get_polymers().at(0).get_atoms();
        Vector cur_dir;
        Vector cur_pos, fin_pos;
        double R_norm_sq = 0.0;
        for (size_t j = 0; j < natoms; ++j){
            cur_pos = cur_atoms.at(j).position;
            fin_pos = fin_atoms.at(j).position;
            cur_dir = subtract(fin_pos, cur_pos);
            deltaR.at(j) = cur_dir;
            R_norm_sq += normsq(cur_dir);
        }
        double R_norm = sqrt(R_norm_sq);
        if (dR > R_norm) {
            dR = R_norm;
        }
        double weight = dR / R_norm;
        for (size_t j = 0; j < natoms; ++j){
            deltaR.at(j) = multiply(deltaR.at(j), weight);
            //fprintf(stderr, "%s:%s\n",
            //    "Delta R",
            //    vector_to_string(deltaR.at(j)).c_str());
        }
    }
    Vector PLERP::_rb(Vector r1, Vector r2){
        Vector r12 = subtract(r2, r1);
        return divide(r12, norm(r12));
    }
    void PLERP::_project_dR(size_t j, Vector rb){
        Vector dr = subtract(deltaR.at(j+1), deltaR.at(j));
        double x, y, z;
        x = rb.x * rb.x * dr.x + rb.x * rb.y * dr.y + rb.x * rb.z * dr.z;
        y = rb.y * rb.x * dr.x + rb.y * rb.y * dr.y + rb.y * rb.z * dr.z;
        z = rb.z * rb.x * dr.x + rb.z * rb.y * dr.y + rb.z * rb.z * dr.z;
        dr = vector(x, y, z);
        deltaR.at(j)   += multiply(dr, 0.5);
        deltaR.at(j+1) -= multiply(dr, 0.5);
    }
    void PLERP::_correct(size_t j, Vector r1, Vector r2, const Record &cur){
        const std::vector<simple::Atom> cur_atoms
            = cur.atom_state().get_polymers().at(0).get_atoms();
        int nbonds = simple::BasePolymer::nb();
        double dsq = pow(simple::BasePolymer::d(), 2.0);
        double cr = 0.0;
        double c_last, c_this, c_next;
        Vector del_cr_1 = vector(0.0, 0.0, 0.0);
        Vector del_cr_2 = vector(0.0, 0.0, 0.0);
        Vector r0, r3;
        Vector r_b;
        r_b = subtract(r2, r1);
        c_this = normsq(r_b) - dsq;
        cr += pow(c_this, 2.0);
        del_cr_1 += multiply(r_b, - 4.0 * c_this);
        del_cr_2 += multiply(r_b, 4.0 * c_this);
        if (j != nbonds - 1){
            r3 = cur_atoms.at(j+2).position;
            r_b = subtract(r3, r2);
            c_next = normsq(r_b) - dsq;
            cr += pow(c_next, 2.0);
            del_cr_2 += multiply(r_b, - 4.0 * c_next);
        }
        if (j != 0){
            r0 = cur_atoms.at(j-1).position;
            r_b = subtract(r1, r0);
            c_last = normsq(r_b) - dsq;
            cr += pow(c_last, 2.0);
            del_cr_1 += multiply(r_b, 4.0 * c_last);
        }
        double norm_del_cr_1 = norm(del_cr_1);
        double norm_del_cr_2 = norm(del_cr_2);
        double a = - cr / (norm_del_cr_1 + norm_del_cr_2);
        deltaR.at(j)   += multiply(divide(del_cr_1, norm_del_cr_1), a);
        deltaR.at(j+1) += multiply(divide(del_cr_2, norm_del_cr_2), a);
    }
    void PLERP::_move1(Record &cur){
        size_t natoms = simple::BasePolymer::nb() + 1;
        std::vector<simple::Atom>& cur_atoms
            = cur.atom_state().polymers.at(0).atoms;
        Vector rb;      // vector whose outerpoduct w itself gives the matrix
        Vector dr;      // difference of deltas of j+1 atom and jth atom
        double x, y, z;
        for (size_t j = 0; j < natoms - 1; ++j){
            // get current dr = delta r_{j+1} - delta r_{j}
            // will be continuously updated
            dr = subtract(deltaR.at(j+1), deltaR.at(j));
            // this unit bond vector is fixed in time
            rb = subtract(cur_atoms.at(j+1).position, cur_atoms.at(j).position);
            rb = divide(rb, norm(rb));
            x = rb.x * rb.x * dr.x + rb.x * rb.y * dr.y + rb.x * rb.z * dr.z;
            y = rb.y * rb.x * dr.x + rb.y * rb.y * dr.y + rb.y * rb.z * dr.z;
            z = rb.z * rb.x * dr.x + rb.z * rb.y * dr.y + rb.z * rb.z * dr.z;
            dr = vector(x, y, z);
            deltaR.at(j)   += multiply(dr, 0.5);
            deltaR.at(j+1) -= multiply(dr, 0.5);
            //fprintf(stderr, "%s:%s\n",
            //    "Delta R",
            //    vector_to_string(deltaR.at(j)).c_str());
        }
        for (size_t j = 0; j < natoms; ++j){
            cur_atoms.at(j).position += deltaR.at(j);
        }
    }
    void PLERP::_move2(Record &cur){
        size_t natoms = simple::BasePolymer::nb() + 1;
        double dsq = pow(simple::BasePolymer::d(), 2.0);
        std::vector<simple::Atom>& cur_atoms
            = cur.atom_state().polymers.at(0).atoms;
        Vector rb;
        double cr = 0.0;
        double constraint = 0.0;
        Vector del_crj;
        double delcr_norm_sq = 0.0;
        // calculate first element
        del_crj = vector(0.0, 0.0, 0.0);
        deltaR.at(0) = vector(0.0, 0.0, 0.0);
        rb = subtract(cur_atoms.at(1).position, cur_atoms.at(0).position);
        constraint = normsq(rb) - dsq;
        cr += pow(constraint, 2.0);
        del_crj += multiply(rb, -4.0 * constraint);
        deltaR.at(0) += del_crj;
        delcr_norm_sq += normsq(del_crj);
        for (size_t j = 1; j < natoms - 1; ++j){
            deltaR.at(j) = vector(0.0, 0.0, 0.0);
            del_crj = vector(0.0, 0.0, 0.0);
            del_crj += multiply(rb, 4.0 * constraint);
            rb = subtract(cur_atoms.at(j+1).position, cur_atoms.at(j).position);
            constraint = normsq(rb) - dsq;
            cr += pow(constraint, 2.0);
            del_crj += multiply(rb, -4.0 * constraint);
            deltaR.at(j) += del_crj;
            delcr_norm_sq += normsq(del_crj);
        }
        //fprintf(stderr, "C(r): %f\n", cr);
        //fprintf(stderr, "|del C(r)|^2: %f\n", delcr_norm_sq);
        double weight = -1.0 * cr / delcr_norm_sq;
        // calculate last element
        deltaR.at(natoms-1) = vector(0.0, 0.0, 0.0);
        del_crj = vector(0.0, 0.0, 0.0);
        del_crj += multiply(rb, 4.0 * constraint);
        deltaR.at(natoms-1) += del_crj;
        // calculate and apply the steps
        for (size_t j = 0; j < natoms; ++j){
            deltaR.at(j) = multiply(deltaR.at(j), weight);
            //fprintf(stderr, "%s:%s\n",
            //    "Delta R",
            //    vector_to_string(deltaR.at(j)).c_str());
            cur_atoms.at(j).position += deltaR.at(j);
        }
    }
    void PLERP::escape(Record &cur, const Record &fin, double param){
        fprintf(stderr, "%s\n", "geodesic::PLERP::escape(): not yet implemented");
    }
} // namespace geodesic
