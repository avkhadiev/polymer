// 2018 Artur Avkhadiev
/*! \file geodesic_observables.cpp
*/
#include "../include/geodesic_observables.h"
namespace geodesic{
    /**************************************************************************
    * LENGTH
    **************************************************************************/
    Length::Length(bool print_inst_val, bool e_format) :
        Observable({.full = "Kinetic Length",
                    .abridged = "ell",
                    .latex = "\\ell",
                    .units = "\\sigma"},
                    observable::MAIN_LOOP,
                    {.mean = false,
                     .meansq = false},
                     print_inst_val, e_format),
        increment_sq(0.0){}
    Length::Length() : Length(true, false){}
    void Length::update(const geodesic::Record &last,
                        const geodesic::Record &next){
        simple::BondState last_state = last.state();
        simple::BondState next_state = next.state();
        increment_sq = 0;
        size_t nb = simple::BasePolymer::nb();
        std::vector<simple::Bond> last_bonds
            = last_state.get_polymers().at(0).get_bonds();
        std::vector<simple::Bond> next_bonds
            = next_state.get_polymers().at(0).get_bonds();
        Vector omega1, omega2;
        for (size_t i = 0; i < nb; ++i){
            omega1 = last_bonds.at(i).position;
            omega2 = next_bonds.at(i).position;
            double delta;
            if (omega1 == omega2) {
                delta = 0.0;
            }
            else {
                delta = acos(dot(omega1, omega2)/(norm(omega1) * norm(omega2)));
            }
            if (std::isnan(delta)){
                fprintf(stderr, "%s\n%s: %s\n%s: %s\n",
                    "Length increment is NaN",
                    "omega1", vector_to_string(omega1).c_str(),
                    "omega2", vector_to_string(omega2).c_str());
                fprintf(stderr, "vectors equal: %d\n", omega1 == omega2);
            }
            increment_sq += pow(delta, 2.0);
        }
        value += sqrt(increment_sq);
    }
    /**************************************************************************
    * LINK OBSERVABLE
    **************************************************************************/
    LinkObservable::LinkObservable(observable::name_t name,
        size_t link_number,
        observable::update_time_t update_time,
        bool calculate_mean,
        bool calculate_error,
        bool print_inst_val,
        bool e_format) :
        Observable( name,
                    update_time,
                    {.mean = calculate_mean || calculate_error,
                     .meansq = calculate_error},
                     print_inst_val, e_format),
        _link_number(link_number)
    {
        amend_names();
    }
    void LinkObservable::amend_names() {
        std::string link_number = " (k=" + std::to_string(_link_number) + ")";
        _name.full += link_number;
        _name.abridged += link_number;
        _name.latex += link_number;
    }
    /**************************************************************************
    * OMEGA PROJ
    **************************************************************************/
    OmegaProj::OmegaProj(
         Vector ini,
         Vector fin,
         size_t link_number,
         bool calculate_mean,
         bool calculate_error,
         bool print_inst_val,
         bool e_format
    ) :
        LinkObservable({.full = "Link Projection",
                    .abridged = "omega_proj",
                    .latex = "\\hat{\\Omega}_{k} \\cdot \\hat{n}_{k}",
                    .units = "1"},
                    link_number,
                    observable::MAIN_LOOP,
                    calculate_mean, calculate_error,
                    print_inst_val, e_format),
        _ini(ini),
        _fin(fin),
        _nhat(vector(0.0,0.0,0.0))
        {
            _update_nhat();
        }
    void OmegaProj::update(Vector current){
        value = dot(current, _nhat);
    }
    void OmegaProj::_update_nhat(){
        if (normsq(cross(_ini, _fin)) < 0.00001){
            _nhat = vector(0.0, 0.0, 0.0);
        }
        else {
            _nhat = divide(cross(_ini, _fin), norm(cross(_ini, _fin)));
        }
    }
    void OmegaProj::set_ini(Vector ini){
        _ini = ini;
        _update_nhat();
    }
    void OmegaProj::set_fin(Vector fin){
        _fin = fin;
        _update_nhat();
    }
    OmegaProj::OmegaProj() :
        OmegaProj(vector(0.0,0.0,0.0), vector(0.0,0.0,0.0),
                0, false, false, false, false){}
    void OmegaProj::update(const geodesic::Record &current){
        Vector vec = current.state().get_polymers().at(0).get_bonds().at(_link_number - 1).position;
        update(vec);
    }
    /**************************************************************************
    * PSI
    **************************************************************************/
    Psi::Psi(size_t link_number, Vector fin,
            bool print_inst_val, bool e_format) :
        LinkObservable({.full = "Remaining angle",
                    .abridged = "psi",
                    .latex = "\\psi_{k}",
                    .units = "\\mathrm{rad}"},
                    link_number,
                    observable::FORCE_LOOP,
                    // calculate_mean, calculate_error
                    false, false,
                    print_inst_val, e_format),
        _fin(fin),
        _finnorm(norm(_fin)){}
    Psi::Psi() : Psi(0, vector(0.0, 0.0, 0.0), false, false){}
    void Psi::set_fin(Vector fin){
        _fin = fin;
        _finnorm = norm(_fin);
    }
    void Psi::update(const geodesic::Record &current){
         simple::BondState state = current.state();
         std::vector<simple::Bond> bonds
             = state.get_polymers().at(0).get_bonds();
         Vector omega = bonds.at(_link_number-1).position;
         double sprod = dot(omega, _fin);
         double norml = norm(omega) * _finnorm;
         value = acos(sprod/norml);
         if (std::isnan(value)){
             fprintf(stderr, "%s %f\n", "bad cosine", sprod/norml);
         }
    }
    /**************************************************************************
    * DELTA PSI
    **************************************************************************/
    DeltaPsi::DeltaPsi(size_t link_number,
        bool calculate_mean, bool calculate_error,
        bool print_inst_val, bool e_format) :
        LinkObservable({.full = "In-plane angular step",
                    .abridged = "delta psi",
                    .latex = "\\Delta \\psi_{k}",
                    .units = "\\mathrm{rad}"},
                    link_number,
                    observable::FORCE_LOOP,
                    calculate_mean, calculate_error,
                    print_inst_val, e_format){}
    DeltaPsi::DeltaPsi() : DeltaPsi(0, false, false, false, false){}
    /**************************************************************************
    * THETA
    **************************************************************************/
    Theta::Theta(size_t link_number, bool print_inst_val, bool e_format) :
        LinkObservable({.full = "Out-of-plane parameter",
                    .abridged = "theta",
                    .latex = "\\theta_{k}",
                    .units = "\\mathrm{rad}"},
                    link_number,
                    observable::FORCE_LOOP,
                    // calculate_mean, calculate_error
                    false, false,
                    print_inst_val, e_format){}
    Theta::Theta() : Theta(0, false, false){}
    /**************************************************************************
    * DELTA THETA
    **************************************************************************/
    DeltaTheta::DeltaTheta(size_t link_number,
        Vector fin,
        bool calculate_mean, bool calculate_error,
        bool print_inst_val, bool e_format) :
        LinkObservable({.full = "Out-of-plane angular step",
                    .abridged = "delta theta",
                    .latex = "\\delta \\theta_{k}",
                    .units = "\\mathrm{rad}"},
                    link_number,
                    observable::FORCE_LOOP,
                    calculate_mean, calculate_error,
                    print_inst_val, e_format),
        _fin(fin),
        _finnorm(norm(_fin)){}
    DeltaTheta::DeltaTheta() :
        DeltaTheta(0, vector(0.0, 0.0, 0.0), false, false, false, false){}
    void DeltaTheta::set_fin(Vector fin){
        _fin = fin;
        _finnorm = norm(_fin);
    }
    void DeltaTheta::update(const geodesic::Record &last,
                            const geodesic::Record &next){
        simple::BondState last_state = last.state();
        simple::BondState next_state = next.state();
        std::vector<simple::Bond> last_bonds
        = last_state.get_polymers().at(0).get_bonds();
        std::vector<simple::Bond> next_bonds
        = next_state.get_polymers().at(0).get_bonds();
        Vector last_omega = last_bonds.at(_link_number-1).position;
        Vector next_omega = next_bonds.at(_link_number-1).position;
        Vector vprod = cross(last_omega, _fin);
        double norml = norm(last_omega) * norm(next_omega) * _finnorm;
        value = asin(dot(next_omega, vprod) / norml);
    }
    /**************************************************************************
    * DELTA PHI
    **************************************************************************/
    DeltaPhi::DeltaPhi(size_t link_number,
        bool calculate_mean, bool calculate_error,
        bool print_inst_val, bool e_format) :
        LinkObservable({.full = "Total angular step",
                    .abridged = "delta phi",
                    .latex = "\\delta \\phi_{k}",
                    .units = "\\mathrm{rad}"},
                    link_number,
                    observable::FORCE_LOOP,
                    calculate_mean, calculate_error,
                    print_inst_val, e_format){}
    DeltaPhi::DeltaPhi() : DeltaPhi(0, false, false, false, false){}
}  // namespace geodesic
