// 2018 Artur Avkhadiev
/*! \file geodesic_observables.h
*/
#ifndef GEODESIC_OBSERVABLES_H
#define GEODESIC_OBSERVABLES_H
#include "observable.h"
#include "geodesic_record.h"
namespace geodesic{
    /**
    * base class for lengh, updated by the main loop
    */
    class Length: public Observable {
    protected:
    public:
        virtual bool update_method_specified() const {return false;};
        virtual void update(const simple::AtomState& state){};
        void update(const geodesic::Record& last, const geodesic::Record& next);
        double increment_sq;        // accumulator for the square of increment
        ~Length(){};
        Length(
            bool print_inst_val,
            bool e_format = false   // use precision notation when writing out?
        );
        Length();
    };
    /**
    * abstract class for observables measured per link
    */
    class LinkObservable: public Observable {
    protected:
        size_t _link_number;
        void amend_names();             // adds link number to observable names
    public:
        size_t link_number() const {return _link_number;};
        virtual bool update_method_specified() const {return false;};
        virtual void update(const simple::AtomState& state) = 0;
        ~LinkObservable(){};
        LinkObservable(observable::name_t name,
            size_t link_number,
            observable::update_time_t update_time,
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,
            bool e_format = false);
    };
    /**
    * Projection of the current link vector onto nhat --- the normal to the
    * plane connecting the endpoints of the link. Computed in the main loop
    */
    class OmegaProj : public LinkObservable {
    protected:
        Vector _ini;
        Vector _fin;
        Vector _nhat;
        void _update_nhat();
    public:
        Vector ini() const {return _ini;};
        Vector fin() const {return _fin;};
        Vector nhat() const {return _nhat;};
        void set_ini(Vector ini);
        void set_fin(Vector fin);
        virtual bool update_method_specified() const {return false;};
        virtual void update(const simple::AtomState& state){};
        void update(const geodesic::Record& current);
        void update(Vector current);
        ~OmegaProj(){};
        OmegaProj();
        OmegaProj(
             Vector ini,
             Vector fin,
             size_t link_number,
             bool calculate_mean,
             bool calculate_error,
             bool print_inst_val,   // print out inst value in log & state?
             bool e_format = false  // use precision notation when writing out?
        );
    };
    /**
    * Remaining angle Psi, measured for a given link,
    *  updated by the force_loop(computer)
    */
    class Psi: public LinkObservable {
    private:
        Vector _fin;
        double _finnorm;
    public:
        Vector fin() const {return _fin;};
        void set_fin(Vector fin);
        virtual bool update_method_specified() const {return false;};
        virtual void update(const simple::AtomState& state){};
        void update(const geodesic::Record& current);
        ~Psi(){};
        Psi();
        Psi( size_t link_number,
             Vector fin,
             bool print_inst_val,   // print out inst value in log & state?
             bool e_format = false);// use precision notation when writing out?
    };
    /**
    * In-plane step deltaPsi, updated by the force_loop(computer)
    */
    class DeltaPsi: public LinkObservable {
    private:
    public:
        virtual bool update_method_specified() const {return false;};
        virtual void update(const simple::AtomState& state){};
        ~DeltaPsi(){};
        DeltaPsi(
            size_t link_number,
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing out?
        );
        DeltaPsi();
    };
    /**
    * Parameter theta, measured for a given link,
    *  updated by the force_loop(computer)
    */
    class Theta: public LinkObservable {
    public:
        virtual bool update_method_specified() const {return false;};
        virtual void update(const simple::AtomState& state){};
        ~Theta(){};
        Theta(
            size_t link_number,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing out?
        );
        Theta();
    };
    /**
    * Out-of-plane step deltaTheta, updated by the force_loop(computer)
    */
    class DeltaTheta: public LinkObservable {
    private:
        Vector _fin;
        double _finnorm;
    public:
        Vector fin() const {return _fin;};
        void set_fin(Vector fin);
        virtual bool update_method_specified() const {return false;};
        virtual void update(const simple::AtomState& state){};
        void update(const geodesic::Record& last, const geodesic::Record& next);
        ~DeltaTheta(){};
        DeltaTheta(
            size_t link_number,
            Vector fin,
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing out?
        );
        DeltaTheta();
    };
    /**
    * angular step deltaPhi, updated by the force_loop(computer)
    */
    class DeltaPhi: public LinkObservable {
    public:
        virtual bool update_method_specified() const {return false;};
        virtual void update(const simple::AtomState& state){};
        ~DeltaPhi(){};
        DeltaPhi(
            size_t link_number,
            bool calculate_mean,
            bool calculate_error,
            bool print_inst_val,   // print out inst value in log & state?
            bool e_format = false  // use precision notation when writing out?
        );
        DeltaPhi();
    };
}   // namespace geodesic
#endif
