// 2018 Artur Avkhadiev
/*! \file geodesic_path_computer.h
*/
#ifndef GEODESIC_PATH_COMPUTER_H
#define GEODESIC_PATH_COMPUTER_H
#include <vector>
#include <sys/stat.h>
#include <string>
#include <cassert>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <cmath>
#include "geodesic_record.h"
#include "geodesic_path.h"
#include "geodesic_manager.h"
#include "observable.h"
#include "potential.h"
#include "force_updater.h"
#include "shove_integrator.h"               /** For the new algorithm */
#include "observable_container.h"
#include "general_observables.h"            /** Potential Energy    */
#include "geodesic_observables.h"           /** angles, steps, etc. */
namespace geodesic{
    /***************************************************************************
    *                 GEODESIC PATH COMPUTER (BASE CLASS)
    ***************************************************************************/
    /**
    * An interface class providing methods for incremental propagation of a path
    * Given a reference to a Record and a constant reference to the final Record
    * these methods help to propagate the Record toward the final Record,
    * one step at a time
    */
    class PathComputer {
    public:
        PathComputer(Potential* polymer_potential,
            Potential* solvent_potential,
            Potential* inter_potential,
            double epsilon = 0.000001);
        ~PathComputer();
        double epsilon() const {return _epsilon;};
        void set_epsilon(double epsilon) {_epsilon = epsilon;};
        void update_PE(Record& rec);
        /** if current record is close enough to final record, return false
        * if not, return true and move current record dtau towards final record
        */
        virtual bool move(Record &cur, const Record &fin, double dr) = 0;
        /**> TODO */
        virtual void escape(Record &cur, const Record &fin, double param) = 0;
    protected:
        bool _is_path_complete(Record &cur, const Record &fin, double epsilon);
        ForceUpdater _fupd;                /** to compute PE    */
        double _epsilon;                   /** measure how close records are  */
    };
    /***************************************************************************
    *                       GEODESIC SLERP PATH COMPUTER
    ***************************************************************************/
    class SLERP : public PathComputer {
    public:
        SLERP(Potential* polymer_potential,
            Potential* solvent_potential,
            Potential* inter_potential,
            double epsilon = 0.000001);
        ~SLERP();
        class Omega {
        public:
            Omega();
            Omega(const simple::Bond& cur, const simple::Bond& fin);
            ~Omega();
            Vector cur() const {return _cur;};
            Vector fin() const {return _fin;};
            double cosPsi() const {return _cosPsi;};
            double sinPsi() const {return _sinPsi;};
            double psi() const {return _psi;};
            double a(double dtau) const {return _a(dtau);};
            double b(double dtau) const {return _b(dtau);};
            virtual void update(const simple::Bond& bond);
            virtual void update(Vector cur);
            void set_should_move(bool should_move){_should_move = should_move;};
            bool should_move() const {return _should_move;};
            bool is_close(double epsilon) const;
            virtual Vector at(double dtau) const;
        protected:
            bool _should_move;
            Vector _cur;    /**> current orientation                        */
            Vector _fin;    /**> final orientation                          */
            double _cosPsi;
            double _sinPsi;
            double _psi;    /**> SLERP angle between ini and fin            */
            void _recompute_angles();
        private:
            //*> sin (Psi(1 - tau))/sinPsi */
            double _a(double dtau) const;
            //*> sin (Psi tau) / sinPsi    */
            double _b(double dtau) const;
        };
        virtual bool move(Record &cur, const Record &fin, double dtau);
        /**> TODO */
        virtual void escape(Record &cur, const Record &fin, double param);
    public:
        Psi psi1, psi2;                  /** remaining angle */
        DeltaPsi delta_psi1, delta_psi2; /** in-plane step   */
    private:
        bool _is_path_complete(std::vector< Omega>& links, double epsilon);
        void _update_bond_from_link(simple::Bond& bond, const Omega& omega);
        void _update_record_from_links(Record& record,
            const std::vector<Omega>& links);
        void _update_links_from_record(std::vector<Omega>& links,
            const Record& record);
    };
    /***************************************************************************
    *             GEODESIC SHORT-STEP PATH COMPUTER (Enhanced SLERP)
    ***************************************************************************/
    class ShortStep : public SLERP {
    public:
        ShortStep(Potential* polymer_potential,
            Potential* solvent_potential,
            Potential* inter_potential,
            double epsilon = 0.000001);
        ~ShortStep();
        class Omega : public SLERP::Omega {
        public:
            Omega();
            Omega(const simple::Bond& cur, const simple::Bond& fin);
            ~Omega();
            void set_theta(double theta){
                _theta = theta;
                _theta_computed = true;
            };
            bool theta_computed() const {return _theta_computed;};
            virtual void update(const simple::Bond& bond);
            virtual void update(Vector cur);
            Vector nhat() const {return _nhat;};
            Vector uhat() const {return _uhat;};
            double theta() const {return _theta;};
            virtual Vector at(double dtau) const;
        protected:
            bool _theta_computed;
            Vector _SLERP(double dtau) const {return SLERP::Omega::at(dtau);};
            Vector _nhat;   /**> normal to the plane defined by cur and fin */
            Vector _uhat;   /**> lin velocity unit vector when omega = cur  */
            double _theta;  /**> ~ to step out of plane of motion           */
            void _recompute_nhat();
            void _recompute_uhat();
        private:
            double _dtheta(double dtau) const {return _theta * dtau;};
        };
        virtual bool move(Record &cur, const Record &fin, double dtau);
        /**> TODO */
        virtual void escape(Record &cur, const Record &fin, double param);
    public:
        Theta theta1, theta2;                    /** out-of-plane parameter */
        DeltaTheta delta_theta1, delta_theta2;   /** out-of-plane step      */
        DeltaPhi delta_phi1, delta_phi2;         /** totals step            */
    private:
        bool _is_path_complete(std::vector< Omega>& links, double epsilon);
        void _update_bond_from_link(simple::Bond& bond, const Omega& omega);
        void _update_record_from_links(Record& record,
            const std::vector< Omega>& links);
        void _update_links_from_record(std::vector< Omega>& links,
            const Record& record);
        /** that's the meat of the algorithm! */
        void _compute_thetas(std::vector< Omega>& links);
    };
    /***************************************************************************
    *                          GEODESIC SHOVE COMPUTER
    ***************************************************************************/
    class SHOVE : public PathComputer {
    public:
        SHOVE(Potential* polymer_potential,
            Potential* solvent_potential,
            Potential* inter_potential,
            double tol = pow(10, -8.0),
            double epsilon = 0.000001);
        ~SHOVE();
        virtual bool move(Record &cur, const Record &fin, double dr);
        /**> TODO */
        virtual void escape(Record &cur, const Record &fin, double param);
        Psi psi1, psi2;                             /** remaining angle */
    private:
        ShoveIntegrator _integrator;
        void _assign_velocities(Record &cur, const Record &fin, double dr);
    };
    /***************************************************************************
    *                          GEODESIC PLERP COMPUTER
    ***************************************************************************/
    class PLERP : public PathComputer {
    public:
        PLERP(Potential* polymer_potential,
            Potential* solvent_potential,
            Potential* inter_potential,
            double epsilon = 0.000001);
        ~PLERP();
        virtual bool move(Record &cur, const Record &fin, double dr);
        /**> TODO */
        virtual void escape(Record &cur, const Record &fin, double param);
    private:
        // stores the change in configuration space vector to be applied
        std::vector<Vector> deltaR;
        void move_at_once(Record &cur, const Record &fin, double dr);
        // initial step: ignores constraints, modifies deltaR
        void _move0(const Record &cur, const Record &fin, double dR);
        // projections (1 - del C del C), modifies deltaR and &cur
        // violates constraints to second order
        void _move1(Record &cur);
        // modifies deltaR and &cur to restore constraints to second order
        void _move2(Record &cur);
        void move_by_bond(Record &cur, const Record &fin, double dr);
        Vector _rb(Vector r1, Vector r2);
        void _project_dR(size_t jbond, Vector rb);
        void _correct(size_t jbond, Vector r1, Vector r2, const Record &cur);
    };
} // namespace geodesic
#endif
