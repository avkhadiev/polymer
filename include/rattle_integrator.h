// 2017 Artur Avkhadiev
/*! \file rattle_integrator.h
*/
#ifndef POLYMER_RATTLE_INTEGRATOR_H
#define POLYMER_RATTLE_INTEGRATOR_H
#include <vector>
#include <cmath>              /* pow */
#include "ljpotential.h"
#include "simulation.h"
#include "molecule.h"
#include "simple_atom.h"
#include "simple_polymer.h"
#include "observable_container.h"
#include "verlet_integrator.h"
class RattleIntegrator :
    public VerletIntegrator {
private:
    int _maxiter;
    double _tol;
    double _rvtol;
    double _tiny;
    double _tol2;     /**> _tol * 2 */
protected:
    int _nb;
    double _dabsq; /**> constant square of bond length for simple polymers */
    double _rm;   /**> constant inverse atomic mass for simple polymers */
    double _inv_timestep;
    // stores inv_timestep in addition to half step and full step
    virtual void _set_timestep(double timestep);
    // pointers to accumulators of observables
    // RATTLE gets its own kinetic energy accumulator, because now that
    // corrections have to be performed on velocities after Verlet full step,
    // it doesn't make sense to calculate kinetic energy from (unconstrained)
    // velocities. Thus these two data members shadow (hide) the ones in the
    // base class
    bool _is_kinetic_energy_acc_set;
    double *_kinetic_energy_acc;
    bool _is_neg_constraint_virial_acc_set;
    double *_neg_constraint_virial_acc;
    // checks whether | r^{i}_{AB}(t + dt) - d_{AB} |^2 < 2 * tol * d^2_{AB}
    bool _is_constraint_within_tol(double dabsq, double difference_of_squares);
    // checks whether r_{AB}(t) * r^{i}_{AB}(t+dt) > tiny * d_{AB}
    // to avoid division by 0
    bool _is_angle_okay(double dabsq, double rr_dot);
    // checks whether | r_{AB}(t+dt) * v_{AB}(t+dt) |^2 < 2 * rvtol * d^2_{AB}
    bool _is_constraint_derivative_within_rvtol(double dabsq, double rv_dot);
    /**
    *   before calling,
    *      call unconstrained verlet: r(t) -> r^0(t + dt), v(t) -> v^0(t + dt)
    *      and save the output molecule as molecule_current_step
    * Then: for each constrained bond r_{AB}, iterative correction:
    *      1. r^i_{AB}(t + dt) -> r^{i+1}_{AB}(t + dt),
    *         v^i_{AB}(t+0.5dt) -> v^{i+1}_{AB}(t+0.5dt), until for all bonds
    *          | r_{i}_{AB}(t + dt) - d_{AB} |^2 < 2 * tol * d^2_{AB},
    *              where d_{AB} = r_{AB}(t) is the constrained bond length,
    *              a factor 2 in front of tol is from Taylor expansion
    *      2. When all bonds are satisfied to within given tolerance, save
    *          the corrected versions of position and velocity
    *          r_{AB}(t + dt) = r^m_{AB}(t + dt),
    *          v_{AB}(t + 0.5dt) = v^m_{AB}(t + 0.5dt)
    */
    // ::state legacy version
    Molecule _move_correct_half_step(Molecule molecule_last_step,
        Molecule molecule_half_step_to_correct);
    // ::simple_state version
    simple::AtomPolymer _move_correct_half_step(
        simple::AtomPolymer molecule_last_step,
        simple::AtomPolymer molecule_half_step_to_correct);
    /**
    * before calling,
    *      call unconstrained verlet: v(t + 0.5dt) -> v^0(t + dt),
    *      and save the output molecule as molecule_full_step
    * Then: for each constrained bond r_{AB}(t+dt), iterative correction:
    *      1. v^i_{AB}(t + dt) -> v^{i+1}_{AB}(t + dt), until for all bonds
    *          | r_{AB}(t + dt) * v_{AB} |^2 < 2 * rvtol * d^2_{AB},
    *              where d_{AB} = r_{AB}(t+dt) is the constrained bond length,
    *              a factor 2 in front of rvtol is from Taylor expansion
    *      2. When all bonds are satisfied to within given tolerance, save
    *          the corrected versions of velocity
    *          v_{AB}(t + dt) = v^m_{AB}(t + dt),
    */
    // ::state legacy version
    Molecule _move_correct_full_step(Molecule molecule_full_step_to_correct,
        bool calculate_observables);
    // ::simple_state version
    simple::AtomPolymer _move_correct_full_step(
        simple::AtomPolymer molecule_full_step_to_correct,
        bool calculate_observables);
    virtual void _zero_accumulators();
    virtual void _correct_accumulators();
    /**
    * The two vectors below manage bookkeeping and cut down on work ---
    * Allen and Tildesley codes at CCL or RATTLE
    */
    // keeps track of which bonds are being corrected at current step i
    // (the bond is corrected if its difference with set length exceeds the
    // specified tolerance)
    std::vector<bool> _moving;
    // keeps track of which bonds were corrected at previous step i - 1
    // (the bond is corrected if its difference with set length exceeds the
    // specified tolerance)
    std::vector<bool> _moved;
    // after verlet steps for a given molecule, resizes the moved-moving vectors
    // number of elements is equal to the number of constrained bonds
    // all bonds are marked as moved; none are marked as moving
    // legacy version
    void _set_up_correction_bookkeeping(Molecule& molecule);
    // simple::state version
    void _set_up_correction_bookkeeping();
public:
    // getters
    double get_tol() const;
    double get_rvtol() const;
    double get_tiny() const;
    virtual double *get_kinetic_energy_acc();
    double *get_neg_constraint_virial_acc();
    // setters
    void set_tol(double tol);
    void set_tiny(double tiny);
    virtual void set_kinetic_energy_acc(double *kinetic_energy_acc);
    void set_neg_constraint_virial_acc(double *neg_constraint_virial_acc);
    // main function
    // check for consistency of the state
    // zero observable accumulators, if necessary
    // call verlet half step
    // perform iterative correction of positions and half-step velocities,
    //      initate the force loop to calculate f(t+dt),
    //      calculate potential and non-constraint virial if necessary
    // call verlet full step (*w/o* observable calculation)
    // perform iterative correction of full-step velocities
    //      calculate constraint virial and kinetic energy if necessary
    virtual void move(double timestep,
        State& state,
        bool calculate_observables = false);
    virtual void move(double timestep,
        simple::AtomState& state,
        bool calculate_observables = false);
    // constructors and a destructor
    RattleIntegrator();
    RattleIntegrator(double tol, double rvtol, double tiny = pow(10, -7.0));
    RattleIntegrator(ObservableContainer *observables);
    RattleIntegrator(ForceUpdater force_updater,
        double tol,
        double tiny = pow(10, -7.0),
        int maxiter = pow(10, 3),
        double *kinetic_energy_acc = NULL,
        double *neg_constraint_virial_acc = NULL);
    ~RattleIntegrator();
};
#endif
