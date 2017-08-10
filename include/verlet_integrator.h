// 2017 Artur Avkhadiev
/*! \file verlet_integrator.h
*/
#ifndef POLYMER_VERLET_INTEGRATOR_H
#define POLYMER_VERLET_INTEGRATOR_H
#include <vector>
#include <cmath>              /* pow */
#include "ljpotential.h"
#include "simulation.h"
#include "molecule.h"
#include "simple_atom.h"
#include "simple_polymer.h"
#include "force_updater.h"
#include "observable_container.h"
#include "integrator.h"
class VerletIntegrator :
    public Integrator {
protected:
    // boolean has_observable_container indicates whether a pointer to the
    // observable container is valid
    ForceUpdater _force_updater;
    double _timestep;
    double _halfstep;
    virtual void _set_timestep(double timestep);
    // pointers to accumulators of observables
    bool _is_kinetic_energy_acc_set;
    double *_kinetic_energy_acc;
    /**
    * The following two helper functions take in the molecule by value:
    * One feeds in the value of the molecule at time t, and receives the value
    * of the molecule at time t + dt. To do so by value rather than by reference
    * is important because integrator classes that deriver from Verlet may
    * require both Verlet output of molecule at (t + dt) and the initial value
    * at t. For example, RATTLE needs both r_{AB}(t) and r^{0}(t+dt)_{AB},
    * where the latter is the (unconstrained) value of the bond vector at t+dt.
    */
    simple::AtomPolymer _move_verlet_half_step(simple::AtomPolymer molecule);
    simple::AtomPolymer _move_verlet_full_step(simple::AtomPolymer molecule,
        bool calculate_observables);
    Molecule _move_verlet_half_step(Molecule molecule);
    Molecule _move_verlet_full_step(Molecule molecule, bool calculate_observables);
    virtual void _zero_accumulators();
    virtual void _correct_accumulators();
public:
    virtual double *get_kinetic_energy_acc();
    virtual ForceUpdater& get_force_updater();
    virtual void set_force_updater(ForceUpdater force_updater);
    virtual void set_kinetic_energy_acc(double *kinetic_energy_acc);
    virtual void move(double timestep,
        State& state,
        bool calculate_observables = false);
    virtual void move(double timestep,
        simple::AtomState& state,
        bool calculate_observables = false);
    // constructors and a destructor
    VerletIntegrator(ForceUpdater force_updater,
        double* kinetic_energy_acc = NULL);
    virtual ~VerletIntegrator();
};
#endif
