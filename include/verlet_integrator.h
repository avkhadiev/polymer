// 2017 Artur Avkhadiev
/*! \file verlet_integrator.h
*/
#ifndef POLYMER_VERLET_INTEGRATOR_H
#define POLYMER_VERLET_INTEGRATOR_H
#include <vector>
#include <cmath>              /* pow */
#include "ljpotential.h"
#include "simple_atom.h"
#include "simple_polymer.h"
#include "force_updater.h"
#include "observable.h"
#include "integrator.h"
class VerletIntegrator :
    public Integrator {
protected:
    // boolean has_observable_container indicates whether a pointer to the
    // observable container is valid
    ForceUpdater& _force_updater;
    double _timestep;
    double _halfstep;
    virtual void _set_timestep(double timestep);
    // pointers to accumulators of observables
    struct obs_t {
        bool is_set;
        Observable *ptr;
    } _ke;                                          /**> kinetic energy*/
    struct mirror_image_t {
        bool is_set;
        double box;
    } _mirror_image;
    // if mirror image is set, uses box to apply the mirror image convention
    // outputs the adjusted position.
    // uses function round() from math.h instead of ANINT (same behavior)
    Vector _get_mirror_image(Vector r);
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
    simple::Solvent _move_verlet_half_step(simple::Solvent molecule);
    simple::Solvent _move_verlet_full_step(simple::Solvent molecule, bool calculate_observables);
    virtual void _zero_accumulators();
public:
    virtual ForceUpdater& get_force_updater();
    virtual void set_ke(Observable *ptr);
    virtual void move(double timestep,
        simple::AtomState& state,
        bool calculate_observables = false);
    // constructors and a destructor
    // if *observable = NULL, that observable will not be updated
    // if box = 0.0, no mirror image convention will be applied to solvent pos's
    VerletIntegrator(ForceUpdater& force_updater,
        Observable* ke = NULL,
        double box = 0.0);
    virtual ~VerletIntegrator();
};
#endif
