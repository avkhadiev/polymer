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
#include "force_updater.h"
#include "observable_container.h"
#include "integrator.h"
class VerletIntegrator :
    public Integrator {
protected:
    // boolean has_observable_container indicates whether a pointer to the
    // observable container is valid
    ForceUpdater _force_updater;
    // pointers to accumulators of observables
    bool _is_kinetic_energy_set;
    double *_kinetic_energy_acc;
    void _move_verlet_half_step(double timestep, State& state);
    void _move_verlet_full_step(double timestep,
        State& state,
        bool calculate_observables);
public:
    virtual ForceUpdater& get_force_updater();
    virtual void set_force_updater(ForceUpdater force_updater);
    virtual void move(double timestep,
        State& state,
        bool calculate_observables = false);
    // constructors and a destructor
    VerletIntegrator(ForceUpdater force_updater,
        double* kinetic_energy_acc = NULL);
    virtual ~VerletIntegrator();
};
#endif
