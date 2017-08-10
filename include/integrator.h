// 2017 Artur Avkhadiev
/*! \file integrator.h
*/
#ifndef POLYMER_INTEGRATOR_H
#define POLYMER_INTEGRATOR_H
#include <vector>
#include "ljpotential.h"
#include "simple_state.h"
#include "force_updater.h"
#include "observable_container.h"
class Integrator {
public:
    virtual ForceUpdater& get_force_updater() = 0;
    virtual void set_force_updater(ForceUpdater force_updater) = 0;
    virtual void move(double timestep,
        simple::AtomState& state,
        bool calculate_observables = false) = 0;
    virtual ~Integrator();
protected:
    Integrator();
};
#endif
