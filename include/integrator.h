// 2017 Artur Avkhadiev
/*! \file integrator.h
*/
#ifndef POLYMER_INTEGRATOR_H
#define POLYMER_INTEGRATOR_H
#include <vector>
#include <ljpotential.h>
#include <state.h>
#include "observable_container.h"
class Integrator {
public:
    virtual void move(double timestep,
        LJPotential& potential,
        State& state,
        bool calculate_observables = false);
    virtual ~Integrator();
};
#endif
