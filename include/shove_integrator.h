// 2017 Artur Avkhadiev
/*! \file shove_integrator.h
*/
#ifndef POLYMER_SHOVE_INTEGRATOR_H
#define POLYMER_SHOVE_INTEGRATOR_H
#include <vector>
#include <cmath>              /* pow */
#include "ljpotential.h"
#include "simple_atom.h"
#include "simple_polymer.h"
#include "force_updater.h"
#include "observable.h"
#include "rattle_integrator.h"
class ShoveIntegrator :
    public RattleIntegrator {
public:
    virtual void move(double timestep,
        simple::AtomState& state,
        bool calculate_observables = false);
    // constructors and a destructor
    ShoveIntegrator(
        double tol,
        Observable *wc = NULL,
        double tiny = pow(10, -7.0),
        int maxiter = pow(10, 3));
    ~ShoveIntegrator();
};
#endif
