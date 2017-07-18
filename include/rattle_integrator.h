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
#include "observable_container.h"
#include "verlet_integrator.h"
class RattleIntegrator :
    public VerletIntegrator {
private:
    typedef struct parameters_t {
        double tol;
        double rvtol;
        double tiny;
    } _Parameters;
    void _move_a(double timestep,
        LJPotential& potential,
        State& state);
    void _move_b(double timestep,
        LJPotential& potential,
        State& state,
        bool calculate_observables);
public:
    // getters
    double get_tol();
    double get_rvtol();
    double get_tiny();
    // setters
    void set_tol(double tol);
    void set_rvtol(double rvtol);
    void set_tiny(double tiny);
    virtual void move(double timestep,
        LJPotential &potential,
        State& state,
        bool calculate_observables = false);
    // constructors and a destructor
    RattleIntegrator();
    RattleIntegrator(double tol, double rvtol, double tiny = pow(10, -7.0));
    RattleIntegrator(ObservableContainer *observables);
    RattleIntegrator(ObservableContainer *observables,
        double tol,
        double rvtol,
        double tiny = pow(10, -7.0));
    ~RattleIntegrator();
};
#endif
