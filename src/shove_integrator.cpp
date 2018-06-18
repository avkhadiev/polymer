// 2017 Artur Avkhadiev
/*! \file shove_integrator.cpp
*/
#include <vector>
#include <string>
#include <algorithm>                                /* std::fill */
#include <cmath>                                    /* pow */
#include <cassert>
#include "../include/shove_integrator.h"
/**
* initialization: setup force updater with zero potential
* no kinetic energy observables
* no periodic boundary conditions
*/
ShoveIntegrator::ShoveIntegrator(
    double tol,
    Observable *wc,
    double tiny,
    int maxiter):
    RattleIntegrator(ForceUpdater(), tol, NULL, NULL, wc, 0.0, tiny, maxiter){}
ShoveIntegrator::~ShoveIntegrator(){}

/**
* move: utilizes functions from RATTLE
*/
void ShoveIntegrator::move(double timestep,
    simple::AtomState &state,
    bool calculate_observables){
    // really there is no time in configuration space.
    // speeds of the atoms in the state will already be such that the
    // the appropriate "timestep" is 1.0
    _set_timestep(1.0);
    if(calculate_observables){
        _zero_accumulators();
    }
    simple::AtomPolymer molecule_last_step;
    simple::AtomPolymer molecule_half_step_to_correct;
    // loop over polymer molecules, perform verlet half-step and correction
    for (int im = 0; im < state.nm(); ++im){
        // molecule contains positions and velocities at last step
        molecule_last_step = state.polymers.at(im);
        // pass molecule at last step to verlet integrator by value
        // the integrator returns the molecule
        // with uncorrected full-step positions
        // and uncorrected half-step velocities
        molecule_half_step_to_correct = _move_verlet_half_step(molecule_last_step);
        // pass both last-step molecule and molecule from verlet integrator
        // to RATTLE for iterative correction
        // molecule with corrected full-step positions and
        // corrected half-step velocities is stored in the state
        try
        {
            state.polymers.at(im) = _move_correct_half_step(molecule_last_step, molecule_half_step_to_correct);
        }
        catch (std::invalid_argument)
        {
            state.write_to_file("/Users/Arthur/stratt/polymer/test/", "log", true, true);
            fprintf(stderr, "%s\n", "RATTLE: wrote out crashing configuration to log");
            throw;
        }
    }
    // loop over solvent molecules, perform verlet half-step
    simple::Solvent solvent;
    for (int is = 0; is < state.nsolvents(); ++is){
        solvent = state.solvents.at(is);
        state.solvents.at(is) = _move_verlet_half_step(solvent);
    }
    // r(t + dt) has been calculated and corrected
    // ?????????????????????????????????????????????
    // SECOND STAGE --- VELOCITY CORRECTION
    // ?????????????????????????????????????????????
    // multiply accumulators by whatever factors necessary
    if(calculate_observables) {
        _correct_accumulators();
    }
}
