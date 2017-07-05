// 2017 Artur Avkhadiev
/*! \file integrator.h
*/
#ifndef POLYMER_INTEGRATOR_H
#define POLYMER_INTEGRATOR_H
#include <vector>
#include <ljpotential.h>
#include <molecule.h>
#include "observable_container.h"
class Integrator {
public:
    virtual void move(double timestep,
        LJPotential &potential,
        std::vector<Molecule> molecules,
        bool calculate_observables = false);
    virtual ~Integrator();
};
#endif
