// 2017 Artur Avkhadiev
/*! \file integrator.h
*/
#ifndef POLYMER_INTEGRATOR_H
#define POLYMER_INTEGRATOR_H
#include <vector>
#include <ljpotential.h>
#include <molecule.h>
class Integrator {
public:
    virtual void get_timestep() = 0;
    virtual void set_timestep(double timestep) = 0;
    virtual void move(LJPotential &potential, std::vector<Molecule> molecules) = 0;
    virtual ~Integrator();
};
#endif
