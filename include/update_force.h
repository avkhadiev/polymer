// 2017 Artur Avkhadiev
/*! \file update_force.h
*/
#ifndef POLYMER_UPDATE_FORCE_H
#define POLYMER_UPDATE_FORCE_H
#include "vector.h"
#include "ljpotential.h"
#include "atom.h"
#include "molecule.h"
void update_forces(Molecule *molecules, LJPotential &potential);
#endif
