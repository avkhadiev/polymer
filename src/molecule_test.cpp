// 2017 Artur Avkhadiev
/*! \file molecule_test.cpp
*/
#include <string>
#include <utility>              /* std::pair, std::make_pair */
#include <stdexcept>
#include <gtest/gtest.h>
#include "../include/vector.h"
#include "../include/atom.h"
#include "../include/bond.h"
#include "../include/molecule.h"

TEST(MoleculeTest, Initialization) {
    Vector r1 = {.x = 1.0, .y = 1.0, .z = 1.0};
    Vector r2 = {.x = -1.0, .y = -1.0, .z = -1.0};
    Vector v = {.x = 0.0, .y = 0.0, .z = 0.0};
    Vector f = {.x = 0.0, .y = 0.0, .z = 0.0};
    double m = 1.0;
    double t = 3.0;
    Atom a1 = initialize_atom(m, r1, v, t);
    Atom a2 = initialize_atom(m, r2, v, t);
    double fixed_length = 1.0;
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    atoms.push_back(a1);
    atoms.push_back(a2);
    Bond b = initialize_bond(&(atoms.at(0)), &(atoms.at(1)), fixed_length);
    bonds.push_back(b);
    ASSERT_NO_THROW(initialize_molecule(atoms, bonds));
    Molecule mol = initialize_molecule(atoms, bonds);
    std::cout << mol << std::endl;
    EXPECT_EQ(mol.na, atoms.size());
    EXPECT_EQ(mol.nb, bonds.size());
    mol.atoms.at(0).position.second = 2 * t;
    mol.atoms.at(0).velocity.second = 2 * t;
    printf("molecule atoms\n");
    std::cout << mol.atoms.at(0) << std::endl;
    std::cout << mol.atoms.at(1) << std::endl;
    printf("bond atoms\n");
    std::cout << *(b.atom1) << std::endl;
    std::cout << *(b.atom2) << std::endl;
    std::cout << mol << std::endl;
    EXPECT_EQ(false, is_time_consistent(mol));
    set_time(&mol, t);
    EXPECT_EQ(true, is_time_consistent(mol, t));
    std::cout << mol << std::endl;
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
