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

class MoleculeTest : public ::testing::Test {
 protected:
     Vector r1;
     Vector r2;
     Vector v;
     Vector f;
     double m;
     double t;
     double fixed_length;
     Atom a1;
     Atom a2;
     std::vector<Atom> atoms;
     Bond b;
     std::vector<Bond> bonds;
     Molecule mol;
     virtual void SetUp() {
         r1.x = 1.0; r1.y = 1.0; r1.z = 1.0;
         r2.x = -1.0; r2.y = -1.0; r2.z = -1.0;
         v.x = 0.0; v.y = 0.0; f.z = 0.0;
         f.x = 0.0; f.y = 0.0; f.z = 0.0;
         m = 1.0;
         t = 3.0;
         fixed_length = 1.0;
         a1 = initialize_atom(m, r1, v, t);
         a2 = initialize_atom(m, r2, v, t);
         atoms.push_back(a1);
         atoms.push_back(a2);
         b = initialize_bond(0, 1, fixed_length);
         bonds.push_back(b);
     }
  // virtual void TearDown() {}
};

TEST_F(MoleculeTest, Initialization) {
    ASSERT_NO_THROW(initialize_molecule(atoms, bonds));
    Molecule mol = initialize_molecule(atoms, bonds);
    EXPECT_EQ(mol.na, atoms.size());
    EXPECT_EQ(mol.nb, bonds.size());
    mol.atoms.at(1).position.second = -1 * t;
    mol.atoms.at(1).velocity.second = 2 * t;
    mol.atoms.at(0).position.second = 3 * t;
    mol.atoms.at(0).velocity.second = 4 * t;
    EXPECT_EQ(false, is_time_consistent(mol));
    std::cout << mol << std::endl;
    set_time(&mol, t);
    std::cout << mol << std::endl;
    EXPECT_EQ(true, is_time_consistent(mol, t));
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
