// 2018 Artur Avkhadiev
/*! \file geodesic_record_test.cpp
*/
#include <string>
#include <stdexcept>
#include "../include/geodesic_record.h"

#include <gtest/gtest.h>

class GeodesicRecordTest : public ::testing::Test {
    protected:
        // state set up
        // IO settings
        bool overwrite;
        bool output_header;
        std::string file;
        // Base State Variables
        int nsolvents;
        double solvent_m;
        int nm;
        int nb;
        double m;
        double d;
        Vector zero_vector;
        Vector cm_displacement;
        double t;
        std::vector<simple::Solvent> solvents;
        // Bond State Variables
        Vector d1, d2, ddot1, ddot2;
        std::vector<simple::Bond> bonds;
        std::vector<simple::BondPolymer> bond_polymers;
        simple::BondState bond_state;
        // other stuff
        double pe;
        geodesic::Record record;
    virtual void SetUp() {
        overwrite = true;
        output_header = true;
        file = "/Users/Arthur/stratt/polymer/test/geodesic/record_check.cfg";
        // SOLVENT molecules
        Vector solvent_r;
        Vector solvent_v;
        for(int i = 0; i < nsolvents; ++i){
            double k = i;
            solvent_v = solvent_r = vector(k, k, k);
            solvents.push_back(simple::Solvent(solvent_r, solvent_v));
        }
        /*
        * STATE:
        *  number of solvent atoms = 2
        *  mass of solvent atoms = 0.5
        *  number of bonds = 2
        *  mass of atoms = 1.0
        *  bond length = 3.0
        * SOLVENT MOLECULES
        * r1, v1 = (1, 1, 1), (1, 1, 1)
        * r2, v2 = (2, 2, 2), (2, 2, 2)
        * TRIATOMIC MOLECULE:      (1)                     (2)
        *   RCM                    origin            origin + cm_displacement
        *   VCM                    (0.0, 0.0, 0.0)   (0.0, 0.0, 0.0)
        */
        nsolvents = 2;
        solvent_m = 0.5;
        simple::BaseState::set_nsolvents(nsolvents);
        simple::Solvent::set_m(solvent_m);
        nm = 2;
        nb = 2;
        m = 1.0;
        d = 3.0;
        t = 0.0;
        zero_vector = vector(0.0, 0.0, 0.0);
        cm_displacement = vector(0.0, 0.0, 25.0);
        Vector rcm[nm];
        Vector vcm[nm];
        simple::BaseState::set_nm(nm);
        simple::BondPolymer::set_nb(nb);
        simple::BondPolymer::set_m(m);
        simple::BondPolymer::set_d(d);
        for(int i = 0; i < nm; ++i){
            rcm[i] = add(zero_vector, multiply(cm_displacement, i));
            vcm[i] = zero_vector;
        }
        for(int i = 0; i < nsolvents; ++i){
            rcm[i] = add(zero_vector, multiply(cm_displacement, i));
            vcm[i] = zero_vector;
        }
        /*
        * BOND REPRESENTATION IN CM FRAME:
        *  BOND1_R = BOND2_R = (0, 1, 0)
        *  BOND1_V = (0, 0, 1), BOND2_V = (0, 0, -1)
        */
        d1 = vector(0.0, 1.0, 0.0);
        d2 = vector(0.0, 1.0, 0.0);
        ddot1 = vector(0, 0, 1.0);
        ddot2 = vector(0, 0, -1.0);
        simple::Bond bond[nm][nb];
        for(int i = 0; i < nm; ++i){
           bond[i][0] = simple::Bond(d1, ddot1);
           bond[i][1] = simple::Bond(d2, ddot2);
           bonds.push_back(bond[i][0]);
           bonds.push_back(bond[i][1]);
           bond_polymers.push_back(simple::BondPolymer(bonds, rcm[i], vcm[i]));
           bonds.clear();
        }
        bond_state = simple::BondState(bond_polymers, solvents, t);
        pe = 42.0;
        record = geodesic::Record(bond_state, pe);
    }
     // virtual void TearDown() {}
};

TEST_F(GeodesicRecordTest, DataMembers) {
    EXPECT_EQ(bond_state, record.state());
    EXPECT_EQ(pe, record.pe());
}

TEST_F(GeodesicRecordTest, BinaryOperators) {
    //fprintf(stderr, "%s\n", "checking equality");
    // FIXME this comparison is not being done yet; see source for
    // == operator in geodesic record
    EXPECT_EQ(record == record, true);
    EXPECT_EQ(record != record, false);
}

TEST_F(GeodesicRecordTest, IO) {
    //fprintf(stderr, "%s\n", "creating new record");
    geodesic::Record r_expect = record;
    //fprintf(stderr, "%s\n", "writing to file");
    r_expect.write(file, overwrite);
    //fprintf(stderr, "%s\n", "assigning from file");
    geodesic::Record r_check = geodesic::Record(file);
    //fprintf(stderr, "%s\n", "checking equality");
    EXPECT_EQ(r_check, r_expect);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
