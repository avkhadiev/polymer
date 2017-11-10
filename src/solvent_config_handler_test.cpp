// 2017 Artur Avkhadiev
/*! \file solvent_config_handler_test.cpp
*/
#include <string>
#include <stdexcept>
#include "../include/vector.h"
#include "../include/simple_state.h"
#include "../include/solvent_config_handler.h"

#include <gtest/gtest.h>

class SolventConfigTest : public ::testing::Test {
protected:
    int nc;
    double density;
    double box;
    double temperature;
    // IO settings
    bool verbose;
    bool overwrite;
    std::string fname;
    std::string outdir;
    std::string indir;
    virtual void SetUp() {
        // 4 * 8 = 32 atoms
        // L = (32/4)^(1/3) = 2 sigma
        nc = 2;
        density = 4.0;
        box = 2;
        temperature = 1.0;
        verbose = true;
        overwrite = true;
        fname ="solvent_config_check";
        outdir = "/Users/Arthur/stratt/polymer/test/";
        indir = "/Users/Arthur/stratt/polymer/test/";
    }
    // virtual void TearDown() {}
};

TEST_F(SolventConfigTest, Init) {
    SolventConfigHandler test = SolventConfigHandler(density, nc);
    EXPECT_EQ(test.nc(), nc);
    EXPECT_EQ(test.box(), box);
}

TEST_F(SolventConfigTest, FCCPos) {
    SolventConfigHandler test = SolventConfigHandler(density, nc);
    test.fcc_positions();
    test.ran_velocities(temperature);
    test.atom_state().write_to_file(outdir, fname, verbose, overwrite);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
