// 2018 Artur Avkhadiev
/*! \file geodesic_manager_test.cpp
*/
#include <string>
#include <stdexcept>
#include <cmath>
#include <list>
#include "../include/ljpotential.h"
#include "../include/geodesic_manager.h"

#include <gtest/gtest.h>

class GeodesicManagerTest : public ::testing::Test {
    protected:
        // IO
        bool overwrite;
        std::string cndir; std::string outdir; std::string sim_name;
        std::string ini_file; std::string fin_file; std::string md_path_file;
        // ForceUpdater
        LJPotential pp;
        AdjustedLJPotential ss;
        AdjustedLJPotential ps;
        ForceUpdater fupd;
        // The manager itself
        geodesic::Manager manager_w_p;
        geodesic::Manager manager_wo_p;
    virtual void SetUp() {
        // setup IO
        overwrite = true;
        cndir = "/Users/Arthur/stratt/polymer/test/";
        outdir = "/Users/Arthur/stratt/polymer/test/geodesic/";
        sim_name = "test";
        ini_file = outdir + sim_name + "_ini_rec.cfg";
        fin_file = outdir + sim_name + "_fin_rec.cfg";
        md_path_file = outdir + sim_name + "_md_path.cfg";
        // setup ForceUpdater
        //               epp   sp
        pp = LJPotential(1.0, 1.0);
        //                       ess   ss  rc   box
        ss = AdjustedLJPotential(1.0, 1.0, 5.0, 10.0);
        //                       eps  sps  rc   box
        ps = AdjustedLJPotential(1.0, 1.0, 5.0, 10.0);
        fupd = ForceUpdater(&pp, &ss, &ps);
        // setup managers
        manager_w_p = geodesic::Manager(&fupd);
        manager_wo_p = geodesic::Manager();
        manager_w_p.read_states(cndir, sim_name);
        manager_wo_p.read_states(cndir, sim_name);
    }
    // virtual void TearDown() {}
};

TEST_F(GeodesicManagerTest, IO) {
    manager_w_p.write_geodesic_inputs(outdir, sim_name + "_w_p");
    manager_wo_p.write_geodesic_inputs(outdir, sim_name + "_wo_p");
    manager_w_p.initial().write(ini_file, overwrite);
    manager_w_p.final().write(fin_file, overwrite);
}

TEST_F(GeodesicManagerTest, MDPath) {
    geodesic::Path path = manager_w_p.MD_path();
    path.write(md_path_file, overwrite);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
