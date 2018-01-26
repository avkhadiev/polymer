// 2018 Artur Avkhadiev
/*! \file geodesic_path_test.cpp
*/
#include <string>
#include <stdexcept>
#include <cmath>
#include <list>
#include "../include/simple_state.h"
#include "../include/geodesic_record.h"
#include "../include/geodesic_path.h"

#include <gtest/gtest.h>

class GeodesicPathTest : public ::testing::Test {
    protected:
        bool overwrite;
        std::string file;
        double ini_time; double fin_time; double inter_time;
        double ini_pe; double fin_pe; double inter_pe;
        geodesic::Record ini_record;
        geodesic::Record fin_record;
        geodesic::Record inter_record;
        std::list<geodesic::Record> records;
        geodesic::Path path;
    virtual void SetUp() {
        overwrite = true;
        file = "/Users/Arthur/stratt/polymer/test/geodesic/path_check.cfg";
        ini_time = 0.0; fin_time = 1.0; inter_time = 0.5;
        ini_pe = 42.0; fin_pe = 84.0; inter_pe = 63.0;
        simple::BondState ini_state = simple::BondState(ini_time);
        simple::BondState fin_state = simple::BondState(fin_time);
        simple::BondState inter_state = simple::BondState(inter_time);
        ini_record = geodesic::Record(ini_state, ini_pe);
        fin_record = geodesic::Record(fin_state, fin_pe);
        inter_record = geodesic::Record(inter_state, inter_pe);
        records.push_back(ini_record);
        records.push_back(inter_record);
        records.push_back(fin_record);
        path = geodesic::Path(ini_record, fin_record);
    }
    // virtual void TearDown() {}
};

TEST_F(GeodesicPathTest, DataMembers) {
    EXPECT_EQ(path.initial(), ini_record);
    EXPECT_EQ(path.final(), fin_record);
}

TEST_F(GeodesicPathTest, BasicGetters) {
    EXPECT_EQ(path.nrecords(), 2);
    EXPECT_EQ(path.time(), abs(fin_time - ini_time));
}

TEST_F(GeodesicPathTest, ListInitialization) {
    geodesic::Path new_path = geodesic::Path(records);
    EXPECT_EQ(new_path.nrecords(), 3);
    EXPECT_EQ(path.time(), new_path.time());
}

TEST_F(GeodesicPathTest, IO) {
    geodesic::Path new_path = geodesic::Path(records);
    geodesic::Path p_expect = new_path;
    p_expect.write(file, overwrite);
    geodesic::Path p_check = geodesic::Path(file);
    EXPECT_EQ(p_check.initial(), p_expect.initial());
    EXPECT_EQ(p_check.final(), p_expect.final());
    EXPECT_EQ(p_check.time(), p_expect.time());
    EXPECT_EQ(p_check.nrecords(), p_expect.nrecords());
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
