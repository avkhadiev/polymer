// 2017 Artur Avkhadiev
/*! \file settings_parser_test.cpp
*/
#include <string>
#include <iostream>
#include <gtest/gtest.h>
#include "../include/settings_parser.h"
#include "../include/parsing.h"

class SettingsParserTest : public ::testing::Test {
 protected:
     std::string outdir;
     std::string name1;
     std::string name2;
     virtual void SetUp() {
        outdir = "/Users/Arthur/stratt/polymer/test/";
        name1 = "/Users/Arthur/stratt/polymer/test/settings_parser_test1.cfg";
        name2 = "/Users/Arthur/stratt/polymer/test/settings_parser_test2.cfg";
     }
  // virtual void TearDown() {}
};

TEST_F(SettingsParserTest, Initalization) {
    SettingsParser test_out = SettingsParser();
    test_out.write(outdir, name1);
    SettingsParser test_in = SettingsParser(outdir + name1);
    test_in.write(outdir, name2);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
