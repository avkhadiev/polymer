// 2017 Artur Avkhadiev
/*! \file parsing_test.cpp
*/
#include <string>
#include <stdexcept>
#include <../include/parsing.h>
#include <gtest/gtest.h>

class ParsingTest : public ::testing::Test {
 protected:
     std::string s_check;
     std::string s_expect;
     char replace_what;
     char replace_with;
     virtual void SetUp() {
     }
    // virtual void TearDown() {}
};

TEST_F(ParsingTest, ReplaceCharacters){
    s_check = "this is a trial string";
    s_expect = "this_is_a_trial_string";
    replace_what = ' ';
    replace_with = '_';
    s_check = replace_characters(s_check, replace_what, replace_with);
    EXPECT_EQ(s_expect, s_check);
}

TEST_F(ParsingTest, NothingToReplace){
    s_check = "this_is_a_trial_string";
    s_expect = "this_is_a_trial_string";
    replace_what = ' ';
    replace_with = '_';
    s_check = replace_characters(s_check, replace_what, replace_with);
    EXPECT_EQ(s_expect, s_check);
}

TEST_F(ParsingTest, ParseString){
    s_check = "Test Run Epsilon=1.0";
    s_expect = "Test_Run_Epsilon=1p0";
    s_check = parse_string(s_check);
    EXPECT_EQ(s_expect, s_check);
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
