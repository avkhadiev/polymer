// 2017 Artur Avkhadiev
/*! \file parsing.h
*/
#ifndef PARSING_H
#define PARSING_H
#include <string>
#include <algorithm>
#include <map>
// parses a string and replaces all replace_what characters with replace_with
std::string replace_characters(std::string s,
    char replace_what,
    char replace_with);
// uses replace_characters to replace a series of characters with pre-defined
// substitutes; makes all characters lowercase before replacement
std::string parse_string(std::string s);
#endif
