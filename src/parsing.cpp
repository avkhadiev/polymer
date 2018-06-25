// 2017 Artur Avkhadiev
/*! \file parsing.cpp
*/
#include <../include/parsing.h>
std::string replace_characters(std::string s,
    char replace_what,
    char replace_with){
    std::replace(s.begin(), s.end(), replace_what, replace_with); 
    return s;
}
std::string parse_string(std::string s){
    std::map<char, char> character_replacement_dict;
    character_replacement_dict[' '] = '_';
    character_replacement_dict['.'] = 'p';
    char replace_what;
    char replace_with;
    // make all characters lowercase first
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    for (std::map<char,char>::iterator it=character_replacement_dict.begin();
        it!=character_replacement_dict.end();
        ++it){
        replace_what = it->first;
        replace_with = it->second;
        s = replace_characters(s, replace_what, replace_with);
    }
    return s;
}
