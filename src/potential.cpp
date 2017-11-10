// 2017 Artur Avkhadiev
/*! \file potential.cpp
*/
#include "../include/potential.h"
Potential::Potential(){}
Potential::~Potential(){};
::std::ostream& operator<<(::std::ostream& os, const Potential& potential){
    return os << potential.get_str().c_str();
}
