// 2017 Artur Avkhadiev
/*! \file vector.cpp
\brief Contains implementations of functions on a Vector struct.

    Functions add, subtract, multiply and divide take in 2 arguments and output
    a vector. They are overaloaded to take in either 2 vectors or a vector and
    a scalar.
    Additionally, dot product and cross product functions are defined.
*/
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include "../include/vector.h"
// constructors
Vector vector(double x, double y, double z) {
    Vector v = {.x = x, .y = y, .z = z};
    return v;
}
// Operations on vectors
// Addition
Vector add(Vector v1,  Vector v2) {
    Vector v;
    v.x = v1.x + v2.x;
    v.y = v1.y + v2.y;
    v.z = v1.z + v2.z;
    return v;
};
Vector add(Vector v1,  double s) {
    Vector v;
    v.x = v1.x + s;
    v.y = v1.y + s;
    v.z = v1.z + s;
    return v;
};
// Subtraction
Vector subtract(Vector v1,  Vector v2) {
    Vector v;
    v.x = v1.x - v2.x;
    v.y = v1.y - v2.y;
    v.z = v1.z - v2.z;
    return v;
};
Vector subtract(Vector v1,  double s) {
    Vector v;
    v.x = v1.x - s;
    v.y = v1.y - s;
    v.z = v1.z - s;
    return v;
};
// Multiplication
Vector multiply(Vector v1, Vector v2) {
    Vector v;
    v.x = v1.x * v2.x;
    v.y = v1.y * v2.y;
    v.z = v1.z * v2.z;
    return v;
};
Vector multiply(Vector v1, double s) {
    Vector v;
    v.x = v1.x * s;
    v.y = v1.y * s;
    v.z = v1.z * s;
    return v;
};
Vector pow(Vector v1, double s) {
    Vector v;
    v.x = pow(v1.x, s);
    v.y = pow(v1.y, s);
    v.z = pow(v1.z, s);
    return v;
};
double dot(Vector v1, Vector v2) {
    double s;
    s = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    return s;
};
Vector cross(Vector v1, Vector v2) {
    Vector v;
    v.x = v1.y * v2.z - v1.z * v2.y;
    v.y = v1.z * v2.x - v1.x * v2.z;
    v.z = v1.x * v2.y - v1.y * v2.x;
    return v;
};
// Division
Vector divide(Vector v1, Vector v2) {
    Vector v;
    v.x = v1.x / v2.x;
    v.y = v1.y / v2.y;
    v.z = v1.z / v2.z;
    return v;
};
Vector divide(Vector v1, double s) {
    Vector v;
    v.x = v1.x / s;
    v.y = v1.y / s;
    v.z = v1.z / s;
    return v;
};
double normsq(Vector v) {
    return dot(v, v);
};
double norm(Vector v) {
    return sqrt(normsq(v));
}
std::string vector_to_string(Vector v1) {
    std::string v_str_x = std::to_string(v1.x);
    std::string v_str_y = std::to_string(v1.y);
    std::string v_str_z = std::to_string(v1.z);
    std::string v_str = v_str_x + " " + v_str_y + " " + v_str_z;
    return v_str;
};
::std::ostream& operator<<(::std::ostream& os, const Vector& v) {
    return os << vector_to_string(v).c_str();
};
Vector string_to_vector(std::string v_str) {
    char delimeter = ' ';
    std::stringstream ss(v_str);
    std::string item;
    std::vector<std::string> v_str_split;
    while (std::getline(ss, item, delimeter))
    {
       v_str_split.push_back(item);
    }
    if (v_str_split.size() != 3) {
        std::string err_msg = "string_to_vector: string does not contain 3 elements.";
        throw std::invalid_argument( err_msg );
    }
    else{
        Vector v;
        v.x = atof(v_str_split[0].c_str());
        v.y = atof(v_str_split[1].c_str());
        v.z = atof(v_str_split[2].c_str());
        return v;
    }
};
