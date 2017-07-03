// 2017 Artur Avkhadiev
/*! \file vector.cpp
    \brief Contains implementations of functions on a Vector struct.

    Functions add, subtract, multiply and divide take in 2 arguments and output
    a vector. They are overaloaded to take in either 2 vectors or a vector and
    a scalar.
    Additionally, dot product and cross product functions are defined.
*/
#include <iostream>
#include <cmath>
#include "../include/vector.h"
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
Vector divide(double s, Vector v1) {
    Vector v;
    v.x = s / v1.x;
    v.y = s / v1.y;
    v.z = s / v1.z;
    return v;
};
double normsq(Vector v) {
    return dot(v, v);
};
double norm(Vector v) {
    return sqrt(normsq(v));
}
std::string get_string(Vector v1) {
    std::string v_str_x = std::to_string(v1.x);
    std::string v_str_y = std::to_string(v1.y);
    std::string v_str_z = std::to_string(v1.z);
    std::string v_str = "(" + v_str_x + ", " + v_str_y + ", " + v_str_z + ")";
    return v_str;
};
::std::ostream& operator<<(::std::ostream& os, const Vector& v) {
    return os << get_string(v).c_str();
};
