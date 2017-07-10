// 2017 Artur Avkhadiev
/*! \file vector.h
    \brief Contains definitions of a vector and related functions.

    Vector is a struct of 3 doubles.
    Functions add, subtract, multiply and divide take in 2 arguments and output
    a vector. They are overaloaded to take in either 2 vectors or a vector and
    a scalar.
    Additionally, dot product and cross product functions are defined.
    Function get_string converts a vector into a string.
*/
#include <string>
#include <iostream>
#ifndef POLYMER_VECTOR_H
#define POLYMER_VECTOR_H
// struct vector contains x, y, z coordinates
// is used to represent all 3d vectors in the simulation
// e.g. position, velocity, acceleration, force, angular momentum
struct vector_t {
    double x;
    double y;
    double z;
    vector_t& operator=(const vector_t& v)
    {
        x=v.x;
        y=v.y;
        z=v.z;
        return *this;
    }
    vector_t& operator-()
    {
        x = -x;
        y = -y;
        z = -z;
        return *this;
    }
    bool operator==(const vector_t& v) const
    {
        return (x == v.x && y == v.y && z == v.z);
    }
    bool operator!=(const vector_t& v) const
    {
        return (x != v.x || y != v.y || z != v.z);
    }
};
typedef struct vector_t Vector;
// Constructors of vectors
Vector vector(double x, double y, double z);
// Operations on vectors
// Addition
Vector add(Vector v1,  Vector v2);
Vector add(Vector v1,  double s);
// Subtraction
Vector subtract(Vector v1,  Vector v2);
Vector subtract(Vector v1,  double s);
// Multiplication
/**
* Takes two vectors v1, v2 and outputs a vector
*   v = (v1.x * v2.x, v1.y * v2.y, v1.z * v2.z)
*/
Vector multiply(Vector v1, Vector v2);
Vector multiply(Vector v1, double s);
double dot(Vector v1, Vector v2);
Vector cross(Vector v1, Vector v2);
// Division
Vector divide(Vector v1, Vector v2);
// Exponentiaion
Vector pow(Vector v1, double s);
/**
* Takes a vector v1, and a scalar s outputs a vector v = v1 / s
*/
Vector divide(Vector v1, double s);
double normsq(Vector v);
double norm(Vector v);
/**
* Takes a vector and ouptuts its std::string representation
*/
std::string vector_to_string(Vector v1);
::std::ostream& operator<<(::std::ostream& os, const Vector& v);
/**
* Takes an std::string representation of a vector and outputs the vector
*/
Vector string_to_vector(std::string v_str);
#endif
