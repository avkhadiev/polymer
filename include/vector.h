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
#ifndef POLYMER_VECTOR_H
#define POLYMER_VECTOR_H
namespace math {
    // struct vector contains x, y, z coordinates
    // is used to represent all 3d vectors in the simulation
    // e.g. position, velocity, acceleration, force, angular momentum
    typedef struct vector_t {
        double x;
        double y;
        double z;
    } Vector;
    // Operations on vectors
    // Addition
    inline Vector add(Vector v1,  Vector v2);
    inline Vector add(Vector v1,  double s);
    // Subtraction
    inline Vector subtract(Vector v1,  Vector v2);
    inline Vector subtract(Vector v1,  double s);
    // Multiplication
    /**
    * Takes two vectors v1, v2 and outputs a vector
    *   v = (v1.x * v2.x, v1.y * v2.y, v1.z * v2.z)
    */
    inline Vector multiply(Vector v1, Vector v2);
    inline Vector multiply(Vector v1, double s);
    inline double dot(Vector v1, Vector v2);
    inline Vector cross(Vector v1, Vector v2);
    // Division
    inline Vector divide(Vector v1, Vector v2);
    /**
    * Takes a vector v1, and a scalar s outputs a vector v = v1 / s
    */
    inline Vector divide(Vector v1, double s);
    /**
    * Takes a scalar s and a vector v1 and outputs a vector v = s / v1
    */
    inline Vector divide(double s, Vector v1);
    /**
    * Takes a vector and ouptuts its std::string representation
    */
    inline std::string get_string(Vector v1);
}   // namespace math
#endif
