// 2017 Artur Avkhadiev
/*! \file vector.cpp
    \brief Contains implementations of functions on a Vector struct.

    Functions add, subtract, multiply and divide take in 2 arguments and output
    a vector. They are overaloaded to take in either 2 vectors or a vector and
    a scalar.
    Additionally, dot product and cross product functions are defined.
*/
#include <iostream>
#include "../include/vector.h"
namespace math {
    // Operations on vectors
    // Addition
    inline Vector add(Vector v1,  Vector v2) {
        Vector v;
        v.x = v1.x + v2.x;
        v.y = v1.y + v2.y;
        v.z = v1.z + v2.z;
        return v;
    };
    inline Vector add(Vector v1,  double s) {
        Vector v;
        v.x = v1.x + s;
        v.y = v1.y + s;
        v.z = v1.z + s;
        return v;
    };
    // Subtraction
    inline Vector subtract(Vector v1,  Vector v2) {
        Vector v;
        v.x = v1.x - v2.x;
        v.y = v1.y - v2.y;
        v.z = v1.z - v2.z;
        return v;
    };
    inline Vector subtract(Vector v1,  double s) {
        Vector v;
        v.x = v1.x - s;
        v.y = v1.y - s;
        v.z = v1.z - s;
        return v;
    };
    // Multiplication
    inline Vector multiply(Vector v1, Vector v2) {
        Vector v;
        v.x = v1.x * v2.x;
        v.y = v1.y * v2.y;
        v.z = v1.z * v2.z;
        return v;
    };
    inline Vector multiply(Vector v1, double s) {
        Vector v;
        v.x = v1.x * s;
        v.y = v1.y * s;
        v.z = v1.z * s;
        return v;
    };
    inline double dot(Vector v1, Vector v2) {
        double s;
        s = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
        return s;
    };
    inline Vector cross(Vector v1, Vector v2) {
        Vector v;
        v.x = v1.y * v2.z - v1.z * v2.y;
        v.y = v1.z * v2.x - v1.x * v2.z;
        v.z = v1.x * v2.y - v1.y * v2.x;
        return v;
    };
    // Division
    inline Vector divide(Vector v1, Vector v2) {
        Vector v;
        v.x = v1.x / v2.x;
        v.y = v1.y / v2.y;
        v.z = v1.z / v2.z;
        return v;
    };
    inline Vector divide(Vector v1, double s) {
        Vector v;
        v.x = v1.x / s;
        v.y = v1.y / s;
        v.z = v1.z / s;
        return v;
    };
    inline Vector divide(double s, Vector v1) {
        Vector v;
        v.x = s / v1.x;
        v.y = s / v1.y;
        v.z = s / v1.z;
        return v;
    };
    inline std::string get_string(Vector v1) {
        std::string v_str_x = std::to_string(v1.x);
        std::string v_str_y = std::to_string(v1.y);
        std::string v_str_z = std::to_string(v1.z);
        std::string v_str = "(" + v_str_x + ", " + v_str_y + ", " + v_str_z + ")";
        return v_str;
    };
} // namespace math
int main(int argc, char **argv) {
    if (argc != 4) {
        std::string err_msg = "USAGE: vector x y z";
        fprintf(stderr, "%s\n", err_msg.c_str());
        return 1;
    }
    else{
        // Print the vector back
        double x = atof(argv[1]);
        double y = atof(argv[2]);
        double z = atof(argv[3]);
        math::Vector v = {.x = x, .y = y, .z = z};
        std::string v_str = math::get_string(v);
        fprintf(stdout, "-30%s%-15s\n", "Vector v:", v_str.c_str());
        return 0;
    }
};
