#ifndef POLYMER_INCLUDE_BASICS_H
#define POLYMER_INCLUDE_BASICS_H
namespace basics {
    typedef struct units_t {
        char *units;
    } Units;
    typedef struct time_t {
        double t;
        const Units *units;
    } Time;
}   // namespace basics
#endif
