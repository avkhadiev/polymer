// 2017 Artur Avkhadiev
/*! \file study_ljpotential.cpp
*/
#include <iostream>
#include <cmath>                            /* pow */
#include <string>
#include "../include/ljpotential.h"
#include "../include/ljpotential_observable_container.h"

// default output directory
const std::string default_outdir = "/Users/Arthur/stratt/polymer/test/";

// number of observations to store in memory before writeout
const int observations_before_writeout = 100;

// distances in units of sigma
const double distance_over_sigma_min = 0.01;
const double distance_over_sigma_max = 2.0;
const double distance_over_sigma_step = 0.01;

const std::string pair_potential_str = "Potential";
const std::string pair_virial_str = "Virial";
const std::string pair_fstrength_str = "Force";
const std::string distance_str = "Distance";

void make_observation(double rij,
    LJPotential *potential,
    LJPotentialObservableContainer *container) {
    // get observable accumulators
    double *distance
        = container->get_scalar_observable_accumulator(distance_str);
    double *pair_potential
        = container->get_scalar_observable_accumulator(pair_potential_str);
    double *pair_virial
        = container->get_scalar_observable_accumulator(pair_virial_str);
    double *pair_fstrength
        = container->get_scalar_observable_accumulator(pair_fstrength_str);
    // calculate inverse rij squarred
    double inv_rijsq = 1 / (pow(rij, 2.0));
    // update accumulators;
    *distance = rij;
    *pair_potential = potential->calculate_pair_potential(inv_rijsq);
    *pair_virial = -1.0 * potential->calculate_neg_pair_virial(inv_rijsq);
    *pair_fstrength = potential->calculate_fstrength_over_r(inv_rijsq) * rij;
    // record observations
    // no time is associated with studying the form of the potential
    double time = 0.0;
    container->update_observable_through_accumulator(distance_str, time);
    container->update_observable_through_accumulator(pair_potential_str, time);
    container->update_observable_through_accumulator(pair_virial_str, time);
    container->update_observable_through_accumulator(pair_fstrength_str, time);
};

// takes in arguments: epsilon & sigma of the potential to study
int main(int argc, char **argv){
    std::string fname;
    if (argc > 2){
        fprintf(stderr,
            "%s\n",
            "Usage: ./study_ljpotential optional <outdir>");
    }
    else {
        // declare the title of the study
        std::string study_name = "ljpot";
        // save output directory
        std::string outdir;
        if (argc == 2){
            outdir = std::string(argv[1]);
        }
        else {
            outdir = default_outdir;
        }
        LJPotential ljpot = LJPotential();
        // create the observable container to store data
        LJPotentialObservableContainer container
            = LJPotentialObservableContainer();
        double dist_min = distance_over_sigma_min;
        double dist_max = distance_over_sigma_max;
        double dist_step = distance_over_sigma_step;
        // for a range of distances, calculate and store the observables
        // number of observations stored in container
        int nobservations = 0;
        // overwrite observables record only on first write
        bool overwrite = true;
        for (double dist = dist_min; dist < dist_max; dist += dist_step) {
            make_observation(dist, &ljpot, &container);
            nobservations += 1;
            if (nobservations > observations_before_writeout) {
                container.writeout_observables_to_file({},
                    outdir,
                    study_name,
                    overwrite);
                // do not overwrite after the first time
                overwrite = false;
                // reset the counter
                nobservations = 0;
                // empty observable_container
                container.clear_observables_records();
            }
        }
        // writeout the remaining observations
        container.writeout_observables_to_file({},
            outdir,
            study_name,
            overwrite);
        fprintf(stdout, "%s %s\n",
            "wrote out csv data to",
            outdir.c_str());
    }
}
