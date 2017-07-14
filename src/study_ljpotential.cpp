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
const double distance_over_sigma_max = 5.0;
const double distance_over_sigma_step = 0.01;

const std::string pair_potential_str = "Pair Potential";
const std::string pair_virial_str = "Pair Virial";
const std::string pair_fstrength_str = "Force Strength";
const std::string distance_str = "Distance";

void make_observation(double rij,
    LJPotential *potential,
    LJPotentialObservableContainer *container) {
    // calculate all the prefactors for observables
    potential->set_calculate_prefactors(true);
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
    if (argc < 3 || argc > 4){
        fprintf(stderr,
            "%s\n",
            "Usage: ./study_ljpotential <epsilon> <sigma> optional <outdir>");
    }
    else {
        // declare the title of the study
        std::string study_name = "ljpot_e_" + std::string(argv[1])
            + "_s_" + std::string(argv[2]);
        // retrieve potential parameters
        double epsilon = atof(argv[1]);
        double sigma = atof(argv[2]);
        // save output directory
        std::string outdir;
        if (argc == 4){
            outdir = std::string(argv[3]);
        }
        else {
            outdir = default_outdir;
        }
        // create potential with indicated parameters
        // always calculate  prefactors (proportional to epsilon) for all
        // observables
        bool calculate_prefactors = true;
        LJPotential ljpot = LJPotential(epsilon, sigma, calculate_prefactors);
        // create the observable container to store data
        LJPotentialObservableContainer container
            = LJPotentialObservableContainer();
        // writeout potential parameters in a csv file
        ljpot.writeout_parameters_to_file(outdir, study_name);
        // convert rij coordinates using the given value of sigma
        double dist_min = distance_over_sigma_min * sigma;
        double dist_max = distance_over_sigma_max * sigma;
        double dist_step = distance_over_sigma_step * sigma;
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
