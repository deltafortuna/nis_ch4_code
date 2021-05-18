#ifndef PARAMS_H
#define PARAMS_H

// must be constant at compile time
extern const int bitlength = 1000;

// the following are defined by values in the parameters file
extern int popsize;
extern double mutrate;
extern int seqlength;
extern int sampsize;
extern int sampfreq;
extern vector<int> demography;
extern vector<double> dem_parameter;
extern vector<int> dem_start_gen;
extern vector<int> dem_end_gen;
extern vector<int> carrying_cap;
extern vector<int> pop_schedule;

#endif
