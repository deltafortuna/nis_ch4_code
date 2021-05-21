#ifndef SUMMARYSTATS_H
#define SUMMARYSTATS_H

#include "params.h"

// declarations
double get_pi (vector<bitset<bitlength>> &sample);
double get_watterson (vector<bitset<bitlength>> &sample, int S);
double get_tajimas_d(double pi, double watterson, int S);

// definitions
double get_pi (vector<bitset<bitlength>> &sample) {
	double sumdiffs = 0.;
	double numcomp = 0.;

	#pragma omp parallel for collapse(2)
	for (int i  = 0; i< sample.size() - 1; ++i) {
		for (int j = i+1; j<sample.size(); ++j) {
			sumdiffs += (sample[i] ^ sample[j]).count();
			numcomp+=1.;
		}
	}
	return (sumdiffs/numcomp);
}

double get_watterson (vector<bitset<bitlength>> &sample, int S) {
	double denominator = 0.;
	for (double i=1.; i<sampsize; ++i)
		denominator += 1./i;
	return (S/denominator);
}

double get_tajimas_d (double pi, double watterson, int S) {
	double d = pi - watterson; // numerator
	double a1 = 0.;
	double a2 = 0.;
	for (double i=1.; i < sampsize; ++i) {
		a1 += 1./i;
		a2 += 1./(i*i);
	}
	double n = sampsize; // for easier expression
	double var = watterson * ( (n+1) / (3*(n-1)) - 1/a1 ) +
	    ( S * (S-1) * ( 1 / (a1*a1+a2)) *
		( (2*(n*n + n+3))/(9*n*(n-1)) - (n+2)/(a1*n) + a2/(a1*a1) ) );
	return(d / sqrt(var));
}

#endif
