#include <algorithm>
#include <string>
#include <iterator>
#include <bitset>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <regex>

using namespace std;

#include "params.h"
#include "params.cc" // inclusion causes parameter values to be read
#include "allele.h"
#include "individual.h"
#include "population.h"

int main(int argc, char *argv[]) {

	int gens = atoi(argv[1]);
	mt19937 engine(time(0));  //initialize the random engine
	Population::e  = engine;
	Population pop;

cout << popsize << " is popsize" << endl;
cout << dem_parameter[1] << endl;

	// simulate for gens generations
	for (int i =0; i < gens; i++) {
		pop.reproduce(i);
		if (i % 100 == 0)
			cout << "gen " << i << endl;
	}
	pop.close_output_files();

	return 0;
}

// static variables for population class
mt19937 Population::e;
