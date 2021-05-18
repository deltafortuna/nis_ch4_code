#ifndef POPULATION_H
#define POPULATION_H

#include "params.h"  // provides access to global parameter values (extern)
#include "summarystats.h" // provides access to  summary statistic calculations/functions

class Population {

private:
double mu_sequence;
vector<Individual*> individuals;
map<int, Allele*> alleles; /// int key is position of the allele
uniform_int_distribution<int> randompos;
uniform_int_distribution<int> randomind;
uniform_real_distribution<double> randomnum;
poisson_distribution<int> randommut;
ofstream pi_file;
ofstream watterson_file;
ofstream allele_file;

vector<vector<int> > mutate(const vector<int> &parents, const int &gen) {
	vector<vector<int> > mutation_results;

	// determine which, if any, positions are mutated
	vector<int> mutnum{randommut(e)};
	mutnum.push_back(randommut(e));

	mutation_results.push_back({mutnum[0]});
	mutation_results.push_back({mutnum[1]});

	// determine which of the two homologs is transmitted by each parent
	mutation_results.push_back({(randomnum(e)<0.5) ? 0 : 1});
	mutation_results.push_back({(randomnum(e)<0.5) ? 0 : 1});

	// resolve any mutation(s) that did occur
	for (int i=0; i<2; ++i)  {
		for (int j = 0; j < mutnum[i]; ++j)  { // loop not entered if no mutation (i.e., mutnum[i] == 0)
			int position = randompos(e);
			if (alleles.find(position) == alleles.end()) { // new mutation to a derived allele in the population
				alleles.insert({position, new Allele(position, gen)});
				mutation_results[i].push_back(position);
			} else { // mutation present in POPULATION; determine if derived allele found in the considered sequence
				vector<int> seq = (*(individuals[parents[i]])).get_sequence(mutation_results[i+2][0]);
				vector<int>::iterator p = find(seq.begin(), seq.end(), position);
				if (p != seq.end())   // back mutation
					mutation_results[i].push_back(position * -1);	 // negative position signals removal of allele by back mutation
				else
					mutation_results[i].push_back(position);
			}
 		}
	}
	return mutation_results;
}

void update_alleles(const int &gen) {
	// reset all counts to zero
	for (auto iter = alleles.begin(); iter != alleles.end(); ++iter)
		(*(iter->second)).set_count(0);
	map<int, int> new_allele_counts;
	for (auto iter = individuals.begin(); iter != individuals.end(); ++iter) {
		for (int i=0; i<2; ++i) {
			vector<int> a = (**iter).get_sequence(i);
			for (auto iter2 = a.begin(); iter2 != a.end(); ++iter2) {
				++new_allele_counts[*iter2];
				(*alleles[*iter2]).increment_count();  // NOTE: alleles[*iter] returns a reference (mapped_type) to the pointer to the Allele object at position *iter2
																							// leading * dereferences the pointer, giving us access
			}
		}
	}

	for (auto iter3 = new_allele_counts.begin(); iter3 != new_allele_counts.end(); ++iter3)
		(*alleles[iter3->first]).set_count(iter3->second);

	// identify lost and fixed alleles and print to allele history file
	vector<int> to_remove;
	for (auto iter = alleles.begin(); iter != alleles.end(); ++iter) {
		int current_count = (*(iter->second)).get_count();
		if (current_count == 0) { // allele LOST from population
			to_remove.push_back(iter->first);  // first is position
			int birthgen = (*(iter->second)).get_birthgen();
			allele_file << iter->first << " " << birthgen << " " << gen - birthgen << " 0" << endl;
		}
		if (current_count == popsize*2) {  // derived allele FIXED in population
			to_remove.push_back(iter->first);
			int birthgen = (*(iter->second)).get_birthgen();
			for (auto iter2 = individuals.begin(); iter2 != individuals.end(); ++iter2)
				(**iter2).remove_fixed_allele(iter->first); // removed fixed allele's position from all individuals' sequences (currently stored in nextgen)
			allele_file << iter->first << " " << birthgen << " " << gen - birthgen << " 1" << endl;
		}
	}

	// free memory associated w/ lost/fixed alleles and remove entry from alleles container
	for (auto iter = to_remove.begin(); iter != to_remove.end(); ++iter) {
		delete alleles[*iter];  // free memory from Allele object itself
		alleles.erase(*iter);  // erase alleles map entry corresponding to the deleted allele
	}
}

void get_sample(int gen) {
	vector<bitset<bitlength>> sample;
	map<int, int> allele_counts; // note that a map is always sorted by keys
	int count = 0;
	for (auto iter = individuals.begin(); iter != individuals.begin()+sampsize; ++iter) { // determines which alleles are present in sample
		vector<int> haplotype = (**iter).get_sequence(0);
		for (auto iter2 = haplotype.begin(); iter2 != haplotype.end(); ++iter2)
			++allele_counts[*iter2];
	}
	for (auto iter = individuals.begin(); iter != individuals.begin()+sampsize; ++iter) {  // creates haplotypes for each sequence in the sample and populates bitset
		vector<int> haplotype = (**iter).get_sequence(0);
		sort(haplotype.begin(), haplotype.end());
		string hap;
		for (auto iter = allele_counts.begin(); iter != allele_counts.end(); ++iter)
			if ( binary_search (haplotype.begin(), haplotype.end(), iter->first))
				hap += "1";
			else
				hap += "0";
		sample.push_back(bitset<bitlength> (hap));
	}
	int S = allele_counts.size();
    pi_file << gen << " " << get_pi(sample) << endl;
    watterson_file << gen << " " << get_watterson(sample, S) << endl;
}

public:
void reproduce(int gen) {

	for (int i=0; i< popsize; ++i) {
		vector<int> parents;
		parents.push_back(randomind(e));
		parents.push_back(randomind(e));

		// create descendant of individuals parents[0] and parents[1]
		individuals.push_back( new Individual(individuals[parents[0]], individuals[parents[1]], mutate(parents, gen)) );
	}

	// delete dynamically allocated individuasl of the last generation
	for (auto iter = individuals.begin(); iter != individuals.end() - popsize; ++iter)
		delete *iter;
	// remove orphaned pointers from individuals
	individuals.erase(individuals.begin(), individuals.end()-popsize);

	// update allele counts on sample generations
	if (gen % sampfreq == 0 )
		update_alleles(gen);

	if ( gen != 0 && gen % sampfreq == 0 ) {
		random_shuffle(individuals.begin(), individuals.end() ) ;
		get_sample(gen);
	}
}

void close_output_files () {
	watterson_file.close();
	pi_file.close();
	allele_file.close();
}

Population () {
	// initialize random number distributions
	mu_sequence = seqlength * mutrate;
	randompos.param(uniform_int_distribution<int>::param_type(1,seqlength));
	randomind.param(uniform_int_distribution<int>::param_type(0,popsize-1));
	randomnum.param(uniform_real_distribution<double>::param_type(0.,1.));
 	randommut.param(poisson_distribution<int>::param_type(mu_sequence));

	individuals.reserve(popsize*2);
	for (int i=0; i<popsize; ++i) {
		vector<int> s1; vector<int> s2;
		vector<vector<int>> ses{s1,s2};
		individuals.push_back(  new Individual(ses)  );
	}

	string fname = "nucleotide_diversity";
	pi_file.open(fname.c_str());
	pi_file << "gen summary.stat" << endl;

	fname = "watterson_estimator";
	watterson_file.open(fname.c_str());
	watterson_file << "gen summary.stat" << endl;

	fname = "allele_info";
	allele_file.open(fname.c_str());
	allele_file << "position birthgen lifespan extinct.fixed" << endl;
}

static mt19937 e;

};

#endif
