#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "allele.h"

class Individual {

private:
	vector<vector<int>> sequences;
	void remove_allele_by_position (int seqnum, int position) {
		auto pos = find(sequences[seqnum].begin(), sequences[seqnum].end(), position);
		if (pos != sequences[seqnum].end())  // ensures position in vector
			sequences[seqnum].erase(pos);
	}

public:
	inline vector<int> get_sequence(int whichseq) { return sequences[whichseq]; }
	void remove_fixed_allele(int to_remove) {
		for (int i = 0; i<2; ++i) {
			vector<int>::iterator p = find(sequences[i].begin(), sequences[i].end(), to_remove);
			if (p != sequences[i].end()) // i.e., element not found
		    		sequences[i].erase(p);
		 }
	}

/*	vector<int> get_alleles() {
		vector<int> a = sequences[0];
		a.insert(a.end(), sequences[1].begin(), sequences[1].end());
		return(a);
	}
*/

	// generation 0 constructor
	Individual (vector<vector<int>> seqs): sequences(seqs) {
			;
	}

	// intra-simulation constructor
	Individual (Individual *p1, Individual *p2, vector<vector<int> > mutation_results) {
		sequences.push_back((*p1).get_sequence(mutation_results[2][0]));
		sequences.push_back((*p2).get_sequence(mutation_results[3][0]));
		for (int i=0; i<2; ++i) {
			for (int j=1; j<mutation_results[i].size(); ++j)  {
				if (mutation_results[i][j] > 0)
					sequences[i].push_back(mutation_results[i][j]);
				else
					remove_allele_by_position(i, -1 * mutation_results[i][j]);
			}
		}
	}

	~Individual() {} // destructor
};

#endif 
