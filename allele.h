#ifndef ALLELE_H
#define ALLELE_H

class Allele {

private:
	int position;
	int birthgen;
	int count;

public:
	inline int get_count() { return count; }
	inline void set_count(int ccount) { count = ccount; }
	inline int get_position() { return position; }
	inline void set_position(int pposition) { position = pposition; cout << pposition << ": " << position << endl; }
	inline void increment_count() { ++count; }
	inline int get_birthgen() { return birthgen; }

	// constructor
	Allele (int pos, int gen): position(pos), birthgen(gen) {
		count = 0;
	}
};

#endif
