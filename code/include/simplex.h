#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <vector>
#include <unordered_set>

using namespace std;

enum class Mode{
	MIN, MAX
};

class Tableau {
public:
	Mode obj;
	size_t varnum;
	vector<vector<double>> mat;

	Tableau();

	// Mode m
	// MIN = MINIMIZE
	// MAX = MAXIMIZE
	Tableau( Mode m, vector<double> vars );

	void add_constraint( vector<double> vars );

	size_t row_size();

	size_t col_size();
};

vector<double> get_solution( Tableau &tab );

bool simplex( Tableau &tab );

#endif
