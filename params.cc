#include "params.h"  // access to declarations of global parameter values

vector<int> get_multi_int_param(const string &parameters_fn, const string &key)
{
	vector<int> vec;
	istringstream iss(key);
	string param;
	while(getline(iss, param, ' '))
		vec.push_back(atoi(param.c_str()));
	return vec;
}

vector<double> get_multi_double_param(const string &key)
{
	vector<double> vec;
	istringstream iss(key.c_str());
	string param;
	while(getline(iss, param, ' '))
		vec.push_back(atof(param.c_str()));
	return vec;
}

vector<int> create_pop_schedule()
{
	vector<int> ps;
	int i=0;
	int cursize = popsize;
	for (int step = 0; step < demography.size(); ++step) {
		for (; i<dem_start_gen[step]; ++i)
			ps.push_back(cursize);
		for (; i <= dem_end_gen[step]; ++i) {
			switch(demography[step]) {
				case 0: ps.push_back(cursize);  // no size change
					    break;
				case 1: if (i == dem_start_gen[step])
							cursize += dem_parameter[step];
						ps.push_back(cursize); // instantaneous
						break;
				case 2: cursize += dem_parameter[step];
						ps.push_back(cursize);  // linear
						break;
				case 3: cursize *= exp(dem_parameter[step]);  // exponential
						ps.push_back(cursize);
						break;
				case 4: cursize = (carrying_cap[step] * cursize) /
									(cursize + (carrying_cap[step] - cursize)*exp(-1*dem_parameter[step]));  // logistic
						ps.push_back(cursize);
			}
		}
	}
	return ps;
}

// returns value of type int for requested parameter key
int getParameter_int(const string &parameters_fn, const string &key)
{
	string item;
	ifstream paramfile;
	int to_return;
	paramfile.open(parameters_fn.c_str());

	while (paramfile >> item)
	{
	    if (item == key) {
	    	paramfile >> item;
	    	paramfile.close();
				to_return = atoi(item.c_str());
	 //   	return(atoi(item.c_str()));
	    } else
	    	paramfile >> item; // read value, which is ignored
	}
	paramfile.close();
	return(to_return);
}

// returns value of type double for requested parameter key
double getParameter_double(const string &parameters_fn, const string &key)
{
	string item;
	ifstream paramfile;
	double to_return;
	paramfile.open(parameters_fn.c_str());
	// Read in a line
	while (paramfile >> item)
	{
	    if (item == key) {
	    	paramfile >> item;
	    	paramfile.close();
				to_return = atof(item.c_str());
	    	//return (atof(item.c_str()));
	    } else
	    	paramfile >> item; // read value, which is ignored
	}
	paramfile.close();
	return(to_return);
}

string parameters_fn("parameters"); // the name of parameters file

// read in the parameter values from the parameters file
int popsize = getParameter_int (parameters_fn, "popsize");
double mutrate = getParameter_double (parameters_fn, "mutrate");
int sampsize = getParameter_int (parameters_fn, "sampsize");
int seqlength = getParameter_double (parameters_fn, "seqlength");
int sampfreq = getParameter_int (parameters_fn, "sampfreq");
vector<int> demography =  get_multi_int_param("demography");
vector<double> dem_parameter = get_multi_double_param("dem_parameter");
vector<int> dem_start_gen =  get_multi_int_param("dem_start_gen");
vector<int> dem_end_gen = get_multi_int_param("dem_end_gen");
vector<int> carrying_cap = get_multi_int_param("carrying_cap");

vector<int> pop_schedule = create_pop_schedule();
