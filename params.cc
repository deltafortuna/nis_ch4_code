#include "params.h"  // access to declarations of global parameter values

map<string, string> read_parameters_file(const string &parameters_fn)
{
	//map<int, map<string, string> > params_by_block;
	map<string,string> params;
	ifstream paramfile(parameters_fn.c_str());
	string line;
	while(getline(paramfile, line)) {
		istringstream iss(line.c_str());
		string key, nextone, value;
		iss >> key;
		while (iss >> nextone)
			value += nextone + " ";
		params[key] = value;
	}
	return params;
}

map<string, string> parameters = read_parameters_file("parameters"); // change name if parameters file not named parameters

vector<int> get_multi_int_param(const string &key)
{
	vector<int> vec;
	istringstream iss(parameters[key].c_str());
	string param;
	while(getline(iss, param, ' '))
		vec.push_back(atoi(param.c_str()));
	return vec;
}

vector<double> get_multi_double_param(const string &key)
{
	vector<double> vec;
	istringstream iss(parameters[key].c_str());
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
		{
			ps.push_back(cursize);
		}
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

// variable names
int popsize, sampsize, seqlength, sampfreq;
double mutrate;
bool useMS;
string mscommand;
vector<int> demography;
vector<double> dem_parameter;
vector<int> dem_start_gen;
vector<int> dem_end_gen;
vector<int> carrying_cap;

int process_parameters() {
	popsize = atoi(parameters["popsize"].c_str());
  mutrate = atof(parameters["mutrate"].c_str());
	sampsize = atoi(parameters["sampsize"].c_str());
	seqlength = atof(parameters["seqlength"].c_str()); // covernsion using atof() enables use of e notation in parameters file
	sampfreq = atoi(parameters["sampfreq"].c_str());
  useMS = atoi(parameters["useMS"].c_str());
  mscommand = parameters["mscommand"];
	demography = get_multi_int_param("demography");
	dem_parameter = get_multi_double_param("dem_parameter");
	dem_start_gen = get_multi_int_param("dem_start_gen");
	dem_end_gen = get_multi_int_param("dem_end_gen");
	carrying_cap = get_multi_int_param("carrying_cap");
	return 1;
}

int good_parameters = process_parameters();

vector<int> pop_schedule = create_pop_schedule();
