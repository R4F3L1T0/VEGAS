#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <numeric>
#include <cmath>

using namespace std;

clock_t start = clock(); 

#define TURNTIME float(clock()-start)/CLOCKS_PER_SEC

/////////////////
/// FUNCTIONS ///
/////////////////

// Integrating function
double f(double x) {
	return sin(x);
}


int rnd(int a) {
	return rand()%a;
}

void save(string file_name, double a) {
	ofstream f;
	f.open(file_name);
	f<<a<<endl;
}

void save(string file_name,vector<double> vec) {

	ofstream f;
	f.open(file_name);
	for (int i = 0; i < vec.size(); ++i)	f<<vec[i]<<endl;
}

// Calculates the standard deviation given a vector of measurements
double calculate_st_dev(int term, vector<double> v) {
	double st_dev = 0;
	double mean_integral_value = accumulate(v.begin() + term, v.end(),0.0)/( (double) v.size());
	for (int i = term; i < v.size(); ++i) st_dev += (v[i] - mean_integral_value)*(v[i] - mean_integral_value);
	st_dev = sqrt(st_dev/(v.size() -1));

	return st_dev;
}

///////////////
/// CLASSES ///
///////////////
class Vegas{
public:

	int N, K;

	double a,b;

	vector<double> lattice;
	vector<double> p;

	double integral_value = 0;

	Vegas(int N, int K, double a, double b)
	{
		this->N = N;
		this->K = K;
		this->a = a;
		this->b = b;
	}

	// Initialize the lattice with equally spaced intervals
	void initialization()
	{
		if (a>b) swap(a,b);

		for (int i = 0; i < N; ++i) lattice.push_back(a + (b-a)*i/N);
		lattice.push_back(b);
	}

	// Heat the lattice and update it 
	void update_lattice()
	{
		vector<int> m; 							// number of subintevals for each interval
		vector<double> f_(N,0.0);   			// not normalized weight of each interval

		vector<double> temp_lattice;			// heated lattice with K subintervals
		temp_lattice.push_back(lattice[0]);
		
		// Calculating the ot normalized weight of each interval
		for (int i = 0; i < N*100; ++i)
		{
			int r = rnd(N);
			double x = lattice[r] + ((double) rand() / RAND_MAX )*(lattice[r+1] - lattice[r]);
			f_[r] += abs(f(x))*(lattice[r+1] - lattice[r]);
		}
		
		double sum_f_ = accumulate(f_.begin(),f_.end(),0.0);  // f_ normalization

		// Calculation of subintervals for each interval
		m.clear();
		for (int i = 0; i<N; ++i)	m.push_back( (int) (double) K*f_[i]/sum_f_); 
		
		int sum = accumulate(m.begin(),m.end(),0);
		
		// Adding the missing subintervals. This lack is due to the approximation of m as int
		for (int i = 0; i < K-sum; ++i)
		{
			int tmp = rnd(N);
			m[tmp] += 1; 		
		}
		
		// Calculating the temporary lattice with K subintervals
		for (int i = 0; i < lattice.size() - 1; ++i)
		{
			if (m[i] == 0) ++m[i];
			for (int j = 0; j < m[i]; ++j)	temp_lattice.push_back(temp_lattice.back() + (lattice[i+1] - lattice[i])/m[i]);
		}

		// Final lattice and probability update
		lattice.clear();
		for (int i = 0; i < N; ++i) lattice.push_back(temp_lattice[i*K/N]);
		lattice.push_back(b);

		p.clear();
		for (int i = 0; i < lattice.size()-1; ++i) p.push_back(1/(N*(lattice[i+1] - lattice[i])));
	}

	void calculate_integral()
	{
		integral_value = 0;
		for (int i = 0; i < N*100; ++i)
		{
			int r = rnd(N);
			double x = lattice[r] + ((double) rand() / RAND_MAX )*(lattice[r+1] - lattice[r]);
			integral_value +=f(x)/p[r];  
		}
		integral_value/=(N*100);
	}

};	

////////////
/// MAIN ///
////////////
int main()
{
	srand((unsigned)time(0));

	int lattice_bin_number = 100;						// Number of intervals inside the integration domain. The higher the more accurate
	int sub_interval_number = 100000;				    	// Number of subintervals, MUST be a multiple of lattice_bin_number
	int epochs=1000;						 			// Number of lattice warm ups
	int term=100;										// Termalizations steps

	vector<double> integral_values;						// Vector with integral calculations

	double interval_start = 0, interval_end = 1*M_PI; 	// Integration domain

 	Vegas vegas(lattice_bin_number, sub_interval_number, interval_start, interval_end);

 	vegas.initialization();

 	for (int i = 0; i < epochs; ++i)
 	{	
 		vegas.update_lattice();
 		vegas.calculate_integral();
 		integral_values.push_back(vegas.integral_value);
 		
 	}

 	double mean_integral_value = accumulate(integral_values.begin() + term, integral_values.end(),0.0)/( (double) integral_values.size());
 	double st_dev = calculate_st_dev(term,integral_values);

 	cout<<"integral value "<<mean_integral_value<<"		Standard Error "<<st_dev<<endl;

 	save("results/probability.txt",vegas.p);
 }