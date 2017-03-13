
#ifndef generic_SDE
#define generic_SDE

#include <vector>
#include <fstream>

class Generic_SDE{
	
	private:
		
		static constexpr double Pi=3.141592653589793238462643383279502884197l;
		
		int no_of_runs,no_of_points,seeder;
		double beta,dt,T,X_start,**normal_dist,*time;
		double alpha(double X,double t),theta(double X,double t),sigma(double X,double t),forcing(double X,double t);
		double alpha_X,average_X,amp_X,freq_X,phase_X,sigma_X;
		void initial_conditions(double _start);
		void gas_dev(int seed,int runs);
		void generate_dist(void);
		void print_dist(void);
		double mean_dist(double *distribution,int n);
		double variance_dist(double *distribution,int n);
		double Power_Stack(double yy);
		
	public:
		
		~Generic_SDE();
		Generic_SDE();
		double **X_t;
		Generic_SDE(int _points,int _seed,double _T,double _start,double _beta);
		Generic_SDE(int _runs,int _points,int _seed,double _T,double _start,double _beta,double _alpha_X,double _average_X,double _amp_X,double _phase_X,double _sigma_X);
		void solve(void);
		
};

#endif
