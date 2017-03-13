#ifndef _the_generic_solver_h_included_
#define _the_generic_solver_h_included_

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>
#include <fixed_grid.h>
#include <param_maps.h>

namespace PJLib {
  
  class GenericSolver
  {
    std::string solver_name;
    ParameterMaps *solver_params;
    void link_parameters(ParameterMaps *params);
    void erase_parameters();
    
    protected:
      
      bool countdown;
      // solve over some time period
      Fixed_Grid t;
      // parameters in time
      int timesteps;
      double t_min,t_max,theta;
      // solve at a timestep
      virtual void solve_at_timestep(double tau);
      
    public:
      // constructors
      GenericSolver(ParameterMaps *params,std::string filename,bool show_countdown);
      virtual ~GenericSolver();
      // initial conditions
      virtual void update_time();
      virtual void update_storage();
      virtual void initial_conditions();
      virtual void end_conditions();
      virtual void solve();
      virtual double value();
      std::string name(){return solver_name;};
      virtual std::string store_name(){return solver_name;};
  };
  
  class GenericSolverInt: public GenericSolver
  {
    protected:
      
      virtual void solve_at_timestep(int k);
      
    public:
      
      GenericSolverInt(ParameterMaps *params,std::string filename,bool show_countdown);
      virtual ~GenericSolverInt(){};
      // solve with ints instead of doubles
      virtual void solve();
      
  };
  
  class MultiCommitSolver
  {
    std::string solver_name;
    void link_parameters(ParameterMaps *params);
    void erase_parameters();
    void clear_commits();
    // keep a check on whether storage has been created
    bool storage_created;
    
    protected:
      
      int commits;
      ParameterMaps *solver_params;
      // storage
      ParameterMaps **param_array;
      GenericSolver **solver_array;
      // solve some contract time period
      Fixed_Grid t,C;
      // parameters in time
      int no_of_contracts,no_of_commits;
      double contract_min,contract_max,commit_min,commit_max;
      // create storage for the solver
      virtual void choose_solver(int n);
      virtual void create_storage(int n);
      // copy parameter choices across to all contracts
      void update_parameters();
      // set new parameters for the run, note that parameters are reset after updating the storage
      virtual void set_new_parameters();
      // solve for one contract
      virtual void solve_a_contract();
      
    public:
      
      // constructors
      MultiCommitSolver(ParameterMaps *params,std::string filename);
      virtual ~MultiCommitSolver();
      // initial conditions
      virtual void initial_conditions();
      // conditions applied at start of contract
      virtual void begin_conditions();
      // conditions applied at end of contract
      virtual void end_conditions();
      // update storage
      virtual void update_storage();
      // solve for multiple contracts
      virtual void solve();
      // return some value in the state space
      virtual double value();
      // return the name of the solver
      std::string name(){return solver_name;};
      std::string solve_name(){if(storage_created)return solver_array[0]->name();else return " ";};
      
  };
  
  class MultiBoostSolver: public MultiCommitSolver
  {
    
    protected:
      
      // solve for one contract using threads
      void solve_a_contract();
      void solve_all_threads();
      void solve_a_thread(int istart,int ifinish);
    
    public:
      
      MultiBoostSolver(ParameterMaps *params,std::string filename);
      virtual ~MultiBoostSolver();
        
  };
  
}

#endif
