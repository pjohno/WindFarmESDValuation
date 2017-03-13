#include <generic_solver.h>
#include <countdown.h>

using namespace PJLib;

void GenericSolver::link_parameters(ParameterMaps *params)
{
  params->assignParameter(&timesteps,solver_name+"-timesteps"," No of timesteps : "+solver_name+" :: ",101);
  params->assignParameter(&t_min,solver_name+"-t_min"," Option start time : "+solver_name+" :: ",0.);
  params->assignParameter(&t_max,solver_name+"-t_max"," Option end time: "+solver_name+" :: ",24.);
  params->assignParameter(&theta,solver_name+"-theta"," Difference Scheme : "+solver_name+" :: ",0.5);
  params->reset();
}

void GenericSolver::erase_parameters()
{
  solver_params->eraseParameter(solver_name+"-timesteps");
  solver_params->eraseParameter(solver_name+"-t_min");
  solver_params->eraseParameter(solver_name+"-t_max");
  solver_params->eraseParameter(solver_name+"-theta");
}

void GenericSolver::solve_at_timestep(double tau)
{
}

// constructors
GenericSolver::GenericSolver(ParameterMaps *params,std::string filename,bool show_countdown)
{
  solver_name = filename;
  solver_params = params;
  countdown = show_countdown;
  link_parameters(params);
}

GenericSolver::~GenericSolver()
{
//   std::cout << " delete GenericSolver \n";
  erase_parameters();
//   std::cout << " GenericSolver deleted \n";
}
// test
void GenericSolver::update_storage(){}
// initial conditions
void GenericSolver::initial_conditions(){}
void GenericSolver::end_conditions(){}
// generic solve algorithm runs through 
void GenericSolver::solve(void)
{
  Countdown timer;
  // 
  update_time();
  // solve through time
  double tau;
  if(countdown){std::cout << "\n";timer.update(0,t.size()-2);}
  for(int k=0;k<timesteps-1;k++)
  {
    tau = t[k]+theta*t.dx();
//     std::cout << "solving at time " << t[k+1]  << "\n";
    solve_at_timestep(tau);
    if(countdown)timer.update(k+1,t.size()-1);
  }
  if(countdown)std::cout << "\n";
//   end_conditions();
}

void GenericSolver::update_time()
{
  t.setup_x(timesteps,t_min,t_max); 
}

double GenericSolver::value(){return 0.;}

// ------------------------------------------------
//                                                 
// ------------------------------------------------
// Generic solver with ints
// ------------------------------------------------
void GenericSolverInt::solve_at_timestep(int k)
{
}

GenericSolverInt::GenericSolverInt(ParameterMaps *params,std::string filename,bool show_countdown):GenericSolver(params,filename,show_countdown)
{
}

void GenericSolverInt::solve(void)
{
  Countdown timer;
  // 
  update_time();
  // solve through time
  if(countdown){std::cout << "\n";timer.update(0,t.size()-2);}
  for(int k=0;k<timesteps-1;k++)
  {
    solve_at_timestep(k);
    if(countdown)timer.update(k+1,t.size()-1);
  }
  if(countdown)std::cout << "\n";
  //   end_conditions();
}


// ------------------------------------------------
//                                                 
// ------------------------------------------------
// Multi commit solve
// ------------------------------------------------

void MultiCommitSolver::link_parameters(ParameterMaps *params)
{
  params->assignParameter(&no_of_commits,solver_name+"-no_of_commits"," No of commits : "+solver_name+" :: ");
  no_of_commits=11;
  params->assignParameter(&no_of_contracts,solver_name+"-no_of_contracts"," No of contracts : "+solver_name+" :: ");
  no_of_contracts=25;
  params->assignParameter(&contract_min,solver_name+"-contract_min"," Contract start time : "+solver_name+" :: ");
  contract_min=0.;
  params->assignParameter(&contract_max,solver_name+"-contract_max"," Contract end time: "+solver_name+" :: ");
  contract_max=24.;
  params->assignParameter(&commit_min,solver_name+"-commit_min"," Minimum commitment : "+solver_name+" :: ");
  commit_min=0.;
  params->assignParameter(&commit_max,solver_name+"-commit_max"," Maximum commitment : "+solver_name+" :: ");
  commit_max=1.;
  params->reset();
}

void MultiCommitSolver::erase_parameters()
{
  solver_params->eraseParameter(solver_name+"-no_of_contracts");
  solver_params->eraseParameter(solver_name+"-contract_min");
  solver_params->eraseParameter(solver_name+"-contract_max");
  solver_params->eraseParameter(solver_name+"-no_of_commits");
  solver_params->eraseParameter(solver_name+"-commit_min");
  solver_params->eraseParameter(solver_name+"-commit_max");
}

MultiCommitSolver::MultiCommitSolver(ParameterMaps *params,std::string filename)
{
  // set object name and parameters
  solver_params = params;
  solver_name = filename;
  // link parameters
  link_parameters(params);
  storage_created = false;
}

MultiCommitSolver::~MultiCommitSolver()
{
  clear_commits();
  erase_parameters();
}

void MultiCommitSolver::choose_solver(int n)
{
}

void MultiCommitSolver::create_storage(int n)
{
//   std::cout << " create storage \n";
  param_array = new ParameterMaps*[n];
  solver_array = new GenericSolver*[n];
  // link first parameter class to global parameters
  param_array[0] = solver_params;
  // create new instances of the parameter class for the rest of the array
  for(int i=1;i<n;i++)param_array[i] = new ParameterMaps;
  for(int i=0;i<n;i++)solver_array[i] = new GenericSolver(param_array[i],solver_name,false);
  // store number of commits
  commits = n;
}

void MultiCommitSolver::update_storage()
{
//   std::cout << " update storage " << storage_created << "\n";
  if(!storage_created)
  {
//     std::cout << " new storage \n";
    create_storage(no_of_commits);
    storage_created=true;
  }
  else if(no_of_commits!=commits)
  {
//     std::cout << " renew storage \n";
    clear_commits();
    create_storage(no_of_commits);
  }
  C.setup_x(no_of_commits,commit_min,commit_max);
  t.setup_x(no_of_contracts,contract_min,contract_max);
//   std::cout << " storage updated \n";
}

void MultiCommitSolver::clear_commits()
{
//   std::cout << " clear commits\n ";
  // delete parameter classes (apart from global parameter class)
  for(int i=0;i<commits;i++)delete solver_array[i];
  for(int i=1;i<commits;i++)delete param_array[i];
  delete [] param_array;
  delete [] solver_array;
//   std::cout << " cleared commits\n ";
}

void MultiCommitSolver::update_parameters()
{
//   std::cout << " sync all parameter values to the global parameter set \n";
  for(int i=1;i<commits;i++)*param_array[i] = *param_array[i-1];
//   std::cout << " set up individual commits for each solver \n";
  for(int i=0;i<commits;i++)param_array[i]->changeParameter(C[i],solver_array[i]->store_name()+"-C_commit");
}

void MultiCommitSolver::set_new_parameters()
{
  param_array[0]->changeParameter(11,solver_array[0]->name()+"-timesteps");
}

// solve a contract 
void MultiCommitSolver::solve_a_contract()
{
  update_parameters();
  begin_conditions();
  for(int i=0;i<commits;i++)
  {
    //std :: cout << " \n ### C="<< C[i] << "\n";
    //param_array[i]->print_all_params(&std::cout);
    solver_array[i]->solve();
  }
  end_conditions();
}

void MultiCommitSolver::begin_conditions()
{
  
}

void MultiCommitSolver::initial_conditions()
{
  
}

void MultiCommitSolver::solve()
{
//   std::cout << " update storage \n";
  update_storage();
//   std::cout << " set new parameters \n";
  set_new_parameters();
  update_parameters();
//   std::cout << " initial conditions \n";
  initial_conditions();
  
  #ifdef DEBUG
  bool countdown=true;
  #else
  bool countdown=false;
  #endif
  Countdown timer;
  // 
  t.setup_x(no_of_contracts,contract_min,contract_max);
  // solve through time
  if(countdown){std::cout << "\n";timer.update(0,t.size()-1);}
  for(int k=1;k<t.size();k++)
  {
    param_array[0]->changeParameter(t[k-1],solver_array[0]->name()+"-t_min");
    param_array[0]->changeParameter(t[k],solver_array[0]->name()+"-t_max");
    solve_a_contract();
    if(countdown)timer.update(k,t.size()-1);
  }
  if(countdown)std::cout << "\n";
  //end_conditions();
}



void MultiCommitSolver::end_conditions()
{
}

double MultiCommitSolver::value()
{
  return 0.;
}

// boost solver 


MultiBoostSolver::MultiBoostSolver(ParameterMaps *params,std::string filename):MultiCommitSolver(params,filename){};

MultiBoostSolver::~MultiBoostSolver(){};

void MultiBoostSolver::solve_a_thread(int istart,int ifinish)
{
  for(int id=istart;id<ifinish;id++)
  {
    solver_array[id]->solve();
  }
}
// solving for all threads not worrying about how they are solved due to way code is parallizable
void MultiBoostSolver::solve_all_threads()
{
  int numberOfThreads = commits,remainder,base,current,istart,ifinish;
  
  boost::thread_group worker_threads;
  
  current = 0;
  base = commits/numberOfThreads;
  remainder = commits % numberOfThreads;
//   std :: cout << " base " << base << " remainder " << remainder << "\n";
  
  for(int i=0;i<numberOfThreads;i++)
  {
    // work out how to portion the work between threads
    istart = current;
    ifinish = current + base + (i<remainder);
//     std::cout << " create a thread  " << i << " ... \n";
//     std::cout << " working from commit " << istart << " to " << ifinish << "\n";
    worker_threads.create_thread(boost::bind(&MultiBoostSolver::solve_a_thread,this,istart,ifinish));
    current = ifinish;
  }
  
  worker_threads.join_all();
  
  // test
  
}

// here we solve for a single hour contract
void MultiBoostSolver::solve_a_contract()
{
  update_parameters();
  begin_conditions();
  solve_all_threads();
  end_conditions();
}
