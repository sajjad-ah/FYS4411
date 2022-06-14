using namespace std ;

#include "mpi.h"
#include "system.h"
#include "energy/kinetic.h"
#include "energy/harosc.h"
#include "energy/haroscjastrow.h"
#include "readfile/read_input.h"
#include "setstate/initstate.h"
#include "setstate/unirand.h"
#include "wavefunc/wavefunc.h"
#include "wavefunc/gauss.h"
#include "wavefunc/gaussjastrow.h"
#include "montec/montec.h"
#include "montec/metrop.h"
#include "montec/hastings.h"
#include "lin.h"
#include <iostream>
#include <ctime>

int main ( int argc , char ** argv )
{
  //initializing MPI
  int my_rank ; int numprocs ;

  MPI_Init ( &argc , &argv ) ;
  MPI_Comm_size ( MPI_COMM_WORLD , &numprocs ) ;
  MPI_Comm_rank ( MPI_COMM_WORLD , &my_rank) ;

  //reading input from file with Read class
  Read input ( argv[1] ) ;

  //defining values of variables from input file
  double  alpha           = atof ( input.value ( "alpha"             ) .c_str() ) ;
  double  beta            = atof ( input.value ( "beta"              ) .c_str() ) ;
  double  step_length     = atof ( input.value ( "step_length"       ) .c_str() ) ;
  double  equilibrate     = atof ( input.value ( "equilibration"     ) .c_str() ) ;
  double  learning_rate   = atof ( input.value ( "learning_rate"     ) .c_str() ) ;
  double  eps             = atof ( input.value ( "eps"               ) .c_str() ) ;
  double  num_of_opt_cyc  = atof ( input.value ( "num_of_opt_cycles" ) .c_str() ) ;
  double  max_iter        = atof ( input.value ( "max_iterations"    ) .c_str() ) ;
  double  a               = atof ( input.value ( "jastrow_parameter" ) .c_str() ) ;
  double  D               = atof ( input.value ( "diffusion_coeff"   ) .c_str() ) ;
  double  x_min           = atof ( input.value ( "onebody_min"       ) .c_str() ) ;
  double  x_max           = atof ( input.value ( "onebody_max"       ) .c_str() ) ;
  double  dr              = atof ( input.value ( "onebody_dr"        ) .c_str() ) ;
  int     num_of_part     = atof ( input.value ( "num_of_particles"  ) .c_str() ) ;
  int     num_of_dim      = atof ( input.value ( "num_of_dimensions" ) .c_str() ) ;
  int     num_of_cyc      = atof ( input.value ( "num_of_cycles"     ) .c_str() ) ;
  int     printout_freq   = atof ( input.value ( "printout_freq"     ) .c_str() ) ;

  //initializing the system class
  System * system = new System ;
  system -> set_state ( new Unirand ( system , num_of_dim , num_of_part , beta ) ) ;
  system -> set_mpi           ( my_rank , numprocs  ) ;
  system -> set_equilibrate   ( equilibrate         ) ;
  system -> set_num_of_dim    ( num_of_dim          ) ;
  system -> set_num_of_part   ( num_of_part         ) ;
  system -> set_num_of_cyc    ( num_of_cyc          ) ;
  system -> set_step_length   ( step_length         ) ;

/*************************************************/

  //choosing (from input file) whether or not to use importance sampling
  if ( input.value( "hastings" ) == "no" )
  {
    system -> set_montec  ( new Metrop ( system , my_rank , num_of_dim , num_of_part , step_length ) ) ;
  }
  else if ( input.value( "hastings" ) == "yes" )
  {
    system -> set_montec  ( new Hastings ( system , my_rank , num_of_dim , num_of_part , step_length , D ) ) ;
  }
  else
  {
    cout << "Input parameter hastings takes values yes or no." << endl ;
    return EXIT_FAILURE ;
  }

/*************************************************/

  //choosing (from input file) whether or not to use hard-sphere potential
  if ( input.value( "jastrow" ) == "no" )
  {
    system -> set_hamiltonian ( new Harosc ( system , alpha , beta , num_of_dim , num_of_part ) ) ;
    system -> set_wave_func ( new Gauss ( system , num_of_dim , num_of_part , alpha ) ) ;
  }
  else if ( input.value( "jastrow" ) == "yes" )
  {
    system -> set_hamiltonian ( new HaroscJastrow ( system , alpha , beta , a , num_of_dim , num_of_part ) ) ;
    system -> set_wave_func ( new GaussJastrow ( system , num_of_dim , num_of_part , alpha , a ) ) ;
  }
  else
  {
    cout << "Input parameter jastrow takes values yes or no." << endl ;
    return EXIT_FAILURE ;
  }

/*************************************************/

  //choosing run-type, either onebody density, optimization of alpha, or standard MC run.
  if  ( input.value( "run_type" ) == "onebody_density" )
  {
    int direc = 0 ;
    if ( input.value( "direction" ) == "x" ){ direc = 0 ; }
    else if ( input.value( "direction" ) == "z" ){ direc = 1 ; }
    system -> onebody_rho_run ( x_min , x_max , dr , beta , direc ) ;
    MPI_Finalize () ;
    return 0 ;
  }

  if ( input.value( "run_type" ) == "parameter_opt" )
  {
    system -> parameter_opt_run ( alpha , learning_rate , eps , num_of_opt_cyc , max_iter , equilibrate ) ;
    MPI_Finalize () ;
    return 0 ;
  }

  if ( input.value( "run_type" ) == "standard" )
  {
    double time_i ;
    double time_f ;
    if ( my_rank == 0 ) { time_i = clock () ; }
    system -> montec_run ( printout_freq ) ;
    if ( my_rank == 0 ) { time_f = clock () ; cout << "CPU time: " << time_f - time_i << endl ; }
    MPI_Finalize () ;
    return 0 ;
  }
}
