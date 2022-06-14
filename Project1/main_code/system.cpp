
#include "mpi.h"
#include "system.h"
#include "energy/kinetic.h"
#include "energy/harosc.h"
#include "energy/haroscjastrow.h"
#include "wavefunc/wavefunc.h"
#include "wavefunc/gauss.h"
#include "montec/montec.h"
#include "montec/metrop.h"
#include "montec/hastings.h"
#include "setstate/initstate.h"
#include "sampler.h"
#include <iostream>
#include <random>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <iomanip>
#include <string>

System::System (){ m_optimization = 0 ; }

void System::set_hamiltonian ( Kinetic * hamiltonian )
{
  m_hamiltonian = hamiltonian ;
}

void System::set_wave_func ( WaveFunc * wavefunc )
{
  m_wavefunc = wavefunc ;
}

void System::set_state ( InitState * initstate )
{
  m_initstate = initstate ;
  m_initstate -> state_generat () ;
}

void System::set_montec ( Montec * montec )
{
  m_montec = montec ;
  generat.seed ( clock()*1.4*m_my_rank + 0.2*m_my_rank*m_my_rank ) ;
}

//one-body density run
void System::onebody_rho_run ( double x_min , double x_max , double dx , double beta , int direc )
{
  //fixing the x-value of the position at x_min
  double x = x_min ;
  if ( direc == 0 )
  {
    *pos[0] = Lin(x,0,0,beta) ;
    *pos_old[0] = Lin(x,0,0,beta) ;
  }
  else if ( direc == 1 )
  {
    *pos[0] = Lin(0,0,x,beta) ;
    *pos_old[0] = Lin(0,0,x,beta) ;
  }

  m_montec -> set_start ( 1 ) ; //do not move particle 0 in the run
  m_montec -> set_seed () ;
  m_montec -> set_wf_old ( m_wavefunc -> get_wavefunc_sq () ) ;


  //equilibration
  for ( int i=0 ; i<m_equilibrate ; i++ )
  {
    m_montec -> test ( ) ;
  }

  //creating a Sampler object
  m_sampler = new Sampler ( m_num_of_cyc ) ;


  int num_of_cyc_per_proc = (int) m_num_of_cyc / m_numprocs ;
  while ( x <= x_max )
  {
    for ( int i=0 ; i<num_of_cyc_per_proc ; i++ )
    {
      m_montec -> test ( ) ;
      m_sampler -> sample_wf ( m_wavefunc -> get_wavefunc_sq () ) ;
    }

    m_sampler -> average () ;

    double tot_wf_mean    = 0 ;

    double wf_mean    = m_sampler -> get_average_wf    () ;
    MPI_Reduce ( & wf_mean    , & tot_wf_mean     , 1 , MPI_DOUBLE , MPI_SUM , 0 ,MPI_COMM_WORLD ) ;
    if (m_my_rank == 0)
    {
      tot_wf_mean /= m_numprocs ;
      std::cout << std::setprecision(5) << x << "       " << std::setprecision(15) << tot_wf_mean << std::endl ;
    }
    x += dx ;

    if ( direc == 0 )
    {
      *pos[0] = Lin(x,0,0,beta) ;
      *pos_old[0] = Lin(x,0,0,beta) ;
    }
    else if ( direc == 1 )
    {
      *pos[0] = Lin(0,0,x,beta) ;
      *pos_old[0] = Lin(0,0,x,beta) ;
    }
    m_sampler -> set_average_wf ( 0 ) ;
    m_sampler -> set_sum_wf     ( 0 ) ;
  }
  delete m_sampler ;
}

//run for optimization of alpha
void System::parameter_opt_run ( double alpha , double learning_rate , double eps , int num_of_opt_cyc , int max_iter , int equilibrate )
{
  m_optimization = 1 ;
  m_alpha = alpha ;
  double alpha_vals [max_iter+1] ;
  double temp_der_wf_energy ;
  double temp_der_wf ;
  double temp_energy ;

  m_montec -> set_seed () ;
  m_montec -> set_wf_old ( m_wavefunc -> get_wavefunc_sq () ) ;

  int num_of_cyc_per_proc = (int) num_of_opt_cyc / m_numprocs ;
  //m_montec -> set_energy ( m_hamiltonian -> local_energy () ) ;
  m_sampler = new Sampler ( num_of_cyc_per_proc*m_numprocs ) ;

  //equlibrating the system
  //std::cout<<std::setprecision(9) <<pos[1]->X()<<std::endl;
  for ( int i=0 ; i<equilibrate ; i++ )
  {
    m_montec -> test ( ) ;
  }

  for ( int i=0 ; i<max_iter ; i++ )
  {
    alpha_vals[i] = m_alpha ;
    for ( int j=0 ; j<num_of_cyc_per_proc ; j++ )
    {
      m_montec -> test ( ) ;
      m_sampler -> sample_der ( m_hamiltonian -> local_energy () , m_wavefunc -> get_der_alpha () ) ;
    }
    m_sampler -> average () ;

    temp_der_wf_energy = m_sampler -> get_average_der_wf_energy () ;
    temp_der_wf = m_sampler -> get_average_der_wf () ;
    temp_energy = m_sampler -> get_average_energy () ;

    double tot_temp_der_wf_energy = 0 ;
    double tot_temp_der_wf        = 0 ;
    double tot_temp_energy        = 0 ;

    MPI_Reduce ( & temp_der_wf_energy , & tot_temp_der_wf_energy  , 1 , MPI_DOUBLE , MPI_SUM , 0 ,MPI_COMM_WORLD ) ;
    MPI_Reduce ( & temp_der_wf        , & tot_temp_der_wf  , 1 , MPI_DOUBLE , MPI_SUM , 0 ,MPI_COMM_WORLD ) ;
    MPI_Reduce ( & temp_energy        , & tot_temp_energy  , 1 , MPI_DOUBLE , MPI_SUM , 0 ,MPI_COMM_WORLD ) ;
    if ( m_my_rank == 0 )
    {
      m_alpha -= learning_rate*2*( tot_temp_der_wf_energy - tot_temp_energy*tot_temp_der_wf ) ;
    }

    MPI_Bcast( & m_alpha , 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD) ;

    m_wavefunc -> set_alpha ( m_alpha ) ;
    m_hamiltonian -> set_alpha ( m_alpha ) ;
    if ( m_my_rank == 0 )
    {
      std::cout << "alpha: " << std::setprecision(15) << m_alpha << std::endl ;
    }
    if ( (std::abs(m_alpha - alpha_vals[i]) ) < eps )
    {
      break ;
    }
    m_sampler -> set_average ( 0 , 0 , 0 , 0 ) ;
    m_sampler -> set_sum     ( 0 , 0 , 0 , 0 ) ;
  }
  delete m_sampler ;
}

//standard MC run
void System::montec_run ( int printout_freq )
{
  //std::ofstream outfile;
  int num_of_cyc_per_proc = (int) m_num_of_cyc / m_numprocs ;
  std::string base ;
  std::string base_r ;
  base = ".dat" ;
  base_r = ".pos" ;
  for(int i=0;i<m_numprocs;++i){
    if ( m_my_rank == i )
    {
      std::ofstream(std::to_string(i)+base) ;
      std::ofstream(std::to_string(i)+base_r) ;
    }
  }
  m_montec -> set_seed () ;
  m_montec -> set_wf_old ( m_wavefunc -> get_wavefunc_sq () ) ;
  //m_montec -> set_energy ( m_hamiltonian -> local_energy () ) ;
  m_sampler = new Sampler ( num_of_cyc_per_proc*m_numprocs ) ;

  //equlibrating the system
  if ( m_optimization != 1 )
  {
    for ( int i=0 ; i<m_equilibrate ; i++ )
      {
        m_montec -> test ( ) ;
      }
  }

  //running montecarlo cycles
  //outfile.open("energy.dat") ;
  int counter_print = 0 ;
  for ( int i=0 ; i<num_of_cyc_per_proc ; i++ )
  {
    m_montec -> test ( ) ;
    m_sampler -> sample ( m_hamiltonian -> local_energy () ) ;
    m_sampler -> sample_acceptance ( m_montec -> get_acceptance() ) ;
    counter_print += 1 ;
    if ( counter_print == printout_freq )
    {
      std::ofstream(std::to_string(m_my_rank)+base, std::ios::app)
      << std::setprecision(15) << m_hamiltonian -> local_energy () << std::endl ;
      std::ofstream(std::to_string(m_my_rank)+base_r, std::ios::app)
      << std::setprecision(15) << pos[0] -> abs () << std::endl ;
    //  outfile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    //  outfile << std::setprecision(15) << m_hamiltonian -> local_energy () << std::endl ;
    counter_print = 0 ;
    }
  }
  m_sampler -> average () ;
  double acceptance_rate = m_sampler -> acceptance_rate ( m_num_of_part ) ;
  double tot_acceptance_rate = 0 ;

  double tot_energy_mean    = 0 ;
  double tot_energy_mean_sq = 0 ;

  double energy_mean    = m_sampler -> get_average_energy    () ;
  double energy_mean_sq = m_sampler -> get_average_energy_sq () ;

  MPI_Reduce ( & energy_mean    , & tot_energy_mean     , 1 , MPI_DOUBLE , MPI_SUM , 0 ,MPI_COMM_WORLD ) ;
  MPI_Reduce ( & energy_mean_sq , & tot_energy_mean_sq  , 1 , MPI_DOUBLE , MPI_SUM , 0 ,MPI_COMM_WORLD ) ;
  MPI_Reduce ( & acceptance_rate , & tot_acceptance_rate  , 1 , MPI_DOUBLE , MPI_SUM , 0 ,MPI_COMM_WORLD ) ;

  if ( m_my_rank == 0 )
  {
    m_sampler -> set_average ( tot_energy_mean , tot_energy_mean_sq , zero , zero ) ;
    m_sampler -> print () ;
    std::cout << "Acceptance rate: " << std::setprecision(7) << tot_acceptance_rate << std::endl ;
  }
  //outfile.close() ;
  for(int i=0;i<m_numprocs;++i)
    if ( m_my_rank == i )
    {
      std::ofstream(std::to_string(i)+base,std::ios::app).close();
      std::ofstream(std::to_string(i)+base_r,std::ios::app).close();
    }
  delete m_sampler ;
}

double System::get_wavefunc_sq  () { return (double) (m_wavefunc -> get_wavefunc_sq ()) ; }
double System::get_local_energy () { return (double) (m_hamiltonian -> local_energy ()) ; }
Lin    System::get_drift_force  ( int k ) { return m_hamiltonian -> drift_force ( k )   ; }


void System::set_equilibrate  ( double  equilibrate     ) { m_equilibrate = equilibrate ; }
void System::set_num_of_dim   ( int     num_of_dim      ) { m_num_of_dim  = num_of_dim  ; }
void System::set_num_of_cyc   ( int     num_of_cyc      ) { m_num_of_cyc  = num_of_cyc  ; }
void System::set_step_length  ( double  step_length     ) { m_step_length = step_length ; }
void System::set_mpi          ( int my_rank , int numprocs  ) { m_my_rank     = my_rank     ; m_numprocs = numprocs           ; }
void System::set_num_of_part  ( int     num_of_part         ) { m_num_of_part = num_of_part ; num_of_particles = num_of_part  ; }
