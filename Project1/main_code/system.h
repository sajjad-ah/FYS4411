#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include <random>
#include "lin.h"
//#include "lin.h"
//#include "particles.h"

class System
{
private:
  class Kinetic *         m_hamiltonian = nullptr ;
  class WaveFunc *        m_wavefunc    = nullptr ;
  class Montec *          m_montec      = nullptr ;
  class InitState *       m_initstate   = nullptr ;
  class Sampler *         m_sampler     = nullptr ;

  int m_num_of_dim  ;
  int m_num_of_part ;
  int m_equilibrate ;
  int m_num_of_cyc ;
  int m_my_rank ;
  int m_numprocs ;
  int m_optimization ;
  double m_alpha ;
  double m_energy ;
  double m_step_length ;
  double m_wf ;
  double m_wf_old ;
  double m_wf_rel ;
  double rando ;
  double zero = 0.0 ;


  std::mt19937_64 generat ;
public:
  System ( ) ;
  void set_hamiltonian    ( Kinetic *               ) ;
  void set_wave_func      ( WaveFunc *              ) ;
  void set_state          ( InitState *             ) ;
  void set_montec         ( Montec *                ) ;
  void set_num_of_dim     ( int                     ) ;
  void set_num_of_part    ( int                     ) ;
  void set_num_of_cyc     ( int                     ) ;
  void set_mpi            ( int , int               ) ;
  void set_equilibrate    ( double                  ) ;
  void set_step_length    ( double                  ) ;
  void montec_run         ( int                     ) ;
  void metropolis_test    (                         ) ;
  void parameter_opt_run  ( double , double , double , int , int , int ) ;
  void onebody_rho_run    ( double , double , double , double , int ) ;
  double get_wavefunc_sq  ()  ;
  double get_local_energy ()  ;
  Lin get_drift_force ( int ) ;



  int num_of_particles ;

  std::vector < class Lin * > pos ;
  std::vector < class Lin * > pos_old ;
};
#endif
