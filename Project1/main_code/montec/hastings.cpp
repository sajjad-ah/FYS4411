
#include "hastings.h"
#include <iostream>
#include <cmath>
#include <iomanip>

Hastings::Hastings ( System * system , int my_rank , int num_of_dim , int num_of_part , double step_length , double D ) : Montec ( system , my_rank , num_of_dim , num_of_part , step_length )
{
  m_step_length_sqrt = sqrt ( step_length ) ;
  m_D = D ;
}

void Hastings::test ( )
{
  m_acceptance = 0 ;
  for ( int k=0 ; k<m_num_of_part ; k++ ){ *(m_system -> pos_old[k]) = *(m_system -> pos[k]) ; }
  std::normal_distribution<double> RNG(0.0,1.0) ;
  std::uniform_real_distribution<double> RNG_positive(0,1) ;
  for ( int i=m_start ; i<m_num_of_part ; i++ )
  {
    m_df = m_system -> get_drift_force ( i ) ;
    for ( int j=0 ; j<m_num_of_dim ; j++ )
    {
      m_rand_temp = RNG ( generat ) ;
      m_df_comp = m_df.get ( j ) ;
      //std::cout << m_rand_temp*m_step_length_sqrt + m_df_comp*m_step_length*m_D << std::endl ;
      m_adjust = m_rand_temp*m_step_length_sqrt + m_df_comp*m_step_length*m_D ;
      //if ( aaa == 0 ) { std::cout << m_adjust << " " << (*m_system -> pos[i]).X() << std::endl ; }
      m_system -> pos[i] -> tweak ( j , m_adjust ) ;
      //if ( aaa == 0 ) { std::cout << (*m_system -> pos[i]).X() << std::endl ; aaa=1 ; }
    }
    m_wf = m_system -> get_wavefunc_sq () ;
    m_df_trial = m_system -> get_drift_force ( i ) ;
    m_greens_ratio = greens ( *m_system -> pos[i] , *m_system -> pos_old[i] ) ;
    //std::cout << m_greens_ratio << std::endl ;
    if ( (m_greens_ratio*(m_wf/m_wf_old)) > RNG_positive ( generat ) )
    {
      m_wf_old = m_wf ;
      m_acceptance += 1 ;
    }
    else
    {
      *(m_system -> pos[i]) = *(m_system -> pos_old[i]) ;
    }
  }
  //if (aaa <10){std::cout<<"aaa"<<std::setprecision(15) <<m_system->pos[1]->X()<<std::endl; }
  //m_energy = m_system -> get_local_energy () ;
}

double Hastings::greens ( Lin pos_temp , Lin pos_old_temp )
{
  double greens_sum = 0 ;
  for ( int i=0 ; i<m_num_of_dim ; i++ )
  {
    greens_sum += 0.5 * ( m_df_trial.get(i) + m_df.get(i) )
    * (
        ( m_D*m_step_length*0.5 ) * ( m_df_trial.get(i) - m_df.get(i) )
        - pos_temp.get(i) + pos_old_temp.get(i)
      ) ;
  }
  return exp(greens_sum) ;
}
