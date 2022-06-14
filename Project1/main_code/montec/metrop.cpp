
#include "metrop.h"
#include "../lin.h"
#include <iostream>

Metrop::Metrop ( System * system , int my_rank , int num_of_dim , int num_of_part , double step_length ) : Montec ( system , my_rank , num_of_dim , num_of_part , step_length )
{
}

void Metrop::test ( )
{
  m_acceptance = 0 ;
  for ( int k=0 ; k<m_num_of_part ; k++ ){ *(m_system -> pos_old[k]) = *(m_system -> pos[k]) ; }
  std::uniform_real_distribution<double> RNG_positive(0,1) ;
  std::uniform_real_distribution<double> RNG(-0.5,0.5) ;
  for ( int i=m_start ; i<m_num_of_part ; i++ )
  {
    for ( int j=0 ; j<m_num_of_dim ; j++ )
    {
      m_system -> pos[i] -> tweak ( j , m_step_length * RNG ( generat ) ) ;
    }
    m_wf = m_system -> get_wavefunc_sq () ;
    m_wf_rel = m_wf / m_wf_old ;
    if ( m_wf_rel >= RNG_positive ( generat ) )
    {
      m_wf_old = m_wf ;
      m_acceptance += 1 ;
    }
    else
    {
      *(m_system -> pos[i]) = *(m_system -> pos_old[i]) ;
    }
  }
}
