
#include "harosc.h"
#include <iostream>
#include "../lin.h"

Harosc::Harosc ( System * system , double alpha , double beta , int num_of_dim , int num_of_part ) : Kinetic ( system , alpha , beta , num_of_dim , num_of_part )
{
}

double Harosc::local_energy ()
{
  kinetic_energy () ;
  m_V = 0 ;
  for ( int i=0 ; i<m_num_of_part ; i++ )
  {
    m_V += m_system -> pos[i] -> R_sq_gamma_sq () ;
  }
  return 0.5 * m_V + m_K ;
}

Lin Harosc::drift_force ( int k )
{
  m_temp_pos_beta = (*m_system -> pos[k]) ;
  m_temp_pos_beta.tweak_star ( 2 , m_beta ) ;
  return -4*m_alpha * (m_temp_pos_beta) ;
}
