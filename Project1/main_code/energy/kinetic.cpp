
#include "kinetic.h"
#include "../lin.h"
#include <iostream>

Kinetic::Kinetic ( System * system , double alpha , double beta , int num_of_dim , int num_of_part )
{
  m_system = system ;
  m_alpha = alpha ;
  m_alpha_sq = alpha*alpha ;
  m_beta = beta ;
  m_num_of_dim = num_of_dim ;
  m_num_of_part = num_of_part ;
}

void Kinetic::kinetic_energy (  )
{
  m_K = 0 ;
  for ( int i=0 ; i<m_num_of_part ; i++ )
  {
    m_K +=  m_system -> pos[i] -> R_sq_gamma_sq () ;
  }
  //std::cout << m_alpha_sq << std::endl ;
  //std::cout << m_K << std::endl ;
  m_K *= -2*m_alpha_sq ;
  m_K += m_alpha*(2+m_beta)* m_num_of_part ;
}

void Kinetic::set_alpha ( double alpha ) { m_alpha = alpha ; m_alpha_sq = alpha*alpha ; }
