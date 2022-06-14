
#include "gauss.h"
#include "../lin.h"
#include <cmath>
#include <iostream>


Gauss::Gauss ( System * system , int num_of_dim , int num_of_part , double alpha ) : WaveFunc (system)
{
  m_num_of_param = 1 ;
  m_num_of_dim = num_of_dim ;
  m_num_of_part = num_of_part ;
  m_alpha = alpha ;
  //param.resize ( m_num_of_param ) ;
  //param[0] = m_system -> m_alpha ;
}

double Gauss::get_wavefunc ()
{
  double value = 0 ;
  for ( int i=0 ; i<m_num_of_part ; i++)
  {
    value -= m_alpha * ( m_system -> pos[i] -> R_sq_gamma () )  ;
  }
  return exp( value ) ;
}

double Gauss::get_wavefunc_sq ()
{
  double value = get_wavefunc () ;
  return value*value ;
}
