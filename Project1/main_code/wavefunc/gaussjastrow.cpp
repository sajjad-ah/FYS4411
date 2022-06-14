
#include "gaussjastrow.h"
#include "../lin.h"
#include <cmath>
#include <iostream>


GaussJastrow::GaussJastrow ( System * system , int num_of_dim , int num_of_part , double alpha , double a ) : WaveFunc (system)
{
  m_num_of_param = 1 ;
  m_num_of_dim = num_of_dim ;
  m_num_of_part = num_of_part ;
  m_alpha = alpha ;
  m_a = a ;
  //param.resize ( m_num_of_param ) ;
  //param[0] = m_system -> m_alpha ;
}

double GaussJastrow::get_wavefunc ()
{
  double value = 0 ;
  Lin r_temp ;
  double r_temp_abs ;
  for ( int i=0 ; i<m_num_of_part ; i++)
  {
    value -= m_alpha * ( m_system -> pos[i] -> R_sq_gamma () )  ;
  }
  value = exp( value ) ;
  for ( int i=0 ; i<m_num_of_part-1 ; i++ )
  {
    for ( int j=i+1 ; j<m_num_of_part ; j++ )
    {
      r_temp = ( *m_system -> pos[i] - *m_system -> pos[j] ) ;
      r_temp_abs = r_temp.abs() ;
      //std::cout << r_temp_abs  << std::endl ;
      if ( r_temp_abs >= m_a )
      {
        value *= ( 1 - m_a / r_temp_abs ) ;
      }
      else
      {
        return 0 ;
      }
    }
  }
  return value ;
}

double GaussJastrow::get_wavefunc_sq ()
{
  double value = get_wavefunc () ;
  //std::cout << value*value << std::endl ;
  return value*value ;
}
