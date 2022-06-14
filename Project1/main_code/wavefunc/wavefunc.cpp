
#include "wavefunc.h"
#include <iostream>


WaveFunc::WaveFunc ( System * system )
{
  m_system = system ;
}

double WaveFunc::get_der_alpha ()
{
  m_temp_r2 = 0 ;
  Lin a ;
  for ( int i=0 ; i<m_num_of_part ; i++ )
  {
    a = *(m_system -> pos[i]) ;
    m_temp_r2 += ( a . R_sq_gamma () ) ;
  }
  return -m_temp_r2 ;
}

void WaveFunc::set_alpha ( double alpha ) { m_alpha = alpha ; }
