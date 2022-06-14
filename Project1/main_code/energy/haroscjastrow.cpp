
#include "haroscjastrow.h"
#include "../lin.h"
#include <iostream>

HaroscJastrow::HaroscJastrow ( System * system , double alpha , double beta , double a , int num_of_dim , int num_of_part ) : Kinetic ( system , alpha , beta , num_of_dim , num_of_part )
{
  m_a = a ;
  kinetic_energy () ;
}

double HaroscJastrow::local_energy ()
{
  //1st term
  kinetic_energy () ;
  m_V = 0 ;
  m_K2 = 0 ;

  Lin r_temp ;
  Lin r_kj ;
  Lin r_ki ;
  Lin r_beta ;
  double r_kj_abs ;
  double r_ki_abs ;
  double r_kj_sq ;
  double r_ki_sq ;
  double r_kj_ki ;
  double u_kj ;

  for ( int k=0 ; k<m_num_of_part ; k++ )
  {
    r_temp = Lin ( 0 , 0 , 0 ) ;
    for ( int j=0 ; j<m_num_of_part ; j++ )
    {
      if ( j != k )
      {
        r_kj = ( *m_system -> pos[k] - *m_system -> pos[j] ) ;
        r_kj_abs = r_kj.abs() ;
        r_kj_sq = r_kj_abs*r_kj_abs ;
        u_kj = ( m_a / ( r_kj_abs*(r_kj_abs - m_a) ) ) ;
        //2nd term
        r_temp = ( u_kj / r_kj_abs ) * r_kj + r_temp ;
        //4th term
        m_K2 += ( ( 1 - 2*(r_kj_abs/m_a) )*u_kj + 2/r_kj_abs )*u_kj ;
        /*  3rd term (double sum)  */
        for ( int i=0 ; i<m_num_of_part ; i++ )
        {
          if ( k != i )
          {
            r_ki = ( *m_system -> pos[k] - *m_system -> pos[i] ) ;
            r_kj_ki = ( r_kj , r_ki ) ;
            r_ki_abs = r_ki.abs() ;
            r_ki_sq = r_ki_abs*r_ki_abs ;
            m_K2 += ( r_kj_ki / ( r_kj_abs*r_ki_abs ) ) * ( m_a / ( r_kj_sq - m_a*r_kj_abs ) ) * ( m_a / ( r_ki_sq - m_a*r_ki_abs ) ) ;
          }
        }
      }
    }
    r_beta = (*m_system -> pos[k]) ;
    r_beta.tweak_star( 2 , m_beta ) ;
    // "," is overloaded as a dot product
    m_K2 += -2*m_alpha* ( r_beta , r_temp ) ;
  }
  m_K2 *= ( -0.5 ) ;

  //computing potential energy contribution
  for ( int i=0 ; i<m_num_of_part ; i++ )
  {
    m_V += m_system -> pos[i] -> R_sq_gamma_sq() ;
  }
  m_V = 0.5 * m_V ;
  return m_K + m_K2 + m_V ;
}

Lin HaroscJastrow::drift_force ( int k )
{
  Lin df ( 0 , 0 , 0 ) ;
  Lin r_temp ;
  double r_abs_temp ;
  for ( int i=0 ; i<k ; i++ )
  {
    r_temp = *m_system -> pos[i] - *m_system -> pos[k] ;
    r_abs_temp = r_temp.abs() ;
    df = df +  ( 1 / ( r_abs_temp*(r_abs_temp - m_a) ) ) * r_temp.unit() ;
  }

  for ( int i=k+1 ; i<m_num_of_part ; i++ )
  {
    r_temp = *m_system -> pos[k] - *m_system -> pos[i] ;
    r_abs_temp = r_temp.abs() ;
    df = df - ( 1 / ( r_abs_temp*(r_abs_temp - m_a) ) ) * r_temp.unit() ;
  }
  df = m_a*df ;
  df = df - (2*m_alpha)* (*m_system -> pos[k]) ;
  return df ;
}


//3rd term
/*for ( int j=0 ; j<m_system -> num_of_particles ; j++ )
{
  for ( int i=0 ; i<m_system -> num_of_particles ; i++ )
  {
    if ( k != i and k != j )
    {
      r_kj = ( *m_system -> pos[k] - *m_system -> pos[j] ) ;
      r_ki = ( *m_system -> pos[k] - *m_system -> pos[i] ) ;
      r_kj_ki = ( r_kj , r_ki ) ;
      r_kj_abs = r_kj.abs() ;
      r_ki_abs = r_ki.abs() ;
      r_kj_sq = r_kj_abs*r_kj_abs ;
      r_ki_sq = r_ki_abs*r_ki_abs ;
      m_K2 += ( r_kj_ki / ( r_kj_abs*r_ki_abs ) ) * ( m_a / ( r_kj_sq - m_a*r_kj_abs ) ) * ( m_a / ( r_ki_sq - m_a*r_ki_abs ) ) ;
    }
  }
}*/
