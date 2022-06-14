#include <cmath>
#include <iomanip>
#include <iostream>
#include "lin.h"

Lin::Lin ()
{
  x = 0; y = 0; z = 0 ;
  pos[0] = & x ; pos[1] = & y ; pos[2] = & z ;
}

Lin::Lin ( double x0 , double y0 , double z0 )
{
  x = x0; y = y0; z = z0 ;
  pos[0] = & x ; pos[1] = & y ; pos[2] = & z ;
}

Lin::Lin ( double x0 , double y0 , double z0 , double gamma )
{
  x = x0; y = y0; z = z0 ;
  pos[0] = & x ; pos[1] = & y ; pos[2] = & z ;
  m_gamma = gamma ;
  m_gamma_sq = gamma*gamma ;
}

void Lin::set ( int i , double xi )
{
  * pos[i] = xi ;
}

void Lin::set_gamma ( double gamma )
{
  m_gamma = gamma ;
  m_gamma_sq = gamma*gamma ;
}

double Lin::get ( int i ) const
{
  return * pos[i] ;
}

void Lin::tweak ( int i , double xi )
{
  * pos[i] += xi ;
}

void Lin::tweak_star ( int i , double xi )
{
  * pos[i] *= xi ;
}

double Lin::X () const
{
  return x ;
}
double Lin::Y () const
{
  return y ;
}
double Lin::Z () const
{
  return z ;
}
double Lin::GAMMA () const
{
  return m_gamma ;
}
double Lin::GAMMA_SQ () const
{
  return m_gamma_sq ;
}

double Lin::R_sq () const
{
  return x*x + y*y + z*z ;
}

double Lin::R_sq_gamma () const
{
  return x*x + y*y + m_gamma*z*z ;
}

double Lin::R_sq_gamma_sq () const
{
  return x*x + y*y + m_gamma_sq*z*z ;
}

double Lin::abs () const
{
  return sqrt( x*x + y*y + z*z ) ;
}

Lin Lin::unit ()
{
  norm = sqrt ( x*x + y*y + z*z ) ;
  return Lin ( x/norm , y/norm, z/norm ) ;
}

Lin& Lin::operator= ( const Lin & c )
{
   x = c.X() ;
   y = c.Y() ;
   z = c.Z() ;
   m_gamma = c.GAMMA() ;
   m_gamma_sq = c.GAMMA_SQ() ;
   return *this ;
}

void Lin::read() const
{
  std::cout << std::setprecision(9) << x << " " << y << " " << z << std::endl ;
}
