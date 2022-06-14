#ifndef KINETIC_H
#define KINETIC_H

#include "../system.h"

class Kinetic
{
protected:
  double m_K ;
  double m_alpha ;
  double m_alpha_sq ;
  double m_beta ;
  int m_num_of_dim ;
  int m_num_of_part ;
  System * m_system ;
public:
  Kinetic ( System * , double , double , int , int ) ;
  void kinetic_energy ( ) ;
  void set_alpha ( double ) ;
  virtual double local_energy () = 0 ;
  virtual Lin drift_force  ( int ) = 0 ;
};

#endif
