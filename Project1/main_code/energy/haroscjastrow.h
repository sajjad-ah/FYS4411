#ifndef HAROSCJASTROW_H
#define HAROSCJASTROW_H
#include "kinetic.h"

class HaroscJastrow : public Kinetic
{
private:
  class System * _system = nullptr ;
  double m_V ;
  double m_K2 ;
  double m_omega ;
  double m_omega_sq ;
  double m_a ;

public:
  HaroscJastrow ( System * , double , double , double , int , int ) ;
  double local_energy () ;
  Lin drift_force ( int ) ;
};

#endif
