#ifndef HAROSC_H
#define HAROSC_H

#include "kinetic.h"

class Harosc : public Kinetic
{
private:
  double m_V ;
  Lin m_temp_pos_beta ;
public:
  Harosc ( System * , double , double , int , int ) ;
  double local_energy () ;
  Lin drift_force  ( int ) ;
};

#endif
