#ifndef InitState_H
#define InitState_H

#include "../lin.h"
#include "../system.h"

class InitState
{
protected:
  System * m_system ;
  double m_beta ;
public:
  InitState ( System * , double ) ;
  void create_state () ;
  void add_pos ( double * ) ;
  int m_num_of_dim ;
  int m_num_of_part ;
  virtual void state_generat () = 0;
};

#endif
