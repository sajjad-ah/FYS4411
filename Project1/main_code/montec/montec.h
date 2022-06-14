#ifndef MONTEC_H
#define MONTEC_H

#include "../system.h"
#include <random>

class Montec
{
protected:
  std::mt19937_64 generat ;
  int m_my_rank ;
  int m_num_of_dim ;
  int m_num_of_part ;
  int m_start ;
  int m_acceptance ;
  double m_wf ;
  double m_wf_old ;
  double m_wf_rel ;
  double m_energy ;
  double m_step_length ;
  System * m_system ;

public:
  Montec ( System * , int , int , int , double ) ;
  void set_seed () ;
  void set_wf_old ( double ) ;
  void set_energy ( double ) ;
  void set_start ( int ) ;
  double get_energy () ;
  double get_acceptance () ;
  virtual void test ( ) = 0 ;

};


#endif
