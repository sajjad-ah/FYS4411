#ifndef HASTINGS_H
#define HASTINGS_H

#include "../system.h"
#include "montec.h"

class Hastings : public Montec
{
protected:
  int aaa;
  double m_step_length_sqrt ;
  double m_rand_temp ;
  double m_D ;
  double m_df_comp ;
  double m_adjust ;
  //double m_wf_old ;
  double m_greens_ratio ;
  //std::vector < class Lin * > pos ;
  //std::vector < class Lin * > pos ;
  Lin m_df ;
  Lin m_df_trial ;
public:
  Hastings ( System * , int , int , int , double , double ) ;
  void test ( ) ;
  double greens ( Lin , Lin ) ;

};


#endif
