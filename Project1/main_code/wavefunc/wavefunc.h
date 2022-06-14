#ifndef WAVEFUNC_H
#define WAVEFUNC_H

#include "../system.h"
#include <vector>

class WaveFunc
{
protected:
  System * m_system   ;
  int m_num_of_param  ;
  int m_num_of_dim    ;
  int m_num_of_part   ;
  //double * m_param ;
  //std::vector <double> param ;
  double m_alpha ;
  double m_beta ;
  double m_temp_r2 ;
public:
  WaveFunc ( System * ) ;
  double get_der_alpha () ;
  void set_alpha ( double ) ;
  virtual double get_wavefunc () = 0 ;
  virtual double get_wavefunc_sq () = 0 ;
};

#endif
