#ifndef GAUSSJASTROW_H
#define GAUSSJASTROW_H

#include "wavefunc.h"

class GaussJastrow : public WaveFunc
{
protected:
  double x_sq ; double y_sq ; double z_sq ;
  double m_a ;
public:
  GaussJastrow ( System * , int , int , double , double ) ;
  double get_wavefunc () ;
  double get_wavefunc_sq () ;
};

#endif
