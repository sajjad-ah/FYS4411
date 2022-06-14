#ifndef GAUSS_H
#define GAUSS_H

#include "wavefunc.h"

class Gauss : public WaveFunc
{
protected:
  double x_sq ; double y_sq ; double z_sq ;
public:
  Gauss ( System * , int , int , double ) ;
  double get_wavefunc () ;
  double get_wavefunc_sq () ;
};

#endif
