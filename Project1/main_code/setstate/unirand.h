#ifndef UNIRANd_H
#define UNIRAND_H

#include "initstate.h"
#include "../lin.h"
#include "../system.h"
#include <random>
#include <cstdlib>

class Unirand : public InitState
{
private:
  std::mt19937_64 generat ;
public:
  Unirand ( System * , int , int , double ) ;
  void state_generat () ;
};

#endif
