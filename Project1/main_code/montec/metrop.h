#ifndef METROP_H
#define METROP_H

#include "../system.h"
#include "montec.h"

class Metrop : public Montec
{
protected:
public:
  Metrop ( System * , int , int , int , double ) ;
  void test ( ) ;

};


#endif
