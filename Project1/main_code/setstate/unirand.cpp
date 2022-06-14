
#include "unirand.h"
#include "../lin.h"
#include "../system.h"
#include <random>
#include <iostream>
#include <ctime>

Unirand::Unirand ( System * system , int num_of_dim , int num_of_part , double beta ) : InitState ( system , beta )
{
  generat.seed ( clock() ) ;
  InitState::m_num_of_dim = num_of_dim ;
  InitState::m_num_of_part = num_of_part ;
  create_state () ;
}

void Unirand::state_generat ()
{
  std::uniform_real_distribution<double> RNG(-0.1,0.1) ;
  int _arraylength = m_num_of_dim * m_num_of_part ;
  double * _temp_randarray = new double [ _arraylength ] ;
  for ( int i=0 ; i<_arraylength ; i++ )
  {
    _temp_randarray [i] = RNG ( generat ) ;
  }
  InitState::add_pos ( _temp_randarray ) ;
  delete [] _temp_randarray ;
}
