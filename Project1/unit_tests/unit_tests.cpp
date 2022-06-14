using namespace std;

#define CATCH_CONFIG_MAIN
#include <iostream>
#include <cmath>

#include "catch.hpp"

#include "mpi.h"
#include "../main_code/system.h"
#include "../main_code/energy/kinetic.h"
#include "../main_code/energy/harosc.h"
#include "../main_code/energy/haroscjastrow.h"
#include "../main_code/readfile/read_input.h"
#include "../main_code/setstate/initstate.h"
#include "../main_code/setstate/unirand.h"
#include "../main_code/wavefunc/wavefunc.h"
#include "../main_code/wavefunc/gauss.h"
#include "../main_code/wavefunc/gaussjastrow.h"
#include "../main_code/montec/montec.h"
#include "../main_code/montec/metrop.h"
#include "../main_code/montec/hastings.h"
#include "../main_code/lin.h"
double eps = 0.1 ;

TEST_CASE( "" )
{
  double alpha = 0.5 ;
  double beta = 1.0 ;
  int num_of_dim = 3 ;
  int num_of_part = 2 ;
  System * system = new System ;
  system -> pos.resize ( 2 ) ;
  Lin * new_particle1 = new Lin ;
  Lin * new_particle2 = new Lin ;
  system -> pos[0] = new_particle1 ;
  system -> pos[1] = new_particle2 ;
  (*system -> pos[0]) = Lin (1,0,0) ;
  (*system -> pos[1]) = Lin (0,0,0) ;
  system -> set_hamiltonian ( new Harosc ( system , alpha , beta , num_of_dim , num_of_part ) ) ;
  system -> set_wave_func ( new Gauss ( system , num_of_dim , num_of_part , alpha ) ) ;
  double loc_energy ;
  double wavefunc ;
  wavefunc = system -> get_wavefunc_sq () ;
  loc_energy = system -> get_local_energy () ;

  REQUIRE( 3.0 == loc_energy ) ;
  REQUIRE( 0.36787944117144233 == wavefunc ) ;
  //REQUIRE(0.07709975801630 == Approx(vec1[5]).epsilon(eps));
  //REQUIRE(0.28884436449600 == Approx(vec1[7]).epsilon(eps));
  //REQUIRE(0.07709975801630 == Approx(vec2[5]).epsilon(eps));
  //REQUIRE(0.28884436449600 == Approx(vec2[7]).epsilon(eps));
  //REQUIRE(0.07709975801630 == Approx(vec3[5]).epsilon(eps));
  //REQUIRE(0.28884436449600 == Approx(vec3[7]).epsilon(eps));
}
