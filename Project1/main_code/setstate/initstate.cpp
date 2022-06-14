
#include "initstate.h"
#include "../lin.h"
#include "../system.h"
#include <iostream>

InitState::InitState ( System * system , double beta )
{
  m_system = system ;
  m_beta = beta ;
}

void InitState::create_state ()
{
  m_system -> pos.resize ( m_num_of_part ) ;
  m_system -> pos_old.resize ( m_num_of_part ) ;
}

void InitState::add_pos ( double * _temp_randarray )
{
  int counter = 0 ;
  for ( int i=0 ; i<m_num_of_part ; i++ )
  {
    Lin * new_particle = new Lin ;
    Lin * new_particle_old = new Lin ;
    m_system -> pos[i] = new_particle ;
    m_system -> pos_old[i] = new_particle_old ;
    m_system -> pos[i] -> set_gamma ( m_beta ) ;
    m_system -> pos_old[i] -> set_gamma ( m_beta ) ;

    for ( int j=0 ; j<m_num_of_dim ; j++ )
    {
      m_system -> pos[i] -> set ( j , _temp_randarray[ counter ] ) ;
      m_system -> pos_old[i] -> set ( j , 0 ) ;
      counter ++ ;
    }
  }
}
