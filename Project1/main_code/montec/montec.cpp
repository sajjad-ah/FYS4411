
#include "montec.h"
#include <random>
#include <iostream>

Montec::Montec ( System * system , int my_rank , int num_of_dim , int num_of_part , double step_length )
{
  m_system = system ;
  m_my_rank = my_rank ;
  m_num_of_dim = num_of_dim ;
  m_num_of_part = num_of_part ;
  m_step_length = step_length ;
  m_start = 0 ;
  m_acceptance = 0 ;
}

void Montec::set_seed ()
{
  generat.seed ( clock() + 3.14*m_my_rank ) ;
}

void Montec::set_wf_old ( double wf_old )
{
  m_wf_old = wf_old ;
}

void Montec::set_energy ( double energy )
{
  m_energy = energy ;
}

double Montec::get_energy ( )
{
  return m_energy ;
}

//used for one-body density: keeps the positions
//of particle number < start fixed
void Montec::set_start ( int start )
{
  m_start = start ;
}

double Montec::get_acceptance ()
{
  return m_acceptance ;
}
