
#include "sampler.h"
#include <iomanip>
#include <iostream>

Sampler::Sampler ( int num_of_cyc ){
  m_num_of_cyc = num_of_cyc ;
  m_energy = 0 ;
  m_energy_sq = 0 ;
  m_der_wf = 0 ;
  m_der_wf_energy = 0 ;
  m_wf_sq = 0 ;
  m_acceptance = 0 ;
}

void Sampler::sample ( double energy )
{
  //std::cout << energy << std::endl ;
  m_energy += energy ;
  m_energy_sq += energy*energy ;
}

//used for parameter opt
void Sampler::sample_der ( double energy , double der_wf )
{
  m_energy += energy ;
  m_energy_sq += energy*energy ;
  m_der_wf += der_wf ;
  m_der_wf_energy += der_wf * energy ;
}

void Sampler::sample_wf ( double wf_sq )
{
  m_wf_sq += wf_sq ;
}

void Sampler::sample_acceptance ( int acceptance )
{
  m_acceptance += acceptance ;
}

void Sampler::average ()
{
  m_energy_mean = m_energy / m_num_of_cyc ;
  m_energy_sq_mean = m_energy_sq / m_num_of_cyc ;
  m_der_wf_mean = m_der_wf / m_num_of_cyc ;
  m_der_wf_energy_mean = m_der_wf_energy / m_num_of_cyc ;
  m_wf_sq_mean = m_wf_sq / m_num_of_cyc ;
}

double Sampler::acceptance_rate ( int num_of_part )
{
  double double_acceptance = (double) m_acceptance ;
  return double_acceptance / (m_num_of_cyc*num_of_part) ;
}

void Sampler::set_average ( double energy_mean , double energy_sq_mean , double der_wf_mean , double der_wf_energy_mean )
{
  m_energy_mean = energy_mean ;
  m_energy_sq_mean = energy_sq_mean ;
  m_der_wf_mean = der_wf_mean ;
  m_der_wf_energy_mean = der_wf_energy_mean ;
}

void Sampler::set_average_wf ( double wf_sq_mean )
{
  m_wf_sq_mean = wf_sq_mean ;
}

void Sampler::set_sum ( double energy , double energy_sq , double der_wf , double der_wf_energy )
{
  m_energy = energy ;
  m_energy_sq = energy_sq ;
  m_der_wf = der_wf ;
  m_der_wf_energy = der_wf_energy ;
}

void Sampler::set_sum_wf ( double wf_sq )
{
  m_wf_sq = wf_sq ;
}

double Sampler::get_average_energy ()
{
  return m_energy_mean ;
}

double Sampler::get_average_energy_sq ()
{
  return m_energy_sq_mean ;
}

double Sampler::get_average_der_wf ()
{
  return m_der_wf_mean ;
}

double Sampler::get_average_der_wf_energy ()
{
  return m_der_wf_energy_mean ;
}

double Sampler::get_average_wf ()
{
  return m_wf_sq_mean ;
}

double Sampler::get_wf ()
{
  return m_wf_sq ;
}

void Sampler::print ()
{
  std::cout << "Mean energy: " << std::setprecision(9) << m_energy_mean << std::endl ;
  std::cout << "Mean squared energy: " << std::setprecision(9) << m_energy_sq_mean << std::endl ;
  std::cout << "Standard deviation: " << std::setprecision(9) << m_energy_sq_mean - m_energy_mean*m_energy_mean << std::endl ;
}

void Sampler::print_energies ( double energy )
{
  std::cout << energy << std::endl ;
}
