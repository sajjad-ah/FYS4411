#ifndef SAMPLER_H
#define SAMPLER_H

class Sampler
{
private:
  int m_num_of_cyc ;
  double m_energy ;
  double m_energy_sq ;
  double m_der_wf ;
  double m_der_wf_energy ;
  double m_wf_sq ;
  int m_acceptance ;

  double m_energy_mean ;
  double m_energy_sq_mean ;
  double m_der_wf_mean ;
  double m_der_wf_energy_mean ;
  double m_wf_sq_mean ;


public:
  Sampler ( int ) ;
  void sample ( double ) ;
  void sample_der ( double , double ) ;
  void sample_wf ( double ) ;
  void sample_acceptance ( int ) ;
  void average () ;
  void average_wf () ;
  void set_average  ( double , double , double , double ) ;
  void set_sum      ( double , double , double , double ) ;
  void set_average_wf  ( double ) ;
  void set_sum_wf      ( double ) ;
  double acceptance_rate ( int ) ;
  double get_average_energy () ;
  double get_average_energy_sq () ;
  double get_average_der_wf () ;
  double get_average_der_wf_energy () ;
  double get_average_wf () ;
  double get_wf () ;
  void print () ;
  void print_energies ( double ) ;
};


#endif
