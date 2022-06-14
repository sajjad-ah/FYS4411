
#include "read_input.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

Read::Read ( char const * filename )
{
  _num_of_elements = 0 ;
  std::ifstream input_data ( filename ) ;
  std::string line ;
  while ( getline ( input_data , line ) )
  {
    if ( line.substr( 0 , 1 ) == "#" ) { continue ; }
    //removing whitespace
    std::string::iterator end_pos = std::remove ( line.begin() , line.end() , ' ' ) ;
    line.erase( end_pos , line.end() ) ;
    //separate by ":"
    std::stringstream ss ( line ) ;
    while ( getline ( ss , line , ':' ) )
    {
      _input.push_back ( line ) ;
      _num_of_elements += 1 ;
    }
  }
  input_data.close () ;
}

std::string Read::value ( std::string name )
{
  std::string value ;
  for ( int i=0 ; i<_num_of_elements ; i += 2 )
  {
    if ( name == _input[i] )
    {
      value = _input[i+1] ;
      break ;
    }
  }
  return value ;
}

void Read::print ()
{
  for ( int i=0 ; i<_num_of_elements ; i++ ){ std::cout << _input[i] << std::endl ; }
}
