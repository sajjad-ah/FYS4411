#ifndef READ_INPUT_H
#define READ_INPUT_H

#include <vector>
#include <string>

class Read
{
private:
  std::vector < std::string > _input ;
  int _num_of_elements ;
public:
  Read ( char const * ) ;
  std::string value ( std::string ) ;
  void print () ;
  double element ( char * ) ;
};


#endif
