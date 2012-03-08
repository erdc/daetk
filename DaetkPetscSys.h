#ifndef DAETKPETSCSYS_H
#define DAETKPETSCSYS_H

#include "Definitions.h"

namespace Daetk 
{
namespace Petsc 
{
class Err
{
public:
  Err();
  Err(int i);
  int operator=(int i);
};


class Sys
{  
public:
  Sys();
  Sys(int& argc, char **argv,char* help=0, char* file=0);
  virtual ~Sys();
  void barrier();
  bool master();
  bool isInitialized();
  void beginSequential(int ng=1);
  void endSequential(int ng=1);
  bool catchError(bool error);
  int getRank();
  int getSize();
  static bool initialized;
protected:
  Err ierr;
  bool commCreator;
  static int rank,size;
};
}//Petsc
}//Daetk
#endif
