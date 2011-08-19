#ifndef TRACER_H
#define TRACER_H

#include "Definitions.h"
#include IOSTREAM
#include CSTRING

namespace Daetk 
{
struct Tracer
{
  char msg[200];
#ifdef DEBUG_TRACE
  inline Tracer (const char* message)
  {
    std::strcpy(msg,message);
    std::cout << "---> Entering function: " << std::msg << std::endl << std::flush;
  }
#else
  inline Tracer (const char*) 
  {}
#endif
  
#ifdef DEBUG_TRACE
  inline ~Tracer ()
  {
    std::cout << "<--- Leaving function : " << std::msg << std::endl << std::flush; 
  }
#else
  inline ~Tracer () 
  {}
#endif
};
}//Daetk
#endif

