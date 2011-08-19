#ifndef TRACER_H
#define TRACER_H
#include <string>
struct Tracer
{
  std::string msg;
#ifdef DEBUG_TRACE
    Tracer (const char* message)
    {
      msg = message;
      cout << "---> Entering function: " << msg.c_str() << endl << flush;
    }
   #else
   Tracer (const char*) { ; }
   #endif

   #ifdef DEBUG_TRACE
    ~Tracer ()
    {cout << "<--- Leaving function : " << msg << endl << flush; }
   #else
    ~Tracer () { ; }
   #endif
  };

#endif

