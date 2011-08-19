// chrono.h - A Chronograph class
//            for program timing.
// Greg Messer

/**
 This class provides a Chronograph for an
 application.  It works like a stopwatch.
 It provides timing for laps (time since
 last lap) and splits (same as elapsed time).
 It can be stopped, restarted and reset to
 zero.

 This class uses the ANSI C clock() function
 to get the time.  The clock() function
 returns a clock_t that is the number of
 clock ticks since a program started.
 The lap, split and elapsed times are
 calculated as the difference between two
 clock_t values.  If the program is a
 long-running one, the clock_t may overflow.
 So, don't use this for timing longer than
 the value of MAX_LONG divided by the ANSI C
 macro CLOCKS_PER_SECOND.

 The compiler-supplied copy constructor and
 assignment operator are OK for this class,
 so they are not implemented here.

 Use the istream function precision(2) to
 format the returned double values from the
 member functions.

 Use "%.2lf" as the printf() format string
 to format the returned double values from
 the member functions.

 */

#ifndef CHRONO_H
#define CHRONO_H

#include "Definitions.h"
#include <ctime>

namespace Daetk 
{

class Chronograph
{
private:
  bool _isStopped;

//   clock_t _start,
//     _lapstart;
  std::clock_t _start_ct, _lapstart_ct;
  std::time_t _start_tt, _lapstart_tt;
  double _elapsed,_lap;
  
  double diff(std::clock_t start_ct, std::clock_t end_ct,std::time_t start_tt, std::time_t end_tt);
  //double diff(timeval* start, timeval* end);

public:

  Chronograph();
   virtual ~Chronograph(){}
   
  //updates and returns elapsed time, resets lap timer
  double start(void);

  //updates and returns elapsed time, sets lap time
  double stop(void);
  
  double reset(void);
  
  double lap(void);

  double split(void);
  
  int isStopped(void);
  
  double elapsed(void);
  
  double elapsedHMS(double &Hours,
                    double &Mins,
                    double &Secs);
};
}//Daetk
#endif


