// chrono.cpp - Chronograph class member
//              function.
// Greg Messer

#include "Chronograph.h"

namespace Daetk 
{
  using std::clock_t;
  using std::time_t;
  using std::clock;
  using std::time;
  using std::difftime;

#define ULONG unsigned long

// // double Chronograph::diff(clock_t start, clock_t end)
// double Chronograph::diff(timeval* start,timeval* end)
// {
// //   return double((end - start)) / CLOCKS_PER_SEC;
//   return double(end->tv_sec - start->tv_sec) + double(end->tv_usec - start->tv_usec)/1.0e6); 
// }
  
double Chronograph::diff(std::clock_t start_ct, std::clock_t end_ct,std::time_t start_tt, std::time_t end_tt)
{
  double seconds,frac_of_second,short_time,long_time;
  short_time = double((end_ct - start_ct)) / CLOCKS_PER_SEC;
  frac_of_second = std::modf(short_time,&seconds);
  long_time = difftime(end_tt,start_tt);
  if (frac_of_second > 0.0)
    return long_time + frac_of_second;
  else
    return long_time;
}

Chronograph::Chronograph():
  _isStopped(true),
  _elapsed(0.0),
  _lap(0.0)
{}

double Chronograph::start(void)
{
//   timezone tz;
//   timeval ct;
//   gettimeofday(&ct,&tz);
  clock_t ct=clock();
  time_t tt=time(0);

  if (!_isStopped)
    {
      double split = diff(_start_ct,ct,_start_tt,tt);
      _elapsed+=split;
      _start_ct=ct;
      _lapstart_ct=ct;
      _start_tt=tt;
      _lapstart_tt=tt;
    }
  else
    {
      _isStopped=false;
      _start_ct=ct;
      _lapstart_ct=ct;
      _start_tt=tt;
      _lapstart_tt=tt;
     }
  
  return _elapsed;
}

double Chronograph::stop(void)
{
//   timezone tz;
//   timeval ct;
//   gettimeofday(&ct,&tz);
  clock_t ct=clock();
  time_t tt=time(0);
  if (!_isStopped)
    {
      _elapsed+=diff(_start_ct,ct,_start_tt,tt);
      _lap=diff(_lapstart_ct,ct,_lapstart_tt,tt);
      _isStopped=true;
    }
  
  return _elapsed;
}

double Chronograph::reset(void)
{
//   timezone tz;
//   timeval ct;
//   gettimeofday(&ct,&tz);
  clock_t ct=clock();
  time_t tt=time(0);
  double oldElapsed=_elapsed;
  if (!_isStopped)
    oldElapsed = _elapsed + diff(_start_ct,ct,_start_tt,tt);
  _isStopped=true;
//   timerclear(&_start);
//   timerclear(&_lapstart);
  _start_ct=0;_start_tt=0;
  _lapstart_ct=0;_lapstart_tt=0;
  _elapsed=0.0;
  _lap=0.0;
  return oldElapsed;
}

double Chronograph::lap(void)
{
  //if running calculates lap time and resets lap timer, otherwise returns last lap time
//   timezone tz;
//   timeval ct;
//   gettimeofday(&ct,&tz);
 clock_t ct=clock();
 time_t tt=time(0);
  if (!_isStopped)
    {
      _lap = diff(_lapstart_ct,ct,_lapstart_tt,tt);
      _lapstart_ct=ct;
      _lapstart_tt=tt;
    }
  
  return _lap;
}

double Chronograph::split(void)
{
  return elapsed();
}

int Chronograph::isStopped(void)
{
  return _isStopped;
}
  
double Chronograph::elapsed(void)
{
 //  timezone tz;
//   timeval ct;
//   gettimeofday(&ct,&tz);
  clock_t ct=clock();     
  time_t tt=time(0);
  if (_isStopped)
    return _elapsed;
  else
    {
      return _elapsed + diff(_start_ct,ct,_start_tt,tt);
    }
}


double Chronograph::elapsedHMS(double &Hours,
                               double &Mins,
                               double &Secs)
{
  double dHour   = 60.0 * 60.0;
  double dMinute = 60.0;
    double dElapsed;
    double dWork;
    double dHours;
    double dMins;
    double dSecs;

    dElapsed = elapsed();

    if(dElapsed < dMinute)
    {
        dHours = 0.0;
        dMins = 0.0;
        dSecs = dElapsed;
    }
    else if(dElapsed < dHour)
    {
        dWork = dElapsed;

        dHours = 0.0;

        dMins = (ULONG) (dElapsed / dMinute);
        dWork -= dMins * dMinute;

        dSecs = dWork;
    }
    else
    {
        dWork = dElapsed;

        dHours = (ULONG) (dElapsed / dHour);
        dWork -= dHours * dHour;

        dMins = (ULONG) (dWork / dMinute);
        dWork -= dMins * dMinute;

        dSecs = dWork;
    }

    Hours = dHours;
    Mins  = dMins;
    Secs  = dSecs;

    return dElapsed;
}

}//Daetk

