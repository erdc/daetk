// chrono.cpp - Chronograph class member
//              function.
// Greg Messer

#include "Chronograph.h"

namespace Daetk 
{
#define ULONG unsigned long

double Chronograph::diff(clock_t start, clock_t end)
{
  return double((end - start)) / CLOCKS_PER_SEC;
}
  
Chronograph::Chronograph():
  _isStopped(true),
  _start(0),
  _lapstart(0),
  _elapsed(0.0),
  _lap(0.0)
{}

double Chronograph::start(void)
{
  clock_t ct=clock();
  
  if (!_isStopped)
    {
      double split = diff(_start,ct);
      _elapsed+=split;
      _start=ct;
      _lapstart=ct;
    }
  else
    {
      _isStopped=false;
      _start=ct;
      _lapstart=ct;
    }
  
  return _elapsed;
}

double Chronograph::stop(void)
{
  clock_t ct=clock();
  
  if (!_isStopped)
    {
      _elapsed+=diff(_start,ct);
      _lap=diff(_lapstart,ct);
      _isStopped=true;
    }
  
  return _elapsed;
}

double Chronograph::reset(void)
{
  clock_t ct=clock();
  double oldElapsed=_elapsed;
  if (!_isStopped)
    oldElapsed = _elapsed + diff(_start,ct);
  _isStopped=true;
  _start=0;
  _lapstart=0;
  _elapsed=0.0;
  _lap=0.0;
  return oldElapsed;
}

double Chronograph::lap(void)
{
  //if running calculates lap time and resets lap timer, otherwise returns last lap time
  clock_t ct=clock();
  
  if (!_isStopped)
    {
      _lap = diff(_lapstart,ct);
      _lapstart=ct;
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
  clock_t ct=clock();     
  
  if (_isStopped)
    return _elapsed;
  else
    {
      return _elapsed + diff(_start,ct);
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

