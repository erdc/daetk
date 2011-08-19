#ifndef UTILITIES_H
#define UTILITIES_H

#include "Definitions.h"
#include "Chronograph.h"
#include "Tracer.h"

namespace Daetk 
{

#ifdef CRAYCC
inline real copysign(const real& a,const real& b)
{ return (b<0) ? -a:a;}
inline real sqr(const real& a){return a*a;}
#endif
#ifdef __GNUC__
inline real sqr(const real& a){return a*a;}
#endif

inline const real& max2(const real& a,const real& b)
{
  return (b<a) ? a:b;
}
inline const real& min2(const real& a,const real& b)
{
  return (b>a) ? a:b;
}

inline const real& max3(const real& a,const real& b,const real& c)
{
  if (b<a)
    return (c<a) ? a:c;
  else
    return (c<b) ? b:c;
}

inline const real& min3(const real& a,const real& b,const real& c)
{
  if (b>a)
    return (c>a) ? a:c;
  else
    return (c>b) ? b:c;
}

extern "C"
{
  double ccsqrt(const double& x);
  double ccfabs(const double& x);
  double ccpow(const double& x, const double& id);
  double f77powdd(const double& x, const double& id);
  double  f77powdi(const double& x, const int& i);
  int  f77powii(const int& x,const int& i);
}

inline double mypow(const double& x,const int& k)
{
  double temp=1.0/x;
  for(int i=1;i<k;i++)
    temp*=temp;
  return temp;
}
      
}//Daetk
#endif
