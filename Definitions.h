#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#ifndef CSTDDEF
#define CSTDDEF <cstddef>
#endif

#ifndef CSTRING
#define CSTRING <cstring>
#endif

#ifndef CFLOAT
#define CFLOAT <cfloat>
#endif

#ifndef CMATH
#define CMATH  <cmath>
#endif

#ifndef ALGORITHM
#define ALGORITHM <algorithm>
#endif

#ifndef STRING
#define STRING    <string>
#endif

#ifndef CASSERT
#define CASSERT <cassert>
#endif

#ifndef IOSTREAM
#define IOSTREAM <iostream>
#endif

#ifndef IOMANIP
#define IOMANIP <iomanip>
#endif

#include <bitset>
#include IOSTREAM
#include IOMANIP
#include STRING
#include ALGORITHM

#ifdef __sgi
extern "C" 
{
#endif
#include CSTDDEF
#include CFLOAT
#include CMATH
#include CSTRING
#include CASSERT
#ifdef __sgi
}
#endif
#include "mpi++.h"
//#define DEBUG_TRACE

namespace Daetk 
{
#ifdef __sgi
  using ::size_t;
  using ::fabs;
  using ::sqrt;
  using ::pow;
  using ::log;
  using ::log10;
  using ::exp;
#else
  using std::size_t;
  using std::fabs;
  using std::sqrt;
  using std::pow;
  using std::log;
  using std::log10;
  using std::exp;
#endif
typedef double real;

//#define USE_SINGLE_PRECISION

#ifdef F77_POST_UNDERSCORE
#define F77NAME(name) name ## _
#else
#define F77NAME(name) name
#endif

#ifndef USE_SINGLE_PRECISION
const double MACHINE_EPSILON=DBL_EPSILON;
const double SQRT_MACHINE_EPSILON=sqrt(MACHINE_EPSILON);
#else
const float MACHINE_EPSILON=FLT_EPSILON;
const float SQRT_MACHINE_EPSILON=sqrt(MACHINE_EPSILON);
#endif
}//Daetk
#endif
