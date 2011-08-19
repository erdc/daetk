#ifndef CMRVECOPERATORS_H
#define CMRVECOPERATORS_H

#include "CMRVec.h"
#include "CMRVecIterator.h"
#include "CMRVecConstIterator.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <iomanip.h>


template<class T>
CMRVec<T>& operator*=(CMRVec<T> &x, const T &a)
{
  T* ix(x.begin());
  while (ix < x.end())
    {
      (*ix)*=a;
      ix+=stride_;
    }
  return x;
}

template<class T>
CMRVec<T>& operator/=(CMRVec<T> &x, const T &a)
{
  T* ix(x.begin());
  while (ix < x.end())
    {
      (*ix)/=a;
      ix+=stride_;
    }
  return x;
}

CMRVec<T>& operator+=(CMRVec<T> &x, const CMRVec<T> &y)
{
  assert(x.size() == y.size());
  const T* iy(y.begin());
  T* ix(x.begin());
  while (ix <x.end() )
    {
      (*ix)+=(*iy);
      ix+=stride_;
      iy+=stride_;
    }
  return x;
}
          
      
template<class T>
CMRVec<T>& operator-=(CMRVec<T> &x, const CMRVec<T> &y)
{
  assert(x.size() == y.size());
  /*
  if (x.size() != y.size())
    {
      cout << "Incompatible vector lengths in -." << endl;
      exit(1);
    }
    */
  const T* iy(y.begin());
  T* ix(x.begin());
  while (ix < x.end())
    {
      (*ix)-=(*iy);
      ix+=stride_;
      iy+=stride_;
    }
  return x;
}
          
      
template<class T>
T dot(const CMRVec<T> &x, const CMRVec<T> &y)
{
        
  //  Check for compatible dimensions:
  assert(x.size() == y.size());
  /*
  if (x.size() != y.size())
      {
         cout << "Incompatible dimensions in dot(). " << endl;
         exit(1);
      }
      */
  double temp =  0.0;
  const T* ix(x.begin()),iy(y.begin());
  while (ix < x.end())
    {
      temp+=(*ix)*(*iy);
      ix+=stride_; 
      iy+=stride_;
    }
  return temp;
}

template<class T>
T norm(const CMRVec<T> &x)
{
      T temp = dot(x,x);
      return sqrt(temp);
}

template<class T>
ostream&  operator<<(ostream& s, const CMRVec<T>& V)
{
  s.setf(ios::scientific);
  s.precision(17);
  for (int i=0; i< V.dim() ; i++)
    s<<V(i)<<endl;
  return s;
}

template<class T>
istream&  operator>>(istream& s, CMRVec<T>& V)
{
  int length;
  s>>length;
  V.newsize(length);
  for (int i=0; i<length ; i++)
    s>>V(i);
  return s;
}

#endif








