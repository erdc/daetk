#ifndef DAETKPETSCVECCONSTITERATOR_H
#define DAETKPETSCVECCONSTITERATOR_H

#include "Definitions.h"
#include "DaetkPetscVec.h"  

///
namespace Daetk 
{
namespace Petsc 
{
class VecConstIterator
{
public:
  ///
  inline VecConstIterator(const Vec& V);
  ///
  inline VecConstIterator(const VecConstIterator& I);
  ///
  inline const real& operator*();
  ///
  inline VecConstIterator& operator=(const real* rhs);
  ///
  inline VecConstIterator& operator++();
  ///
  inline VecConstIterator operator++(int);
  ///
  inline VecConstIterator& operator--();
  ///
  inline VecConstIterator operator--(int);
  ///
  inline VecConstIterator& operator+=(const int& n);
  ///
  inline VecConstIterator& operator-=(const int& n);
  ///
  inline VecConstIterator operator+(const int& n);
  ///
  inline VecConstIterator operator-(const int& n);
  ///
  inline const real& operator[](const int& n);
  ///
  inline bool operator==(const real* rhs);
  ///
  inline bool operator!=(const real* rhs);
  ///
  inline bool operator<(const real* rhs);
private:
  const real* p_;
  const unsigned int& stride_;
};



inline VecConstIterator::VecConstIterator(const Vec& V):
  p_(V.begin()),
  stride_(V.stride_)
{}


inline VecConstIterator::VecConstIterator(const VecConstIterator& I):
  stride_(I.stride_),
  p_(I.p_)
{}


inline const real& VecConstIterator::operator*()
{
  return *p_;
}


inline VecConstIterator& VecConstIterator::operator=(const real* rhs)
{
  p_=rhs;
  return *this;
}


inline VecConstIterator& VecConstIterator::operator++()
{
  p_+=stride_;
  return *this;
}


inline VecConstIterator VecConstIterator::operator++(int)
{
  VecConstIterator temporary(*this);
  p_+=stride_;
  return temporary;
}


inline VecConstIterator& VecConstIterator::operator--()
{
  p_-=stride_;
  return *this;
}


inline VecConstIterator VecConstIterator::operator--(int)
{
  VecConstIterator temporary(*this);
  p_-=stride_;
  return temporary;
}
  

inline VecConstIterator& VecConstIterator::operator+=(const int& n)
{
  p_+=n*stride_;
  return *this;
}


inline VecConstIterator& VecConstIterator::operator-=(const int& n)
{
  p_-=n*stride_;
  return *this;
}


inline VecConstIterator VecConstIterator::operator+(const int& n)
{
  VecConstIterator temporary(*this);
  temporary+=n;
  return temporary;
}


inline VecConstIterator VecConstIterator::operator-(const int& n)
{
  VecConstIterator temporary(*this);
  temporary-=n;
  return temporary;
}


inline const real& VecConstIterator::operator[](const int& n)
{
  return *(p_+n*stride_);
}


inline bool VecConstIterator::operator==(const real* rhs)
{
  return p_ == rhs;
}
 

inline bool VecConstIterator::operator!=(const real* rhs)
{
  return p_ != rhs;
}
  

inline bool VecConstIterator::operator<(const real* rhs)
{
  return p_ < rhs;
}
}//Petsc
}//Daetk
#endif
