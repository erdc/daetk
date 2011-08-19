#ifndef DAETKPETSCVECITERATOR_H
#define DAETKPETSCVECITERATOR_H

#include "Definitions.h"
#include "DaetkPetscVec.h"
namespace Daetk 
{
namespace Petsc 
{
///
class VecIterator
{
public:
  ///
  inline VecIterator(Vec& V);
  ///
  inline VecIterator(const VecIterator& I);
  ///
  inline real& operator*();
  ///
  inline VecIterator& operator=(real* rhs);
  ///
  inline VecIterator& operator++();
  ///
  inline VecIterator operator++(int);
  ///
  inline VecIterator& operator--();
  ///
  inline VecIterator operator--(int);
  ///
  inline VecIterator& operator+=(const int& n);
  ///
  inline VecIterator& operator-=(const int& n);
  ///
  inline VecIterator operator+(const int& n);
  ///
  inline VecIterator operator-(const int& n);
  ///
  inline real& operator[](const int& n);
  ///
  inline bool operator==(const real* rhs);
  ///
  inline bool operator!=(const real* rhs);
  ///
  inline bool operator<(const real* rhs);
private:
  real* p_;
  const unsigned int& stride_;
};

 
inline VecIterator::VecIterator(Vec& V):
  p_(V.begin()),
  stride_(V.stride_)
{}


 
inline VecIterator::VecIterator(const VecIterator& I):
  p_(I.p_),
  stride_(I.stride_)
{}

 
inline real& VecIterator::operator*()
{
  return *p_;
}

 
inline VecIterator& VecIterator::operator=(real* rhs)
{
  p_=rhs;
  return *this;
}

 
inline VecIterator& VecIterator::operator++()
{
  p_+=stride_;
  return *this;
}

 
inline VecIterator VecIterator::operator++(int)
{
  VecIterator temporary(*this);
  p_+=stride_;
  return temporary;
}

  
 
inline VecIterator& VecIterator::operator--()
{
  p_-=stride_;
  return *this;
}

 
inline VecIterator VecIterator::operator--(int)
{
  VecIterator temporary(*this);
  p_-=stride_;
  return temporary;
}

 
inline VecIterator& VecIterator::operator+=(const int& n)
{
  p_+=n*stride_;
  return *this;
}

 
inline VecIterator& VecIterator::operator-=(const int& n)
{
  p_-=n*stride_;
  return *this;
}

 
inline VecIterator VecIterator::operator+(const int& n)
{
  VecIterator temporary(*this);
  return temporary+=n;
}
 
inline VecIterator VecIterator::operator-(const int& n)
{
  VecIterator temporary(*this);
  return temporary-=n;
}

 
inline real& VecIterator::operator[](const int& n)
{
  return *(p_+n*stride_);
}

 
inline bool VecIterator::operator==(const real* rhs)
{
  return p_ == rhs;
}


 
inline bool VecIterator::operator!=(const real* rhs)
{
  return p_ != rhs;
}

 
inline bool VecIterator::operator<(const real* rhs)
{
  return p_ < rhs;
}

}//Petsc
}//Daetk
#endif
