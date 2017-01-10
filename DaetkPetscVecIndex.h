#ifndef DAETKPETSCVECINDEX_H
#define DAETKPETSCVECINDEX_H

#include "Definitions.h"
///
namespace Daetk 
{
namespace Petsc 
{
class VecIndex

/** This class is used to index unit stride subvectors of Vec's.  */

{
public:
  ///
  inline VecIndex();
  /** Create a VecIndex that accesses the entire vector */

  ///
  inline VecIndex( int i1);
  /** Create a VecIndex that accesses the element i1 */

  ///
  inline VecIndex( int i1,  int i2);
  /** Create a VecIndex that accesses a subvector starting at i1
      and ending at i2 */

  ///
  inline VecIndex(const VecIndex &s);
  /** Create a copy of s */

  ///
  inline int start() const;
  /** Return the first index of subvectors accessed by this VecIndex */
  
  ///
  inline int end() const;
  /** Return the last index of subvectors accessed by this VecIndex */

  ///
  inline int length() const;
  /** Return the length of a subvector accessed by this VecIndex */

  ///
  inline bool all() const;
  /** Return true if the VecIndex accesses an entire vector */

  ///
  inline VecIndex& operator=(const VecIndex& VI);
  /** Copy VI */

  ///
  inline VecIndex operator+(int i);
  /** Return a VecIndex with start=this->start_+i and
      end=this->end_+i */

  ///
  inline VecIndex& operator+=(int i);
  /** Add i to start and end */

  ///
  inline VecIndex operator-(int i);
  /** Return a VecIndex with start= this->start_-i and
      end=this->end_-i */
  
  ///
  inline VecIndex& operator-=(int i);
  /** Subtract i from start and end */

private:
  int start_;
  int end_;       
  bool all_;
};

VecIndex::VecIndex():
  start_(0),
  end_(0),
  all_(true) 
{}

VecIndex::VecIndex( int i1):
  start_(i1),
  end_(i1),
  all_(false)
{}

VecIndex::VecIndex( int i1,  int i2): 
  start_(i1), 
  end_(i2),
  all_(false)
{
  assert(i1 <= i2);
}

VecIndex::VecIndex(const VecIndex &s):
  start_(s.start_),
  end_(s.end_), 
  all_(s.all_)
{}

inline int VecIndex::start() const 
{ 
  return start_;
}

int VecIndex::end() const 
{ 
  return end_;
}

int VecIndex::length() const 
{ 
  return (all_==true) ? 0 : (end_-start_+1);
}

bool VecIndex::all() const 
{ 
  return all_; 
}

VecIndex& VecIndex::operator=(const VecIndex& VI)
{ 
  start_=VI.start_; 
  end_ = VI.end_; 
  return *this;
}

VecIndex VecIndex::operator+(int i)
{ 
  return VecIndex(start_ +i, end_ +i); 
}

VecIndex& VecIndex::operator+=(int i)
{ 
  start_ += i; 
  end_ += i; 
  return *this; 
}

VecIndex VecIndex::operator-(int i)
{ 
  return VecIndex(start_ -i, end_ -i); 
}

VecIndex& VecIndex::operator-=(int i)
{ 
  start_ -= i; 
  end_ -= i; 
  return *this; 
}
}//Petsc
}//Daetk
#endif  
