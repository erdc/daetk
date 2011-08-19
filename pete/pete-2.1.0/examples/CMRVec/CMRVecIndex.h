#ifndef CMRVECINDEX_H
#define CMRVECINDEX_H

#include <assert.h>

///
class CMRVecIndex

/** This class is used to index unit stride subvectors of CMRVec's.  */

{
public:
  ///
  inline CMRVecIndex();
  /** Create a CMRVecIndex that accesses the entire vector */

  ///
  inline CMRVecIndex(unsigned int i1);
  /** Create a CMRVecIndex that accesses the element i1 */

  ///
  inline CMRVecIndex(unsigned int i1, unsigned int i2);
  /** Create a CMRVecIndex that accesses a subvector starting at i1
      and ending at i2 */

  ///
  inline CMRVecIndex(const CMRVecIndex &s);
  /** Create a copy of s */

  ///
  inline int start() const;
  /** Return the first index of subvectors accessed by this CMRVecIndex */
  
  ///
  inline int end() const;
  /** Return the last index of subvectors accessed by this CMRVecIndex */

  ///
  inline int length() const;
  /** Return the length of a subvector accessed by this CMRVecIndex */

  ///
  inline bool all() const;
  /** Return true if the CMRVecIndex accesses an entire vector */

  ///
  inline CMRVecIndex& operator=(const CMRVecIndex& I);
  /** Copy I */

  ///
  inline CMRVecIndex operator+(int i);
  /** Return a CMRVecIndex with start=this->start_+i and
      end=this->end_+i */

  ///
  inline CMRVecIndex& operator+=(int i);
  /** Add i to start and end */

  ///
  inline CMRVecIndex operator-(int i);
  /** Return a CMRVecIndex with start= this->start_-i and
      end=this->end_-i */
  
  ///
  inline CMRVecIndex& operator-=(int i);
  /** Subtract i from start and end */

private:
  unsigned int start_;
  unsigned int end_;       
  bool all_;
};

CMRVecIndex::CMRVecIndex():
  start_(0),
  end_(0),
  all_(true) 
{}

CMRVecIndex::CMRVecIndex(unsigned int i1):
  start_(i1),
  end_(i1),
  all_(false)
{}

CMRVecIndex::CMRVecIndex(unsigned int i1, unsigned int i2): 
  start_(i1), 
  end_(i2),
  all_(false)
{
  assert(i1 <= i2);
}

CMRVecIndex::CMRVecIndex(const CMRVecIndex &s):
  start_(s.start_),
  end_(s.end_), 
  all_(s.all_)
{}

inline int CMRVecIndex::start() const 
{ 
  return (all_==true) ? 0 : start_;
}

int CMRVecIndex::end() const 
{ 
  return (all_ ==true) ? 0 : end_;
}

int CMRVecIndex::length() const 
{ 
  return (all_==true) ? 0 : (end_-start_+1);
}

bool CMRVecIndex::all() const 
{ 
  return all_; 
}

CMRVecIndex& CMRVecIndex::operator=(const CMRVecIndex& I)
{ 
  start_=I.start_; 
  end_ = I.end_; 
  return *this;
}

CMRVecIndex CMRVecIndex::operator+(int i)
{ 
  return CMRVecIndex(start_ +i, end_ +i); 
}

CMRVecIndex& CMRVecIndex::operator+=(int i)
{ 
  start_ += i; 
  end_ += i; 
  return *this; 
}

CMRVecIndex CMRVecIndex::operator-(int i)
{ 
  return CMRVecIndex(start_ -i, end_ -i); 
}

CMRVecIndex& CMRVecIndex::operator-=(int i)
{ 
  start_ -= i; 
  end_ -= i; 
  return *this; 
}

#endif  





