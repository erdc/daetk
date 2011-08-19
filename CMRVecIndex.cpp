#include "CMRVecIndex.h"

namespace Daetk
{

CMRVecIndex::CMRVecIndex():
  start_(0),
  end_(0),
  all_(true) 
{}

CMRVecIndex::CMRVecIndex( int i1):
  start_(i1),
  end_(i1),
  all_(false)
{}

CMRVecIndex::CMRVecIndex( int i1,  int i2): 
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

int CMRVecIndex::start() const 
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

}//Daetk
