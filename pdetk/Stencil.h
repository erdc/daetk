#ifndef STENCIL_H
#define STENCIL_H

#include <vector>
#include <list>
#include <functional>
#include "Definitions.h"

namespace Daetk 
{
class Stencil
{
public:
  Stencil(int i=0);
  virtual ~Stencil();
  class Location
  {
  public:
    int i,j,k,
      globalNodeNumber,
      point;
    Location(int I, int J, int K,int GN,int P);
    Location();
    Location(const Location& rhs);
  };
  
  class PointListIterator;

  class PointList : public std::vector<Location>
  {
  public:
    PointList();
    void reserve(size_type n);
    void push_back(const Location& v);
    friend class Daetk::Stencil::PointListIterator;
    typedef Daetk::Stencil::PointListIterator iterator;
    inline iterator end();
    inline iterator begin();
  private:
     std::list<int> points;
  };
  
  class PointListIterator
  {
  public:
    inline PointListIterator(const PointListIterator& rhs);
    inline PointListIterator(PointList& pl);
    inline PointListIterator(PointList& pl,Location* l, std::list<int>::iterator i);
    inline ~PointListIterator();
    inline Location& operator*();
    inline Location* operator->();
    inline bool operator!=(const PointListIterator& r);
    inline PointListIterator& operator++();
    inline PointListIterator operator++(int);
  private:
    Location* l_;
    std::list<int>::iterator i_;
    PointList& pl_;
  };

  Location  *anchor;
  PointList *pointList;
  typedef PointList::iterator iterator;
  inline iterator begin();
  inline iterator end();
};
  

  inline Stencil::PointListIterator::PointListIterator(const PointListIterator& rhs):
    l_(rhs.l_),
    i_(rhs.i_),
    pl_(rhs.pl_)
  {}
  
  inline Stencil::PointListIterator::PointListIterator(PointList& pl):
    l_(0),
    i_(pl.points.begin()),
    pl_(pl)
  {
    l_=&pl_[*i_];
  }
  
  inline Stencil::PointListIterator::PointListIterator(PointList& pl,Location* l, std::list<int>::iterator i):
    l_(l),
    i_(i),
    pl_(pl)
  {}
  
  inline Stencil::PointListIterator::~PointListIterator(){}
  
  inline Stencil::Location& Stencil::PointListIterator::operator*(){return *l_;}
  
  inline Stencil::Location* Stencil::PointListIterator::operator->(){return l_;}
  
  inline bool Stencil::PointListIterator::operator!=(const Stencil::PointListIterator& r)
  {
    return i_ != r.i_;
  }
  
  inline Stencil::PointListIterator& Stencil::PointListIterator::operator++()
  {
    ++i_; 
    l_=&pl_[*i_]; 
    return *this;
  }
  
  inline Stencil::PointListIterator Stencil::PointListIterator::operator++(int)
  {
    PointListIterator tmp(pl_);
    tmp.l_ = l_;
    tmp.i_ = i_;
    ++i_; 
    l_=&pl_[*i_];
    return tmp;
  }
  
  inline Stencil::PointList::iterator Stencil::PointList::end()
  {
    return iterator(*this,0,points.end());
  }
  
  inline Stencil::PointList::iterator Stencil::PointList::begin()
  {
    return iterator(*this,&(*this)[*points.begin()],points.begin());
  }

  inline Stencil::iterator Stencil::begin()
  {
    return pointList->begin();
  }
  
  inline Stencil::iterator Stencil::end()
  {
    return pointList->end();
  }

}//Daetk
#endif

