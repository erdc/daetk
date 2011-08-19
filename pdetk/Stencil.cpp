#include "Stencil.h"

namespace Daetk
{
  Stencil::Stencil(int i){}

  Stencil::~Stencil(){}
    
  Stencil::Location::Location(int I, int J, int K,int GN,int P):
    i(I),j(J),k(K),
    globalNodeNumber(GN),
    point(P)
  {}

  Stencil::Location::Location():i(0),j(0),k(0),globalNodeNumber(0),point(0){}

  Stencil::PointList::PointList(){}
  
  Stencil::Location::Location(const Stencil::Location& rhs):
    i(rhs.i),
    j(rhs.j),
    k(rhs.k),
    globalNodeNumber(rhs.globalNodeNumber),
    point(rhs.point)
  {}

  void Stencil::PointList::reserve(size_type n)
  {
    std::vector<Location>::resize(n,Location(-1,-1,-1,-1,-1));
  }
  
  void Stencil::PointList::push_back(const Location& v)
  {
    //mwf added v arg to resize
    if (v.point >=  signed(size()))
      std::vector<Location>::resize(v.point+1,v);
    else
      (*this)[v.point] = v;
    points.push_back(v.point);
    points.sort(std::less<int>());
  }
  
  
}//Daetk

