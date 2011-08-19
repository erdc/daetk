#ifndef SECONDORDERFD_H
#define SECONDORDERFD_H

#include "Definitions.h"
#include "Stencil.h"

namespace Daetk 
{
class SecondOrderFd : public Stencil
{
public:
  int nxNodes,
    nyNodes,
    nzNodes,
    nxyNodes,
    nxyzNodes,
    center,
    left,
    right,
    bottom,
    top,
    front,
    back,
    interLeft,
    interRight,
    interBottom,
    interTop,
    interFront,
    interBack;
  
  inline int nNodes() { return nxyzNodes;}
  
  SecondOrderFd(int nxNodesIn=1, int nyNodesIn=1,int nzNodesIn=1);
 
  virtual ~SecondOrderFd();
  
  void setDomainSize(int nxNodesIn=1, int nyNodesIn=1,int nzNodesIn=1);
  
  inline SecondOrderFd& operator()(int k)
    {
      center     = k;
      left       = k-1;
      right      = k+1;
      interRight = k+1;
      interLeft  = k;
      pointList = &pointListField[k];
      anchor = &pointListField[k].back();
      return *this;
    }

  inline SecondOrderFd& operator()(int j,int k)
    {
      center      = j*nxNodes+k;
      back        = (j+1)*nxNodes+k;
      front       = (j-1)*nxNodes+k;
      right       = j*nxNodes + k+1;
      left        = j*nxNodes + k-1;
      interBack   = (j+1)*nxNodes+k;
      interFront  = j*nxNodes+k;
      interRight  = j*(nxNodes+1)+k+1;
      interLeft   = j*(nxNodes+1)+k;
      pointList = &pointListField[j*nxNodes+k];
      anchor = &pointListField[j*nxNodes+k].back();
      return *this;
    }
  
  inline SecondOrderFd& operator()(int i, int j, int k)
    {
      center      = i*nxyNodes + j*nxNodes+k;
      top         = (i+1)*nxyNodes + j*nxNodes + k;
      bottom      = (i-1)*nxyNodes + j*nxNodes + k;
      back        = i*nxyNodes + (j+1)*nxNodes+k;
      front       = i*nxyNodes + (j-1)*nxNodes+k;
      right       = i*nxyNodes + j*nxNodes + k+1;
      left        = i*nxyNodes + j*nxNodes + k-1;
      interTop    = (i+1)*nxyNodes + j*nxNodes + k;
      interBottom = i*nxyNodes + j*nxNodes + k;
      interBack   = i*(nxyNodes+nxNodes) + (j+1)*nxNodes+k;
      interFront  = i*(nxyNodes+nxNodes) + j*nxNodes+k;
      interRight  = i*(nxyNodes+nyNodes) + j*(nxNodes+1)+k+1;
      interLeft   = i*(nxyNodes+nyNodes) + j*(nxNodes+1)+k;
      pointList = &pointListField[i*nxyNodes + j*nxNodes+k];
      anchor =&pointListField[i*nxyNodes + j*nxNodes+k].back();
      return *this;
    }
  
  inline SecondOrderFd& globalNode(int nodeNumber)
    {
      pointList = &pointListField[nodeNumber];
      anchor =&pointListField[nodeNumber].back();
      //cout<<anchor->globalNodeNumber<<'\t'<<nodeNumber<<endl;
      assert(anchor->globalNodeNumber == nodeNumber);
      center      = nodeNumber;
      top         = center + nxyNodes;
      bottom      = center - nxyNodes;
      back        = center + nxNodes;
      front       = center - nxNodes;
      right       = center+1;
      left        = center-1;

      interTop    = (anchor->i+1)*nxyNodes + anchor->j*nxNodes + anchor->k;
      interBottom = anchor->i*nxyNodes + anchor->j*nxNodes + anchor->k;
      interBack   = anchor->i*(nxyNodes+nxNodes) + (anchor->j+1)*nxNodes+anchor->k;
      interFront  = anchor->i*(nxyNodes+nxNodes) + anchor->j*nxNodes+anchor->k;
      interRight  = anchor->i*(nxyNodes+nyNodes) + anchor->j*(nxNodes+1)+anchor->k+1;
      interLeft   = anchor->i*(nxyNodes+nyNodes) + anchor->j*(nxNodes+1)+anchor->k;

      return (*this);
    }
protected:
  typedef std::vector< std::list<Location> > ListField; 
  ListField pointListField;
};

class SecondOrderFd1d : public SecondOrderFd
{
public:
  enum Point {LEFT,CENTER,RIGHT};
  
  SecondOrderFd1d(int nxNodesIn, int nyNodes=1, int nzNodes =1);
  inline int size() { return 3;}
};

class SecondOrderFd2d : public SecondOrderFd
{
public:
  enum Point {FRONT,LEFT,CENTER,RIGHT,BACK};

  SecondOrderFd2d(int nxNodesIn, int nyNodesIn, int nzNodes=1);
  inline int size() { return 5;}
};

class SecondOrderFd3d : public SecondOrderFd
{
public:
  enum Point {BOTTOM,FRONT,LEFT,CENTER,RIGHT,BACK,TOP};
  SecondOrderFd3d(int nxNodesIn, int nyNodesIn, int nzNodesIn);
  
  inline int size() { return 7;}
};
}//Daetk
#endif
