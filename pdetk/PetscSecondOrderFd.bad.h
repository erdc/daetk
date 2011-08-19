#ifndef PETSCSECONDORDERFD_H
#define PETSCSECONDORDERFD_H

#include <vector>
#include "Definitions.h"
#include "Stencil.h"
#include "PetscSys.h"

//notice: global interface numbers do not work as all interface arrays are local
// a local array with ghost points is large enough to hold the interface arrays so I simply
// use a local array with local node numbering done as for the old stencil 

namespace Daetk 
{
namespace Petsc
{
  namespace cc
    {
  struct _p_DA;
  struct _p_AO;
    }
class SecondOrderFd : public Stencil
{
public:
//    enum Point {BOTTOM,FRONT,LEFT,CENTER,RIGHT,BACK,TOP};
  enum Point {CENTER,LEFT,RIGHT,FRONT,BACK,BOTTOM,TOP};
  bool isLocal;

  Petsc::cc::_p_DA* da_,*dadof_;
  Petsc::cc::_p_AO* ao_;
  int *ltog;
  int nxNodes,
    nyNodes,
    nzNodes,
    nxyNodes,
    nxyzNodes,
    local_x0,
    local_y0,
    local_z0,
    globalLow,
    globalHigh,
    local_nxNodes,
    local_nyNodes,
    local_nzNodes,
    local_nxyNodes,
    local_nxyzNodes,
    ghost_x0,
    ghost_y0,
    ghost_z0,
    ghost_nxNodes,
    ghost_nyNodes,
    ghost_nzNodes,
    ghost_nxyNodes,
    ghost_nxyzNodes,
    ghost_xOffSet,
    ghost_yOffSet,
    ghost_zOffSet,
    ghost_globalLow,
    ghost_globalHigh,
    center,
    center_noGhost,
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
    interBack,
    npx,
    npy,
    npz;
  
  inline int nNodes();
  
  SecondOrderFd(int nxNodesIn=1, int nyNodesIn=1,int nzNodesIn=1);
  
  virtual ~SecondOrderFd();
  
  void setDomainSize(int nxNodesIn=1, int nyNodesIn=1,int nzNodesIn=1);
  
  inline SecondOrderFd& operator()(int k);
  
  inline SecondOrderFd& operator()(int j,int k);
  
  inline SecondOrderFd& operator()(int i, int j, int k);
  
  inline int globalToLocal(int nodeNumber);

  inline SecondOrderFd& globalNode(int nodeNumber);
  
  inline void setGlobalPointList(int i, int j, int k);

  inline void setGlobalPointList(int nodeNumber);

  int petscToGlobal(int nodeNumber);

  inline void setInterfaces(int i, int j, int k);

  inline void setPoints();
  
  inline SecondOrderFd& setLocalPointList();

  inline SecondOrderFd& localIndex(int k);
  
  inline SecondOrderFd& localIndex(int j,int k);
  
  inline SecondOrderFd& localIndex(int i, int j, int k);
  
//    inline SecondOrderFd& localNode(int nodeNumber);
  
//    inline int localToGlobal(int localNodeNumber) 
protected:

//    typedef std::vector< std::map<int, Location> > MapField; 
//    MapField pointMapField, localPointMapField;
//    typedef std::vector< std::list<Location> > ListField; 
//    ListField pointListField, localPointListField;
  bool ownsVectorField;
  typedef std::vector< std::vector<Location> > VectorField; 
  VectorField pointVectorField;
};

class SecondOrderFd1d : public SecondOrderFd
{
public:
  SecondOrderFd1d(int nxNodesIn, int nyNodes=1, int nzNodes =1, int dof=1);
  SecondOrderFd1d(SecondOrderFd1d& fineStencil, int nxNodesIn, int nyNodes=1, int nzNodes =1, int dof=1);
  inline int size() { return 3;}
};

class SecondOrderFd2d : public SecondOrderFd
{
public:
  SecondOrderFd2d(int nxNodesIn, int nyNodesIn, int nzNodes=1, int dof=1);
  SecondOrderFd2d(SecondOrderFd2d& fineStencil, int nxNodesIn, int nyNodesIn, int nzNodes=1, int dof=1);
  inline int size() { return 5;}
};

class SecondOrderFd3d : public SecondOrderFd
{
public:
  SecondOrderFd3d(int nxNodesIn, int nyNodesIn, int nzNodesIn, int dof=1);
  SecondOrderFd3d(SecondOrderFd3d& fineStencil, int nxNodesIn, int nyNodesIn, int nzNodesIn, int dof=1);
  int size() { return 7;}
};

  inline int SecondOrderFd::nNodes() { return nxyzNodes;}

  inline SecondOrderFd& SecondOrderFd::operator()(int k)
  {
    int i=0,j=0;
    if (k>= local_x0 && k < (local_x0 + local_nxNodes))
      {
        setInterfaces(i,j,k);
        isLocal=true;
      }
    else
      isLocal=false;
//      pointMap = &pointMapField[k];
//      anchor = &pointMapField[k][CENTER];
//      pointList = &pointListField[k];
//      anchor = &pointListField[k].back();
    pointVector = &pointVectorField[k];
    anchor = &pointVectorField[k][CENTER];
    
    setPoints();
    return *this;
  }
  
  inline SecondOrderFd& SecondOrderFd::operator()(int j,int k)
  {
    int i=0;
    if (k>= local_x0 && k < (local_x0 + local_nxNodes) &&
        j>= local_y0 && j < (local_y0 + local_nyNodes))
      {
        setInterfaces(i,j,k);
        isLocal=true;
      }
    else
      isLocal=false;
//      pointMap = &pointMapField[j*nxNodes+k];
//      anchor = &pointMapField[j*nxNodes+k][CENTER];
//      pointList = &pointListField[j*nxNodes+k];
//      anchor = &pointListField[j*nxNodes+k].back();
    pointVector = &pointVectorField[j*nxNodes+k];
    anchor = &pointVectorField[j*nxNodes+k][CENTER];
    
    setPoints();
    return *this;
  }
  
  inline SecondOrderFd& SecondOrderFd::operator()(int i, int j, int k)
  {
    if ((k>= local_x0) && (k < (local_x0 + local_nxNodes)) &&
        (j>= local_y0) && (j < (local_y0 + local_nyNodes)) &&
        (i>= local_z0) && (i < (local_z0 + local_nzNodes)))
      {
        setInterfaces(i,j,k);
        isLocal=true;
      }
    else
      isLocal=false;
//      pointMap = &pointMapField[i*nxyNodes + j*nxNodes+k];
//      anchor =&pointMapField[i*nxyNodes + j*nxNodes+k][CENTER];
//      pointList = &pointListField[i*nxyNodes + j*nxNodes+k];
//      anchor =&pointListField[i*nxyNodes + j*nxNodes+k].back();
    pointVector = &pointVectorField[i*nxyNodes + j*nxNodes+k];
    anchor =&pointVectorField[i*nxyNodes + j*nxNodes+k][CENTER];
    
    setPoints();
    return *this;
  }
  
  inline int SecondOrderFd::globalToLocal(int nodeNumber)
  {
//      Location& n=pointListField[nodeNumber].back();
    Location& n=pointVectorField[nodeNumber][CENTER];
    return (n.i-ghost_z0)*ghost_nxyNodes 
      + (n.j - ghost_y0)*ghost_nxNodes
      + (n.k - ghost_x0);
  }

  inline SecondOrderFd& SecondOrderFd::globalNode(int nodeNumber)
  {
//      pointMap = &pointMapField[nodeNumber];
//      anchor =&pointMapField[nodeNumber][CENTER];
//      pointList = &pointListField[nodeNumber];
//      anchor =&pointListField[nodeNumber].back();
    pointVector = &pointVectorField[nodeNumber];
    anchor =&pointVectorField[nodeNumber][CENTER];
    int i = anchor->i;
    int j = anchor->j;
    int k = anchor->k;
    
    if ((k>= local_x0) && (k < (local_x0 + local_nxNodes)) &&
        (j>= local_y0) && (j < (local_y0 + local_nyNodes)) &&
        (i>= local_z0) && (i < (local_z0 + local_nzNodes)))
      {
        setInterfaces(i,j,k);
        isLocal=true;
      }
    else
      isLocal=false;
    setPoints();
    return (*this);
  }
  
  inline void SecondOrderFd::setGlobalPointList(int i, int j, int k)
  {
//      pointList = &pointListField[i*nxyNodes + j*nxNodes+k];
//      anchor =&pointListField[i*nxyNodes + j*nxNodes+k].back();
    pointVector = &pointVectorField[i*nxyNodes + j*nxNodes+k];
    anchor =&pointVectorField[i*nxyNodes + j*nxNodes+k][CENTER];
  }

  inline void SecondOrderFd::setGlobalPointList(int nodeNumber)
  {
//      pointList = &pointListField[nodeNumber];
//      anchor =&pointListField[nodeNumber].back();
    pointVector = &pointVectorField[nodeNumber];
    anchor =&pointVectorField[nodeNumber][CENTER];
  }

  inline void SecondOrderFd::setInterfaces(int i, int j, int k)
  {
    i-= local_z0;
    j-= local_y0;
    k-= local_x0;
    interTop    = (i+1)*local_nxyNodes + j*local_nxNodes + k;
    interBottom = i*local_nxyNodes + j*local_nxNodes + k;
    interBack   = i*(local_nxyNodes+local_nxNodes) + (j+1)*local_nxNodes+k;
    interFront  = i*(local_nxyNodes+local_nxNodes) + j*local_nxNodes+k;
    interRight  = i*(local_nxyNodes+local_nyNodes) + j*(local_nxNodes+1)+k+1;
    interLeft   = i*(local_nxyNodes+local_nyNodes) + j*(local_nxNodes+1)+k;
  }

  inline void SecondOrderFd::setPoints()
  {
//      this won't work
//      center = (*pointVector)[CENTER].globalNodeNumber;
//      left = (*pointVector)[LEFT].globalNodeNumber;
//      right = (*pointVector)[RIGHT].globalNodeNumber;
//      front = (*pointVector)[FRONT].globalNodeNumber;
//      back = (*pointVector)[BACK].globalNodeNumber;
//      bottom = (*pointVector)[BOTTOM].globalNodeNumber;
//      top = (*pointVector)[TOP].globalNodeNumber;
    
//      PointMap::iterator sit=pointMap->begin(),end=pointMap->end();
//      PointList::iterator sit=pointList->begin(),end=pointList->end();
    PointVector::iterator sit=pointVector->begin(),end=pointVector->end();
    while (sit!=end)
      {
//          switch (sit->first)
        switch (sit->point)
          {
          case CENTER:
//              center      = sit->second.globalNodeNumber;
            center      = sit->globalNodeNumber;
            break;
          case RIGHT:
//              right       = sit->second.globalNodeNumber;
            right       = sit->globalNodeNumber;
            break;
          case LEFT:
//              left        = sit->second.globalNodeNumber;
            left        = sit->globalNodeNumber;
            break;
          case BACK:
//              back        = sit->second.globalNodeNumber;
            back        = sit->globalNodeNumber;
            break;
          case FRONT:
//              front       = sit->second.globalNodeNumber;
            front       = sit->globalNodeNumber;
            break;
          case TOP:
//              top         = sit->second.globalNodeNumber;
            top         = sit->globalNodeNumber;
            break;
          case BOTTOM:
//              bottom      = sit->second.globalNodeNumber;
            bottom      = sit->globalNodeNumber;
            break;
          default:
            std::cerr<<"enum value out of range in SecondOrderFD::setPoints()"<<std::endl; 
            break;
          }
        ++sit;
      }
  }
  
//    inline SecondOrderFd& SecondOrderFd::setLocalPointList()
//    {
//      pointList = &localPointListField[center];
//      anchor =&localPointListField[center].back();
//      return *this;
//    }

  inline SecondOrderFd& SecondOrderFd::localIndex(int k)
  {
    isLocal = true;

    interRight = k+1;
    interLeft  = k;

    center_noGhost = k;

    k+=ghost_xOffSet;
    
    center     = k;
    left       = k-1;
    right      = k+1;
    
    return *this;
  }
  
  inline SecondOrderFd& SecondOrderFd::localIndex(int j,int k)
  {
    isLocal = true;

    interBack   = (j+1)*local_nxNodes+k;
    interFront  = j*local_nxNodes+k;
    interRight  = j*(local_nxNodes+1)+k+1;
    interLeft   = j*(local_nxNodes+1)+k;

    center_noGhost = interFront;

    j+=ghost_yOffSet; 
    k+=ghost_xOffSet;
    
    center      = j*ghost_nxNodes+k;
    back        = (j+1)*ghost_nxNodes+k;
    front       = (j-1)*ghost_nxNodes+k;
    right       = j*ghost_nxNodes + k+1;
    left        = j*ghost_nxNodes + k-1;
    
    return *this;
  }
  
  inline SecondOrderFd& SecondOrderFd::localIndex(int i, int j, int k)
  {
//      std::cout<<i<<'\t'<<j<<'\t'<<k<<std::endl;
//      std::cout<<local_nxNodes<<'\t'<<local_nyNodes<<'\t'<<local_nzNodes<<std::endl;
//      std::cout<<local_nxyNodes<<std::endl;
    isLocal = true;
    interTop    = (i+1)*local_nxyNodes + j*local_nxNodes + k;
    interBottom = i*local_nxyNodes + j*local_nxNodes + k;
    interBack   = i*(local_nxyNodes+local_nxNodes) + (j+1)*local_nxNodes+k;
    interFront  = i*(local_nxyNodes+local_nxNodes) + j*local_nxNodes+k;
    interRight  = i*(local_nxyNodes+local_nyNodes) + j*(local_nxNodes+1)+k+1;
    interLeft   = i*(local_nxyNodes+local_nyNodes) + j*(local_nxNodes+1)+k;

    center_noGhost = interBottom;

    i+=ghost_zOffSet;
    j+=ghost_yOffSet;
    k+=ghost_xOffSet;
    
    center      = i*ghost_nxyNodes + j*ghost_nxNodes+k;
    top         = (i+1)*ghost_nxyNodes + j*ghost_nxNodes + k;
    bottom      = (i-1)*ghost_nxyNodes + j*ghost_nxNodes + k;
    back        = i*ghost_nxyNodes + (j+1)*ghost_nxNodes+k;
    front       = i*ghost_nxyNodes + (j-1)*ghost_nxNodes+k;
    right       = i*ghost_nxyNodes + j*ghost_nxNodes + k+1;
    left        = i*ghost_nxyNodes + j*ghost_nxNodes + k-1;
    
    return *this;
  }
  
//    inline SecondOrderFd& SecondOrderFd::localNode(int nodeNumber)
//    {
//      isLocal = true;
//  //      pointMap = &localPointMapField[nodeNumber];
//  //      anchor =&localPointMapField[nodeNumber][CENTER];
//      pointVector = &localPointVectorField[nodeNumber];
//      anchor =&localPointVectorField[nodeNumber].back();
//      center      = anchor->globalNodeNumber;
//      top         = center + ghost_nxyNodes;
//      bottom      = center - ghost_nxyNodes;
//      back        = center + ghost_nxNodes;
//      front       = center - ghost_nxNodes;
//      right       = center+1;
//      left        = center-1;
    
//  //      interTop    = center;
//  //      interBottom = bottom;
//  //      interBack   = center;
//  //      interFront  = front;
//  //      interRight  = center;
//  //      interLeft   = left;
    
//      return (*this);
//    }
  
//    inline int SecondOrderFd::localToGlobal(int localNodeNumber) 
//    {
//      return globalLow + localNodeNumber;
//    }

}//Petsc
}//Daetk
#endif
