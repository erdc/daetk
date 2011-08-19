#ifndef PETSCSECONDORDERFD_H
#define PETSCSECONDORDERFD_H

#include <vector>
#include "Definitions.h"
#include "Stencil.h"
#include "DaetkPetscSys.h"

// using namespace Daetk;
// using namespace Petsc;
// using Daetk::Petsc::Err;

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
    npz,
    NOT_ONED,
    NOT_TWOD,
    stencilSize;
  
  inline int nNodes();

  SecondOrderFd(int nxNodesIn, int nyNodes=1, int nzNodes =1, int dof=1, int stencilWidth=1);
  SecondOrderFd(SecondOrderFd& fineStencil, int nxNodesIn, int nyNodes=1, int nzNodes =1, int dof=1, int stencilWidth=1);
  virtual ~SecondOrderFd();
  
  void setDomainSize(int nxNodesIn=1, int nyNodesIn=1,int nzNodesIn=1);
  
  inline SecondOrderFd& operator()(int k);
  
  inline SecondOrderFd& operator()(int j,int k);
  
  inline SecondOrderFd& operator()(int i, int j, int k);

  inline void operator++();
  //special fast global index for jacobian--must be used locally paying attention to boundaries 
  inline SecondOrderFd& globalIndex(int i, int j, int k);
  
  inline int globalToLocal(int nodeNumber);

  inline SecondOrderFd& globalNode(int nodeNumber);
  
  inline void setGlobalPointList(int i, int j, int k);

  inline void setGlobalPointList(int nodeNumber);

  int petscToGlobal(int nodeNumber);

  inline void setInterfaces(int i, int j, int k);

  inline void setPoints();
  
  inline SecondOrderFd& localIndex(int k);
  
  inline SecondOrderFd& localIndex(int j,int k);
  
  inline SecondOrderFd& localIndex(int i, int j, int k);
  
  inline SecondOrderFd& localNode(int nodeNumber);
  
//    inline int localToGlobal(int localNodeNumber) 
  inline int size() { return stencilSize;}

protected:

  typedef std::vector< PointList > ListField; 
  ListField pointListField;

  void ctor1d(int nxNodesIn, int nyNodes=1, int nzNodes =1, int dof=1, int stencilWidth=1);
  void ctor1d(SecondOrderFd& fineStencil, int nxNodesIn, int nyNodes=1, int nzNodes =1, int dof=1, int stencilWidth=1);
  void ctor2d(int nxNodesIn, int nyNodesIn, int nzNodes=1, int dof=1, int stencilWidth=1);
  void ctor2d(SecondOrderFd& fineStencil, int nxNodesIn, int nyNodesIn, int nzNodes=1, int dof=1, int stencilWidth=1);
  void ctor3d(int nxNodesIn, int nyNodesIn, int nzNodesIn, int dof=1, int stencilWidth=1);
  void ctor3d(SecondOrderFd& fineStencil, int nxNodesIn, int nyNodesIn, int nzNodesIn, int dof=1, int stencilWidth=1);
};

typedef SecondOrderFd SecondOrderFd1d;
typedef SecondOrderFd SecondOrderFd2d;
typedef SecondOrderFd SecondOrderFd3d;

  inline int SecondOrderFd::nNodes() { return nxyzNodes;}

  inline SecondOrderFd& SecondOrderFd::operator()(int k)
  {
    int i=0,j=0;
    if (k>= local_x0 && k < (local_x0 + local_nxNodes))
      {
        setInterfaces(i-local_z0,j-local_y0,k-local_x0);
        isLocal=true;
      }
    else
      isLocal=false;

    pointList = &pointListField[k];
    anchor = &pointListField[k][CENTER];
    
    setPoints();

    return *this;
  }
  
  inline SecondOrderFd& SecondOrderFd::operator()(int j,int k)
  {
    int i=0;
    if (k>= local_x0 && k < (local_x0 + local_nxNodes) &&
        j>= local_y0 && j < (local_y0 + local_nyNodes))
      {
        setInterfaces(i-local_z0,j-local_y0,k-local_x0);
        isLocal=true;
      }
    else
      isLocal=false;

    pointList = &pointListField[j*nxNodes+k];
    anchor = &pointListField[j*nxNodes+k][CENTER];
    
    setPoints();

    return *this;
  }
  
  inline SecondOrderFd& SecondOrderFd::operator()(int i, int j, int k)
  {
    if ((k>= local_x0) && (k < (local_x0 + local_nxNodes)) &&
        (j>= local_y0) && (j < (local_y0 + local_nyNodes)) &&
        (i>= local_z0) && (i < (local_z0 + local_nzNodes)))
      {
        setInterfaces(i-local_z0,j-local_y0,k-local_x0);
        isLocal=true;
      }
    else
      isLocal=false;

    pointList = &pointListField[i*nxyNodes + j*nxNodes+k];
    anchor =&pointListField[i*nxyNodes + j*nxNodes+k][CENTER];
    
    setPoints();

    return *this;
  }
  

  inline SecondOrderFd& SecondOrderFd::globalIndex(int i, int j, int k)
  {
    setInterfaces(i-local_z0,j-local_y0,k-local_x0);
    pointList = &pointListField[i*nxyNodes + j*nxNodes+k];
    anchor =&pointListField[i*nxyNodes + j*nxNodes+k][CENTER];
    setPoints();
    return *this;
  }
  
  inline int SecondOrderFd::globalToLocal(int nodeNumber)
  {
    Location& n=pointListField[nodeNumber][CENTER];
    return (n.i-ghost_z0)*ghost_nxyNodes 
      + (n.j - ghost_y0)*ghost_nxNodes
      + (n.k - ghost_x0);
  }

  inline SecondOrderFd& SecondOrderFd::globalNode(int nodeNumber)
  {
    pointList = &pointListField[nodeNumber];
    anchor =&pointListField[nodeNumber][CENTER];

    int i = anchor->i;
    int j = anchor->j;
    int k = anchor->k;
    
    if ((k>= local_x0) && (k < (local_x0 + local_nxNodes)) &&
        (j>= local_y0) && (j < (local_y0 + local_nyNodes)) &&
        (i>= local_z0) && (i < (local_z0 + local_nzNodes)))
      {
        setInterfaces(i-local_z0,j-local_y0,k-local_x0);
        isLocal=true;
      }
    else
      isLocal=false;

    setPoints();

    return (*this);
  }
  
  inline void SecondOrderFd::setGlobalPointList(int i, int j, int k)
  {
    pointList = &pointListField[i*nxyNodes + j*nxNodes+k];
    anchor =&pointListField[i*nxyNodes + j*nxNodes+k][CENTER];
  }

  inline void SecondOrderFd::setGlobalPointList(int nodeNumber)
  {
    pointList = &pointListField[nodeNumber];
    anchor =&pointListField[nodeNumber][CENTER];
  }

inline void SecondOrderFd::setInterfaces(int i, int j, int k)
{
  interTop    = NOT_TWOD*((i+1)*local_nxyNodes + j*local_nxNodes + k);
  interBottom = NOT_TWOD*(i*local_nxyNodes + j*local_nxNodes + k);
  interBack   = NOT_ONED*(i*(local_nxyNodes+local_nxNodes) + (j+1)*local_nxNodes+k);
  interFront  = NOT_ONED*(i*(local_nxyNodes+local_nxNodes) + j*local_nxNodes+k);
  interRight  = i*(local_nxyNodes+local_nyNodes) + j*(local_nxNodes+1)+k+1;
  interLeft   = i*(local_nxyNodes+local_nyNodes) + j*(local_nxNodes+1)+k;
}

inline void SecondOrderFd::setPoints()
{
  center = (*pointList)[CENTER].globalNodeNumber;
  left   = (*pointList)[LEFT].globalNodeNumber;
  right  = (*pointList)[RIGHT].globalNodeNumber;
  front  = (*pointList)[FRONT].globalNodeNumber;
  back   = (*pointList)[BACK].globalNodeNumber;
  bottom = (*pointList)[BOTTOM].globalNodeNumber;
  top    = (*pointList)[TOP].globalNodeNumber;
}
  
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

  inline SecondOrderFd& SecondOrderFd::localNode(int nodeNumber)
  {
    center_noGhost = nodeNumber;
    int i,j,k;
    i = nodeNumber/local_nxyNodes;
    j = (nodeNumber%local_nxyNodes)/local_nxNodes;
    k = nodeNumber%local_nxNodes;
    
    if(center_noGhost != i*local_nxyNodes + j*local_nxNodes + k)
      std::cerr<<"incorrect calculation"<<std::endl;
    
    setInterfaces(i,j,k);
    
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
    
    return (*this);
  }

inline void SecondOrderFd::operator++()
{
    ++center;
    ++center_noGhost;
    ++left;
    ++right;
    ++bottom;
    ++top;
    ++front;
    ++back;
    ++interLeft;
    ++interRight;
    ++interBottom;
    ++interTop;
    ++interFront;
    ++interBack;
  
}
}//Petsc
}//Daetk
#endif
