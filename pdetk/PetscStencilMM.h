#ifndef PETSC_STENCILMM_H
#define PETSC_STENCILMM_H
#include "Stencil.h"
#include "DaetkPetscSys.h"
#include "Definitions.h"
#include <vector>

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
class StencilMM : public Stencil
{
public:
  enum Point {BOTTOM,FRONT,LEFT,CENTER,RIGHT,BACK,TOP,LEFT_FRONT,
	      LEFT_BACK,RIGHT_FRONT,RIGHT_BACK};
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
    //mwf add for corner terms 2d only for now
    leftBack,   //(i-1,j+1)
    leftFront,  //(i-1,j-1)
    rightBack,  //(i+1,j+1)
    rightFront, //(i+1,j-1)
    // 
    interLeft,  
    interRight,
    interBottom,
    interTop,
    interFront,
    interBack,
    //mwf add for corner interface terms over ghosted region
    //mwf 2d only for now
    interLeftBackGhost,   //(i-1,j+1/2)
    interLeftFrontGhost,  //(i-1,j-1/2)
    interRightBackGhost,  //(i+1,j+1/2)
    interRightFrontGhost, //(i+1,j-1/2)
    interBackRightGhost,  //(i+1/2,j+1)
    interBackLeftGhost,   //(i-1/2,j+1)
    interFrontRightGhost, //(i+1/2,j-1)
    interFrontLeftGhost,  //(i-1/2,j-1)
    //mwf these also account for ghosted region interfaces
    interLeftGhost,   //(i-1/2,j,k)
    interRightGhost,  //(i+1/2,j,k)
    interFrontGhost,  //(i,j-1/2,k)
    interBackGhost,   //(i,j+1/2,k)
    interTopGhost,    //(i,j,k+1/2)
    interBottomGhost, //(i,j,k-1/2)
    npx,
    npy,
    npz,
    stencilSize;

  inline int size(){return stencilSize;}
  //mwf added for boundary conditions in transverse direction?
  bool isInGhostRegion;
  
  int nNodes() { return nxyzNodes;}
  
  StencilMM(int nxNodesIn, int nyNodes=1, int nzNodesIn =1, int dof=1, int stencilWidth=1);
  StencilMM(StencilMM& fineStencil, int nxNodesIn, int nyNodesIn=1, int nzNodesIn =1, int dof=1, int stencilWidth=1);
  
  virtual ~StencilMM();
  
  void setDomainSize(int nxNodesIn=1, int nyNodesIn=1,int nzNodesIn=1)
  {
    nxNodes = nxNodesIn;
    nyNodes = nyNodesIn;
    nzNodes = nzNodesIn;
    nxyNodes =nxNodesIn*nyNodesIn;
    nxyzNodes = nxyNodes*nzNodes;
  }
  
  StencilMM& operator()(int k)
  {
    int i=0,j=0;
    if (k>= local_x0 && k < (local_x0 + local_nxNodes))
      {
        setInterfaces(i,j,k);
        isLocal=true;
        isInGhostRegion=true;
      }
    else
      {
	isLocal=false;
	//mwf now determine if in ghostRegion as well?
	if (k>= ghost_x0 && k < (ghost_x0 + ghost_nxNodes))
	  {
	    isInGhostRegion=true;
	  }
	else
	  isInGhostRegion=false;
      }
//      pointMap = &pointMapField[k];
//      anchor = &pointMapField[k][CENTER];
    pointList = &pointListField[k];
    anchor = &pointListField[k][CENTER];
    
    setPoints();
    return *this;
  }
  
  StencilMM& operator()(int j,int k)
  {
    int i=0;
    if (k>= local_x0 && k < (local_x0 + local_nxNodes) &&
        j>= local_y0 && j < (local_y0 + local_nyNodes))
      {
        setInterfaces(i,j,k);
        isLocal=true;
        isInGhostRegion=true;
      }
    else
      {
	isLocal=false;
	//mwf now determine if in ghostRegion as well?
	if (k>= ghost_x0 && k < (ghost_x0 + ghost_nxNodes) &&
	    j>= ghost_y0 && j < (ghost_y0 + ghost_nyNodes))
	  {
	    isInGhostRegion=true;
	  }
	else
	  isInGhostRegion=false;
      }

//      pointMap = &pointMapField[j*nxNodes+k];
//      anchor = &pointMapField[j*nxNodes+k][CENTER];
    pointList = &pointListField[j*nxNodes+k];
    anchor = &pointListField[j*nxNodes+k][CENTER];
    
    setPoints();
    return *this;
  }
  
  StencilMM& operator()(int i, int j, int k)
  {
    if ((k>= local_x0) && (k < (local_x0 + local_nxNodes)) &&
        (j>= local_y0) && (j < (local_y0 + local_nyNodes)) &&
        (i>= local_z0) && (i < (local_z0 + local_nzNodes)))
      {
        setInterfaces(i,j,k);
        isLocal=true;
        isInGhostRegion=true;
      }
    else
      {
	isLocal=false;

	//mwf now determine if in ghostRegion as well?
	if ((k>= ghost_x0) && (k < (ghost_x0 + ghost_nxNodes)) &&
	    (j>= ghost_y0) && (j < (ghost_y0 + ghost_nyNodes)) &&
	    (i>= ghost_z0) && (i < (ghost_z0 + ghost_nzNodes)))
	  {
	    isInGhostRegion=true;
	  }
	else
	  isInGhostRegion=false;
      }
//      pointMap = &pointMapField[i*nxyNodes + j*nxNodes+k];
//      anchor =&pointMapField[i*nxyNodes + j*nxNodes+k][CENTER];
    pointList = &pointListField[i*nxyNodes + j*nxNodes+k];
    anchor =&pointListField[i*nxyNodes + j*nxNodes+k][CENTER];
    
    setPoints();
    return *this;
  }
  
  int globalToLocal(int nodeNumber)
  {
    Location& n=pointListField[nodeNumber][CENTER];
    return (n.i-ghost_z0)*ghost_nxyNodes 
      + (n.j - ghost_y0)*ghost_nxNodes
      + (n.k - ghost_x0);
  }

  StencilMM& globalNode(int nodeNumber)
  {
//      pointMap = &pointMapField[nodeNumber];
//      anchor =&pointMapField[nodeNumber][CENTER];
    pointList = &pointListField[nodeNumber];
    anchor =&pointListField[nodeNumber][CENTER];
    int i = anchor->i;
    int j = anchor->j;
    int k = anchor->k;
    
    if ((k>= local_x0) && (k < (local_x0 + local_nxNodes)) &&
        (j>= local_y0) && (j < (local_y0 + local_nyNodes)) &&
        (i>= local_z0) && (i < (local_z0 + local_nzNodes)))
      {
        setInterfaces(i,j,k);
        isLocal=true;
        isInGhostRegion=true;
      }
    else
      {
	isLocal=false;
	//mwf try and see if it's in the ghost region for bc's
	if ((k>= ghost_x0) && (k < (ghost_x0 + ghost_nxNodes)) &&
	    (j>= ghost_y0) && (j < (ghost_y0 + ghost_nyNodes)) &&
	    (i>= ghost_z0) && (i < (ghost_z0 + ghost_nzNodes)))
	  {
	    isInGhostRegion=true;
	  }
	else
	  isInGhostRegion=false;
      }
    setPoints();
    return (*this);
  }
  
  void setGlobalPointList(int i, int j, int k)
  {
    pointList = &pointListField[i*nxyNodes + j*nxNodes+k];
    anchor =&pointListField[i*nxyNodes + j*nxNodes+k][CENTER];
  }

  void setGlobalPointList(int nodeNumber)
  {
    pointList = &pointListField[nodeNumber];
    anchor =&pointListField[nodeNumber][CENTER];
  }

  int petscToGlobal(int nodeNumber);

  void setInterfaces(int i, int j, int k)
  {
    //mwf go ahead and try to calculate global interfaces?
    //now do local
    i-= local_z0;
    j-= local_y0;
    k-= local_x0;

    //(i-1,j+1/2)
    interLeftBackGhost = i*(ghost_nxyNodes+ghost_nxNodes) 
      + (j+1)*ghost_nxNodes+k-1;
    //(i-1,j-1/2)
    interLeftFrontGhost= i*(ghost_nxyNodes+ghost_nxNodes) 
      + j*ghost_nxNodes+k-1;
    //(i+1,j+1/2)
    interRightBackGhost = i*(ghost_nxyNodes+ghost_nxNodes) 
      + (j+1)*ghost_nxNodes+k+1;
    //(i+1,j-1/2)
    interRightFrontGhost= i*(ghost_nxyNodes+ghost_nxNodes) 
      + j*ghost_nxNodes+k+1;
    //(i+1/2,j+1)
    interBackRightGhost = i*(ghost_nxyNodes+ghost_nyNodes) 
      + (j+1)*(ghost_nxNodes+1)+k+1;
    //(i-1/2,j+1)
    interBackLeftGhost  = i*(ghost_nxyNodes+ghost_nyNodes) 
      + (j+1)*(ghost_nxNodes+1)+k;
    //(i+1/2,j-1)
    interFrontRightGhost= i*(ghost_nxyNodes+ghost_nyNodes) 
      + (j-1)*(ghost_nxNodes+1)+k+1;
    //(i-1/2,j-1)
    interFrontLeftGhost = i*(ghost_nxyNodes+ghost_nyNodes) 
      + (j-1)*(ghost_nxNodes+1)+k;

    interTopGhost    = (i+1)*ghost_nxyNodes + j*ghost_nxNodes + k;
    interBottomGhost = i*ghost_nxyNodes + j*ghost_nxNodes + k;
    interBackGhost   = i*(ghost_nxyNodes+ghost_nxNodes) 
      + (j+1)*ghost_nxNodes+k;
    interFrontGhost  = i*(ghost_nxyNodes+ghost_nxNodes) 
      + j*ghost_nxNodes+k;
    interRightGhost  = i*(ghost_nxyNodes+ghost_nyNodes) 
      + j*(ghost_nxNodes+1)+k+1;
    interLeftGhost   = i*(ghost_nxyNodes+ghost_nyNodes) 
      + j*(ghost_nxNodes+1)+k;
    //
    //i-= local_z0;
    //j-= local_y0;
    //k-= local_x0;
    //
    interTop    = (i+1)*local_nxyNodes + j*local_nxNodes + k;
    interBottom = i*local_nxyNodes + j*local_nxNodes + k;
    interBack   = i*(local_nxyNodes+local_nxNodes) + (j+1)*local_nxNodes+k;
    interFront  = i*(local_nxyNodes+local_nxNodes) + j*local_nxNodes+k;
    interRight  = i*(local_nxyNodes+local_nyNodes) + j*(local_nxNodes+1)+k+1;
    interLeft   = i*(local_nxyNodes+local_nyNodes) + j*(local_nxNodes+1)+k;

  }

  void setPoints()
  {
//      PointMap::iterator sit=pointMap->begin(),end=pointMap->end();
    PointList::iterator sit=pointList->begin(),end=pointList->end();
    
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
	    {
	      right       = sit->globalNodeNumber;
	      //mwf debug this
	      //std::cerr<<"setPoints sit->point==RIGHT ("<<sit->i<<","<<sit->j
	      //  <<","<<sit->k<<")= "<<right<<std::endl;
	    }
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
	  case LEFT_FRONT:
	    leftFront   = sit->globalNodeNumber;
	    break;
	  case LEFT_BACK:
	    leftBack    = sit->globalNodeNumber;
	    break;
	  case RIGHT_FRONT:
	    rightFront  = sit->globalNodeNumber;
	    break;
	  case RIGHT_BACK:
	    rightBack   = sit->globalNodeNumber;
	    break;
          default:
            std::cerr<<"enum value out of range in StencilMM::setPoints()"<<std::endl; 
            break;
          }
        ++sit;
      }
  }
  
  StencilMM& setLocalPointList()
  {
    pointList = &localPointListField[center];
    anchor =&localPointListField[center][CENTER];
    return *this;
  }

  StencilMM& localIndex(int k)
  {
    isLocal = true;

    interRight = k+1;
    interLeft  = k;

    center_noGhost = k;


    k+=ghost_xOffSet;
    
    center     = k;
    left       = k-1;
    right      = k+1;

    //mwf for flux arrays with dimensions over whole ghosted region?
    interLeftGhost = k;
    interRightGhost= k+1;

    return *this;
  }
  
  StencilMM& localIndex(int j,int k)
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
    //mwf add for corner terms 2d only for now
    //(i-1,j+1)
    leftBack    = (j+1)*ghost_nxNodes+k-1;
    //(i-1,j-1)
    leftFront   = (j-1)*ghost_nxNodes+k-1; 
    //(i+1,j+1)
    rightBack   = (j+1)*ghost_nxNodes + k+1; 
    //(i+1,j-1)
    rightFront  = (j-1)*ghost_nxNodes + k+1;

    //mwf for flux arrays with dimensions over whole ghosted region?
    interBackGhost   = (j+1)*ghost_nxNodes+k;
    interFrontGhost  = j*ghost_nxNodes+k;
    interRightGhost  = j*(ghost_nxNodes+1)+k+1;
    interLeftGhost   = j*(ghost_nxNodes+1)+k;

    //mwf now add these for corner points
    //(i-1,j+1/2)
    interLeftBackGhost = (j+1)*ghost_nxNodes+k-1;
    //(i-1,j-1/2)
    interLeftFrontGhost= j*ghost_nxNodes+k-1;
    //(i+1,j+1/2)
    interRightBackGhost = (j+1)*ghost_nxNodes+k+1;
    //(i+1,j-1/2)
    interRightFrontGhost= j*ghost_nxNodes+k+1;
    //(i+1/2,j+1)
    interBackRightGhost = (j+1)*(ghost_nxNodes+1)+k+1;
    //(i-1/2,j+1)
    interBackLeftGhost  = (j+1)*(ghost_nxNodes+1)+k;
    //(i+1/2,j-1)
    interFrontRightGhost= (j-1)*(ghost_nxNodes+1)+k+1;
    //(i-1/2,j-1)
    interFrontLeftGhost = (j-1)*(ghost_nxNodes+1)+k;
    
    return *this;
  }
  
  StencilMM& localIndex(int i, int j, int k)
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
    //mwf also missing corner points for 3d
    //(i-1,j+1)
    leftBack    = i*ghost_nxyNodes + (j+1)*ghost_nxNodes+k-1;
    //(i-1,j-1)
    leftFront   = i*ghost_nxyNodes + (j-1)*ghost_nxNodes+k-1; 
    //(i+1,j+1)
    rightBack   = i*ghost_nxyNodes + (j+1)*ghost_nxNodes + k+1; 
    //(i+1,j-1)
    rightFront  = i*ghost_nxyNodes + (j-1)*ghost_nxNodes + k+1;

    //mwf for flux arrays with dimensions over whole ghosted region?
    interTopGhost    = (i+1)*ghost_nxyNodes + j*ghost_nxNodes + k;
    interBottomGhost = i*ghost_nxyNodes + j*ghost_nxNodes + k;
    interBackGhost   = i*(ghost_nxyNodes+ghost_nxNodes) 
      + (j+1)*ghost_nxNodes+k;
    interFrontGhost  = i*(ghost_nxyNodes+ghost_nxNodes) 
      + j*ghost_nxNodes+k;
    interRightGhost  = i*(ghost_nxyNodes+ghost_nyNodes) 
      + j*(ghost_nxNodes+1)+k+1;
    interLeftGhost   = i*(ghost_nxyNodes+ghost_nyNodes) 
      + j*(ghost_nxNodes+1)+k;


    //mwf missing most of interfaces needed for 3d
    //(i-1,j+1/2)
    interLeftBackGhost = i*(ghost_nxyNodes+ghost_nxNodes) 
      + (j+1)*ghost_nxNodes+k-1;
    //(i-1,j-1/2)
    interLeftFrontGhost= i*(ghost_nxyNodes+ghost_nxNodes) 
      + j*ghost_nxNodes+k-1;
    //(i+1,j+1/2)
    interRightBackGhost = i*(ghost_nxyNodes+ghost_nxNodes) 
      + (j+1)*ghost_nxNodes+k+1;
    //(i+1,j-1/2)
    interRightFrontGhost= i*(ghost_nxyNodes+ghost_nxNodes) 
      + j*ghost_nxNodes+k+1;
    //(i+1/2,j+1)
    interBackRightGhost = i*(ghost_nxyNodes+ghost_nyNodes) 
      + (j+1)*(ghost_nxNodes+1)+k+1;
    //(i-1/2,j+1)
    interBackLeftGhost  = i*(ghost_nxyNodes+ghost_nyNodes) 
      + (j+1)*(ghost_nxNodes+1)+k;
    //(i+1/2,j-1)
    interFrontRightGhost= i*(ghost_nxyNodes+ghost_nyNodes) 
      + (j-1)*(ghost_nxNodes+1)+k+1;
    //(i-1/2,j-1)
    interFrontLeftGhost = i*(ghost_nxyNodes+ghost_nyNodes) 
      + (j-1)*(ghost_nxNodes+1)+k;


   
    return *this;
  }
  
  StencilMM& localNode(int nodeNumber)
  {
    isLocal = true;
//      pointMap = &localPointMapField[nodeNumber];
//      anchor =&localPointMapField[nodeNumber][CENTER];
    pointList = &localPointListField[nodeNumber];
    anchor =&localPointListField[nodeNumber][CENTER];
    center      = anchor->globalNodeNumber;
    top         = center + ghost_nxyNodes;
    bottom      = center - ghost_nxyNodes;
    back        = center + ghost_nxNodes;
    front       = center - ghost_nxNodes;
    right       = center+1;
    left        = center-1;
    //mwf I hope these are right?
    //(i-1,j+1)
    leftBack    = center - 1 + ghost_nxNodes; 
    //(i-1,j-1)
    leftFront   = center - 1 - ghost_nxNodes; 
    //(i+1,j+1)
    rightBack   = center + 1 + ghost_nxNodes;
    //(i+1,j-1)
    rightFront  = center + 1 - ghost_nxNodes; 
    //mwf should I do also the ghost interface points?
    //mwf missing corner points for 3d
    //just go ahead and try to use anchor's i,j,k here?
    //otherwise run into problems with dimension?
    //try and recreate localIndex routine
    int i = anchor->i;
    int j = anchor->j;
    int k = anchor->k;

    //mwf try and do flux indeces as well
    interTop    = (i+1)*local_nxyNodes + j*local_nxNodes + k;
    interBottom = i*local_nxyNodes + j*local_nxNodes + k;
    interBack   = i*(local_nxyNodes+local_nxNodes) + (j+1)*local_nxNodes+k;
    interFront  = i*(local_nxyNodes+local_nxNodes) + j*local_nxNodes+k;
    interRight  = i*(local_nxyNodes+local_nyNodes) + j*(local_nxNodes+1)+k+1;
    interLeft   = i*(local_nxyNodes+local_nyNodes) + j*(local_nxNodes+1)+k;
    center_noGhost = interBottom;

   //mwf I think i,j,k should be relative to "ghosted" region already?
    //mwf no, it looks like the i,j,k that get set can be negatiev
    i+=ghost_zOffSet;
    j+=ghost_yOffSet;
    k+=ghost_xOffSet;
    //hopefully i == 0 for 2d
    //mwf for flux arrays with dimensions over whole ghosted region?
    interTopGhost    = (i+1)*ghost_nxyNodes + j*ghost_nxNodes + k;
    interBottomGhost = i*ghost_nxyNodes + j*ghost_nxNodes + k;
    interBackGhost   = i*(ghost_nxyNodes+ghost_nxNodes) 
      + (j+1)*ghost_nxNodes+k;
    interFrontGhost  = i*(ghost_nxyNodes+ghost_nxNodes) 
      + j*ghost_nxNodes+k;
    interRightGhost  = i*(ghost_nxyNodes+ghost_nyNodes) 
      + j*(ghost_nxNodes+1)+k+1;
    interLeftGhost   = i*(ghost_nxyNodes+ghost_nyNodes) 
      + j*(ghost_nxNodes+1)+k;


    //mwf missing most of interfaces needed for 3d
    //(i-1,j+1/2)
    interLeftBackGhost = i*(ghost_nxyNodes+ghost_nxNodes) 
      + (j+1)*ghost_nxNodes+k-1;
    //(i-1,j-1/2)
    interLeftFrontGhost= i*(ghost_nxyNodes+ghost_nxNodes) 
      + j*ghost_nxNodes+k-1;
    //(i+1,j+1/2)
    interRightBackGhost = i*(ghost_nxyNodes+ghost_nxNodes) 
      + (j+1)*ghost_nxNodes+k+1;
    //(i+1,j-1/2)
    interRightFrontGhost= i*(ghost_nxyNodes+ghost_nxNodes) 
      + j*ghost_nxNodes+k+1;
    //(i+1/2,j+1)
    interBackRightGhost = i*(ghost_nxyNodes+ghost_nyNodes) 
      + (j+1)*(ghost_nxNodes+1)+k+1;
    //(i-1/2,j+1)
    interBackLeftGhost  = i*(ghost_nxyNodes+ghost_nyNodes) 
      + (j+1)*(ghost_nxNodes+1)+k;
    //(i+1/2,j-1)
    interFrontRightGhost= i*(ghost_nxyNodes+ghost_nyNodes) 
      + (j-1)*(ghost_nxNodes+1)+k+1;
    //(i-1/2,j-1)
    interFrontLeftGhost = i*(ghost_nxyNodes+ghost_nyNodes) 
      + (j-1)*(ghost_nxNodes+1)+k;
    
    return (*this);
  }
  
//    int localToGlobal(int localNodeNumber) 
//    {
//      return globalLow + localNodeNumber;
//    }
protected:

//    typedef std::vector< std::map<int, Location> > MapField; 
//    MapField pointMapField, localPointMapField;
  //  typedef std::vector< std::list<Location> > ListField; 
  typedef std::vector< PointList > ListField; 
  ListField pointListField, localPointListField;

 void ctor1d(int nxNodesIn, int nyNodes=1, int nzNodes =1, int dof=1, int stencilWidth=1);
  void ctor1d(StencilMM& fineStencil, int nxNodesIn, int nyNodes=1, int nzNodes =1, int dof=1, int stencilWidth=1);
  void ctor2d(int nxNodesIn, int nyNodesIn, int nzNodes=1, int dof=1, int stencilWidth=1);
  void ctor2d(StencilMM& fineStencil, int nxNodesIn, int nyNodesIn, int nzNodes=1, int dof=1, int stencilWidth=1);
  void ctor3d(int nxNodesIn, int nyNodesIn, int nzNodesIn, int dof=1, int stencilWidth=1);
  void ctor3d(StencilMM& fineStencil, int nxNodesIn, int nyNodesIn, int nzNodesIn, int dof=1, int stencilWidth=1);
};

typedef StencilMM StencilMM1d;
typedef StencilMM StencilMM2d;
typedef StencilMM StencilMM3d;

}//Petsc
}//Daetk
#endif






















