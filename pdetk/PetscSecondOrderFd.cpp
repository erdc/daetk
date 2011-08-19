#include "PetscSecondOrderFd.h"

#include "Tracer.h"
namespace Daetk 
{
  namespace Petsc
  {
    namespace cc
    {
      extern "C"
      {
#include "petsc.h"
#include "petscda.h"
      }
    }

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::istream;
  
void SecondOrderFd::setDomainSize(int nxNodesIn, int nyNodesIn,int nzNodesIn)
  {
    nxNodes = nxNodesIn;
    nyNodes = nyNodesIn;
    nzNodes = nzNodesIn;
    nxyNodes =nxNodesIn*nyNodesIn;
    nxyzNodes = nxyNodes*nzNodes;
    if (nyNodes > 1)
      NOT_ONED=1;
    if (nzNodes > 1)
      NOT_TWOD=1;
  }

  SecondOrderFd::SecondOrderFd(int nxNodesIn, int nyNodesIn, int nzNodesIn, int dof, int stencilWidth):
    nxNodes(nxNodesIn),
    nyNodes(nyNodesIn),
    nzNodes(nzNodesIn),
    nxyNodes(nxNodesIn*nyNodesIn),
    nxyzNodes(nxyNodes*nzNodes),
    local_x0(0),
    local_y0(0),
    local_z0(0),
    globalLow(0),
    globalHigh(nxyzNodes),
    local_nxNodes(nxNodes),
    local_nyNodes(nyNodes),
    local_nzNodes(nzNodes),
    local_nxyNodes(nxyNodes),
    local_nxyzNodes(nxyzNodes),
    ghost_x0(0),
    ghost_y0(0),
    ghost_z0(0),
    ghost_nxNodes(nxNodes),
    ghost_nyNodes(nyNodes),
    ghost_nzNodes(nzNodes),
    ghost_nxyNodes(nxyNodes),
    ghost_nxyzNodes(nxyzNodes),
    ghost_xOffSet(0),
    ghost_yOffSet(0),
    ghost_zOffSet(0),
    ghost_globalLow(0),
    ghost_globalHigh(nxyzNodes),
    center(0),
    center_noGhost(0),
    left(0),
    right(0),
    bottom(0),
    top(0),
    front(0),
    back(0),
    interLeft(0),
    interRight(0),
    interBottom(0),
    interTop(0),
    interFront(0),
    interBack(0),
    npx(1),
    npy(1),
    npz(1),
    NOT_ONED(0),
    NOT_TWOD(0),
    pointListField(nxyzNodes)
  {    
    if (nyNodes > 1)
      NOT_ONED=1;
    if (nzNodes > 1)
      NOT_TWOD=1;

    if (NOT_TWOD)
      ctor3d(nxNodesIn,nyNodesIn,nzNodesIn,dof,stencilWidth);
    else if (NOT_ONED)
      ctor2d(nxNodesIn,nyNodesIn,1,dof,stencilWidth);
    else
      ctor1d(nxNodesIn,1,1,dof,stencilWidth);
  }

SecondOrderFd::SecondOrderFd(SecondOrderFd& fineStencil, int nxNodesIn, int nyNodesIn, int nzNodesIn, int dof, int stencilWidth):
    nxNodes(nxNodesIn),
    nyNodes(nyNodesIn),
    nzNodes(nzNodesIn),
    nxyNodes(nxNodesIn*nyNodesIn),
    nxyzNodes(nxyNodes*nzNodes),
    local_x0(0),
    local_y0(0),
    local_z0(0),
    globalLow(0),
    globalHigh(nxyzNodes),
    local_nxNodes(nxNodes),
    local_nyNodes(nyNodes),
    local_nzNodes(nzNodes),
    local_nxyNodes(nxyNodes),
    local_nxyzNodes(nxyzNodes),
    ghost_x0(0),
    ghost_y0(0),
    ghost_z0(0),
    ghost_nxNodes(nxNodes),
    ghost_nyNodes(nyNodes),
    ghost_nzNodes(nzNodes),
    ghost_nxyNodes(nxyNodes),
    ghost_nxyzNodes(nxyzNodes),
    ghost_xOffSet(0),
    ghost_yOffSet(0),
    ghost_zOffSet(0),
    ghost_globalLow(0),
    ghost_globalHigh(nxyzNodes),
    center(0),
    center_noGhost(0),
    left(0),
    right(0),
    bottom(0),
    top(0),
    front(0),
    back(0),
    interLeft(0),
    interRight(0),
    interBottom(0),
    interTop(0),
    interFront(0),
    interBack(0),
    npx(1),
    npy(1),
    npz(1),
    NOT_ONED(0),
    NOT_TWOD(0),
    pointListField(nxyzNodes)
  {    
    if (nyNodes > 1)
      NOT_ONED=1;
    if (nzNodes > 1)
      NOT_TWOD=1;

    if (NOT_TWOD)
      ctor3d(fineStencil,nxNodesIn,nyNodesIn,nzNodesIn,dof,stencilWidth);
    else if (NOT_ONED)
      ctor2d(fineStencil,nxNodesIn,nyNodesIn,1,dof,stencilWidth);
    else
      ctor1d(fineStencil,nxNodesIn,1,1,dof,stencilWidth);
  }


SecondOrderFd::~SecondOrderFd()
{
  cc::DADestroy(da_);
  cc::DADestroy(dadof_);
}

  
int SecondOrderFd::petscToGlobal(int nodeNumber)
{
  Err ierr;
  ierr = cc::AOPetscToApplication(ao_,1,&nodeNumber);
  return nodeNumber;
}

void SecondOrderFd::ctor1d(int nxNodesIn, int nyNodes, int nzNodes , int dof, int stencilWidth)
{
  stencilSize=3;
  using namespace cc;
  assert(nyNodes ==1);
  assert(nzNodes ==1);
  Err ierr;
  ierr = DACreate1d(PETSC_COMM_WORLD,DA_NONPERIODIC,
                    nxNodes,
                    1,stencilWidth,
                    PETSC_NULL,&da_);   
  ierr = DACreate1d(PETSC_COMM_WORLD,DA_NONPERIODIC,
                    nxNodes,
                    dof,stencilWidth,
                    PETSC_NULL,&dadof_);   
  int dim,M,N,P,doftest,s;
  DAPeriodicType wrap;
  DAStencilType st;
  ierr =  DAGetInfo(da_,&dim,&M,&N,&P,&npx,&npy,&npz,&doftest,&s,&wrap,&st);
  
  ierr =  DAGetCorners(da_, &local_x0, &local_y0, &local_z0,
                       &local_nxNodes, &local_nyNodes, &local_nzNodes);
  local_nxyNodes = local_nxNodes*local_nyNodes;
  local_nxyzNodes = local_nzNodes*local_nxyNodes;

  ierr =  DAGetGhostCorners(da_, &ghost_x0, &ghost_y0, &ghost_z0,
                            &ghost_nxNodes, &ghost_nyNodes, &ghost_nzNodes);
  ghost_nxyNodes = ghost_nxNodes*ghost_nyNodes;
  ghost_nxyzNodes = ghost_nzNodes*ghost_nxyNodes;

   ghost_xOffSet = local_x0 - ghost_x0;
  ghost_yOffSet = 0;
  ghost_zOffSet = 0;

  ierr = DAGetAO(da_,&ao_);

  int i=0,j=0;
  for (int k=0;k<nxNodes;k++)
    {
      center     = k;
      left       = k-1;
      right      = k+1;
      
      pointList = &pointListField[center];
      pointList->reserve(7);

      ierr = AOApplicationToPetsc(ao_,1,&center);
      pointList->push_back(Location(i,j,k,center,CENTER));      

      if (k-1 >=0)
        {
          ierr = AOApplicationToPetsc(ao_,1,&left);
          pointList->push_back(Location(i,j,k-1,left,LEFT));
        }
      if (k+1 < nxNodes)
        {
          ierr = AOApplicationToPetsc(ao_,1,&right);
          pointList->push_back(Location(i,j,k+1,right,RIGHT));
        }    

    }      
  
  (*this)(local_x0);
  globalLow = center;
  globalHigh = center + local_nxNodes;
}

void SecondOrderFd::ctor1d(SecondOrderFd& fineStencil, int nxNodesIn, int nyNodes, int nzNodes , int dof, int stencilWidth)
{
  stencilSize=3;
  using namespace cc;
  assert(nyNodes ==1);
  assert(nzNodes ==1);
  Err ierr;
  int* lc = new int[fineStencil.npx];
  for (int i=0;i<fineStencil.npx;i++)
    lc[i]=nxNodes/fineStencil.npx;;
  ierr = DACreate1d(PETSC_COMM_WORLD,DA_NONPERIODIC,
                    nxNodes,
                    1,stencilWidth,
                    lc,&da_);   
  ierr = DACreate1d(PETSC_COMM_WORLD,DA_NONPERIODIC,
                    nxNodes,
                    dof,stencilWidth,
                    lc,&dadof_);   
  delete [] lc;
  int dim,M,N,P,doftest,s;
  DAPeriodicType wrap;
  DAStencilType st;
  ierr =  DAGetInfo(da_,&dim,&M,&N,&P,&npx,&npy,&npz,&doftest,&s,&wrap,&st);
  assert(npx == fineStencil.npx);
  assert(M == nxNodes);

  ierr =  DAGetCorners(da_, &local_x0, &local_y0, &local_z0,
                       &local_nxNodes, &local_nyNodes, &local_nzNodes);
  local_nxyNodes = local_nxNodes*local_nyNodes;
  local_nxyzNodes = local_nzNodes*local_nxyNodes;

  ierr =  DAGetGhostCorners(da_, &ghost_x0, &ghost_y0, &ghost_z0,
                            &ghost_nxNodes, &ghost_nyNodes, &ghost_nzNodes);
  ghost_nxyNodes = ghost_nxNodes*ghost_nyNodes;
  ghost_nxyzNodes = ghost_nzNodes*ghost_nxyNodes;

  ghost_xOffSet = local_x0 - ghost_x0;
  ghost_yOffSet = 0;
  ghost_zOffSet = 0;

  ierr = DAGetAO(da_,&ao_);

  int i=0,j=0;
  for (int k=0;k<nxNodes;k++)
    {
      center     = k;
      left       = k-1;
      right      = k+1;
      
       pointList = &pointListField[center];
       pointList->reserve(7);

      ierr = AOApplicationToPetsc(ao_,1,&center);
      pointList->push_back(Location(i,j,k,center,CENTER));      

      if (k-1 >=0)
        {
          ierr = AOApplicationToPetsc(ao_,1,&left);
          pointList->push_back(Location(i,j,k-1,left,LEFT));
        }
      if (k+1 < nxNodes)
        {
          ierr = AOApplicationToPetsc(ao_,1,&right);
          pointList->push_back(Location(i,j,k+1,right,RIGHT));
        }    

    }      
  
  (*this)(local_x0);
  globalLow = center;
  globalHigh = center + local_nxNodes;
}

void SecondOrderFd::ctor2d(int nxNodesIn, int nyNodesIn, int nzNodes, int dof, int stencilWidth)
{
  stencilSize=5;
  using namespace cc;
  assert(nzNodes==1);
  Err ierr;
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,
                    nxNodes,nyNodes,
                    PETSC_DECIDE,PETSC_DECIDE,
                    1,stencilWidth,
                    PETSC_NULL,PETSC_NULL,&da_);
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,
                    nxNodes,nyNodes,
                    PETSC_DECIDE,PETSC_DECIDE,
                    dof,stencilWidth,
                    PETSC_NULL,PETSC_NULL,&dadof_);
  int dim,M,N,P,doftest,s;
  DAPeriodicType wrap;
  DAStencilType st;
  ierr =  DAGetInfo(da_,&dim,&M,&N,&P,&npx,&npy,&npz,&doftest,&s,&wrap,&st);
  
  ierr =  DAGetCorners(da_, &local_x0, &local_y0, &local_z0,
                       &local_nxNodes, &local_nyNodes, &local_nzNodes);
  local_nxyNodes = local_nxNodes*local_nyNodes;
  local_nxyzNodes = local_nzNodes*local_nxyNodes;
  
  ierr =  DAGetGhostCorners(da_, &ghost_x0, &ghost_y0, &ghost_z0,
                            &ghost_nxNodes, &ghost_nyNodes, &ghost_nzNodes);
  ghost_nxyNodes = ghost_nxNodes*ghost_nyNodes;
  ghost_nxyzNodes = ghost_nzNodes*ghost_nxyNodes;


  ghost_xOffSet = local_x0 - ghost_x0;
  ghost_yOffSet = local_y0 - ghost_y0;
  ghost_zOffSet = 0;

  ierr = DAGetAO(da_,&ao_);

  int i=0;
  for (int j=0;j<nyNodes;j++)
    for (int k=0;k<nxNodes;k++)
      {
        center      = j*nxNodes+k;
        back        = (j+1)*nxNodes+k;
        front       = (j-1)*nxNodes+k;
        right       = j*nxNodes + k+1;
        left        = j*nxNodes + k-1;
        
        pointList = &pointListField[center];
        pointList->reserve(7);

        ierr = AOApplicationToPetsc(ao_,1,&center);
        pointList->push_back(Location(i,j,k,center,CENTER));

        if (k-1 >=0)
          {
            ierr = AOApplicationToPetsc(ao_,1,&left);
            pointList->push_back(Location(i,j,k-1,left,LEFT));
          }
        if (k+1 < nxNodes)
          {
            ierr = AOApplicationToPetsc(ao_,1,&right);
            pointList->push_back(Location(i,j,k+1,right,RIGHT));
          }
        if (j-1 >=0)
          {
            ierr = AOApplicationToPetsc(ao_,1,&front);
            pointList->push_back(Location(i,j-1,k,front,FRONT));
          }
        if (j+1 < nyNodes)
          {
            ierr = AOApplicationToPetsc(ao_,1,&back);
            pointList->push_back(Location(i,j+1,k,back,BACK));
          }
      }

  (*this)(local_y0,local_x0);
  globalLow = center;
  globalHigh = center + local_nxyNodes;
}

void SecondOrderFd::ctor2d(SecondOrderFd& fineStencil, int nxNodesIn, int nyNodesIn, int nzNodes, int dof, int stencilWidth)
{
  stencilSize=5;
  using namespace cc;
  assert(nzNodes==1);
  Err ierr;
  int* lx = new int[fineStencil.npx];
  for (int i=0;i<fineStencil.npx;i++)
    lx[i]=nxNodes/fineStencil.npx;

  int* ly = new int[fineStencil.npy];
  for (int i=0;i<fineStencil.npy;i++)
    ly[i]=nyNodes/fineStencil.npy;
  
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,
                    nxNodes,nyNodes,
                    fineStencil.npx, fineStencil.npy,
                    1,stencilWidth,
                    lx,ly,&da_);
  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,
                    nxNodes,nyNodes,
                    fineStencil.npx, fineStencil.npy,
                    dof,stencilWidth,
                    lx,ly,&dadof_);
  delete [] lx;
  delete [] ly;

  int dim,M,N,P,doftest,s;
  DAPeriodicType wrap;
  DAStencilType st;
  ierr =  DAGetInfo(da_,&dim,&M,&N,&P,&npx,&npy,&npz,&doftest,&s,&wrap,&st);
  assert( M == nxNodes);
  assert( N == nyNodes);
  assert( npx == fineStencil.npx);
  assert( npy == fineStencil.npy);

  ierr =  DAGetCorners(da_, &local_x0, &local_y0, &local_z0,
                       &local_nxNodes, &local_nyNodes, &local_nzNodes);
  local_nxyNodes = local_nxNodes*local_nyNodes;
  local_nxyzNodes = local_nzNodes*local_nxyNodes;
  
  ierr =  DAGetGhostCorners(da_, &ghost_x0, &ghost_y0, &ghost_z0,
                            &ghost_nxNodes, &ghost_nyNodes, &ghost_nzNodes);
  ghost_nxyNodes = ghost_nxNodes*ghost_nyNodes;
  ghost_nxyzNodes = ghost_nzNodes*ghost_nxyNodes;

  ghost_xOffSet = local_x0 - ghost_x0;
  ghost_yOffSet = local_y0 - ghost_y0;
  ghost_zOffSet = 0;

  ierr = DAGetAO(da_,&ao_);

  int i=0;
  for (int j=0;j<nyNodes;j++)
    for (int k=0;k<nxNodes;k++)
      {
        center      = j*nxNodes+k;
        back        = (j+1)*nxNodes+k;
        front       = (j-1)*nxNodes+k;
        right       = j*nxNodes + k+1;
        left        = j*nxNodes + k-1;
        
        pointList = &pointListField[center];
        pointList->reserve(7);

        ierr = AOApplicationToPetsc(ao_,1,&center);
        pointList->push_back(Location(i,j,k,center,CENTER));
        if (k-1 >=0)
          {
            ierr = AOApplicationToPetsc(ao_,1,&left);
            pointList->push_back(Location(i,j,k-1,left,LEFT));
          }
        if (k+1 < nxNodes)
          {
            ierr = AOApplicationToPetsc(ao_,1,&right);
            pointList->push_back(Location(i,j,k+1,right,RIGHT));
          }
        if (j-1 >=0)
          {
            ierr = AOApplicationToPetsc(ao_,1,&front);
            pointList->push_back(Location(i,j-1,k,front,FRONT));
          }
        if (j+1 < nyNodes)
          {
            ierr = AOApplicationToPetsc(ao_,1,&back);
            pointList->push_back(Location(i,j+1,k,back,BACK));
          }
      }

  (*this)(local_y0,local_x0);
  globalLow = center;
  globalHigh = center + local_nxyNodes;
}

void SecondOrderFd::ctor3d(int nxNodesIn, int nyNodesIn, int nzNodesIn, int dof, int stencilWidth)
{ 
  stencilSize=7;
  using namespace cc;
  Err ierr;
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,
                    nxNodes,nyNodes,nzNodes,
                    PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                    1,stencilWidth,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,&da_);
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,
                    nxNodes,nyNodes,nzNodes,
                    PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                    dof,stencilWidth,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,&dadof_);

  int dim,M,N,P,doftest,s;
  DAPeriodicType wrap;
  DAStencilType st;
  ierr =  DAGetInfo(da_,&dim,&M,&N,&P,&npx,&npy,&npz,&doftest,&s,&wrap,&st);
  ierr =  DAGetCorners(da_, &local_x0, &local_y0, &local_z0,
                       &local_nxNodes, &local_nyNodes, &local_nzNodes);
  local_nxyNodes = local_nxNodes*local_nyNodes;
  local_nxyzNodes = local_nzNodes*local_nxyNodes;
  ierr =  DAGetGhostCorners(da_, &ghost_x0, &ghost_y0, &ghost_z0,
                            &ghost_nxNodes, &ghost_nyNodes, &ghost_nzNodes);
  ghost_nxyNodes = ghost_nxNodes*ghost_nyNodes;
  ghost_nxyzNodes = ghost_nzNodes*ghost_nxyNodes;

  ghost_xOffSet = local_x0 - ghost_x0;
  ghost_yOffSet = local_y0 - ghost_y0;
  ghost_zOffSet = local_z0 - ghost_z0;
  
  ierr = DAGetAO(da_,&ao_);
  for (int i=0;i<nzNodes;i++)
    for (int j=0;j<nyNodes;j++)
      for (int k=0;k<nxNodes;k++)
        {
          center      = i*nxyNodes + j*nxNodes+k;
          top         = (i+1)*nxyNodes + j*nxNodes + k;
          bottom      = (i-1)*nxyNodes + j*nxNodes + k;
          back        = i*nxyNodes + (j+1)*nxNodes+k;
          front       = i*nxyNodes + (j-1)*nxNodes+k;
          right       = i*nxyNodes + j*nxNodes + k+1;
          left        = i*nxyNodes + j*nxNodes + k-1;
                    
          pointList = &pointListField[center];
          pointList->reserve(7);

          ierr = AOApplicationToPetsc(ao_,1,&center);
          pointList->push_back(Location(i,j,k,center,CENTER));
          if (k-1 >=0)
            {
              ierr = AOApplicationToPetsc(ao_,1,&left);
              pointList->push_back(Location(i,j,k-1,left,LEFT));
            }
          if (k+1 < nxNodes)
            {
              ierr = AOApplicationToPetsc(ao_,1,&right);
              pointList->push_back(Location(i,j,k+1,right,RIGHT));
            }
          if (j-1 >=0)
            {
              ierr = AOApplicationToPetsc(ao_,1,&front);
              pointList->push_back(Location(i,j-1,k,front,FRONT));
           }
          if (j+1 < nyNodes)
            {
              ierr = AOApplicationToPetsc(ao_,1,&back);
              pointList->push_back(Location(i,j+1,k,back,BACK));
            }
          if (i-1 >=0)
            {
              ierr = AOApplicationToPetsc(ao_,1,&bottom);
              pointList->push_back(Location(i-1,j,k,bottom,BOTTOM));
            }
          if (i+1 < nzNodes)
            {
              ierr = AOApplicationToPetsc(ao_,1,&top);
              pointList->push_back(Location(i+1,j,k,top,TOP));
            }
        }
  
  (*this)(local_z0,local_y0,local_x0);
  globalLow = center;
  globalHigh = center + local_nxyzNodes;
}

void SecondOrderFd::ctor3d(SecondOrderFd& fineStencil, int nxNodesIn, int nyNodesIn, int nzNodesIn, int dof, int stencilWidth)
{ 
  stencilSize=7;
  using namespace cc;
  Err ierr;
  int* lx = new int[fineStencil.npx];
  for (int i=0;i<fineStencil.npx;i++)
    lx[i]=nxNodes/fineStencil.npx;

  int* ly = new int[fineStencil.npy];
  for (int i=0;i<fineStencil.npy;i++)
    ly[i]=nyNodes/fineStencil.npy;
  
  int* lz = new int[fineStencil.npz];
  for (int i=0;i<fineStencil.npz;i++)
    lz[i]=nzNodes/fineStencil.npz;
  
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
                    nxNodes,nyNodes,nzNodes,
                    fineStencil.npx, fineStencil.npy, fineStencil.npz,
                    1,stencilWidth,
                    lx,ly,lz,&da_);
  ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
                    nxNodes,nyNodes,nzNodes,
                    fineStencil.npx, fineStencil.npy, fineStencil.npz,
                    dof,stencilWidth,
                    lx,ly,lz,&dadof_);

  delete [] lx;
  delete [] ly;
  delete [] lz;

  int dim,M,N,P,doftest,s;
  DAPeriodicType wrap;
  DAStencilType st;
  ierr =  DAGetInfo(da_,&dim,&M,&N,&P,&npx,&npy,&npz,&doftest,&s,&wrap,&st);
  assert( fineStencil.npx == npx);
  assert( fineStencil.npy == npy);
  assert( fineStencil.npz == npz);
  assert( nxNodes == M);
  assert( nyNodes == N);
  assert( nzNodes == P);

  ierr =  DAGetCorners(da_, &local_x0, &local_y0, &local_z0,
                       &local_nxNodes, &local_nyNodes, &local_nzNodes);
  local_nxyNodes = local_nxNodes*local_nyNodes;
  local_nxyzNodes = local_nzNodes*local_nxyNodes;
  ierr =  DAGetGhostCorners(da_, &ghost_x0, &ghost_y0, &ghost_z0,
                            &ghost_nxNodes, &ghost_nyNodes, &ghost_nzNodes);
  ghost_nxyNodes = ghost_nxNodes*ghost_nyNodes;
  ghost_nxyzNodes = ghost_nzNodes*ghost_nxyNodes;

  ghost_xOffSet = local_x0 - ghost_x0;
  ghost_yOffSet = local_y0 - ghost_y0;
  ghost_zOffSet = local_z0 - ghost_z0;
  
  ierr = DAGetAO(da_,&ao_);
  for (int i=0;i<nzNodes;i++)
    for (int j=0;j<nyNodes;j++)
      for (int k=0;k<nxNodes;k++)
        {
          center      = i*nxyNodes + j*nxNodes+k;
          top         = (i+1)*nxyNodes + j*nxNodes + k;
          bottom      = (i-1)*nxyNodes + j*nxNodes + k;
          back        = i*nxyNodes + (j+1)*nxNodes+k;
          front       = i*nxyNodes + (j-1)*nxNodes+k;
          right       = i*nxyNodes + j*nxNodes + k+1;
          left        = i*nxyNodes + j*nxNodes + k-1;
                    
          pointList = &pointListField[center];
          pointList->reserve(7);

          ierr = AOApplicationToPetsc(ao_,1,&center);
          pointList->push_back(Location(i,j,k,center,CENTER));
          if (k-1 >=0)
            {
              ierr = AOApplicationToPetsc(ao_,1,&left);
              pointList->push_back(Location(i,j,k-1,left,LEFT));
            }
          if (k+1 < nxNodes)
            {
              ierr = AOApplicationToPetsc(ao_,1,&right);
              pointList->push_back(Location(i,j,k+1,right,RIGHT));
            }
          if (j-1 >=0)
            {
              ierr = AOApplicationToPetsc(ao_,1,&front);
              pointList->push_back(Location(i,j-1,k,front,FRONT));
           }
          if (j+1 < nyNodes)
            {
              ierr = AOApplicationToPetsc(ao_,1,&back);
              pointList->push_back(Location(i,j+1,k,back,BACK));
            }
          if (i-1 >=0)
            {
              ierr = AOApplicationToPetsc(ao_,1,&bottom);
              pointList->push_back(Location(i-1,j,k,bottom,BOTTOM));
            }
          if (i+1 < nzNodes)
            {
              ierr = AOApplicationToPetsc(ao_,1,&top);
              pointList->push_back(Location(i+1,j,k,top,TOP));
            }
        }
  

  (*this)(local_z0,local_y0,local_x0);
  globalLow = center;
  globalHigh = center + local_nxyzNodes;
}

}//Petsc
}//Daetk
