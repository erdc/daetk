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
  
void SecondOrderFd::setDomainSize(int nxNodesIn=1, int nyNodesIn=1,int nzNodesIn=1)
  {
    nxNodes = nxNodesIn;
    nyNodes = nyNodesIn;
    nzNodes = nzNodesIn;
    nxyNodes =nxNodesIn*nyNodesIn;
    nxyzNodes = nxyNodes*nzNodes;
  }

  SecondOrderFd::SecondOrderFd(int nxNodesIn=1, int nyNodesIn=1,int nzNodesIn=1):
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
    pointVectorField(nxyzNodes)
//      ownsVectorField(false);
  {
//      pointVectorField = new VectorField(nxyzNodes);
//      ownsVectorField(true);
  }

SecondOrderFd::~SecondOrderFd()
{
//    if (ownsVectorField)
//      delete pointVectorField;
  cc::DADestroy(da_);
  cc::DADestroy(dadof_);
}

  
int SecondOrderFd::petscToGlobal(int nodeNumber)
{
  Petsc::Err ierr;
  ierr = cc::AOPetscToApplication(ao_,1,&nodeNumber);
  return nodeNumber;
}

SecondOrderFd1d::SecondOrderFd1d(int nxNodesIn, int nyNodes, int nzNodes , int dof):
  SecondOrderFd(nxNodesIn,1,1)
{
  using namespace cc;
  assert(nyNodes ==1);
  assert(nzNodes ==1);
  Petsc::Err ierr;
  int stencilWidth=1;
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

//    localPointMapField.resize(ghost_nxyzNodes);
//    localPointListField.resize(ghost_nxyzNodes);
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
      
//        pointMap = &pointMapField[center];
//        pointList = &pointListField[center];
      pointVector = &pointVectorField[center];
      pointVector->reserve(3);

      ierr = AOApplicationToPetsc(ao_,1,&center);
//        pointMap->insert(PointMap::value_type(CENTER, Location(i,j,k,center,CENTER)));      
//        pointList->push_back(Location(i,j,k,center,CENTER));      
      (*pointVector)[CENTER]=Location(i,j,k,center,CENTER);      

      if (k-1 >=0)
        {
          ierr = AOApplicationToPetsc(ao_,1,&left);
//            pointMap->insert(PointMap::value_type(LEFT,Location(i,j,k-1,left,LEFT)));
//            pointList->push_back(Location(i,j,k-1,left,LEFT));
          (*pointVector)[LEFT]=Location(i,j,k-1,left,LEFT);
        }
      if (k+1 < nxNodes)
        {
          ierr = AOApplicationToPetsc(ao_,1,&right);
//            pointMap->insert(PointMap::value_type(RIGHT,Location(i,j,k+1,right,RIGHT)));
//            pointList->push_back(Location(i,j,k+1,right,RIGHT));
          (*pointVector)[RIGHT]=Location(i,j,k+1,right,RIGHT);
        }    


    }      
  
//    i=0;j=0;
//    for (int k=-ghost_xOffSet;k < ghost_nxNodes-ghost_xOffSet; k++)
//      {
//        this->localIndex(k);
//        this->setLocalPointList();
//        if (k-1 >=-ghost_xOffSet)
//  //          pointMap->insert(PointMap::value_type(LEFT, Location(i,j,k-1,left,LEFT)));
//          pointList->push_back(Location(i,j,k-1,left,LEFT));
//        if (k+1 < ghost_nxNodes-ghost_xOffSet)
//  //          pointMap->insert(PointMap::value_type(RIGHT, Location(i,j,k+1,right,RIGHT)));
//          pointList->push_back(Location(i,j,k+1,right,RIGHT));
//  //        pointMap->insert(PointMap::value_type(CENTER, Location(i,j,k,center,CENTER)));
//        pointList->push_back(Location(i,j,k,center,CENTER));
//      }          
  (*this)(local_x0);
  globalLow = center;
  globalHigh = center + local_nxNodes;
}

SecondOrderFd1d::SecondOrderFd1d(SecondOrderFd1d& fineStencil, int nxNodesIn, int nyNodes, int nzNodes , int dof):
  SecondOrderFd(nxNodesIn,1,1)
{
  using namespace cc;
  assert(nyNodes ==1);
  assert(nzNodes ==1);
  Petsc::Err ierr;
  int stencilWidth=1;
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

//    localPointMapField.resize(ghost_nxyzNodes);
//    localPointListField.resize(ghost_nxyzNodes);
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
      
//        pointMap = &pointMapField[center];
//        pointList = &pointListField[center];
      pointVector = &pointVectorField[center];
      pointVector->reserve(3);

      ierr = AOApplicationToPetsc(ao_,1,&center);
//        pointMap->insert(PointMap::value_type(CENTER, Location(i,j,k,center,CENTER)));      
//        pointList->push_back(Location(i,j,k,center,CENTER));      
      (*pointVector)[CENTER] = Location(i,j,k,center,CENTER);      

      if (k-1 >=0)
        {
          ierr = AOApplicationToPetsc(ao_,1,&left);
//            pointMap->insert(PointMap::value_type(LEFT,Location(i,j,k-1,left,LEFT)));
//            pointList->push_back(Location(i,j,k-1,left,LEFT));
          (*pointVector)[LEFT] = Location(i,j,k-1,left,LEFT);
        }
      if (k+1 < nxNodes)
        {
          ierr = AOApplicationToPetsc(ao_,1,&right);
//            pointMap->insert(PointMap::value_type(RIGHT,Location(i,j,k+1,right,RIGHT)));
//            pointList->push_back(Location(i,j,k+1,right,RIGHT));
          (*pointVector)[RIGHT] = Location(i,j,k+1,right,RIGHT);
        }    


    }      
  
//    i=0;j=0;
//    for (int k=-ghost_xOffSet;k < ghost_nxNodes-ghost_xOffSet; k++)
//      {
//        this->localIndex(k);
//        this->setLocalPointList();
//        if (k-1 >=-ghost_xOffSet)
//  //          pointMap->insert(PointMap::value_type(LEFT, Location(i,j,k-1,left,LEFT)));
//          pointList->push_back(Location(i,j,k-1,left,LEFT));
//        if (k+1 < ghost_nxNodes-ghost_xOffSet)
//  //          pointMap->insert(PointMap::value_type(RIGHT, Location(i,j,k+1,right,RIGHT)));
//          pointList->push_back(Location(i,j,k+1,right,RIGHT));
//  //        pointMap->insert(PointMap::value_type(CENTER, Location(i,j,k,center,CENTER)));
//        pointList->push_back(Location(i,j,k,center,CENTER));
//      }          
  (*this)(local_x0);
  globalLow = center;
  globalHigh = center + local_nxNodes;
}

SecondOrderFd2d::SecondOrderFd2d(int nxNodesIn, int nyNodesIn, int nzNodes, int dof):
  SecondOrderFd(nxNodesIn,nyNodesIn,1)
{
  using namespace cc;
  assert(nzNodes==1);
  Err ierr;
  int stencilWidth=1;
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
//    localPointMapField.resize(ghost_nxyzNodes);
//    localPointListField.resize(ghost_nxyzNodes);

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
        
//          pointMap = &pointMapField[center];
//          pointList = &pointListField[center];
        pointVector = &pointVectorField[center];
        pointVector->reserve(5);
   
        ierr = AOApplicationToPetsc(ao_,1,&center);
//          pointMap->insert(PointMap::value_type(CENTER, Location(i,j,k,center,CENTER)));
//          pointList->push_back(Location(i,j,k,center,CENTER));
        pointVector->resize(pointVector->size() + 1);
        (*pointVector)[CENTER] =Location(i,j,k,center,CENTER);

        if (k-1 >=0)
          {
            ierr = AOApplicationToPetsc(ao_,1,&left);
//              pointMap->insert(PointMap::value_type(LEFT, Location(i,j,k-1,left,LEFT)));
//              pointList->push_back(Location(i,j,k-1,left,LEFT));
            pointVector->resize(pointVector->size() + 1);
            (*pointVector)[LEFT] =Location(i,j,k-1,left,LEFT);
          }
        if (k+1 < nxNodes)
          {
            ierr = AOApplicationToPetsc(ao_,1,&right);
//              pointMap->insert(PointMap::value_type(RIGHT, Location(i,j,k+1,right,RIGHT)));
//              pointList->push_back(Location(i,j,k+1,right,RIGHT));
            pointVector->resize(pointVector->size() + 1);
            (*pointVector)[RIGHT] =Location(i,j,k+1,right,RIGHT);
          }
        if (j-1 >=0)
          {
            ierr = AOApplicationToPetsc(ao_,1,&front);
//              pointMap->insert(PointMap::value_type(FRONT, Location(i,j-1,k,front,FRONT)));
//              pointList->push_back(Location(i,j-1,k,front,FRONT));
            pointVector->resize(pointVector->size() + 1);
            (*pointVector)[FRONT] =Location(i,j-1,k,front,FRONT);
          }
        if (j+1 < nyNodes)
          {
            ierr = AOApplicationToPetsc(ao_,1,&back);
//              pointMap->insert(PointMap::value_type(BACK, Location(i,j+1,k,back,BACK)));
//              pointList->push_back(Location(i,j+1,k,back,BACK));
            pointVector->resize(pointVector->size() + 1);
            (*pointVector)[BACK] =Location(i,j+1,k,back,BACK);
          }
      }

//    i=0;
//    for (int j=-ghost_yOffSet; j<ghost_nyNodes-ghost_yOffSet; j++)
//      for (int k=-ghost_xOffSet; k<ghost_nxNodes-ghost_xOffSet; k++)
//        {
//          this->localIndex(j,k);
//          this->setLocalPointList();
//          if (k-1 >=-ghost_xOffSet)
//  //            pointMap->insert(PointMap::value_type(LEFT, Location(i,j,k-1,left,LEFT)));
//            pointList->push_back(Location(i,j,k-1,left,LEFT));
//          if (k+1 < ghost_nxNodes-ghost_xOffSet)
//  //            pointMap->insert(PointMap::value_type(RIGHT, Location(i,j,k+1,right,RIGHT)));
//            pointList->push_back(Location(i,j,k+1,right,RIGHT));
//          if (j-1 >=-ghost_yOffSet)
//  //            pointMap->insert(PointMap::value_type(FRONT, Location(i,j-1,k,front,FRONT)));
//            pointList->push_back(Location(i,j-1,k,front,FRONT));
//          if (j+1 < ghost_nyNodes-ghost_yOffSet)
//  //            pointMap->insert(PointMap::value_type(BACK, Location(i,j+1,k,back,BACK)));
//            pointList->push_back(Location(i,j+1,k,back,BACK));
//  //          pointMap->insert(PointMap::value_type(CENTER, Location(i,j,k,center,CENTER)));
//          pointList->push_back(Location(i,j,k,center,CENTER));
//        }
  (*this)(local_y0,local_x0);
  globalLow = center;
  globalHigh = center + local_nxyNodes;
}

SecondOrderFd2d::SecondOrderFd2d(SecondOrderFd2d& fineStencil, int nxNodesIn, int nyNodesIn, int nzNodes, int dof):
  SecondOrderFd(nxNodesIn,nyNodesIn,1)
{
  using namespace cc;
  assert(nzNodes==1);
  Err ierr;
  int stencilWidth=1;
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
//    localPointMapField.resize(ghost_nxyzNodes);
//    localPointListField.resize(ghost_nxyzNodes);

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
        
//          pointMap = &pointMapField[center];
//          pointList = &pointListField[center];
        pointVector = &pointVectorField[center];
        pointVector->reserve(5);

        ierr = AOApplicationToPetsc(ao_,1,&center);
//          pointMap->insert(PointMap::value_type(CENTER, Location(i,j,k,center,CENTER)));
//          pointList->push_back(Location(i,j,k,center,CENTER));
        (*pointVector)[CENTER] =Location(i,j,k,center,CENTER);

        if (k-1 >=0)
          {
            ierr = AOApplicationToPetsc(ao_,1,&left);
//              pointMap->insert(PointMap::value_type(LEFT, Location(i,j,k-1,left,LEFT)));
//              pointList->push_back(Location(i,j,k-1,left,LEFT));
            (*pointVector)[LEFT] =Location(i,j,k-1,left,LEFT);
          }
        if (k+1 < nxNodes)
          {
            ierr = AOApplicationToPetsc(ao_,1,&right);
//              pointMap->insert(PointMap::value_type(RIGHT, Location(i,j,k+1,right,RIGHT)));
//              pointList->push_back(Location(i,j,k+1,right,RIGHT));
            (*pointVector)[RIGHT] =Location(i,j,k+1,right,RIGHT);
          }
        if (j-1 >=0)
          {
            ierr = AOApplicationToPetsc(ao_,1,&front);
//              pointMap->insert(PointMap::value_type(FRONT, Location(i,j-1,k,front,FRONT)));
//              pointList->push_back(Location(i,j-1,k,front,FRONT));
            (*pointVector)[FRONT] =Location(i,j-1,k,front,FRONT);
          }
        if (j+1 < nyNodes)
          {
            ierr = AOApplicationToPetsc(ao_,1,&back);
//              pointMap->insert(PointMap::value_type(BACK, Location(i,j+1,k,back,BACK)));
//              pointList->push_back(Location(i,j+1,k,back,BACK));
            (*pointVector)[BACK] =Location(i,j+1,k,back,BACK);
          }
      }

//    i=0;
//    for (int j=-ghost_yOffSet; j<ghost_nyNodes-ghost_yOffSet; j++)
//      for (int k=-ghost_xOffSet; k<ghost_nxNodes-ghost_xOffSet; k++)
//        {
//          this->localIndex(j,k);
//          this->setLocalPointList();
//          if (k-1 >=-ghost_xOffSet)
//  //            pointMap->insert(PointMap::value_type(LEFT, Location(i,j,k-1,left,LEFT)));
//            pointList->push_back(Location(i,j,k-1,left,LEFT));
//          if (k+1 < ghost_nxNodes-ghost_xOffSet)
//  //            pointMap->insert(PointMap::value_type(RIGHT, Location(i,j,k+1,right,RIGHT)));
//            pointList->push_back(Location(i,j,k+1,right,RIGHT));
//          if (j-1 >=-ghost_yOffSet)
//  //            pointMap->insert(PointMap::value_type(FRONT, Location(i,j-1,k,front,FRONT)));
//            pointList->push_back(Location(i,j-1,k,front,FRONT));
//          if (j+1 < ghost_nyNodes-ghost_yOffSet)
//  //            pointMap->insert(PointMap::value_type(BACK, Location(i,j+1,k,back,BACK)));
//            pointList->push_back(Location(i,j+1,k,back,BACK));
//  //          pointMap->insert(PointMap::value_type(CENTER, Location(i,j,k,center,CENTER)));
//          pointList->push_back(Location(i,j,k,center,CENTER));
//        }
  (*this)(local_y0,local_x0);
  globalLow = center;
  globalHigh = center + local_nxyNodes;
}

SecondOrderFd3d::SecondOrderFd3d(int nxNodesIn, int nyNodesIn, int nzNodesIn, int dof):
  SecondOrderFd(nxNodesIn,nyNodesIn,nzNodesIn)
{ 
  using namespace cc;
 //cout<<nxNodes<<nyNodes<<nzNodes<<endl;
  Err ierr;
  int stencilWidth=1;
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
  //cout<<"da"<<endl;
  ierr =  DAGetCorners(da_, &local_x0, &local_y0, &local_z0,
                       &local_nxNodes, &local_nyNodes, &local_nzNodes);
  local_nxyNodes = local_nxNodes*local_nyNodes;
  local_nxyzNodes = local_nzNodes*local_nxyNodes;
  //cout<<"get corners"<<endl;
  ierr =  DAGetGhostCorners(da_, &ghost_x0, &ghost_y0, &ghost_z0,
                            &ghost_nxNodes, &ghost_nyNodes, &ghost_nzNodes);
  //cout<<"get ghost corners"<<endl;
  ghost_nxyNodes = ghost_nxNodes*ghost_nyNodes;
  ghost_nxyzNodes = ghost_nzNodes*ghost_nxyNodes;
//    localPointMapField.resize(ghost_nxyzNodes);
//    localPointListField.resize(ghost_nxyzNodes);

  ghost_xOffSet = local_x0 - ghost_x0;
  ghost_yOffSet = local_y0 - ghost_y0;
  ghost_zOffSet = local_z0 - ghost_z0;
  
//    cout<<"DA Information"<<endl
//        <<"local "<<local_nxNodes<<'\t'<<local_nyNodes<<'\t'<<local_nzNodes<<'\t'<<local_nxyNodes<<'\t'<<local_nxyzNodes<<endl
//        <<"ghost "<<ghost_nxNodes<<'\t'<<ghost_nyNodes<<'\t'<<ghost_nzNodes<<'\t'<<ghost_nxyNodes<<'\t'<<ghost_nxyzNodes<<endl
//        <<"offSets "<<ghost_xOffSet<<'\t'<<ghost_yOffSet<<'\t'<<ghost_zOffSet<<endl
//        <<"local corner "<<local_x0<<'\t'<<local_y0<<'\t'<<local_z0<<endl
//        <<"ghost corner "<<ghost_x0<<'\t'<<ghost_y0<<'\t'<<ghost_z0<<endl;
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
                    
//            pointMap = &pointMapField[center];
//            pointList = &pointListField[center];
          pointVector = &pointVectorField[center];
          pointVector->reserve(7);
          ierr = AOApplicationToPetsc(ao_,1,&center);
//            pointMap->insert(PointMap::value_type(CENTER, Location(i,j,k,center,CENTER)));
          (*pointVector)[CENTER] =Location(i,j,k,center,CENTER);

          if (k-1 >=0)
            {
              ierr = AOApplicationToPetsc(ao_,1,&left);
//                pointMap->insert(PointMap::value_type(LEFT, Location(i,j,k-1,left,LEFT)));
//                pointList->push_back(Location(i,j,k-1,left,LEFT));
              (*pointVector)[LEFT] =Location(i,j,k-1,left,LEFT);
            }
          if (k+1 < nxNodes)
            {
              ierr = AOApplicationToPetsc(ao_,1,&right);
//                pointMap->insert(PointMap::value_type(RIGHT, Location(i,j,k+1,right,RIGHT)));
//                pointList->push_back(Location(i,j,k+1,right,RIGHT));
              (*pointVector)[RIGHT] =Location(i,j,k+1,right,RIGHT);
            }
          if (j-1 >=0)
            {
              ierr = AOApplicationToPetsc(ao_,1,&front);
//                pointMap->insert(PointMap::value_type(FRONT, Location(i,j-1,k,front,FRONT)));
//                pointList->push_back(Location(i,j-1,k,front,FRONT));
              (*pointVector)[FRONT] =Location(i,j-1,k,front,FRONT);
           }
          if (j+1 < nyNodes)
            {
              ierr = AOApplicationToPetsc(ao_,1,&back);
//                pointMap->insert(PointMap::value_type(BACK, Location(i,j+1,k,back,BACK)));
//                pointList->push_back(Location(i,j+1,k,back,BACK));
              (*pointVector)[BACK] =Location(i,j+1,k,back,BACK);
            }
          if (i-1 >=0)
            {
              ierr = AOApplicationToPetsc(ao_,1,&bottom);
//                pointMap->insert(PointMap::value_type(BOTTOM, Location(i-1,j,k,bottom,BOTTOM)));
//                pointList->push_back(Location(i-1,j,k,bottom,BOTTOM));
              (*pointVector)[BOTTOM] =Location(i-1,j,k,bottom,BOTTOM);
            }
          if (i+1 < nzNodes)
            {
              ierr = AOApplicationToPetsc(ao_,1,&top);
//                pointMap->insert(PointMap::value_type(TOP, Location(i+1,j,k,top,TOP)));
//                pointList->push_back(Location(i+1,j,k,top,TOP));
              (*pointVector)[TOP] =Location(i+1,j,k,top,TOP);
            }
        }
  
  //cout<<"global pointMaps"<<endl;

//    for (int i=-ghost_zOffSet;i<ghost_nzNodes-ghost_zOffSet;i++)
//      for (int j=-ghost_yOffSet;j<ghost_nyNodes-ghost_yOffSet;j++)
//        for (int k=-ghost_xOffSet;k<ghost_nxNodes-ghost_xOffSet;k++)
//          {
//            this->localIndex(i,j,k);
//            this->setLocalPointList();
//            if (k-1 >=-ghost_xOffSet)
//  //              pointMap->insert(PointMap::value_type(LEFT,Location(i,j,k-1,left,LEFT)));
//              pointList->push_back(Location(i,j,k-1,left,LEFT));
//            if (k+1 < ghost_nxNodes-ghost_xOffSet)
//  //              pointMap->insert(PointMap::value_type(RIGHT,Location(i,j,k+1,right,RIGHT)));
//              pointList->push_back(Location(i,j,k+1,right,RIGHT));
//            if (j-1 >=-ghost_yOffSet)
//  //              pointMap->insert(PointMap::value_type(FRONT,Location(i,j-1,k,front,FRONT)));
//              pointList->push_back(Location(i,j-1,k,front,FRONT));
//            if (j+1 < ghost_nyNodes-ghost_yOffSet)
//  //              pointMap->insert(PointMap::value_type(BACK,Location(i,j+1,k,back,BACK)));
//              pointList->push_back(Location(i,j+1,k,back,BACK));
//            if (i-1 >=-ghost_zOffSet)
//  //              pointMap->insert(PointMap::value_type(BOTTOM,Location(i-1,j,k,bottom,BOTTOM)));
//              pointList->push_back(Location(i-1,j,k,bottom,BOTTOM));
//            if (i+1 < ghost_nzNodes-ghost_zOffSet)
//  //              pointMap->insert(PointMap::value_type(TOP,Location(i+1,j,k,top,TOP)));
//              pointList->push_back(Location(i+1,j,k,top,TOP));
//  //            pointMap->insert(PointMap::value_type(CENTER,Location(i,j,k,center,CENTER)));
//            pointList->push_back(Location(i,j,k,center,CENTER));
//          }
  //cout<<"out of stencil"<<endl;
  (*this)(local_z0,local_y0,local_x0);
  globalLow = center;
  globalHigh = center + local_nxyzNodes;
}

SecondOrderFd3d::SecondOrderFd3d(SecondOrderFd3d& fineStencil, int nxNodesIn, int nyNodesIn, int nzNodesIn, int dof):
  SecondOrderFd(nxNodesIn,nyNodesIn,nzNodesIn)
{ 
  using namespace cc;
  //cout<<nxNodes<<nyNodes<<nzNodes<<endl;
  Err ierr;
  int stencilWidth=1;
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

  //cout<<"da"<<endl;
  ierr =  DAGetCorners(da_, &local_x0, &local_y0, &local_z0,
                       &local_nxNodes, &local_nyNodes, &local_nzNodes);
  local_nxyNodes = local_nxNodes*local_nyNodes;
  local_nxyzNodes = local_nzNodes*local_nxyNodes;
  //cout<<"get corners"<<endl;
  ierr =  DAGetGhostCorners(da_, &ghost_x0, &ghost_y0, &ghost_z0,
                            &ghost_nxNodes, &ghost_nyNodes, &ghost_nzNodes);
  //cout<<"get ghost corners"<<endl;
  ghost_nxyNodes = ghost_nxNodes*ghost_nyNodes;
  ghost_nxyzNodes = ghost_nzNodes*ghost_nxyNodes;
//    localPointMapField.resize(ghost_nxyzNodes);
//    localPointListField.resize(ghost_nxyzNodes);

  ghost_xOffSet = local_x0 - ghost_x0;
  ghost_yOffSet = local_y0 - ghost_y0;
  ghost_zOffSet = local_z0 - ghost_z0;
  
//    cout<<"DA Information"<<endl
//        <<"local "<<local_nxNodes<<'\t'<<local_nyNodes<<'\t'<<local_nzNodes<<'\t'<<local_nxyNodes<<'\t'<<local_nxyzNodes<<endl
//        <<"ghost "<<ghost_nxNodes<<'\t'<<ghost_nyNodes<<'\t'<<ghost_nzNodes<<'\t'<<ghost_nxyNodes<<'\t'<<ghost_nxyzNodes<<endl
//        <<"offSets "<<ghost_xOffSet<<'\t'<<ghost_yOffSet<<'\t'<<ghost_zOffSet<<endl
//        <<"local corner "<<local_x0<<'\t'<<local_y0<<'\t'<<local_z0<<endl
//        <<"ghost corner "<<ghost_x0<<'\t'<<ghost_y0<<'\t'<<ghost_z0<<endl;
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
                    
//            pointMap = &pointMapField[center];
//            pointList = &pointListField[center];
          pointVector = &pointVectorField[center];
          pointVector->reserve(7);

          ierr = AOApplicationToPetsc(ao_,1,&center);
//            pointMap->insert(PointMap::value_type(CENTER, Location(i,j,k,center,CENTER)));
//            pointList->push_back(Location(i,j,k,center,CENTER));
          (*pointVector)[CENTER] =Location(i,j,k,center,CENTER);

          if (k-1 >=0)
            {
              ierr = AOApplicationToPetsc(ao_,1,&left);
//                pointMap->insert(PointMap::value_type(LEFT, Location(i,j,k-1,left,LEFT)));
//                pointList->push_back(Location(i,j,k-1,left,LEFT));
                (*pointVector)[LEFT] =Location(i,j,k-1,left,LEFT);
          }
          if (k+1 < nxNodes)
            {
              ierr = AOApplicationToPetsc(ao_,1,&right);
//                pointMap->insert(PointMap::value_type(RIGHT, Location(i,j,k+1,right,RIGHT)));
//                pointList->push_back(Location(i,j,k+1,right,RIGHT));
              (*pointVector)[RIGHT] =Location(i,j,k+1,right,RIGHT);
            }
          if (j-1 >=0)
            {
              ierr = AOApplicationToPetsc(ao_,1,&front);
//                pointMap->insert(PointMap::value_type(FRONT, Location(i,j-1,k,front,FRONT)));
//                pointList->push_back(Location(i,j-1,k,front,FRONT));
              (*pointVector)[FRONT] =Location(i,j-1,k,front,FRONT);
           }
          if (j+1 < nyNodes)
            {
              ierr = AOApplicationToPetsc(ao_,1,&back);
//                pointMap->insert(PointMap::value_type(BACK, Location(i,j+1,k,back,BACK)));
//                pointList->push_back(Location(i,j+1,k,back,BACK));
              (*pointVector)[BACK] =Location(i,j+1,k,back,BACK);
            }
          if (i-1 >=0)
            {
              ierr = AOApplicationToPetsc(ao_,1,&bottom);
//                pointMap->insert(PointMap::value_type(BOTTOM, Location(i-1,j,k,bottom,BOTTOM)));
//                pointList->push_back(Location(i-1,j,k,bottom,BOTTOM));
              (*pointVector)[BOTTOM] =Location(i-1,j,k,bottom,BOTTOM);
            }
          if (i+1 < nzNodes)
            {
              ierr = AOApplicationToPetsc(ao_,1,&top);
//                pointMap->insert(PointMap::value_type(TOP, Location(i+1,j,k,top,TOP)));
//                pointList->push_back(Location(i+1,j,k,top,TOP));
              (*pointVector)[TOP] =Location(i+1,j,k,top,TOP);
            }
        }
  
  //cout<<"global pointMaps"<<endl;

//    for (int i=-ghost_zOffSet;i<ghost_nzNodes-ghost_zOffSet;i++)
//      for (int j=-ghost_yOffSet;j<ghost_nyNodes-ghost_yOffSet;j++)
//        for (int k=-ghost_xOffSet;k<ghost_nxNodes-ghost_xOffSet;k++)
//          {
//            this->localIndex(i,j,k);
//            this->setLocalPointList();
//            if (k-1 >=-ghost_xOffSet)
//  //              pointMap->insert(PointMap::value_type(LEFT,Location(i,j,k-1,left,LEFT)));
//              pointList->push_back(Location(i,j,k-1,left,LEFT));
//            if (k+1 < ghost_nxNodes-ghost_xOffSet)
//  //              pointMap->insert(PointMap::value_type(RIGHT,Location(i,j,k+1,right,RIGHT)));
//              pointList->push_back(Location(i,j,k+1,right,RIGHT));
//            if (j-1 >=-ghost_yOffSet)
//  //              pointMap->insert(PointMap::value_type(FRONT,Location(i,j-1,k,front,FRONT)));
//              pointList->push_back(Location(i,j-1,k,front,FRONT));
//            if (j+1 < ghost_nyNodes-ghost_yOffSet)
//  //              pointMap->insert(PointMap::value_type(BACK,Location(i,j+1,k,back,BACK)));
//              pointList->push_back(Location(i,j+1,k,back,BACK));
//            if (i-1 >=-ghost_zOffSet)
//  //              pointMap->insert(PointMap::value_type(BOTTOM,Location(i-1,j,k,bottom,BOTTOM)));
//              pointList->push_back(Location(i-1,j,k,bottom,BOTTOM));
//            if (i+1 < ghost_nzNodes-ghost_zOffSet)
//  //              pointMap->insert(PointMap::value_type(TOP,Location(i+1,j,k,top,TOP)));
//              pointList->push_back(Location(i+1,j,k,top,TOP));
//  //            pointMap->insert(PointMap::value_type(CENTER,Location(i,j,k,center,CENTER)));
//            pointList->push_back(Location(i,j,k,center,CENTER));
//          }
  //cout<<"out of stencil"<<endl;
  (*this)(local_z0,local_y0,local_x0);
  globalLow = center;
  globalHigh = center + local_nxyzNodes;
}

}//Petsc
}//Daetk
