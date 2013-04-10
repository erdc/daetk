#include "DaetkPetscVec.h"
#include "Tracer.h"
namespace Daetk 
{
namespace Petsc
{
  namespace cc
  {
    extern "C"
    {
#ifndef DAETK_DEF_CPLUSPLUS_FOR_PETSC_H
#undef __cplusplus
#endif
#include "petsc.h"
#include "petscvec.h"
#include "petscdm.h"
#include "petscis.h"
#ifndef DAETK_DEF_CPLUSPLUS_FOR_PETSC_H
#define __cplusplus
#endif
    }
  }

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::istream;
 void Vec::copy(Vec& vIn)
   {
     using namespace cc;
     VecCopy(vIn.castToPetsc(),rep_);
   }
Vec::Registry Vec::vecRegistry;

Vec::Vec(int dim):
  ownsVec_(false),
  ownsScatter_(false),
  ownsIndexSet_(false),  
  hasGhostPoints_(false),
  firstAttachment_(true),
  firstSetStride_(true),
  isRegistered_(false),
  dim_(dim),
  base_(PetscVec_BASE),
  stride_(1),
  start_(0),
  storageLength_(0),    
  lstorageLength_(0),   
  globalHigh_(0),       
  globalLow_(0),        
  localHigh_(0),        
  localLow_(0),         
  offSet_(0),           
  ldim_(0),             
  virtualGlobalHigh_(0),
  virtualGlobalLow_(0), 
  p_(0),
  pSave_(0),
  globalIndexSet_(0),
  myGlobalIndices_(0),
  localIndexSet_(0),
  myLocalIndices_(0),
  virtualGlobalIndexSet_(0),
  myVirtualGlobalIndices_(0),
  dm_(0),
  rep_(0)
{
  Tracer tr("Petsc::Vec::Vec(int length)");
  using namespace  cc;
  mpiComm_ = PETSC_COMM_SELF;
  
  if (petSys_.isInitialized())
    {
      if (dim)
        {
          enterRegistry(dim);
          if (!ownsVec_)
            create(dim); //calls getStorageInfo and createStrideMappings
          if (storageLength_ != dim)
            cerr<<"problem in ctor"<<endl;
        }
    }
  else
    {
      if (dim!=0)
        {
          cerr<<"Construct a PetscSys object to initialize MPI before before attempting to construct vectors"<<endl;
          exit(1);
        }
    }
 
}

void Vec::getStorageInfo()
{
  using namespace cc;
  Petsc::Err ierr;
  ierr = VecGetSize(rep_,&storageLength_);
  ierr = VecGetLocalSize(rep_,&lstorageLength_);
  ierr = VecGetOwnershipRange(rep_,&globalLow_,&globalHigh_);
  assert( (globalHigh_-globalLow_) == lstorageLength_ );
  localLow_ = 0;
  localHigh_ = lstorageLength_;
}

void Vec::create(int length)
{
  using namespace cc;
  Petsc::Err ierr;
  ownsVec_ = true;
  mpiComm_ = PETSC_COMM_WORLD;
  ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,length,&rep_);
  ierr = VecSetFromOptions(rep_);
  getLocal();
  getStorageInfo();
  createStrideMappings();
}

void Vec::createWithLocal(int local_length)
{
  using namespace cc;
  Petsc::Err ierr;
  ownsVec_ = true;
  mpiComm_ = PETSC_COMM_WORLD;
  ierr = VecCreateMPI(PETSC_COMM_WORLD,local_length,PETSC_DETERMINE,&rep_);
  ierr = VecSetFromOptions(rep_);
  getLocal();
  getStorageInfo();
  dim_ = storageLength_/stride_;
  createStrideMappings();
}

Vec::Vec(int dim, const real& initialValue, int base,int stride):
  ownsVec_(false),
  ownsScatter_(false),
  ownsIndexSet_(false),  
  hasGhostPoints_(false),
  firstAttachment_(true),
  firstSetStride_(true),
  isRegistered_(false),
  dim_(dim),
  base_(base),
  stride_(stride),
  start_(0),
  storageLength_(0),    
  lstorageLength_(0),   
  globalHigh_(0),       
  globalLow_(0),        
  localHigh_(0),        
  localLow_(0),         
  offSet_(0),           
  ldim_(0),             
  virtualGlobalHigh_(0),
  virtualGlobalLow_(0),
  p_(0),
  pSave_(0),
  globalIndexSet_(0),
  myGlobalIndices_(0),
  localIndexSet_(0),
  myLocalIndices_(0),
  virtualGlobalIndexSet_(0),
  myVirtualGlobalIndices_(0),
  dm_(0),
  rep_(0)
{
  Tracer tr("Petsc::Vec::Vec(int dim, const real& initialValue, int base,int stride)");
  using namespace cc;
  mpiComm_ = PETSC_COMM_SELF; 
  assert(petSys_.isInitialized());

  int tmp=dim*stride;
  enterRegistry(dim*stride);
  if (!ownsVec_)
    create(dim*stride);
  if (storageLength_ != tmp)
    cerr<<"prob in ctor2"<<endl;
  (*this) = initialValue;
}

Vec::Vec(const Vec& v):
  ownsVec_(false),
  ownsScatter_(false),
  ownsIndexSet_(false),  
  hasGhostPoints_(false),
  firstAttachment_(true),
  firstSetStride_(true),
  isRegistered_(false),
  dim_(v.dim_),
  base_(v.base_),
  stride_(v.stride_),
  start_(v.start_),
  storageLength_(0),    
  lstorageLength_(0),   
  globalHigh_(0),       
  globalLow_(0),        
  localHigh_(0),        
  localLow_(0),         
  offSet_(0),           
  ldim_(0),             
  virtualGlobalHigh_(0),
  virtualGlobalLow_(0), 
    p_(0),
  pSave_(0),
  globalIndexSet_(0),
  myGlobalIndices_(0),
  localIndexSet_(0),
  myLocalIndices_(0),
  virtualGlobalIndexSet_(0),
  myVirtualGlobalIndices_(0),
  dm_(0),
  rep_(0)
{
  using namespace cc;
  Petsc::Err ierr;
  mpiComm_ = PETSC_COMM_SELF;

  Tracer tr("Petsc::Vec::Vec(const Petsc::Vec& v)");

  if (petSys_.isInitialized())
    {
      if (dim_)
        {
          enterRegistry(v.storageLength_);
          if (!ownsVec_)
            {
              ownsVec_ = true;
              mpiComm_ = v.mpiComm_;
              ierr = VecDuplicate(v.rep_,&rep_);
              ierr = VecCopy(v.rep_,rep_);
              getLocal();
              getStorageInfo();
              createStrideMappings();
            }
          else
            (*this) = v;
          if (storageLength_ != v.storageLength_)
            cerr<<"problem in copy ctor"<<endl;
        }
    }
  else
    if (dim_!=0)
      {
        cerr<<"Construct a PetscSys object to initialize MPI before before attempting to construct vectors"<<endl;
        exit(1);
      } 
}

Vec::Vec(Vec::CopyType t, real* a, int localDim, int base, int stride):
  ownsVec_(false),
  ownsScatter_(false),
  ownsIndexSet_(false),  
  hasGhostPoints_(false),
  firstAttachment_(true),
  firstSetStride_(true),
  isRegistered_(false),
  dim_(0),
  base_(PetscVec_BASE),
  stride_(1),
  start_(0),
  storageLength_(0),    
  lstorageLength_(0),   
  globalHigh_(0),       
  globalLow_(0),        
  localHigh_(0),        
  localLow_(0),         
  offSet_(0),           
  ldim_(0),             
  virtualGlobalHigh_(0),
  virtualGlobalLow_(0), 
    p_(0),
  pSave_(0),
  globalIndexSet_(0),
  myGlobalIndices_(0),
  localIndexSet_(0),
  myLocalIndices_(0),
  virtualGlobalIndexSet_(0),
  myVirtualGlobalIndices_(0),
  dm_(0),
  rep_(0)
{
  using namespace cc;
  mpiComm_ = PETSC_COMM_SELF;

  if (t)
    {
      attachToArray(a,localDim,base,stride);
    }
  else
    {
      clear();
      base_ = base;
      stride_ = stride;
      createWithLocal(localDim*stride);
      for (int i=0;i<ldim_;i++)
        set(i) = a[i*stride_];
    }
}


Vec::Vec(Vec::CopyType ref,const Vec& V, const VecIndex& I):
  ownsVec_(false),
  ownsScatter_(false),
  ownsIndexSet_(false), 
  hasGhostPoints_(false),
  firstAttachment_(true),
  firstSetStride_(true),
  isRegistered_(false),
  dim_(I.end() - I.start() + 1),
  base_(V.base_),
  stride_(V.stride_),
  start_(I.start() - V.base_ + V.start_),
  storageLength_(0),    
  lstorageLength_(0),   
  globalHigh_(0),       
  globalLow_(0),        
  localHigh_(0),        
  localLow_(0),         
  offSet_(0),           
  ldim_(0),             
  virtualGlobalHigh_(0),
  virtualGlobalLow_(0), 
    p_(V.p_),
  pSave_(V.pSave_),
  globalIndexSet_(0),
  myGlobalIndices_(0),
  localIndexSet_(0),
  myLocalIndices_(0),
  virtualGlobalIndexSet_(0),
  myVirtualGlobalIndices_(0),
  dm_(V.dm_),
  rep_(V.rep_)
{
  using namespace cc;
  mpiComm_ = PETSC_COMM_SELF;
 if (I.all())
    {
      dim_=V.dim_;
      start_=V.start_;
    }
  Tracer tr("Petsc::Vec::Vec(Petsc::Vec::CopyType t,const Vec& V, const Petsc::VecIndex& I)");
  assert(petSys_.isInitialized());
  
  attachToVec(ref,V,I);
}

Vec::Vec(Vec::StorageType t, cc::_p_DM* dm):
  ownsVec_(true),
  ownsScatter_(false),
  ownsIndexSet_(false),  
  firstAttachment_(true),
  firstSetStride_(true),
  isRegistered_(false),
  dim_(0),
  base_(PetscVec_BASE),
  stride_(1),
  start_(0),
  storageLength_(0),    
  lstorageLength_(0),   
  globalHigh_(0),       
  globalLow_(0),        
  localHigh_(0),        
  localLow_(0),         
  offSet_(0),           
  ldim_(0),             
  virtualGlobalHigh_(0),
  virtualGlobalLow_(0), 
    p_(0),
  pSave_(0),
  globalIndexSet_(0),
  myGlobalIndices_(0),
  localIndexSet_(0),
  myLocalIndices_(0),
  virtualGlobalIndexSet_(0),
  myVirtualGlobalIndices_(0),
  dm_(dm),
  rep_(0)
{

  Tracer tr("Petsc::Vec::StorageType t, Petsc::cc:_p_DM* dm)");

  using namespace cc;
  Petsc::Err ierr;
  mpiComm_=PETSC_COMM_SELF;

  assert(petSys_.isInitialized());
  
  if (t == LOCAL)
    {
      hasGhostPoints_ = true;
      mpiComm_=PETSC_COMM_SELF;
      ierr =  DMCreateLocalVector(dm,&rep_); 
    }
  else
    {
      mpiComm_=PETSC_COMM_WORLD;
      hasGhostPoints_ = false;
      ierr =  DMCreateGlobalVector(dm,&rep_);
    }
  ierr = VecSetFromOptions(rep_);
  getLocal();
  getStorageInfo();
  dim_ = storageLength_;
  createStrideMappings();
}

Vec::Vec(const CMRVec<real>& vin):
  ownsVec_(false),
  ownsScatter_(false),
  ownsIndexSet_(false),  
  hasGhostPoints_(false),
  firstAttachment_(true),
  firstSetStride_(true),
  isRegistered_(false),
  dim_(vin.dim()),
  base_(PetscVec_BASE),
  stride_(1),
  start_(0),
  storageLength_(0),    
  lstorageLength_(0),   
  globalHigh_(0),       
  globalLow_(0),        
  localHigh_(0),        
  localLow_(0),         
  offSet_(0),           
  ldim_(0),             
  virtualGlobalHigh_(0),
  virtualGlobalLow_(0), 
    p_(0),
  pSave_(0),
  globalIndexSet_(0),
  myGlobalIndices_(0),
  localIndexSet_(0),
  myLocalIndices_(0),
  virtualGlobalIndexSet_(0),
  myVirtualGlobalIndices_(0),
  dm_(0),
  rep_(0)
{
  using namespace cc;
  Petsc::Err ierr;
  mpiComm_=PETSC_COMM_SELF;
  Tracer tr("Petsc::Vec::Vec(int length)");
  
  
  if (petSys_.isInitialized())
    {
      int dim=vin.dim();
      ierr = MPI_Bcast(&dim,1,MPI_INT,0,cc::PETSC_COMM_WORLD);
      enterRegistry(dim);
      if (!ownsVec_)
        create(dim); //calls getStorageInfo and createStrideMappings
      cc::Vec gVec;
      ierr =  DMDACreateNaturalVector(dm_,&gVec);
      restoreLocal();
      
      if (petSys_.master())
        {
          int* ix = new int[dim];
          for (int i=0;i<dim;i++)
            ix[i]=i;
          ierr = VecSetValues(gVec,dim,ix,vin.castToConstArray(),INSERT_VALUES);
          ierr = VecAssemblyBegin(gVec);
          ierr = VecAssemblyEnd(gVec);
          delete [] ix;      
        }
      else
        {
          ierr =  VecAssemblyBegin(gVec);
          ierr =  VecAssemblyEnd(gVec);
        }
      ierr =  DMDANaturalToGlobalBegin(dm_,gVec,INSERT_VALUES,rep_);
      ierr =  DMDANaturalToGlobalEnd(dm_,gVec,INSERT_VALUES,rep_);     
      getLocal();
      VecDestroy(&gVec);
      if (storageLength_ != dim)
        cerr<<"trouble in ctor4"<<endl;
    }
  else
    {
      if (vin.dim()!=0)
        {
          cerr<<"Construct a PetscSys object to initialize MPI before before attempting to construct vectors"<<endl;
          exit(1);
        }
    }
}

Vec::~Vec()
{
  clear();
}



Vec& Vec::operator=(const Vec& V)
{
#ifdef PETSCVEC_BOUNDS_CHECK
  assert(checkConformance(V));
#endif
  int i=0;
  for (i=0;i<ldim_;i++)
    set(i) = V.get(i);
  return *this;
}

Vec& Vec::operator=(const real& a)
{
  int i=0;
  for (i=0;i<ldim_;i++)
    set(i) = a;
  return *this;
}

Vec& Vec::operator*=(const real &a)
{
  int i=0;
  for (i=0;i<ldim_;i++)
    set(i)*=a;
  return *this;
}

Vec& Vec::operator/=(const real &a)
{
  int i=0;
  for (i=0;i<ldim_;i++)
    set(i)/=a;
  return *this;
}

Vec& Vec::operator+=(const Vec& y)
{
#ifdef PETSCVEC_BOUNDS_CHECK
  assert(checkConformance(y));
#endif
  int i=0;
  for (i=0;i<ldim_;i++)
    set(i)+=y.get(i);
  return *this;
}

Vec& Vec::operator-=(const Vec& y)
{
#ifdef PETSCVEC_BOUNDS_CHECK
  assert(checkConformance(y));
#endif
  int i=0;
  for (i=0;i<ldim_;i++)
    set(i)-=y.get(i);
  return *this;
}


Vec& Vec::newsize( Vec::StorageType t, cc::_p_DM* dm)
{
  using namespace cc;
  Petsc::Err ierr;

  Tracer tr("Petsc::Vec::newsize(Petsc::Vec::StorageType t, Petsc::cc:_p_DM* dm)");
  assert(petSys_.isInitialized());
  
  if (t == LOCAL)
    {
      clear();
      mpiComm_=PETSC_COMM_SELF;
      hasGhostPoints_ = true;
      ierr =  DMCreateLocalVector(dm,&rep_); 
      ownsVec_ = true;
    }
  else
    {
      cc::_p_Vec* newRep_;
      ierr =  DMCreateGlobalVector(dm,&newRep_);
      int newStorageLength_;
      ierr = VecGetSize(newRep_,&newStorageLength_);
      if (storageLength_ == newStorageLength_)
        {
          restoreLocal();
          //copy old vector into new vector layout
          VecScatter ctx;
          IS is;
          ierr = ISCreateStride(mpiComm_,storageLength_,0,1,&is);
          ierr = VecScatterCreate(rep_,is,newRep_,is,&ctx);
//           ierr = VecScatterBegin(rep_,newRep_,INSERT_VALUES ,SCATTER_FORWARD, ctx);
          ierr = VecScatterBegin(ctx,rep_,newRep_,INSERT_VALUES ,SCATTER_FORWARD);
//           ierr = VecScatterEnd(rep_,newRep_,INSERT_VALUES ,SCATTER_FORWARD, ctx);
          ierr = VecScatterEnd(ctx,rep_,newRep_,INSERT_VALUES ,SCATTER_FORWARD);
          ierr = VecScatterDestroy(&ctx);
          ierr = ISDestroy(&is);
          getLocal();
        }
      clear();
      mpiComm_=PETSC_COMM_WORLD;
      rep_ = newRep_;
      ownsVec_ = true;
    }
  dm_=dm;
  //ierr = VecSetFromOptions(rep_);
  getLocal();
  getStorageInfo();
  dim_ = storageLength_;
  createStrideMappings();
  return *this;
}

Vec& Vec::newsize(int n)
{
  Tracer tr("Petsc::Vec::newsize( int n)");

  clear();
  enterRegistry(n);
  if (!ownsVec_)
    {
      dim_=n;
      stride_=1;
      create(n);
    }
  if (storageLength_ != n)
    cerr<<"st ctor5"<<endl;
  return *this;
}

Vec& Vec::newsizeSerial(int n)
{
  using namespace cc;
  Petsc::Err ierr;
  Tracer tr("Petsc::Vec::newsize( int n)");

  clear();
  dim_=n;
  stride_=1;
  ownsVec_ = true;
  mpiComm_ = PETSC_COMM_SELF;
  ierr = VecCreateSeq(PETSC_COMM_SELF,n,&rep_);
  ierr = VecSetFromOptions(rep_);
  getLocal();
  getStorageInfo();
  createStrideMappings();  
  if (storageLength_ != n)
    cerr<<"st ctor5"<<endl;
  return *this;
}

Vec& Vec::newsize(Vec::StorageType t,int n)
{
  //confusing, the LOCAL vector here would be global
  Tracer tr("Vec::newsize( int n)");

  clear();
  if (t == GLOBAL)
    {
      enterRegistry(n);
      if (!ownsVec_)
        {
          dim_=n;
          stride_=1;
          create(n);
        }
      if (storageLength_ != n)
        cerr<<"probleme"<<endl;
    }
  else
    {
      ldim_=n;
      stride_=1;
      createWithLocal(n);
    }
  
  return *this;
}

Vec& Vec::redim( int n)
{
  newsize(n);
  return *this;
} 

Vec& Vec::setStride( int newStride)
{
  int virtualStorageLength = (dim_*stride_);
  dim_ = virtualStorageLength / newStride + ((virtualStorageLength%newStride > 0) ? 1 : 0);
  stride_ = newStride;
  int storBefore=storageLength_;
  getStorageInfo();
  if (storBefore!=storageLength_)
    {
      cerr<<"storage length changed"<<endl;
      if (isRegistered_)
        cerr<<"vec is registered as length "<<storBefore<<endl;
    }
  createStrideMappings();
  return *this;
}

Vec& Vec::setStrideMulti( int newStride)
{
  if (firstSetStride_)
    {
      int virtualStorageLength = (dim_*stride_);
      dim_ = virtualStorageLength / newStride + ((virtualStorageLength%newStride > 0) ? 1 : 0);
      stride_ = newStride;
      int storBefore=storageLength_;
      getStorageInfo();
      if (storBefore!=storageLength_)
        {
          cerr<<"storage length changed"<<endl;
          if (isRegistered_)
            cerr<<"vec is registered as length "<<storBefore<<endl;
        }
      createStrideMappings();
      firstSetStride_ = false;
    }
  //else do nothing
  return *this;
}

Vec& Vec::setBlockStrideMulti( int newStride, int newBlock)
{
  if (firstSetStride_)
    {
      int virtualStorageLength = (dim_*stride_);
      dim_ = (virtualStorageLength / newStride)*newBlock + ((virtualStorageLength%newStride > 0) ? 1 : 0);
      stride_ = newStride;
      int storBefore=storageLength_;
      getStorageInfo();
      if (storBefore!=storageLength_)
        {
          cerr<<"storage length changed"<<endl;
          if (isRegistered_)
            cerr<<"vec is registered as length "<<storBefore<<endl;
        }
      createBlockStrideMappings(newBlock);
      firstSetStride_ = false;
    }
  //else do nothing
  return *this;
}

Vec& Vec::setBase(int b)
{  
  base_=b;
  createStrideMappings();
  return *this;
}

// bool Vec::isValid() const
// {
//   cc::PetscBool flg;
//   Petsc::Err ierr;
//   ierr = Petsc::VecValid(rep_,&flg);
//   return flg;
// }



real Vec::sum()
{
  real m;
  Petsc::Err ierr;
  restoreLocal();
  ierr =  cc::VecSum(rep_,&m);
  getLocal();
  return m;
}

void Vec::restoreLocal() const
{
  Petsc::Err ierr;
  if (p_ && rep_)
    {
      pSave_ = p_;
      //mwf gcc 3.3 has problem parsing cc:: here because VecRestoreArray
      // is a macro?
      //ierr = cc::VecRestoreArray(rep_,&p_);
      using namespace cc;
      ierr = VecRestoreArray(rep_,&p_);
      p_=0;
      hasLocal_=false;
    }
}

void Vec::getLocal() const
{
  Petsc::Err ierr;
  if (!p_ && rep_)
    {
      //mwf gcc 3.3 has problem parsing cc:: here because VecRestoreArray
      // is a macro?
      //ierr = cc::VecGetArray(rep_,&p_);
      using namespace cc;
      ierr = VecGetArray(rep_,&p_);
      if (pSave_ &&  pSave_ != p_)
        {
          cerr<<"Petsc has changed the location of the local array."<<endl
              <<"References to this vector are now out of date"<<endl;
          pSave_  = p_;
        }
      else
        pSave_ = p_;
      hasLocal_=true;
    }
}

void Vec::startAssembly()
{
  Petsc::Err ierr;
  restoreLocal();
  ierr =  cc::VecAssemblyBegin(rep_);
}

void Vec::endAssembly()
{
  Petsc::Err ierr;
  ierr =  cc::VecAssemblyEnd(rep_);
  getLocal();
}

void Vec::startSetFromGlobalMulti(const Vec& global)
{
  using namespace cc;
  Petsc::Err ierr;
  restoreLocal();
  global.restoreLocal();
  if (!ownsScatter_)
    {
      VecScatter dm_ltog_,dm_gtol_,dm_ltol_;
      ierr =  DMDAGetScatter(dm_,&dm_ltog_,&dm_gtol_,&dm_ltol_);
      ierr =  VecScatterCopy(dm_gtol_,&my_gtol_);
      ownsScatter_=true;
    }
//   ierr =  VecScatterBegin(global.rep_,rep_,
//                           INSERT_VALUES,SCATTER_FORWARD,
//                           my_gtol_);
  ierr =  VecScatterBegin(my_gtol_,global.rep_,rep_,
                          INSERT_VALUES,SCATTER_FORWARD);
  //  ierr = DMGlobalToLocalBegin(dm_,global.rep_,INSERT_VALUES,rep_);
//    global.getLocal();
//    getLocal();
}

void Vec::endSetFromGlobalMulti(const Vec& global)
{
  Petsc::Err ierr;
//    restoreLocal();
//    global.restoreLocal();
  //ierr = DMGlobalToLocalEnd(dm_,global.rep_,INSERT_VALUES,rep_);
  using namespace cc;
//   ierr =  VecScatterEnd(global.rep_,rep_,
//                         INSERT_VALUES,SCATTER_FORWARD,
//                         my_gtol_);
  ierr =  VecScatterEnd(my_gtol_,
                        global.rep_,rep_,
                        INSERT_VALUES,SCATTER_FORWARD);
  
  global.getLocal();
  getLocal();
}
  
void Vec::startSetFromGlobal(const Vec& global)
{
  Petsc::Err ierr;
  restoreLocal();
  global.restoreLocal();
  using namespace cc;
  ierr = DMGlobalToLocalBegin(dm_,global.rep_,INSERT_VALUES,rep_);
}

void Vec::endSetFromGlobal(const Vec& global)
{
  using namespace cc;
  Petsc::Err ierr;
  ierr = DMGlobalToLocalEnd(dm_,global.rep_,INSERT_VALUES,rep_);
  global.getLocal();
  getLocal();
}
  
void Vec::setFromLocal(const Vec& local)
{
  using namespace cc;
  Petsc::Err ierr;
  restoreLocal();
  local.restoreLocal();
  ierr = DMLocalToGlobalBegin(dm_,local.rep_,INSERT_VALUES,rep_);
  ierr = DMLocalToGlobalEnd(dm_,local.rep_,INSERT_VALUES,rep_);
  local.getLocal();
  getLocal();
}
  
Vec& Vec::clear(bool deRegister)
{
  Petsc::Err ierr;
  //take this vec out of the registry
  
  if (isRegistered_ && deRegister)
    {
      std::set<Vec*>::iterator it=vecRegistry[storageLength_].second.find(this);
      if (this == vecRegistry[storageLength_].first)
        vecRegistry[storageLength_].first = (Vec*)(0);
      else if (it!=vecRegistry[storageLength_].second.end())
        vecRegistry[storageLength_].second.erase(this);
      else
        {
          //error. print the registry to help find problem
          cerr<<"couldn't find vector in registry "<<storageLength_<<'\t'<<this<<endl;
          cerr<<"printing registry"<<endl;
          Vec::Registry::iterator ip=Vec::vecRegistry.begin();
          while (ip != Vec::vecRegistry.end())
            {
              cerr<<"first "<<ip->first<<endl;
              cerr<<"first a"<<ip->second.first<<endl;
              std::set<Vec*>::iterator i=ip->second.second.begin();
              while (i !=ip->second.second.end())
                {
                  cerr<<*i<<endl;
                  ++i;
                }
              ++ip;
            }
          
          assert(0);
        }
      isRegistered_=false;
    }
  using namespace cc;
  if (ownsVec_)
    {
      restoreLocal();
      ierr = VecDestroy(&rep_);
      ownsVec_ = false;
    }
  if (ownsScatter_)
    {
      VecScatterDestroy(&my_gtol_);
      ownsScatter_ = false;
    }
  if (ownsIndexSet_)
    {
      ierr = ISRestoreIndices(globalIndexSet_,&myGlobalIndices_);
      ierr = ISDestroy(&globalIndexSet_);
      ierr = ISRestoreIndices(localIndexSet_,&myLocalIndices_);
      ierr = ISDestroy(&localIndexSet_);
      ierr = ISRestoreIndices(virtualGlobalIndexSet_,&myVirtualGlobalIndices_);
      ierr = ISDestroy(&virtualGlobalIndexSet_);
      ownsIndexSet_=false;
    }
  hasGhostPoints_=false;
  firstAttachment_ = true;
  firstSetStride_ = true;
  dim_=0;
  base_=0;
  stride_=1;
  start_=0;
  storageLength_=0;    
  lstorageLength_=0;   
  globalHigh_=0;       
  globalLow_=0;        
  localHigh_=0;        
  localLow_=0;         
  offSet_=0;           
  ldim_=0;             
  virtualGlobalHigh_=0;
  virtualGlobalLow_=0; 
  mpiComm_=PETSC_COMM_SELF;
  pSave_=0;
  return (*this);
}

Vec& Vec::attachToVec(Vec::CopyType ref,const Vec& V,const VecIndex& I)
{
  using namespace  cc;
  Petsc::Err ierr;
  clear();
  if (!V.p_)
    V.getLocal();
  if (I.all())
    {
      dim_=V.dim_;
      start_=V.start_;
    }
  else
    {
      dim_=I.end() - I.start() + 1;
      start_=I.start() - V.base_ + V.start_;
    }
  base_=V.base_;
  stride_=V.stride_;
  p_=V.p_;
  pSave_ = V.pSave_;
  rep_=V.rep_;
  mpiComm_=V.mpiComm_;

  //if this is just a ref this is all the information we need
  
  getStorageInfo();
  createStrideMappings();
  
  if (!ref) 
    {
      if (V.hasGhostPoints_)
        {
          std::cerr<<"attachToVec not implemented for non-reference copy of a local-ghosted vector, exiting"<<std::endl;
          exit(1);
        }
      else
        {
          stride_ = 1;
          start_ = 0;
          createWithLocal(ldim_);
          restoreLocal();
          V.restoreLocal();
          VecScatter ctx;
          IS ix,iy;
          ierr = ISCreateStride(mpiComm_,dim_,start_,V.stride_,&ix);
          ierr = ISCreateStride(mpiComm_,dim_,start_,stride_,&iy);
          ierr = VecScatterCreate(V.rep_,ix,rep_,iy,&ctx);
//           ierr = VecScatterBegin(V.rep_,rep_,INSERT_VALUES ,SCATTER_FORWARD, ctx);
//           ierr = VecScatterEnd(V.rep_,rep_,INSERT_VALUES ,SCATTER_FORWARD, ctx);
          ierr = VecScatterBegin(ctx,V.rep_,rep_,INSERT_VALUES ,SCATTER_FORWARD);
          ierr = VecScatterEnd(ctx,V.rep_,rep_,INSERT_VALUES ,SCATTER_FORWARD);
          ierr = VecScatterDestroy(&ctx);
          ierr = ISDestroy(&ix);
          ierr = ISDestroy(&iy);
          V.getLocal();
          getLocal();
        }
    }
  else
    {
      if (!V.ownsVec_)
        cerr<<"warning: Referencing a reference. Behavior is undefined."<<endl; 
    }
  //  firstAttachment_ = false;
  firstAttachment_=true;
  return *this;
}

Vec& Vec::attachToVecMulti(Vec::CopyType ref,const Vec& V,const VecIndex& I)
{
  using namespace cc;
  Petsc::Err ierr;
  if (firstAttachment_ || !ref)
    {
      clear();
      if (!V.p_)
        V.getLocal();
      if (I.all())
        {
          dim_ = V.dim_;
          start_ = V.start_;
        }
      else
        {
          dim_=I.end() - I.start() + 1;
          start_=I.start() - V.base_ + V.start_;
        }
      base_=V.base_;
      stride_=V.stride_;
      p_=V.p_;
      pSave_ = V.pSave_;
      rep_=V.rep_;
      mpiComm_ = V.mpiComm_;
       //if this is just a ref this is all the information we need
      
      getStorageInfo();
      createStrideMappings();
      
      if (!ref) 
        {
          if (V.hasGhostPoints_)
            {
              std::cerr<<"attachToVec not implemented for non-reference copy of a local-ghosted vector, exiting"<<std::endl;
              exit(1);
            }
          else
            {
              //shrink to stride 1 to save space
              stride_ = 1;
              start_ = 0;
              createWithLocal(ldim_);
              restoreLocal();
              V.restoreLocal();
              VecScatter ctx;
              IS ix,iy;
              ierr = ISCreateStride(mpiComm_,dim_,start_,V.stride_,&ix);
              ierr = ISCreateStride(mpiComm_,dim_,start_,stride_,&iy);
              ierr = VecScatterCreate(V.rep_,ix,rep_,iy,&ctx);
//               ierr = VecScatterBegin(V.rep_,rep_,INSERT_VALUES ,SCATTER_FORWARD, ctx);
//               ierr = VecScatterEnd(V.rep_,rep_,INSERT_VALUES ,SCATTER_FORWARD, ctx);
              ierr = VecScatterBegin(ctx,V.rep_,rep_,INSERT_VALUES ,SCATTER_FORWARD);
              ierr = VecScatterEnd(ctx,V.rep_,rep_,INSERT_VALUES ,SCATTER_FORWARD);
              ierr = VecScatterDestroy(&ctx);
              ierr = ISDestroy(&ix);
              ierr = ISDestroy(&iy);
              V.getLocal();
              getLocal();
            }
        }
      else
        {
          if (!V.ownsVec_)
            cerr<<"warning: Referencing a reference. Behavior is undefined."<<endl; 
        }
      firstAttachment_ = false;
    }
  else
    {
      p_=V.p_;
      pSave_ = V.pSave_;
      rep_=V.rep_;
      
    }
  return *this;
}

Vec& Vec::attachToArray(real *a, int localDim,int base, int stride)
{
  using namespace cc;
  Petsc::Err ierr;
  clear();
  base_=base;
  start_ = 0;
  stride_=stride;
  mpiComm_ = PETSC_COMM_WORLD;
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,stride_,localDim*stride,PETSC_DETERMINE,a,&rep_);
  ownsVec_ = true; //owns the Petsc Vec, not the array
  getLocal();
  assert(a==p_);
  getStorageInfo();
  dim_ = storageLength_ / stride_;
  createStrideMappings();
  return *this;
}

Vec& Vec::attachToPetscRepMulti(cc::_p_Vec* pv)
{
  using namespace cc;
  if (firstAttachment_)
    {
      clear();
      base_=0;
      start_=0;
      stride_=1;
      ownsVec_=false;
      rep_ = pv;
      mpiComm_=PETSC_COMM_WORLD;
      getLocal();
      getStorageInfo();
      dim_=storageLength_;
      createStrideMappings();
      firstAttachment_=false;
    }
  else
    {
      rep_=pv;
      getLocal();
    }
  return *this;
}

Vec& Vec::detachFromPetscRepMulti()
{
  restoreLocal();
  rep_=0;
  p_=0;
  pSave_=0;
  return *this;
}

void Vec::createStrideMappings()
{  
  using namespace cc;
  Petsc::Err ierr;
  //create the index sets needed by a strided vector or vector slice
  //this function sets
  //   ldim_
  //   offSet_
  //   globalLow_
  //   globalHigh_
  //   localLow_
  //   localHigh_
  //   virtualGlobalLow_
  //   virtualGlobalHigh_
  //   globalIndexSet_        
  //   localIndexSet_          
  //   virtualGlobalIndexSet_
  //   myGlobalIndices_
  //   myLocalIndices_
  //   myVirtualGlobalIndices_

  if (ownsIndexSet_)
    {
      ierr = ISRestoreIndices(globalIndexSet_,&myGlobalIndices_);
      ierr = ISDestroy(&globalIndexSet_);
      ierr = ISRestoreIndices(localIndexSet_,&myLocalIndices_);
      ierr = ISDestroy(&localIndexSet_);
      ierr = ISRestoreIndices(virtualGlobalIndexSet_,&myVirtualGlobalIndices_);
      ierr = ISDestroy(&virtualGlobalIndexSet_);
      ownsIndexSet_=false;
    }
  else
    {
      myGlobalIndices_=0;
      globalIndexSet_=0;
      myLocalIndices_=0;
      localIndexSet_=0;
      myVirtualGlobalIndices_=0;
      virtualGlobalIndexSet_=0;
    }

  if (start_ >= globalHigh_)
    {
      ldim_ = 0;
      offSet_ = 0;
      globalLow_ = -1;
      globalHigh_ = -1;
      localLow_ = -1;
      localHigh_ = -1;
      virtualGlobalLow_ = -1;
      virtualGlobalHigh_ = -1;
    }
  else if (start_ < globalLow_)
    {
      offSet_ = (globalLow_ - start_ ) % stride_;
      globalLow_ += offSet_;
      ldim_=(lstorageLength_ - offSet_)/stride_;
      if ((lstorageLength_ - offSet_) % stride_ > 0)
        ldim_+=1;
      globalHigh_ = std::min(start_ + dim_*stride_,globalLow_ + ldim_*stride_);
     
      //reset ldim_ to correct value
      ldim_ = (globalHigh_ - globalLow_) / stride_;

      localLow_ = offSet_;
      localHigh_ = localLow_ + ldim_*stride_;
      virtualGlobalLow_ = (globalLow_ - start_)/stride_ + base_;
      virtualGlobalHigh_ = virtualGlobalLow_ + ldim_;
    }
  else
    {
      offSet_ = start_ - globalLow_;
      globalLow_ = start_;
      ldim_ = (lstorageLength_ - offSet_) / stride_;
      if ((lstorageLength_ - offSet_) % stride_ > 0)
        ldim_+=1;
      globalHigh_ = std::min(start_ + dim_*stride_,globalLow_ + ldim_*stride_);

      //reset ldim_ to correct value
      ldim_ = (globalHigh_ - globalLow_) / stride_;

      localLow_ = offSet_;
      localHigh_ = localLow_ + ldim_*stride_;
      virtualGlobalLow_ = 0 + base_;
      virtualGlobalHigh_ = virtualGlobalLow_ + ldim_;
    }
  if (ldim_ != 0)
    {
      ownsIndexSet_ = true;      
      ierr = ISCreateStride(PETSC_COMM_SELF,ldim_,globalLow_,stride_,&globalIndexSet_);
      ierr = ISCreateStride(PETSC_COMM_SELF,ldim_,localLow_,stride_,&localIndexSet_);
      ierr = ISCreateStride(PETSC_COMM_SELF,ldim_,virtualGlobalLow_,1,&virtualGlobalIndexSet_);

      int ldimtest;

      //check global index set
      ierr = ISGetLocalSize(globalIndexSet_,&ldimtest);
      assert(ldimtest == ldim_);

      ierr = ISGetIndices(globalIndexSet_,&myGlobalIndices_);
      
      assert(globalLow_  == myGlobalIndices_[0]);
      assert(globalHigh_  == myGlobalIndices_[ldim_-1] + stride_);

      //check local index set
      ierr = ISGetLocalSize(localIndexSet_,&ldimtest);
      assert(ldimtest == ldim_);

      ierr = ISGetIndices(localIndexSet_,&myLocalIndices_);
  
      assert(localLow_ == myLocalIndices_[0]);
      assert(localHigh_ == myLocalIndices_[ldim_-1] + stride_);

      ierr = ISGetIndices(virtualGlobalIndexSet_,&myVirtualGlobalIndices_);
  
      assert(virtualGlobalLow_ == myVirtualGlobalIndices_[0]);
      assert(virtualGlobalHigh_ == myVirtualGlobalIndices_[ldim_-1] + 1);

    }
  else
    {
      globalIndexSet_=0;        
      localIndexSet_ =0;          
      virtualGlobalIndexSet_=0; 
      myGlobalIndices_=0; 
      myLocalIndices_=0; 
      myVirtualGlobalIndices_=0; 
    }
}

void Vec::createBlockStrideMappings(int block)
{  
  using namespace cc;
  Petsc::Err ierr;
  //create the index sets needed by a strided vector or vector slice
  //this function sets
  //   ldim_
  //   offSet_
  //   globalLow_
  //   globalHigh_
  //   localLow_
  //   localHigh_
  //   virtualGlobalLow_
  //   virtualGlobalHigh_
  //   globalIndexSet_        
  //   localIndexSet_          
  //   virtualGlobalIndexSet_
  //   myGlobalIndices_
  //   myLocalIndices_
  //   myVirtualGlobalIndices_

  if (ownsIndexSet_)
    {
      ierr = ISRestoreIndices(globalIndexSet_,&myGlobalIndices_);
      ierr = ISDestroy(&globalIndexSet_);
      ierr = ISRestoreIndices(localIndexSet_,&myLocalIndices_);
      ierr = ISDestroy(&localIndexSet_);
      ierr = ISRestoreIndices(virtualGlobalIndexSet_,&myVirtualGlobalIndices_);
      ierr = ISDestroy(&virtualGlobalIndexSet_);
      ownsIndexSet_=false;
    }
  else
    {
      myGlobalIndices_=0;
      globalIndexSet_=0;
      myLocalIndices_=0;
      localIndexSet_=0;
      myVirtualGlobalIndices_=0;
      virtualGlobalIndexSet_=0;
    }

  if (start_ >= globalHigh_)
    {
      ldim_ = 0;
      offSet_ = 0;
      globalLow_ = -1;
      globalHigh_ = -1;
      localLow_ = -1;
      localHigh_ = -1;
      virtualGlobalLow_ = -1;
      virtualGlobalHigh_ = -1;
    }
  else if (start_ < globalLow_)
    {
      offSet_ = (globalLow_ - start_ ) % stride_;
      globalLow_ += offSet_;
      ldim_=((lstorageLength_ - offSet_)/stride_)*block;
      if ((lstorageLength_ - offSet_) % stride_ > 0)
        ldim_+=block;
      globalHigh_ = std::min(start_ + (dim_/block)*stride_+block-1,globalLow_ + (ldim_/block)*stride_+block-1);
     
      //reset ldim_ to correct value
      ldim_ = ((globalHigh_- block + 1 - globalLow_) / stride_)*block;

      localLow_ = offSet_;
      localHigh_ = localLow_ + (ldim_/block)*stride_+block-1;
      virtualGlobalLow_ = (globalLow_ - start_)/stride_ + base_;
      virtualGlobalHigh_ = virtualGlobalLow_ + ldim_;
    }
  else
    {
      offSet_ = start_ - globalLow_;
      globalLow_ = start_;
      ldim_ = ((lstorageLength_ - offSet_) / stride_)*block;
      if ((lstorageLength_ - offSet_) % stride_ > 0)
        ldim_+=block;
      globalHigh_ = std::min(start_ + (dim_/block)*stride_+block-1,globalLow_ + (ldim_/block)*stride_+block-1);

      //reset ldim_ to correct value
      ldim_ = ((globalHigh_ - block + 1 - globalLow_) / stride_)*block;

      localLow_ = offSet_;
      localHigh_ = localLow_ + (ldim_/block)*stride_+block-1;
      virtualGlobalLow_ = 0 + base_;
      virtualGlobalHigh_ = virtualGlobalLow_ + ldim_;
    }
  if (ldim_ != 0)
    {
      //first generate some index sets for the stride sets without blocking
      cc::_p_IS *globalStrideIndexSet_, *localStrideIndexSet_, *virtualGlobalStrideIndexSet_;
      const int *global_idx,*local_idx,*virtualGlobal_idx;
      
      ierr = ISCreateStride(PETSC_COMM_SELF,ldim_,globalLow_,stride_,&globalStrideIndexSet_);
      ierr = ISCreateStride(PETSC_COMM_SELF,ldim_,localLow_,stride_,&localStrideIndexSet_);
      ierr = ISCreateStride(PETSC_COMM_SELF,ldim_,virtualGlobalLow_,block,&virtualGlobalStrideIndexSet_);

      //get the integer index sets and pass them to the block index set ctor

      ierr = ISGetIndices(globalStrideIndexSet_,&global_idx);
      ierr = ISGetIndices(localStrideIndexSet_,&local_idx);
      ierr = ISGetIndices(virtualGlobalStrideIndexSet_,&virtualGlobal_idx);



      ownsIndexSet_ = true;      

      ierr = ISCreateBlock(PETSC_COMM_SELF,block,ldim_/block,global_idx,PETSC_COPY_VALUES,&globalIndexSet_);
      ierr = ISCreateBlock(PETSC_COMM_SELF,block,ldim_/block,local_idx,PETSC_COPY_VALUES,&localIndexSet_);
      ierr = ISCreateBlock(PETSC_COMM_SELF,block,ldim_/block,virtualGlobal_idx,PETSC_COPY_VALUES,&virtualGlobalIndexSet_);

      //get rid of the stride sets

      ierr = ISRestoreIndices(globalStrideIndexSet_,&global_idx);
      ierr = ISRestoreIndices(localStrideIndexSet_,&local_idx);
      ierr = ISRestoreIndices(virtualGlobalStrideIndexSet_,&virtualGlobal_idx);
      ierr = ISDestroy(&globalStrideIndexSet_);
      ierr = ISDestroy(&localStrideIndexSet_);
      ierr = ISDestroy(&virtualGlobalStrideIndexSet_);

      int ldimtest;

      //check global index set
      ierr = ISGetLocalSize(globalIndexSet_,&ldimtest);
      assert(ldimtest == ldim_);

      ierr = ISGetIndices(globalIndexSet_,&myGlobalIndices_);
      
      assert(globalLow_  == myGlobalIndices_[0]);
      assert(globalHigh_  == myGlobalIndices_[ldim_-1] + stride_);

      //check local index set
      ierr = ISGetLocalSize(localIndexSet_,&ldimtest);
      assert(ldimtest == ldim_);

      ierr = ISGetIndices(localIndexSet_,&myLocalIndices_);
  
      assert(localLow_ == myLocalIndices_[0]);
      assert(localHigh_ == myLocalIndices_[ldim_-1] + stride_);

      const Daetk::Petsc::cc::PetscInt* a;
      ierr = ISGetIndices(virtualGlobalIndexSet_,&a);
      ierr = ISGetIndices(virtualGlobalIndexSet_,&myVirtualGlobalIndices_);
  
      assert(virtualGlobalLow_ == myVirtualGlobalIndices_[0]);
      assert(virtualGlobalHigh_ == myVirtualGlobalIndices_[ldim_-1] + 1);

    }
  else
    {
      globalIndexSet_=0;        
      localIndexSet_ =0;          
      virtualGlobalIndexSet_=0; 
      myGlobalIndices_=0; 
      myLocalIndices_=0; 
      myVirtualGlobalIndices_=0; 
    }
}


real* Vec::castToArray()
{
  return &p_[localLow_];
}

const real* Vec::castToConstArray() const
{
  return &p_[localLow_];
}
  
cc::_p_Vec* Vec::castToPetsc()
{
  return rep_;
}

const cc::_p_Vec* Vec::castToConstPetsc() const
{
  return rep_;
}

real* Vec::begin()
{
  return p_ + localLow_;
}

const real* Vec::begin() const
{
  return p_ + localLow_;
}

real* Vec::end()
{
  return p_ + myLocalIndices_[ldim_-1] + 1;
}

const real* Vec::end() const
{
  return p_ + myLocalIndices_[ldim_-1] + 1;
}

bool Vec::checkConformance(const Vec& v) const
{
  if (dim_ == v.dim_ && ldim_ == v.ldim_)
    {
      return true;
    }
  else if (dim_ == v.dim_)
    {
      cerr<<"non-conforming vectors in expression: local dimension mismatch"<<endl
          <<ldim_<<'\t'<<v.ldim_<<endl;
    }
  else
    {
      cerr<<"non-conforming vectors in expression: global dimension mismatch"<<endl
          <<dim_<<'\t'<<v.dim_<<endl;
      //mwf put in for debugging
      //assert(0);
    }
  return false;
}

void Vec::setExample()
{
  Petsc::Err ierr;
  isRegistered_=true;
  if (vecRegistry[storageLength_].first)
    {
      if (!checkConformance(*vecRegistry[storageLength_].first))
        {
          cerr<<"attempting to use non conforming global vectors of the same dimension as example layout"<<endl;
          exit(1);
        }
      else
        vecRegistry[storageLength_].second.insert(this);
    }
  else
    {
      vecRegistry[storageLength_].first = this;
 
      std::set<Vec*>::iterator i=vecRegistry[storageLength_].second.begin();
      while (i != vecRegistry[storageLength_].second.end())
        {
          using namespace cc;

          _p_Vec* newRep_;
          if (this->dm_ != 0)
            ierr =  DMCreateGlobalVector(this->dm_,&newRep_);
          else
            ierr = VecDuplicate(this->rep_,&newRep_);
          (*i)->restoreLocal();
          //copy old vector into new vector layout
          VecScatter ctx;
          IS is;
          ierr = ISCreateStride((*i)->mpiComm_,(*i)->storageLength_,0,1,&is);
          ierr = VecScatterCreate((*i)->rep_,is,newRep_,is,&ctx);
//           ierr = VecScatterBegin((*i)->rep_,newRep_,INSERT_VALUES ,SCATTER_FORWARD, ctx);
//           ierr = VecScatterEnd((*i)->rep_,newRep_,INSERT_VALUES ,SCATTER_FORWARD, ctx);
          ierr = VecScatterBegin(ctx,(*i)->rep_,newRep_,INSERT_VALUES ,SCATTER_FORWARD);
          ierr = VecScatterEnd(ctx,(*i)->rep_,newRep_,INSERT_VALUES ,SCATTER_FORWARD);
          ierr = VecScatterDestroy(&ctx);
          ierr = ISDestroy(&is);
          
          (*i)->clear(false);//don't deRegister
          (*i)->rep_=newRep_;
          (*i)->ownsVec_=true;
          (*i)->mpiComm_=PETSC_COMM_WORLD;
          (*i)->dm_  = dm_;
          
          //from createStrideMappings -- borrow instead of copy
          (*i)->ownsIndexSet_=false;
          (*i)->myGlobalIndices_=myGlobalIndices_;
          (*i)->globalIndexSet_=globalIndexSet_;
          (*i)->myLocalIndices_=myLocalIndices_;
          (*i)->localIndexSet_=localIndexSet_;
          (*i)->myVirtualGlobalIndices_=myVirtualGlobalIndices_;
          (*i)->virtualGlobalIndexSet_=virtualGlobalIndexSet_;      
          (*i)->ldim_ = ldim_;
          (*i)->offSet_ = offSet_;
          (*i)->globalLow_ = globalLow_;
          (*i)->globalHigh_ = globalHigh_;
          (*i)->localLow_ = localLow_;
          (*i)->localHigh_ = localHigh_;
          (*i)->virtualGlobalLow_ = virtualGlobalLow_;
          (*i)->virtualGlobalHigh_ = virtualGlobalHigh_;
          
          (*i)->storageLength_ = storageLength_;    //should not have changed
          (*i)->lstorageLength_ = lstorageLength_;
          (*i)->dim_=storageLength_;
          
          //get new local pointer
          (*i)->hasLocal_=false;
          (*i)->p_=0;
          (*i)->pSave_=0;
          (*i)->getLocal();
          ++i;
        }
    }
}

void Vec::enterRegistry(int n)
{
  using namespace cc;
  Petsc::Err ierr;
  if (vecRegistry[n].first)
    {
      clear();//do deregister if already registered with different storageLength_
      ownsVec_ = true;
      mpiComm_ = PETSC_COMM_WORLD;
      ierr = VecDuplicate(vecRegistry[n].first->rep_,&rep_);
      dm_  = vecRegistry[n].first->dm_;
      dim_ = n;
      storageLength_ = n;
      lstorageLength_ = vecRegistry[n].first->lstorageLength_;
      
      getLocal();

      //from createStrideMappings -- borrow instead of copy
      ownsIndexSet_=false;
      myGlobalIndices_=vecRegistry[n].first->myGlobalIndices_;
      globalIndexSet_=vecRegistry[n].first->globalIndexSet_;
      myLocalIndices_=vecRegistry[n].first->myLocalIndices_;
      localIndexSet_=vecRegistry[n].first->localIndexSet_;
      myVirtualGlobalIndices_=vecRegistry[n].first->myVirtualGlobalIndices_;
      virtualGlobalIndexSet_=vecRegistry[n].first->virtualGlobalIndexSet_;      
      ldim_ = vecRegistry[n].first->ldim_;
      offSet_ = vecRegistry[n].first->offSet_;
      globalLow_ = vecRegistry[n].first->globalLow_;
      globalHigh_ = vecRegistry[n].first->globalHigh_;
      localLow_ = vecRegistry[n].first->localLow_;
      localHigh_ = vecRegistry[n].first->localHigh_;
      virtualGlobalLow_ = vecRegistry[n].first->virtualGlobalLow_;
      virtualGlobalHigh_ = vecRegistry[n].first->virtualGlobalHigh_;
    }
  isRegistered_=true;
  vecRegistry[n].second.insert(this);
}

}//Petsc
}//Daetk

namespace PetscVecOperators
{
//hack only works when s = cout
//apparently no way to view a slice
std::ostream&  operator<<(std::ostream& s, const Daetk::Petsc::Vec& v)
{ 
  using std::endl;
  using namespace Daetk;
  using namespace Petsc;
  using namespace Petsc::cc;
  Petsc::Err ierr;
  Sys pSys;
  int tag=10;
  int               rank,len,work = v.ldim(),n,size;
  MPI_Status        status;
  real              *values;
  real* array=const_cast<real*>(v.castToConstArray());
  /* determine maximum message to arrive */
  rank = pSys.getRank();
  size = pSys.getSize();
  ierr = MPI_Reduce(&work,&len,1,MPI_INT,MPI_MAX,0,v.mpiComm_);

  if (!rank) 
    {
      for (int i=0; i<v.ldim(); i++)         
        s<<v[i]<<endl;

          /* receive and print messages */

      values = new real[len+1];
      assert(values);
      
      for (int j=1; j<size; j++) 
        {
          ierr = MPI_Recv(values,len,MPIU_SCALAR,j,tag,v.mpiComm_,&status);
          ierr = MPI_Get_count(&status,MPIU_SCALAR,&n);         
          for (int i=0; i<n; i++) 
            s<<values[i]<<endl;
        } 
      delete values;
    } 
  else 
    {
    /* send values */
      ierr = MPI_Send(array,v.ldim(),MPIU_SCALAR,0,tag,v.mpiComm_);
    }
  
  return s;
}

std::istream&  operator>>(std::istream& s, Daetk::Petsc::Vec& v)
{  
  using namespace Daetk;
  using namespace Petsc;
  using namespace Petsc::cc;
  Daetk::Petsc::Err ierr;
  int length;
  if (v.petSys_.master())
    s>>length;
  MPI_Bcast(&length, 1, MPI_INT, 0, PETSC_COMM_WORLD);
  v.newsize(length);  
  v.restoreLocal();
  if (v.petSys_.master())
    {
      int*  vi = new int[length];
      real* vv = new real[length];
      for (int i=0; i<length ; i++)
        {
          s>>vv[i];
          vi[i] = i;
        }
      delete vi;
      delete vv;
      ierr =  VecSetValues(v.rep_,length,vi,vv,INSERT_VALUES); 
    }
  ierr =  VecAssemblyBegin(v.rep_);
  ierr =  VecAssemblyEnd(v.rep_);
  v.getLocal();
  return s;
}

Daetk::real norm(const Daetk::Petsc::Vec& v)
{
  return nrm2(v);
//    Daetk::Petsc::Err ierr;
//    Daetk::real n;
//    v.restoreLocal();
//  //    ierr =  VecSetBlockSize(v.rep_,v.stride_);
//    ierr =  VecNorm(v.rep_,v.start_,NORM_2,&n);
//    v.getLocal();
//    return n;
}

Daetk::real dot(const Daetk::Petsc::Vec& w,const Daetk::Petsc::Vec& v)  
{  
  Daetk::Petsc::Err ierr;
  Daetk::real d;
  w.restoreLocal();
  v.restoreLocal();
  ierr =  Daetk::Petsc::cc::VecDot(w.rep_,v.rep_,&d);
  v.getLocal();
  w.getLocal();
  return d;
}

Daetk::real nrm2(const Daetk::Petsc::Vec& v)  
{  
  Daetk::Petsc::Err ierr;
  Daetk::real n;
  v.restoreLocal();
//    ierr =  VecSetBlockSize(v.rep_,v.stride_);
//    ierr =  VecStrideNorm(v.rep_,v.start_,NORM_2,&n);
  ierr =  Daetk::Petsc::cc::VecNorm(v.rep_,Daetk::Petsc::cc::NORM_2,&n);
  v.getLocal();
  return n;
}

Daetk::real max(const Daetk::Petsc::Vec& v)
{
  Daetk::Petsc::Err ierr;
  Daetk::real m;
  Daetk::Petsc::cc::PetscInt loc;
  v.restoreLocal();
//    ierr =  VecSetBlockSize(v.rep_,v.stride_);
//    ierr =  VecStrideMax(v.rep_,v.start_,PETSC_NULL,&m);
  ierr =  Daetk::Petsc::cc::VecMax(v.rep_,&loc,&m);
  v.getLocal();
  return m;
}

//Daetk::real max(const Daetk::Petsc::Vec& v, int& loc)
//{
//  Daetk::Petsc::Err ierr;
//  Daetk::real m;
//  v.restoreLocal();
//  ierr =  Daetk::Petsc::cc::VecSetBlockSize(v.rep_,v.stride_);
//  ierr =  Daetk::Petsc::cc::VecStrideMax(v.rep_,v.start_,&loc,&m);
//  v.getLocal();
//  return m;
//}

Daetk::real min(const Daetk::Petsc::Vec& v)
{
  Daetk::Petsc::Err ierr;
  Daetk::real m;
  Daetk::Petsc::cc::PetscInt loc;
  v.restoreLocal();
//    ierr =  VecSetBlockSize(v.rep_,v.stride_);
//    ierr =  VecStrideMin(v.rep_,v.start_,PETSC_NULL,&m);
  using namespace Daetk::Petsc::cc;
  ierr =  VecMin(v.rep_,&loc,&m);
  v.getLocal();
  return m;
}

//Daetk::real min(const Daetk::Petsc::Vec& v, int& loc)
//{
//  Daetk::Petsc::Err ierr;
//  Daetk::real m;
//  v.restoreLocal();
//  ierr =  Daetk::Petsc::cc::VecSetBlockSize(v.rep_,v.stride_);
//  ierr =  Daetk::Petsc::cc::VecStrideMin(v.rep_,v.start_,&loc,&m);
//  v.getLocal();
//  return m;
//}

}//PetscVecOperators





