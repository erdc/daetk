#ifndef PETSCVEC
#define PETSCVEC
#include "PETE/PETE.h"
#include "PetscSys.h"
#include "PetscVecIndex.h"
#include "Definitions.h"
#include <set>
#include <map>
#include "CMRVec.h"
#include "mpi.h"
#define PetscVec_BASE 0

namespace Petsc
{
  namespace cc
    {
  struct _p_IS;
  struct _p_DA;
  struct _p_Vec;
  struct _p_VecScatter;
    }
class Vec
{
public:
  //combing CMRVec default and dim based ctor
  Vec(int global_dim=0);

  //only initializes the slice defined by the stride
  Vec(int global_dim, const real& initialValue, int base=PetscVec_BASE,int stride=1);

  //deep copy
  Vec(const Vec& v);
  
  enum CopyType { NOT_REF, REF };
  
  //a is the local array of local length localDim*stride
  Vec(Vec::CopyType t, real* a, int localDim, int base=PetscVec_BASE, int stride=1);

  //ref copy or deep copy; deep copy has only enough storage for V(I) but global entries stay on the same pe's
  Vec(Vec::CopyType t,const Vec& V, const VecIndex& I); 

  enum StorageType{ GLOBAL, LOCAL};

  Vec(Vec::StorageType t, Petsc::cc::_p_DA* da);

  Vec(const CMRVec<real>& vin);

  virtual ~Vec();

  //new Petsc stuff
  //does Petsc internal check of the rep_
  bool isValid() const;
  //checks that the global and local dimensions are the same so the vector operations will without communication
  bool checkConformance(const Vec& v) const;

  //really only uses stdout, only prints whole vec regardless of base/stride/start can be fixed easily
  friend std::ostream&  operator<<(std::ostream& s, const Vec& v);
  friend std::istream&  operator>>(std::istream& s, Vec& v);
  friend real nrm2(const Vec& v);
  friend real dot(const Vec& w,const Vec& v);
  friend real norm(const Vec& v);
  friend real max(const Vec& v);
  friend real min(const Vec& v);
  inline real max(){return Petsc::max(*this);}
  inline real min(){return Petsc::min(*this);}
  friend real max(const Vec& v, int& loc);
  friend real min(const Vec& v, int& loc);
  real sum();

  void startAssembly();
  void endAssembly();

  void startSetFromGlobal(const Vec& global);
  void endSetFromGlobal(const Vec& global);
  //use these if you're going to repeatedly scatter to this local
  //and want to interleave scatters 
  void startSetFromGlobalMulti(const Vec& global);
  void endSetFromGlobalMulti(const Vec& global);

  void setFromLocal(const Vec& local);

  //indexing via the local base-zero indices
  const real& get(int i) const {return p_[myLocalIndices_[i]];}
  real& set(int i) {return p_[myLocalIndices_[i]];}

  //old CMRVec interface; note some functions were dropped or the meaning is slightly different in parallel

  //global indeces; i must be on this pe; need not be base 0;
  inline real& operator()(int i);
  inline const real& operator()(int i) const; 

  //changed [] behavior to use the local array index; need not be base 0;
  inline real& operator[](int i);
  inline const real& operator[](int i) const; 
  
  //dropped
  //Vec operator()(const VecIndex &I) const;
  
  Vec& operator=(const Vec& V);
  Vec& operator=(const real& r);
  
  //expression templates 


  template <class RHS>
  Vec& operator+=(const Expression<RHS> &rhs);
  template <class RHS>
  Vec& operator=(const Expression<RHS> &rhs);
  template <class RHS>
  Vec& operator-=(const Expression<RHS> &rhs);
 
  //standard expressions
  Vec& operator*=(const real &a);
  Vec& operator/=(const real &a);
  Vec& operator+=(const Vec &y);
  Vec& operator-=(const Vec &y);
  
  //dropped
  //const Vec& print(ostream& s) const;
 
   //global info
  int dim() const;
  int size() const; 
  int storageLength() const;
  int base() const;
  int stride() const;
  int start() const;
  bool isRef() const;
  bool isNull() const;
  
  //local info
  int ldim() const;
  int lstorageLength() const;

  //these index values are virtual not the ones for the rep_
  int getGlobalHigh() const;
  int getGlobalLow() const;
  int getLocalHigh() const;
  int getLocalLow() const;
  
  //delete the vector
  Vec& clear(bool deRegister=true);

  //dropped
  //Vec& setSize( int n);
  
  Vec& newsize( Vec::StorageType t, Petsc::cc::_p_DA* da);
  //this next is a hack for now, not really a local vector
  Vec& newsize( Vec::StorageType t, int n);
  Vec& newsize( int n);
  Vec& newsizeSerial( int n);
  Vec& redim( int n);
  
  //dropped
  //Vec& redimCopy( int n);
  //Vec& grow( int n);
  //const real& push(const real& val);
  //const real& pop();
  
  Vec& setStride( int newStride);
  Vec& setStrideMulti( int newStride);
  Vec& setBlockStrideMulti( int newStride, int block=1);
  Vec& setBase(int b);
  
  
  //local
  Vec& attachToVec(Vec::CopyType t,const Vec& V,const VecIndex& I);
  Vec& attachToVecMulti(Vec::CopyType t,const Vec& V,const VecIndex& I);
  Vec& attachToArray(real *a, int localDim,int base, int stride=1);
  Vec& attachToPetscRepMulti(Petsc::cc::_p_Vec* pv);
  Vec& detachFromPetscRepMulti();
  real* castToArray();
  const real* castToConstArray() const;
  Petsc::cc::_p_Vec* castToPetsc();
  const Petsc::cc::_p_Vec* castToConstPetsc() const;
  friend class VecIterator;
  friend class VecConstIterator;
  typedef VecIterator Iterator;
  typedef VecConstIterator ConstIterator;
  typedef real* UnitStrideIterator;
  typedef const real* ConstUnitStrideIterator;
  real* begin();
  const real* begin() const;
  real* end();
  const real* end() const;
  void restoreLocal() const;
  void getLocal() const;
  
  void setExample();
  Petsc::cc::_p_DA*  getDA()
  {
    return da_;
  }
  typedef std::map<int,std::pair<Vec*,std::set<Vec*> > > Registry;
  static Registry vecRegistry;
  void enterRegistry(int length);
protected:
  
  
  bool ownsVec_,                     //set if allocating storage,
    ownsScatter_,
    ownsIndexSet_,                   //set by createStrideMappings
    hasGhostPoints_,
    firstAttachment_,                 //set by all constructors
    firstSetStride_,
    isRegistered_;                 //set by all constructors

  mutable bool hasLocal_;            //set by get/restoreLocal()
  
  int dim_,                          //must be set from user input
    base_,                           //must be set from user input 
    stride_,                         //must be set from user input 
    start_,                          //must be set from user input
    storageLength_,                  //set by getStorageInfo
    lstorageLength_,                 //set by getStorageInfo
    globalHigh_,                     //set by getStorageInfo / reset by createStrideMappings
    globalLow_,                      //set by getStorageInfo / reset by createStrideMappings
    localHigh_,                      //set by getStorageInfo / reset by createStrideMappings
    localLow_,                       //set by getStorageInfo / reset by createStrideMappings
    offSet_,                         //set by createStrideMappings
    ldim_,                           //set by createStrideMappings
    virtualGlobalHigh_,              //set by createStrideMappings
    virtualGlobalLow_;               //set by createStrideMappings
  MPI_Comm  mpiComm_;
  mutable real* p_;     //set by getLocal() / restoreLocal()
  mutable real* pSave_; //used to test what happens after we call restoreLocal(); getLocal();

  //this has static data members
  Petsc::Sys petSys_;

  //these are pointers to the Petsc objects
  //with IS and DA we can do much more general things with these vectors

  // the following three index sets and arrays are set in createStrideMappings
  //the global indices corresponding to this vector in the petsc global numbering
  //0 - ldim -> globalLow_ - globalHigh_
  Petsc::cc::_p_IS*  globalIndexSet_;
  mutable int* myGlobalIndices_;
 
  //the local indices in the petsc local numbering
  //0 - ldim -> localLow_ - localHigh_
  Petsc::cc::_p_IS*  localIndexSet_; 
  mutable int* myLocalIndices_;
  
  //the global indices in this vector's global numbering
  //0 - ldim -> virtualGlobalLow_ - virtualGlobalHigh_
  Petsc::cc::_p_IS*  virtualGlobalIndexSet_; 
  mutable int* myVirtualGlobalIndices_;

  //the data
  Petsc::cc::_p_DA*  da_;  //if this vec comes from a distributed array
  Petsc::cc::_p_VecScatter*    my_gtol_;
  //the data
  Petsc::cc::_p_Vec*  rep_;

private:
  mutable Petsc::Err ierr;       //usefull for calling the petsc error checking

  void createStrideMappings(); 
  void createBlockStrideMappings(int block); 
  void getStorageInfo();
  void create(int length);
  void createWithLocal(int local_length);
};





} //Petsc

class SizeLeaf
{
public:
  
  SizeLeaf(const Petsc::Vec& v): 
    vec_m_ref(v) 
  {}
  SizeLeaf(const SizeLeaf &model):
    vec_m_ref(model.vec_m_ref) 
  {}
  bool operator()(const Petsc::Vec& v) const 
  { 
    return vec_m_ref.checkConformance(v);
  }
  
private:
  
  const Petsc::Vec& vec_m_ref;
  
};

template<class T>
struct LeafFunctor<Scalar<T>, SizeLeaf>
{
  typedef bool Type_t;
  inline static
  bool apply(const Scalar<T> &, const SizeLeaf &) 
  {
    // Scalars always conform.
    return true;
  }
};

//the expression templates use [] right now, would be more general to use iterators
template<>
struct LeafFunctor<Petsc::Vec, SizeLeaf>
{
  typedef bool Type_t;
  inline static
  bool apply(const Petsc::Vec &v, const SizeLeaf &s) 
  {
    return s(v);
  }
};

template<>
struct CreateLeaf<Petsc::Vec >
{
  typedef Reference<Petsc::Vec >  Leaf_t;
  inline static
  Leaf_t make(const Petsc::Vec &a) { return Leaf_t(a); }
};

template<>
struct LeafFunctor<Petsc::Vec, EvalLeaf1>
{
  typedef real Type_t;
  static Type_t apply(const Petsc::Vec &a, const EvalLeaf1 &f)
  { 
    return a.get(f.val1()); 
  }
};

namespace Petsc {
template <class RHS>
Vec& Vec::operator+=(const Expression<RHS> &rhs)
{
    if (forEach(rhs, SizeLeaf(*this), AndCombine()))
      {
        int i;
        for(i=0;i<ldim_;i++)
          {
            set(i)  += forEach(rhs,EvalLeaf1( i ),OpCombine());
          }
      }
    else
      std::cerr<<"Vectors in expression don't conform"<<std::endl;
    return *this;
}

template <class RHS>
Vec& Vec::operator=(const Expression<RHS> &rhs)
{
  if (forEach(rhs, SizeLeaf(*this), AndCombine()))
    {
      int i;
      for(i=0;i<ldim_;i++)
        {
          set(i)  = forEach(rhs,EvalLeaf1( i ),OpCombine());
        }
    }
  else
    std::cerr<<"Vectors in expression don't conform"<<std::endl;
  return *this;
}

template <class RHS>
Vec& Vec::operator-=(const Expression<RHS> &rhs)
{
  if (forEach(rhs, SizeLeaf(*this), AndCombine()))
    {
      int i;
      for(i=0;i<ldim_;i++)
        {
          set(i)  -= forEach(rhs,EvalLeaf1( i ),OpCombine());
        }
    }
  else
    std::cerr<<"Vectors in expression don't conform"<<std::endl;
  return *this;
}

real& Vec::operator[](int localIndex)
{
#ifdef PETSCVEC_BOUNDS_CHECK
  int local_0 = localIndex - base_;
  if( !( local_0 >=0 && local_0 < ldim_) )
    {
      std::cerr<<"Indexing error in Vec::operator[]"<<std::endl
          <<"with localIndex = "<<localIndex<<std::endl
          <<"base_ = "<<base_<<std::endl
          <<"ldim_ = "<<ldim_<<std::endl;
    }
#endif
  return p_[ myLocalIndices_[localIndex-base_] ];
}

const real& Vec::operator[](int localIndex) const
{
#ifdef PETSCVEC_BOUNDS_CHECK
  int local_0 = localIndex - base_;
  if( !( local_0 >=0 && local_0 < ldim_) )
    {
      std::cerr<<"Indexing error in Vec::operator[]"<<std::endl
          <<"with localIndex = "<<localIndex<<std::endl
          <<"base_ = "<<base_<<std::endl
          <<"ldim_ = "<<ldim_<<std::endl;
    }
#endif
  return p_[ myLocalIndices_[localIndex-base_] ];
}

real& Vec::operator()(int globalIndex)
{
  static real dummy;
  //indexing for strided vector only
#ifdef PETSCVEC_BOUNDS_CHECK
  int petsc_global = (globalIndex - base_)*stride_ + start_;
  if ( petsc_global >= globalLow_ && petsc_global < globalHigh_)
    {}
  else
    std::cerr<<"global index out of bounds"<<std::endl<<std::flush;
#endif
  if (globalIndex >= virtualGlobalLow_ && globalIndex < virtualGlobalHigh_)
    return p_[ myLocalIndices_[globalIndex - virtualGlobalLow_] ];
  else
    return dummy=7777777;
}

const real& Vec::operator()(int globalIndex) const
{
  static real dummy;
  //indexing for strided vector only
#ifdef PETSCVEC_BOUNDS_CHECK
  int petsc_global = (globalIndex - base_)*stride_ + start_;
  if ( petsc_global >= globalLow_ && petsc_global < globalHigh_)
    {}
  else
    std::cerr<<"global index out of bounds"<<std::endl<<std::flush;
#endif
  if (globalIndex >= virtualGlobalLow_ && globalIndex < virtualGlobalHigh_)
    return p_[ myLocalIndices_[globalIndex - virtualGlobalLow_] ];
  else
    return dummy=7777777;
}


}//Petsc
#endif













