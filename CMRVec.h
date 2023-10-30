#ifndef CMRVEC_H
#define CMRVEC_H    

#include <iostream>
#include "Definitions.h"
#include "PETE/PETE.h"

namespace Pete
{
  using Daetk::real;
}
#include "CMRVecDefinitions.h"
#include "Utilities.h"
#include "CMRVecIndex.h"
#include "Tracer.h"
namespace Daetk 
{
template<class T> class CMRVecIterator;
template<class T> class CMRVecConstIterator;

//  using namespace std;

template<class T>
class CMRVec
{ 
public:
  CMRVec();
  CMRVec( int length);
  CMRVec( int length, const T& initialValue, int base=CMRVec_BASE,
          int stride=1);   
  CMRVec(const CMRVec& V);
  enum CopyType { NOT_REF, REF }; 
  CMRVec(CopyType t,T* a,  int storageLength, 
         int base=CMRVec_BASE,  int stride=1); 
  CMRVec(CopyType t,const CMRVec& V,const CMRVecIndex& VI); 
  virtual ~CMRVec();                              
  inline T& operator()(int i);
  inline const T& operator()(int i) const; 
  inline T& operator[](int i);
  inline const T& operator[](int i) const; 
  CMRVec operator()(const CMRVecIndex &VI) const;
  CMRVec& operator=(const CMRVec& V);
  CMRVec& operator=(const T& r);

  template <class RHS>
  CMRVec& operator+=(const Pete::Expression<RHS> &rhs)
  {
    int index=base_,end=base_+dim_;
    T* lhs=p_;
    for(index=base_;index<end;index++)
      {
        *lhs+=forEach(rhs,Pete::EvalLeaf1(index),Pete::OpCombine());
        lhs+=stride_;
      }
    return *this;
  }

  template <class RHS>
  CMRVec& operator=(const Pete::Expression<RHS> &rhs)
  {
    int index=base_,end=base_+dim_;
    T* lhs=p_;
    for(index=base_;index<end;index++)
      {
        *lhs=forEach(rhs,Pete::EvalLeaf1(index),Pete::OpCombine());
        lhs+=stride_;
      }
    return *this;
  }

  template <class RHS>
  CMRVec& operator-=(const Pete::Expression<RHS> &rhs)
  {
    int index=base_,end=base_+dim_;
    T* lhs=p_;
    for(index=base_;index<end;index++)
      {
        *lhs-=forEach(rhs,Pete::EvalLeaf1(index),Pete::OpCombine());
        lhs+=stride_;
      }
    return *this;
  }

  CMRVec& operator*=(const T &a);
  CMRVec& operator/=(const T &a);
  CMRVec& operator+=(const CMRVec &y);
  CMRVec& operator-=(const CMRVec &y);
  const CMRVec& print(std::ostream& s) const;
  inline  int size() const;
  inline  int dim() const;
  inline  int storageLength() const;
  inline int base() const;
  inline  int stride() const;
  inline bool isRef() const;
  inline bool isNull() const;
  CMRVec& clear();
  CMRVec& setSize( int n);
  CMRVec& redim( int n);
  CMRVec& newsize( int n);
  CMRVec& redimCopy( int n);
  CMRVec& grow( int n);
  CMRVec& setStride( int newStride);
  const T& push(const T& val);
  const T& pop();
  CMRVec& setBase(int b);
  CMRVec& attachToVec(CopyType t,const CMRVec& V,const CMRVecIndex& VI);
  CMRVec& attachToArray(T *a, int length,int base, 
                                int stride=1);
  T* castToArray();
  const T* castToConstArray() const;
  friend class CMRVecIterator<T>;
  friend class CMRVecConstIterator<T>;
  typedef CMRVecIterator<T> Iterator;
  typedef CMRVecConstIterator<T> ConstIterator;
  typedef T* UnitStrideIterator;
  typedef const T* ConstUnitStrideIterator;
  inline T* begin();
  inline const T* begin() const;
  inline T* end();
  inline const T* end() const;
  T max() const;
  T min() const;
  
protected:                                                           
  T* p_;
   int dim_;
  bool ref_;
  int base_;
   int stride_;
   int storageLength_;
};                          


template<class T>
inline T& CMRVec<T>::operator()(int i)
{
#ifdef CMRVEC_BOUNDS_CHECK
 assert((i-base_)*stride_ < storageLength_);
#endif
  return *(p_+(i-base_)*stride_);
}

template<class T>
inline const T& CMRVec<T>::operator()(int i) const 
{
#ifdef CMRVEC_BOUNDS_CHECK
  assert((i-base_)*stride_ < storageLength_);
#endif
  return *(p_+(i-base_)*stride_);
}

template<class T>
inline T& CMRVec<T>::operator[](int i)
{
#ifdef CMRVEC_BOUNDS_CHECK
   assert((i-base_)*stride_ < storageLength_);
#endif
  return *(p_+(i-base_)*stride_);
}

template<class T>
inline const T& CMRVec<T>::operator[](int i) const 
{
#ifdef CMRVEC_BOUNDS_CHECK
   assert((i-base_)*stride_ < storageLength_);
#endif
  return *(p_+(i-base_)*stride_);
}

template<class T>
inline  int CMRVec<T>::size() const 
{ 
  return dim_;
}

template<class T>
inline  int CMRVec<T>::dim() const
{
  return dim_;
}
  
template<class T>
inline  int CMRVec<T>::storageLength() const
{
  return storageLength_;
}

template<class T>
inline bool CMRVec<T>::isRef() const
{
  return ref_;
}

template<class T>
inline bool CMRVec<T>::isNull() const
{
  return dim_== 0;
}

template<class T>
inline int CMRVec<T>::base() const
{
  return base_;
}
  
template<class T>
inline  int CMRVec<T>::stride() const
{
  return stride_;
}


template<class T>
const CMRVec<T>& CMRVec<T>::print(std::ostream& c) const
{
  T* end=p_+storageLength_;
  T* pi=p_;
  while (pi < end)
    std::cout<<*pi++<<std::endl;
  return *this;
}

template<class T>
CMRVec<T>& CMRVec<T>::setBase(int b)
{
  base_=b;
  return *this;
}

template<class T>
CMRVec<T>& CMRVec<T>::attachToArray(T *a, int n, int base, 
                                      int stride)
{
  if (!ref_ &&  p_)
    delete [] p_;

  ref_ = REF;
  p_=a;
  dim_=n;
  base_=base;
  stride_=stride;
  return *this;
}
  
template<class T>
T* CMRVec<T>::castToArray()
{
  return p_;
}

template<class T>
const T* CMRVec<T>::castToConstArray() const
{
  return p_;
}

template<class T>
inline T* CMRVec<T>::begin()
{
  return p_;
}
  

template<class T>
inline const T* CMRVec<T>::begin() const
{
  return p_;
}

template<class T>
inline T* CMRVec<T>::end()
{
  return p_+storageLength_;
}
  
template<class T>
inline const T* CMRVec<T>::end() const
{
  return p_+storageLength_;
}

template<class T>
T
CMRVec<T>::max() const
{
  T tmp=0;
  if (dim_ < 1)  
    return tmp;
  else 
    {
      T* ptr = p_;
      T* pend = p_+storageLength_;
      tmp = *ptr;
      while (++ptr < pend)
        (tmp > *ptr) ? 0 : (tmp = *ptr);
      return tmp;
    }
}

template<class T>
T
CMRVec<T>::min() const
{
  T tmp=0;
  if (dim_ < 1)  return tmp;
  else {
    T* ptr = p_;
    T* pend = p_+storageLength_;
    tmp = *ptr;
    while (++ptr < pend)
      (tmp < *ptr) ? 0 : (tmp = *ptr);
    return tmp;
  }
}


template<class T>
CMRVec<T>::CMRVec():
  p_(0),
  dim_(0),
  ref_(NOT_REF),
  base_(CMRVec_BASE),
  stride_(1),
  storageLength_(0)
{
  Tracer tr("CMRVec<T>::CMRVec()");
}

template<class T>
CMRVec<T>::CMRVec( int n):
  p_(NULL),
  dim_(n),
  ref_(NOT_REF),
  base_(CMRVec_BASE),
  stride_(1),
  storageLength_(n)
{  
  Tracer tr("CMRVec<T>::CMRVec( int n)");
  p_=new T[n];
  if (p_ == NULL)
    {
      std::cerr << "Error: NULL pointer in CMRVec(int) constructor " << std::endl;
      std::cerr << "       Most likely out of memory... " << std::endl;
      exit(1);
    }
}
 
template<class T>
CMRVec<T>::CMRVec( int n, const T& v, int base, int stride): 
  p_(NULL),
  dim_(n),
  ref_(NOT_REF),
  base_(base),
  stride_(stride),
  storageLength_(n*stride)
{
  Tracer tr("CMRVec<T>::CMRVec( int n, const T& v, int base, int stride)");
  p_=new T[n*stride];
  if (p_ == NULL)
    {
      std::cerr << "Error: NULL pointer in CMRVec(int) constructor " << std::endl;
      std::cerr << "       Most likely out of memory... " << std::endl;
      exit(1);
    }
  for ( int i=0; i<storageLength_; i+=stride_)
    p_[i] = v;
}


template<class T>
CMRVec<T>::CMRVec(const CMRVec & m): 
  p_(new T[m.storageLength_]), 
  dim_(m.dim_) , 
  ref_(NOT_REF),
  base_(m.base_),
  stride_(m.stride_),
  storageLength_(m.storageLength_)
{
  Tracer tr("CMRVec<T>::CMRVec(const CMRVec & m)");
  if (p_ == NULL)
    {
      std::cerr << "Error:  Null pointer in CMRVec(const MV_Vector&); " << std::endl;
      exit(1);
    }
  for ( int i=0; i< storageLength_; i+=stride_)
    p_[i] = m.p_[i];
}




template<class T>
CMRVec<T>::CMRVec(CMRVec::CopyType ref,T* a, 
                   int storageLength, int base,  int stride):
  p_(a), 
  dim_(storageLength/stride + ((storageLength%stride > 0) ? 1 : 0)),
  ref_(ref),
  base_(base),
  stride_(stride),
  storageLength_(storageLength)
{
  Tracer tr("CMRVec<T>::CMRVec(CMRVec<T>::CopyType ref,T* a,  int storageLength, int base,  int stride)");
  if(!ref)//Copy the logical vector specified by a,stride, and storageLength  
    {
      //use base_,ref_, and dim_ as they are from initialization
      storageLength_=dim_;//create only enough storage to hold the vector
      p_=new T[storageLength_];
      if (p_ == NULL)
        {
          std::cerr << "Error: NULL pointer in CMRVec(int) constructor " << std::endl;
          std::cerr << "       Most likely out of memory... " << std::endl;
          exit(1);
        }
      for ( int i=0;i<dim_;i++)
        p_[i]=a[i*stride];
    }
}

template<class T>
CMRVec<T>& CMRVec<T>::attachToVec(CMRVec::CopyType ref,
                                  const CMRVec& V,
                                  const CMRVecIndex& VI)
{
  if (p_ && !ref_)
    delete [] p_;
  p_=V.p_+ (VI.start() - V.base_)*V.stride_; 
  dim_=VI.end() - VI.start() + 1;
  ref_=ref;
  base_=V.base_;
  stride_=V.stride_;
  storageLength_=(VI.end() - VI.start())*V.stride_ + 1; //this is the minimum storage length
  if (!ref) //if a copy of V(VI) is desired instead of a reference
    {
      //use base_,ref_, and dim_ as they are from initialization
      storageLength_=dim_;//create only enough storage for the elements 
      stride_=1;
      p_= new T[storageLength_];
      for ( int i=0; i<storageLength_; i++)
        p_[i]=V.p_[i*V.stride_];
    }
  return *this;
}

template<class T>
CMRVec<T>::CMRVec(CMRVec::CopyType ref, const CMRVec& V,
                  const CMRVecIndex& VI):
  p_(V.p_+ (VI.start() - V.base_)*V.stride_ ), 
  dim_(VI.end() - VI.start() + 1),
  ref_(ref),
  base_(V.base_),
  stride_(V.stride_),
  storageLength_((VI.end() - VI.start())*V.stride_ + 1)
{
  Tracer tr("CMRVec<T>::CMRVec(CMRVec<T>::CopyType ref, const CMRVec& V, const CMRVecIndex& VI)");
  if (!ref) //if a copy of V(VI) is desired instead of a reference
    {
      //use base_,ref_, and dim_ as they are from initialization
      storageLength_=dim_;//create only enough storage for the elements 
      stride_=1;
      p_= new T[storageLength_];
      for ( int i=0; i<storageLength_; i++)
        p_[i]=V.p_[i*V.stride_];
    }
}


template<class T>
CMRVec<T>::~CMRVec()
{
  Tracer tr("CMRVec<T>::~CMRVec()");
  if (!ref_ && p_) delete [] p_;
}


template<class T>
CMRVec<T> CMRVec<T>::operator()(const CMRVecIndex &VI) const
{
  if (VI.all())
    return CMRVec( REF,p_, storageLength_,base_,stride_);
  else
    {
      // check that index is not out of bounds
      if ( VI.end() >= dim_)
        {
          std::cerr << "MV_VecIndex: (" << VI.start() << ":" << VI.end() << 
            ") too big for matrix (0:" << dim_ - 1 << ") " << std::endl;
          exit(1);
        }
      return CMRVec(REF, p_+ VI.start()*stride_, 
                    (VI.end() - VI.start())*stride_ + 1, base_, stride_);
    }
}

template<class T>
CMRVec<T>& CMRVec<T>::setStride( int newStride)
{
  stride_ = newStride;
  dim_ = storageLength_ / newStride + ((storageLength_%newStride > 0) ? 1 : 0);
  return *this;
}

template<class T>
CMRVec<T>& CMRVec<T>::clear()
{
  dim_ = 0;
  return *this;
}

template<class T>
CMRVec<T>& CMRVec<T>::setSize( int n)
{
  return newsize(n);
}

template<class T>
CMRVec<T>& CMRVec<T>::redim( int n)
{
  return newsize(n);
}

template<class T>
CMRVec<T>& CMRVec<T>::newsize( int n)
{
#ifdef TRACE_VEC
  std::cout << "> CMRVec<T>::newsize( int n) " << std::endl;
#endif
  if (ref_ )                  // is this structure just a pointer?
    {
      {
        std::cerr << "MV_Vector<T>::newsize can't operator on references.\n";
        exit(1);
      }
    }
  else
    if (dim_ != n )                     // only delete and new if
      {                                 // the size of memory is Tly
        delete [] p_;                   // changing, otherwise just
        storageLength_=n*stride_;       // copy in place
        p_ = new T[storageLength_];
        if (p_ == NULL)
          {
            std::cerr << "Error : NULL pointer in newsize " << std::endl;
            exit(1);
          }
        dim_ = n;
      }
  
#ifdef TRACE_VEC
  std::cout << "< CMRVec<T>::newsize( int n) " << std::endl;
#endif
  
  return *this;
}

template<class T>
CMRVec<T>& CMRVec<T>::redimCopy( int n)
{
#ifdef TRACE_VEC
  std::cout << "> CMRVec<T>::redimCopy( int n) " << std::endl;
#endif
  if (ref_ )                  // is this structure just a pointer?
    {
      std::cerr <<"CMRVec<T>::redimCopy  can't operate on references."
           <<std::endl;
      exit(1);
    }
  else
    if (dim_ != n )
      {
         int newStorageLength = n*stride_;
        T* pnew = new T[newStorageLength];
        assert (pnew);
        int ncopy = (storageLength_ < newStorageLength) ? storageLength_ : newStorageLength;
        for (int i=0; i<ncopy; i++)
          pnew[i] = p_[i];
        delete [] p_;
        p_ = pnew;
        dim_ = n;
        storageLength_ = newStorageLength;
      }
#ifdef TRACE_VEC
  std::cout << "< CMRVec<T>::redimCopy( int n) " << std::endl;
#endif
  return *this;
}


template<class T>
CMRVec<T>& CMRVec<T>::grow( int n)
{
  if (ref_ )
    {
      std::cerr <<"CMRVec<T>::redimCopy  can't operate on references."
           <<std::endl;
      exit(1);
    }
  else
    {
       int newStorageLength = (dim_ + n)*stride_;
      T* pnew = new T[newStorageLength];
      assert(pnew);
      int ncopy = (storageLength_ < newStorageLength) ? storageLength_ : newStorageLength;
      for (int i=0; i<ncopy; i++)
        pnew[i] = p_[i];
      delete [] p_;
      p_ = pnew;
      //  dim_ = n;
      storageLength_ = newStorageLength;
    }
  return *this;
}


template<class T> const T&
CMRVec<T>::push(const T& val)
{
  grow(1);
  p_[(dim_ - 1)*stride_]= val;
  return p_[(dim_ - 1)*stride_];
}

template<class T> const T& 
CMRVec<T>::pop()
{
  assert(dim_!=0);
  static T temp(p_[0]);
  temp=p_[(dim_-1)*stride_];
  redimCopy(dim_-1);
  return temp;
}


template<class T>
CMRVec<T>& CMRVec<T>::operator=(const CMRVec & m) 
{
  T *it(p_),*last(p_ + stride_*(dim_-1)),*mIt(m.p_);
  
  if (ref_ || m.ref_)         // is either structure just a pointer?
    {
      if (dim_ != m.dim_)     // check conformance,
        {
          std::cerr << "CMRVec<T>::operator= used with vector refs of different dims";
          exit(1);
        }
      
      // We will not allow m to be overwritten in the case of storage
      // overlap
      
      T *mLast(m.p_ + m.stride_*(m.dim_-1)),*start(p_),*mStart(m.p_);
      
      // case:  [start       [mStart         last]          mLast]
      if ( mStart >= start && mStart <= last )
        {
          while( it <= last )
            {
              if ( it == mStart )//we would begin to overwrite p_
                {
                  std::cerr<<"CMRVec<T>::operator= used with overlapping vector refs";
                }
              *it=*mIt;
              it+=stride_;
              mIt+=m.stride_;
            }
          return *this;
        }
      // case [mStart     [start          mLast]        last]
      if ( start >= mStart && start <= mLast )
	{
	  // we'll copy backwards in this case, so that we could actually
	  // throw an exception before and overwrite occurs
          it=last;
          mIt=mLast;
          while( it >= start )
            {
              if ( it == mLast )//we would begin to overwrite m.p_
                {
                  std::cerr<<"CMRVec<T>::operator= used with overlapping vector refs";
                }
              *it=*mIt;
              it-=stride_;
              mIt-=m.stride_;
            }
          return *this;
        }
      else
        {
          while ( it <= last )
            {
              *it=*mIt;
              it+=stride_;
              mIt+=m.stride_;
            }
          return *this;
        }
    }
  else
    {
      newsize(m.dim_*stride_);//preserve old stride--should do nothing
      //if m.dim_ == dim_ since storageLength_=dim_*stride_
      it = p_ ;
      last = p_+storageLength_-stride_;
      
      
      // no need to test for overlap, since this region is new
      // and both vectors have internal storage
      while ( it <= last)
        {
          *it=*mIt;
          it+=stride_;
          mIt+=m.stride_;
        }
      return *this;   
    }
}


template<class T>
CMRVec<T>& CMRVec<T>::operator=(const T & m) 
{
#ifdef TRACE_VEC
  std::cout << "> CMRVec<T>::operator=(const T & m)  " << std::endl;
#endif
  
  T* it = p_;
  const T* end = p_ +storageLength_;
  while (it < end)
    {
      *it = m;
      it+=stride_;
    }
  
#ifdef TRACE_VEC
  std::cout << "< CMRVec<T>::operator=(const T & m)  " << std::endl;
#endif
  return *this;
}

template<class T>
CMRVec<T>& CMRVec<T>::operator*=(const T &a)
{
  T* ix(begin());
  while (ix < end())
    {
      (*ix)*=a;
      ix+=stride_;
    }
  return *this;
}

template<class T>
CMRVec<T>& CMRVec<T>::operator/=(const T &a)
{
  T* ix(begin());
  while (ix < end())
    {
      (*ix)/=a;
      ix+=stride_;
    }
  return *this;
}

template<class T>
CMRVec<T>& CMRVec<T>::operator+=(const CMRVec<T> &y)
{
  assert(size() == y.size());
  const T* iy(y.begin());
  T* ix(begin());
  while (ix <end() )
    {
      (*ix)+=(*iy);
      ix+=stride_;
      iy+=y.stride_;
    }
  return *this;
}
          
      
template<class T>
CMRVec<T>& CMRVec<T>::operator-=(const CMRVec<T> &y)
{
  assert(size() == y.size());
  /*
  if (x.size() != y.size())
    {
      std::cout << "Incompatible vector lengths in -." << std::endl;
      exit(1);
    }
    */
  const T* iy(y.begin());
  T* ix(begin());
  while (ix < end())
    {
      (*ix)-=(*iy);
      ix+=stride_;
      iy+=y.stride_;
    }
  return *this;
}
          
      
template<class T>
T dot(const CMRVec<T> &x, const CMRVec<T> &y)
{
        
  //  Check for compatible dimensions:
  assert(x.size() == y.size());
  /*
  if (x.size() != y.size())
      {
         std::cout << "Incompatible dimensions in dot(). " << std::endl;
         exit(1);
      }
      */
  T temp =  0.0;
  const T* ix(x.begin()),*iy(y.begin());
  while (ix < x.end())
    {
      temp+=(*ix)*(*iy);
      ix+=x.stride(); 
      iy+=y.stride();
    }
  return temp;
}

template<class T>
T norm(const CMRVec<T> &x)
{
      T temp = dot(x,x);
      return sqrt(temp);
}

//  inline real norm(const CMRVec<real> &x)
//  {
//        real temp = dot(x,x);
//        return sqrt(temp);
//  }

//using namespace std;
template<class T>
std::ostream&  operator<<(std::ostream& s, const CMRVec<T>& V)
{
  s.setf(std::ios::scientific);
  s.precision(17);
  for ( int i=0; i< V.dim() ; i++)
    s<<V(i)<<std::endl;
  return s;
}

template<class T>
std::istream&  operator>>(std::istream& s, CMRVec<T>& V)
{
  int length;
  s>>length;
  V.newsize(length);
  for (int i=0; i<length ; i++)
    s>>V(i);
  return s;
}
}//Daetk

namespace Pete
{
template<class T>
struct CreateLeaf<Daetk::CMRVec<T> >
{
  typedef Reference<Daetk::CMRVec<T> >  Leaf_t;
  inline static
  Leaf_t make(const Daetk::CMRVec<T> &a) { return Leaf_t(a); }
};

template<class T>
struct LeafFunctor<Daetk::CMRVec<T>, EvalLeaf1>
{
  typedef T Type_t;
  static Type_t apply(const Daetk::CMRVec<T> &a, const EvalLeaf1 &f)
    { return a(f.val1()); }
};
    
/*                                   
template<class T,class OP,class RHS>
inline void evaluate(Daetk::CMRVec<T> &lhs, const OP &op, 
                     const Expression<RHS> &rhs)
{
  for(int i=0;i<lhs.dim();i++)
    {
      int index = lhs.base() + i;
      op(lhs(index), forEach(rhs, EvalLeaf1(index), OpCombine()));
    }
}
*/
}//Daek

#endif  

