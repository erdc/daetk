#ifndef CMRVECCONSTITERATOR_H
#define CMRVECCONSTITERATOR_H
#include "Definitions.h"
#include "CMRVec.h"  
namespace Daetk 
{

///
template<class T>
class CMRVecConstIterator
{
public:
  ///
  inline CMRVecConstIterator(const CMRVec<T>& V);
  ///
  inline CMRVecConstIterator(const CMRVecConstIterator& I);
  ///
  inline const T& operator*();
  ///
  inline CMRVecConstIterator& operator=(const T* rhs);
  ///
  inline CMRVecConstIterator& operator++();
  ///
  inline CMRVecConstIterator operator++(int);
  ///
  inline CMRVecConstIterator& operator--();
  ///
  inline CMRVecConstIterator operator--(int);
  ///
  inline CMRVecConstIterator& operator+=(const int& n);
  ///
  inline CMRVecConstIterator& operator-=(const int& n);
  ///
  inline CMRVecConstIterator operator+(const int& n);
  ///
  inline CMRVecConstIterator operator-(const int& n);
  ///
  inline const T& operator[](const int& n);
  ///
  inline bool operator==(const T* rhs);
  ///
  inline bool operator!=(const T* rhs);
  ///
  inline bool operator<(const T* rhs);
private:
  const T* p_;
  const unsigned int& stride_;
};


template<class T>
inline CMRVecConstIterator<T>::CMRVecConstIterator(const CMRVec<T>& V):
  p_(V.p_),
  stride_(V.stride_)
{}

template<class T>
inline CMRVecConstIterator<T>::CMRVecConstIterator(const CMRVecConstIterator<T>& I):
  stride_(I.stride_),
  p_(I.p_)
{}

template<class T>
inline const T& CMRVecConstIterator<T>::operator*()
{
  return *p_;
}

template<class T>
inline CMRVecConstIterator<T>& CMRVecConstIterator<T>::operator=(const T* rhs)
{
  p_=rhs;
  return *this;
}

template<class T>
inline CMRVecConstIterator<T>& CMRVecConstIterator<T>::operator++()
{
  p_+=stride_;
  return *this;
}

template<class T>
inline CMRVecConstIterator<T> CMRVecConstIterator<T>::operator++(int)
{
  CMRVecConstIterator<T> temporary(*this);
  p_+=stride_;
  return temporary;
}

template<class T>
inline CMRVecConstIterator<T>& CMRVecConstIterator<T>::operator--()
{
  p_-=stride_;
  return *this;
}

template<class T>
inline CMRVecConstIterator<T> CMRVecConstIterator<T>::operator--(int)
{
  CMRVecConstIterator<T> temporary(*this);
  p_-=stride_;
  return temporary;
}
  
template<class T>
inline CMRVecConstIterator<T>& CMRVecConstIterator<T>::operator+=(const int& n)
{
  p_+=n*stride_;
  return *this;
}

template<class T>
inline CMRVecConstIterator<T>& CMRVecConstIterator<T>::operator-=(const int& n)
{
  p_-=n*stride_;
  return *this;
}

template<class T>
inline CMRVecConstIterator<T> CMRVecConstIterator<T>::operator+(const int& n)
{
  CMRVecConstIterator<T> temporary(*this);
  temporary+=n;
  return temporary;
}

template<class T>
inline CMRVecConstIterator<T> CMRVecConstIterator<T>::operator-(const int& n)
{
  CMRVecConstIterator<T> temporary(*this);
  temporary-=n;
  return temporary;
}

template<class T>
inline const T& CMRVecConstIterator<T>::operator[](const int& n)
{
  return *(p_+n*stride_);
}

template<class T>
inline bool CMRVecConstIterator<T>::operator==(const T* rhs)
{
  return p_ == rhs;
}
 
template<class T>
inline bool CMRVecConstIterator<T>::operator!=(const T* rhs)
{
  return p_ != rhs;
}
  
template<class T>
inline bool CMRVecConstIterator<T>::operator<(const T* rhs)
{
  return p_ < rhs;
}
}//Daetk
#endif



