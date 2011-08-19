#ifndef CMRVECITERATOR_H
#define CMRVECITERATOR_H

#include "CMRVec.h"
  
///
template<class T>
class CMRVecIterator
{
public:
  ///
  inline CMRVecIterator(CMRVec<T>& V);
  ///
  inline CMRVecIterator(const CMRVecIterator& I);
  ///
  inline T& operator*();
  ///
  inline CMRVecIterator& operator=(T* rhs);
  ///
  inline CMRVecIterator& operator++();
  ///
  inline CMRVecIterator operator++(int);
  ///
  inline CMRVecIterator& operator--();
  ///
  inline CMRVecIterator operator--(int);
  ///
  inline CMRVecIterator& operator+=(const int& n);
  ///
  inline CMRVecIterator& operator-=(const int& n);
  ///
  inline CMRVecIterator operator+(const int& n);
  ///
  inline CMRVecIterator operator-(const int& n);
  ///
  inline T& operator[](const int& n);
  ///
  inline bool operator==(const T* rhs);
  ///
  inline bool operator!=(const T* rhs);
  ///
  inline bool operator<(const T* rhs);
private:
  T* p_;
  const unsigned int& stride_;
};

template<class T> 
inline CMRVecIterator<T>::CMRVecIterator(CMRVec<T>& V):
  p_(V.p_),
  stride_(V.stride_)
{}


template<class T> 
inline CMRVecIterator<T>::CMRVecIterator(const CMRVecIterator<T>& I):
  p_(I.p_),
  stride_(I.stride_)
{}

template<class T> 
inline T& CMRVecIterator<T>::operator*()
{
  return *p_;
}

template<class T> 
inline CMRVecIterator<T>& CMRVecIterator<T>::operator=(T* rhs)
{
  p_=rhs;
  return *this;
}

template<class T> 
inline CMRVecIterator<T>& CMRVecIterator<T>::operator++()
{
  p_+=stride_;
  return *this;
}

template<class T> 
inline CMRVecIterator<T> CMRVecIterator<T>::operator++(int)
{
  CMRVecIterator temporary(*this);
  p_+=stride_;
  return temporary;
}

  
template<class T> 
inline CMRVecIterator<T>& CMRVecIterator<T>::operator--()
{
  p_-=stride_;
  return *this;
}

template<class T> 
inline CMRVecIterator<T> CMRVecIterator<T>::operator--(int)
{
  CMRVecIterator<T> temporary(*this);
  p_-=stride_;
  return temporary;
}

template<class T> 
inline CMRVecIterator<T>& CMRVecIterator<T>::operator+=(const int& n)
{
  p_+=n*stride_;
  return *this;
}

template<class T> 
inline CMRVecIterator<T>& CMRVecIterator<T>::operator-=(const int& n)
{
  p_-=n*stride_;
  return *this;
}

template<class T> 
inline CMRVecIterator<T> CMRVecIterator<T>::operator+(const int& n)
{
  CMRVecIterator<T> temporary(*this);
  return temporary+=n;
}
template<class T> 
inline CMRVecIterator<T> CMRVecIterator<T>::operator-(const int& n)
{
  CMRVecIterator<T> temporary(*this);
  return temporary-=n;
}

template<class T> 
inline T& CMRVecIterator<T>::operator[](const int& n)
{
  return *(p_+n*stride_);
}

template<class T> 
inline bool CMRVecIterator<T>::operator==(const T* rhs)
{
  return p_ == rhs;
}


template<class T> 
inline bool CMRVecIterator<T>::operator!=(const T* rhs)
{
  return p_ != rhs;
}

template<class T> 
inline bool CMRVecIterator<T>::operator<(const T* rhs)
{
  return p_ < rhs;
}


#endif


