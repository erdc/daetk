#ifndef BANDCOLMAT_H
#define BANDCOLMAT_H

#include "Definitions.h"
#include "Vec.h"
#include "Mat.h"
#include "LinearOperator.h"

namespace Daetk 
{

class BandColMat : public LinearOperator
{
 public:
  BandColMat();
  BandColMat(int upperDiagonals,int lowerDiagonals,int Neq);
  virtual ~BandColMat();
  inline int* pivotPtr(){return pivots;}
  void newsize(int upperDiagonals,int lowerDiagonals,int Neq);
  inline real& operator()(const int i,const int j);
  BandColMat& operator=(const BandColMat& C);
  void beginAssembly();
  void endAssembly();

  BandColMat& operator=(const real& C); 
  void zeroRow(int i);
  void zeroAll();

  virtual Vec operator*(const Vec& x);
  virtual void matVec(const Vec& x, Vec& Ax);
  bool apply(const Vec& x, Vec& Ax);
  friend std::ostream& operator<<(std::ostream&,BandColMat&);
  real* castToArray();
  const real* castToConstArray() const;
  inline int getNeq() const;
  inline int getUpperBandWidth() const;
  inline int getLowerBandWidth() const;
 protected:
  Vec row;
  Mat M;
  int *pivots;
  int ku,kl,neq,rowOffSet;
};

inline int BandColMat::getNeq() const {return neq;}
inline int BandColMat::getUpperBandWidth() const {return ku;}
inline int BandColMat::getLowerBandWidth() const {return kl;}
  
inline real& BandColMat::operator()(const int i,const int j)
{
  return M(rowOffSet+ku+1+i-j,j);
}

inline Vec BandColMat::operator*(const Vec& x)
{
  assert(1);
  Vec result;
  return result;
}

}//Daetk
#endif









