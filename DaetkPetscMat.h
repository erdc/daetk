#ifndef DAETKPETSCMAT_H
#define DAETKPETSCMAT_H

#include "Definitions.h"
#include "LinearOperator.h"
#include "Vec.h"
#include "DaetkPetscSys.h"

namespace Daetk 
{
namespace Petsc
{
  namespace cc
    {
  struct _p_Mat;
    }
  class Slammer;

class Mat : public LinearOperator
{
public:
  Mat(int neq=0);
  virtual ~Mat();
  Mat(int nGlobalRows, int nDiagonals, int nOffDiagonals=0, 
           bool isSymmetric=false,int blockSize=1);
  Mat(const Mat& mIn);
  void newsizeSequential(int nGlobalRows, int nDiagonals, int nOffDiagonals=0, 
               bool isSymmetric=false,int blockSize=1);
  void newsize(int nGlobalRows, int nDiagonals, int nOffDiagonals=0, 
               bool isSymmetric=false,int blockSize=1);
  enum MatrixFormat {aij_format,block_aij_format};
  MatrixFormat format;
  void attachToMat(Petsc::cc::_p_Mat* mIn);
  inline Slammer operator()(const int &i,const int &j);
  void finalizeRow(int i);
  void finalizeBlockRow(int i);
  void finalizeAddRow(int i);
  void finalizeAddBlockRow(int i);
  const real& get(const int& i, const int& j);
  bool apply(const Vec& x, Vec& Ax);
  friend std::ostream& operator<<(std::ostream& s, const Mat& m);
  void beginAssembly();
  void endAssembly();
  void copy(Mat& mIn);
  enum RowsOrColumns {ROWS,COLUMNS};
  int getLocalDim(RowsOrColumns r);
  int getGlobalHigh();
  int getGlobalLow();
  int getBlockSize();
  void getStorageInfo();
  void  scaleDiagonal(const Vec& left, const Vec& right);
  void scale(real a);
  enum Insertion {ROW_ORIENTED,COLUMN_ORIENTED,ROWS_SORTED,
                  ROWS_UNSORTED,COLUMNS_SORTED,COLUMNS_UNSORTED};
  void setInsertionOrder(Insertion s);
  void zeroAll();
  void zeroRow(int i);
  bool isAssembled();
  bool isValid();  
  Petsc::cc::_p_Mat** castToPetscLValue();
  Petsc::cc::_p_Mat* castToPetsc();
  const Petsc::cc::_p_Mat* castToConstPetsc() const; 
  void draw();
  void print();
  Vec referenceVec;
private:
  Err ierr;
  bool isSymmetric_,firstAssembly_;
  int  nGlobalRows_,
    nGlobalColumns_,
    nLocalRows_,
    nLocalColumns_,
    globalLow_,
    globalHigh_,
    blockSize_,
    sqr_blockSize_;
  friend class Slammer;
  real* rowValueCache;
  int*  blockColumns;
  int*  colIndexCache;
  int cacheSize;
  Petsc::cc::_p_Mat* mat;
};
  
class Slammer
{
public:
  inline Slammer(Petsc::Mat *mat,const int& i,const int& j);
  inline void operator=(const real& v);
  inline void operator+=(const real& v);
  inline void operator-=(const real& v);
  inline ~Slammer();
private:
  Err ierr;
  const int &ir,&jr;
  Petsc::Mat *matp;
};

inline Slammer::Slammer(Petsc::Mat *mat,const int& i,const int& j):
  ir(i),
  jr(j),
  matp(mat)
{}

inline Slammer::~Slammer(){}
  
inline Slammer Mat::operator()(const int& i,const int& j)
{
  return Slammer(this,i,j);
}

inline void Slammer::operator=(const real& v)
{
  matp->rowValueCache[matp->cacheSize]=v;
  matp->colIndexCache[matp->cacheSize]=jr;
  ++(matp->cacheSize);
}

inline void Slammer::operator+=(const real& v)
{
  matp->rowValueCache[matp->cacheSize]=v;
  matp->colIndexCache[matp->cacheSize]=jr;
  ++(matp->cacheSize);
}

inline void Slammer::operator-=(const real& v)
{
  matp->rowValueCache[matp->cacheSize]=-v;
  matp->colIndexCache[matp->cacheSize]=jr;
  ++(matp->cacheSize);
}

}//Petsc
}//Daetk
#endif
