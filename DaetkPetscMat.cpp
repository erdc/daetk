#include "DaetkPetscMat.h"
#include <cstdio>
namespace Daetk 
{
namespace Petsc
{
  namespace cc
  {
    extern "C"
    {
       //mwf 090104 seems to be some voodoo with order of includes
       //mwf from petsc worked with petscda.h in between petscvec
       //mwf and petcmat
#ifndef DAETK_DEF_CPLUSPLUS_FOR_PETSC_H
#undef __cplusplus
#endif
#include "petsc.h"
#include "petscvec.h"
#include "petscdm.h"
#include "petscmat.h"
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
  
  cc::_p_Mat** Mat::castToPetscLValue(){return &mat;}
  cc::_p_Mat* Mat::castToPetsc(){return mat;}
  const cc::_p_Mat* Mat::castToConstPetsc() const {return mat;}  

void Mat::copy(Mat& mIn)
{
  using namespace cc;
  MatCopy(mIn.castToPetsc(),mat,SAME_NONZERO_PATTERN);
}
Mat::Mat(int neq):
  LinearOperator(neq,neq),
  referenceVec(neq),
  isSymmetric_(false), 
  firstAssembly_(true),
  nLocalRows_(0),
  nLocalColumns_(0),
  globalLow_(0),
  globalHigh_(0),
  blockSize_(1),
  cacheSize(0),
  mat(0)
{
}

Mat::~Mat()
{
  ierr = cc::MatDestroy(&mat);
}

Mat::Mat(int nGlobalRows, int nDiagonals, int nOffDiagonals, 
                   bool isSymmetric,int blockSize):
  LinearOperator(nGlobalRows,nGlobalRows),
  referenceVec(nGlobalRows),
  isSymmetric_(isSymmetric), 
  firstAssembly_(true),
  nGlobalRows_(nGlobalRows),
  nGlobalColumns_(nGlobalRows),
  nLocalRows_(0),
  nLocalColumns_(0),
  globalLow_(0),
  globalHigh_(0),
  blockSize_(blockSize),
  cacheSize(0),
  mat(0)
{
  newsize(nGlobalRows,nDiagonals,nOffDiagonals,isSymmetric,blockSize);
}

Mat::Mat(const Mat& mIn):
  LinearOperator(mIn.nGlobalRows_,mIn.nGlobalColumns_),
  referenceVec(mIn.nGlobalRows_),
  isSymmetric_(mIn.isSymmetric_), 
  firstAssembly_(mIn.firstAssembly_),
  nGlobalRows_(mIn.nGlobalRows_),
  nGlobalColumns_(mIn.nGlobalRows_),
  nLocalRows_(mIn.nLocalRows_),
  nLocalColumns_(mIn.nLocalColumns_),
  globalLow_(mIn.globalLow_),
  globalHigh_(mIn.globalHigh_),
  blockSize_(mIn.blockSize_),
  cacheSize(0),
  mat(0)
{
  format=mIn.format;
  ierr = cc::MatDuplicate(mIn.mat,cc::MAT_COPY_VALUES,&mat);
}

void Mat::newsize(int nGlobalRows, int nDiagonals, int nOffDiagonals, 
                         bool isSymmetric,int blockSize)
{
  LinearOperator::newsize(nGlobalRows,nGlobalRows);
  referenceVec.newsize(nGlobalRows);
  dimDomain_=dimRange_=nGlobalRows;
  isSymmetric_=isSymmetric; 
  firstAssembly_=true;
  nGlobalRows_=nGlobalRows;
  nGlobalColumns_=nGlobalRows;
  nLocalRows_=0;
  nLocalColumns_=0;
  globalLow_=0;
  globalHigh_=0;
  blockSize_=blockSize;
  cacheSize=0;
  mat=0;
  
  Vec testVec(nGlobalRows);
  //default format
  format=Mat::aij_format;
  if (blockSize==1)
    {
      format=Mat::aij_format;
      ierr =  cc::MatCreate(cc::PETSC_COMM_WORLD,&mat);
      ierr =  cc::MatSetSizes(mat,
			     testVec.ldim(),testVec.ldim(),
			     nGlobalRows,nGlobalRows);
      ierr =  cc::MatSetType(mat,MATMPIAIJ);
      ierr =  cc::MatSetUp(mat);
      // ierr =  cc::MatCreateMPIAIJ(cc::PETSC_COMM_WORLD,testVec.ldim(),testVec.ldim(),
      //                         nGlobalRows,nGlobalRows,
      // 				  nDiagonals,
      //                         nOffDiagonals,PETSC_NULL,,PETSC_NULL,
      //                         &mat);
      rowValueCache = new real[nDiagonals + nOffDiagonals];
      colIndexCache = new int[nDiagonals + nOffDiagonals];
    }
  else
    {
//        ierr =  cc::MatCreateMPIAIJ(cc::PETSC_COMM_WORLD,testVec.ldim(),testVec.ldim(),
//                                nGlobalRows,nGlobalRows,
//                                blockSize*nDiagonals,PETSC_NULL,
//                                blockSize*nOffDiagonals,PETSC_NULL,
//                                &mat);
//        format=Mat::aij_format;
      // ierr =  cc::MatCreateMPIBAIJ(cc::PETSC_COMM_WORLD,blockSize,
      //                          testVec.ldim(),testVec.ldim(),
      //                          nGlobalRows,nGlobalRows,
      //                          nDiagonals,PETSC_NULL,
      //                          nOffDiagonals,PETSC_NULL,
      //                          &mat);
      format=Mat::block_aij_format;
      ierr =  cc::MatCreate(cc::PETSC_COMM_WORLD,&mat);
      ierr =  cc::MatSetSizes(mat,
			     testVec.ldim(),testVec.ldim(),
			     nGlobalRows,nGlobalRows);
      ierr =  cc::MatSetType(mat,MATMPIBAIJ);
      ierr =  cc::MatSetUp(mat);
      sqr_blockSize_=blockSize_*blockSize_;
      blockColumns  = new int[blockSize_*(nDiagonals + nOffDiagonals)];
      rowValueCache = new real[sqr_blockSize_*blockSize*(nDiagonals + nOffDiagonals)];
      colIndexCache = new int[sqr_blockSize_*blockSize*(nDiagonals + nOffDiagonals)];
    }
  //    ierr = MatSetOption(mat,MAT_COLUMNS_SORTED);
  getStorageInfo();
  using namespace cc;
  ierr = MatSetOption(mat,MAT_ROW_ORIENTED,PETSC_FALSE);
}

void Mat::newsizeSequential(int nGlobalRows, int nDiagonals, int nOffDiagonals, 
                    bool isSymmetric,int blockSize)
{
  LinearOperator::newsize(nGlobalRows,nGlobalRows);
  referenceVec.newsizeSerial(nGlobalRows);
  dimDomain_=dimRange_=nGlobalRows;
  isSymmetric_=isSymmetric; 
  firstAssembly_=true;
  nGlobalRows_=nGlobalRows;
  nGlobalColumns_=nGlobalRows;
  nLocalRows_=0;
  nLocalColumns_=0;
  globalLow_=0;
  globalHigh_=0;
  blockSize_=blockSize;
  cacheSize=0;
  mat=0;
  
  Vec testVec;testVec.newsizeSerial(nGlobalRows);
  //default format
  format=Mat::aij_format;   

  int nNonzero=nDiagonals+nOffDiagonals;
  if (blockSize==1)
    {
      //mwf gcc 3.3 has a problem with resolving cc:: and macros
      using namespace cc;
      //ierr = cc::MatCreateSeqAIJ(PETSC_COMM_SELF,nGlobalRows,nGlobalRows,nNonzero,PETSC_NULL,&mat); 
      format=Mat::aij_format;
      ierr =  cc::MatCreate(cc::PETSC_COMM_WORLD,&mat);
      ierr =  cc::MatSetSizes(mat,
			     testVec.ldim(),testVec.ldim(),
			     nGlobalRows,nGlobalRows);
      ierr =  cc::MatSetType(mat,MATMPIAIJ);
      ierr =  cc::MatSetUp(mat);      
      rowValueCache = new real[nDiagonals + nOffDiagonals];
      colIndexCache = new int[nDiagonals + nOffDiagonals];
    }
  else
    {
      //mwf gcc 3.3 has a problem with resolving cc:: and macros
      using namespace cc;
      // ierr =  cc::MatCreateSeqBAIJ(PETSC_COMM_SELF,blockSize,
      //                          nGlobalRows,nGlobalRows,
      //                          nDiagonals+nOffDiagonals,
      //                          PETSC_NULL, &mat);
      format=Mat::block_aij_format;
      ierr =  cc::MatCreate(cc::PETSC_COMM_WORLD,&mat);
      ierr =  cc::MatSetSizes(mat,
			     testVec.ldim(),testVec.ldim(),
			     nGlobalRows,nGlobalRows);
      ierr =  cc::MatSetType(mat,MATMPIBAIJ);
      ierr =  cc::MatSetUp(mat);      
      sqr_blockSize_=blockSize_*blockSize_;
      blockColumns  = new int[blockSize_*(nDiagonals + nOffDiagonals)];
      rowValueCache = new real[sqr_blockSize_*blockSize*(nDiagonals + nOffDiagonals)];
      colIndexCache = new int[sqr_blockSize_*blockSize*(nDiagonals + nOffDiagonals)];

//        ierr = cc::MatCreateSeqAIJ(cc::PETSC_COMM_SELF,nGlobalRows,nGlobalRows,blockSize*nNonzero,PETSC_NULL,&mat); 
//        format=Mat::aij_format;
//        rowValueCache = new real[nDiagonals + nOffDiagonals];
//        colIndexCache = new int[nDiagonals + nOffDiagonals];
//        ierr = cc::MatCreateSeqBAIJ(cc::PETSC_COMM_SELF,blockSize_,nGlobalRows,nGlobalRows,nNonzero,PETSC_NULL,&mat);
//        format=Mat::block_aij_format;
//        sqr_blockSize_ = blockSize_*blockSize_;
//        rowValueCache = new real[sqr_blockSize_*(nDiagonals + nOffDiagonals)];
//        colIndexCache = new int[sqr_blockSize_*(nDiagonals + nOffDiagonals)];
//        bockColumns = new int[blockSize_*(nDiagonals + nOffDiagonals)];
    }
  getStorageInfo();
  using namespace cc;
  ierr = MatSetOption(mat,MAT_ROW_ORIENTED,PETSC_FALSE);
  //    ierr = MatSetOption(mat,MAT_COLUMNS_SORTED);
}

void Mat::attachToMat(cc::_p_Mat* mIn)
{
  mat = mIn;
  ierr =  cc::MatGetLocalSize(mat,&nLocalRows_,&nLocalColumns_);
  ierr =  cc::MatGetSize(mat,&nGlobalRows_,&nGlobalColumns_);
  ierr =  cc::MatGetOwnershipRange(mat,&globalLow_,&globalHigh_);
  LinearOperator::newsize(nGlobalRows_,nGlobalRows_);
  referenceVec.newsize(nGlobalRows_);
}

void Mat::beginAssembly()
{
  if (isSymmetric_)
    ierr = cc::MatSetOption(mat,cc::MAT_SYMMETRIC,cc::PETSC_TRUE);
  ierr = cc::MatAssemblyBegin(mat,cc::MAT_FINAL_ASSEMBLY);
}

void Mat::endAssembly()
{
  ierr = cc::MatAssemblyEnd(mat,cc::MAT_FINAL_ASSEMBLY);

  if (firstAssembly_)
    {
      firstAssembly_=false;
      ierr = cc::MatSetOption(mat,cc::MAT_NEW_NONZERO_LOCATION_ERR,cc::PETSC_TRUE);
      if (format==Mat::block_aij_format)
        ierr = cc::MatSetOption(mat,cc::MAT_USE_HASH_TABLE,cc::PETSC_TRUE);
    }
}

int Mat::getLocalDim(RowsOrColumns r)
{
  if (r==ROWS)
    return nLocalRows_;
  else
    return nLocalColumns_;
}

int Mat::getGlobalHigh()
{
  return globalHigh_;
}

int Mat::getGlobalLow()
{
  return globalLow_;
}

int Mat::getBlockSize()
{
  return blockSize_;
}

void Mat::getStorageInfo()
{  
  ierr = cc::MatGetLocalSize(mat,&nLocalRows_,&nLocalColumns_);
  ierr = cc::MatGetOwnershipRange(mat,&globalLow_,&globalHigh_);
  ierr = cc::MatGetSize(mat,&nGlobalRows_,&nGlobalColumns_);
  ierr = cc::MatGetBlockSize(mat,&blockSize_);
}

bool Mat::apply(const Vec& x, Vec& Ax)
{
  int i;
  x.restoreLocal();
  Ax.restoreLocal();
  i = cc::MatMult(mat,const_cast<cc::_p_Vec*>(x.castToConstPetsc()),Ax.castToPetsc()); 
  x.getLocal();
  Ax.getLocal();
  ierr = i;
  return i;
}

ostream& operator<<(ostream& s, const Mat& m)
{
  using namespace cc;
  Err ierr;
  ierr = MatView(m.mat,PETSC_VIEWER_STDOUT_WORLD);
  return s;
}

void  Mat::scaleDiagonal(const Vec& left, const Vec& right)
{
  ierr = cc::MatDiagonalScale(mat,const_cast<cc::_p_Vec*>(left.castToConstPetsc()),
                          const_cast<cc::_p_Vec*>(right.castToConstPetsc()));
}

void Mat::scale(real a)
{
  //mwf 2.3.1 changed the order of args
  ierr = cc::MatScale(mat,a);
}

void Mat::setInsertionOrder(Insertion s)
{
  switch (s)
    {
    case ROW_ORIENTED:
      ierr = cc::MatSetOption(mat,cc::MAT_ROW_ORIENTED,cc::PETSC_TRUE);
      break;
    case COLUMN_ORIENTED:
      ierr = cc::MatSetOption(mat,cc::MAT_ROW_ORIENTED,cc::PETSC_FALSE);
      break;
    // case ROWS_SORTED:
    //   ierr = cc::MatSetOption(mat,cc::MAT_ROWS_SORTED); 
    //   break;
    // case ROWS_UNSORTED:
    //   ierr = cc::MatSetOption(mat,cc::MAT_ROWS_UNSORTED); 
    //   break;
    // case COLUMNS_SORTED:
    //   ierr = cc::MatSetOption(mat,cc::MAT_COLUMNS_SORTED); 
    //   break;
    // case COLUMNS_UNSORTED:
    //   ierr = cc::MatSetOption(mat,cc::MAT_COLUMNS_UNSORTED); 
    //   break;
    }
}

void Mat::zeroAll()
{
  ierr = cc::MatZeroEntries(mat);
}

void Mat::zeroRow(int i)
{
  real diag=0.0;
  cc::IS is;
  int idx[1];
  idx[0] = i;
  //mwf gcc 3.3 has some problems with resolving cc::
  //for macros
  using namespace cc;
  ierr = cc::ISCreateGeneral(PETSC_COMM_SELF,1,idx,PETSC_COPY_VALUES,&is);
  //mwf this is now MatZeroRowsIS ?
  ierr = cc::MatZeroRowsIS(mat,is,diag,(Daetk::Petsc::cc::Vec)(PETSC_NULL),(Daetk::Petsc::cc::Vec)(PETSC_NULL));  
  ierr = cc::ISDestroy(&is);
}

void Mat::draw()
{
  Err ierr;
  using namespace cc;
  PetscViewer v;
  ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,
                             reinterpret_cast<const char*>(NULL),
                             reinterpret_cast<const char*>(NULL),
                             PETSC_DECIDE,PETSC_DECIDE,
                             PETSC_DECIDE,PETSC_DECIDE,
                             &v);
  ierr = MatView(mat,v);
}
void Mat::print()
{
  Err ierr;
  using namespace cc;
  ierr = MatView(mat,PETSC_VIEWER_STDOUT_WORLD);
}

bool Mat::isAssembled()
{
  using namespace cc;
  PetscBool t;
  ierr = MatAssembled(mat,&t);
  return t;
}

bool Mat::isValid()
{
  using namespace cc;
  PetscBool t;
  //cek hack
  //ierr = MatValid(mat,&t);
  return t;
}

//  ierr = MatBDiagGetData(Mat mat,int *nd,int *bs,int **diag,int **bdlen,Scalar ***diagv)
//  ierr = cc::MatCreateShell(MPI_Comm comm,int m,int n,int M,int N,void *ctx,Mat *A)
//  ierr = MatShellSetOperation(Mat mat,MatOperation op,void *f)
  
//this is going to be really slow. will fix later
const real& Mat::get(const int& i, const int& j)
{
  static real rvalue;
  int row[1],column[1];
  row[0]=i;
  column[0]=j;
  ierr = cc::MatGetValues(mat,1,row,1,column,&rvalue);
  return rvalue;
}

void Mat::finalizeRow(int i)
{
  using namespace cc;
  ierr = MatSetValues(mat,1,&i,cacheSize,colIndexCache,
                      rowValueCache,INSERT_VALUES);
  cacheSize=0;
}

void Mat::finalizeAddRow(int i)
{
  using namespace cc;
  ierr = MatSetValues(mat,1,&i,cacheSize,colIndexCache,
                      rowValueCache,ADD_VALUES);
  cacheSize=0;
}

void Mat::finalizeBlockRow(int i)
{
  using namespace cc;
  int blocks=cacheSize/sqr_blockSize_;
  for(int bc=0;bc<blocks;bc++)
    blockColumns[bc]=colIndexCache[bc*sqr_blockSize_]/blockSize_;
  ierr = MatSetValuesBlocked(mat,1,&i,blocks,blockColumns,
                      rowValueCache,INSERT_VALUES);
  cacheSize=0;
}

void Mat::finalizeAddBlockRow(int i)
{
  using namespace cc;
  int blocks=cacheSize/sqr_blockSize_;
  for(int bc=0;bc<blocks;bc++)
    blockColumns[bc]=colIndexCache[bc*sqr_blockSize_]/blockSize_;
  ierr = MatSetValuesBlocked(mat,1,&i,blocks,blockColumns,
                      rowValueCache,ADD_VALUES);
  cacheSize=0;
}

//  void Slammer::operator=(const real& v)
//  { 
//    using namespace cc;
//    ierr = MatSetValues(matp,1,const_cast<int*>(&ir),1,const_cast<int*>(&jr),
//                         const_cast<real*>(&v),INSERT_VALUES);
//  }

//  void Slammer::operator+=(real v)
//  { 
//    using namespace cc;
//    ierr = MatSetValues(matp->mat,1,const_cast<int*>(&ir),1,const_cast<int*>(&jr),
//                         const_cast<real*>(&v),ADD_VALUES);
//  }

//  void Slammer::operator-=(real v)
//  { 
//    using namespace cc;
//    v*=-1.0;
//    ierr = MatSetValues(matp,1,const_cast<int*>(&ir),1,const_cast<int*>(&jr),
//                         const_cast<real*>(&v),ADD_VALUES);
//  }

}//Petsc
}//Daetk
