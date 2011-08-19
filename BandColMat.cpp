#include "BandColMat.h"

namespace Daetk 
{

using std::ostream;
using std::min;
using std::max;
using std::endl;
  
BandColMat::BandColMat():LinearOperator(0,0),M(),pivots(0),ku(0),kl(0),neq(0),rowOffSet(0){}

BandColMat::~BandColMat(){delete [] pivots;}

BandColMat& BandColMat::operator=(const BandColMat& C) { M=C.M; return *this;}

void BandColMat::beginAssembly(){}
 
void BandColMat::endAssembly(){}

BandColMat& BandColMat::operator=(const real& C) 
{ 
  for(unsigned int i=0;i<M.dim(0);i++)
    {
      for(unsigned int j=0;j<M.dim(1);j++)
        {
          M(i,j) = C; 
        }
    }
  return *this;
}
  
void BandColMat::zeroRow(int i)
  {
    for (int j=std::max(i-kl,0);j<=std::min(neq-1,i+ku);j++)
      (*this)(i,j) = 0.0;
  }

void BandColMat::zeroAll()
{
  for(unsigned int i=0;i<M.dim(0);i++)
    {
      for(unsigned int j=0;j<M.dim(1);j++)
        {
          M(i,j) = 0.0; 
        }
    }
}
  
bool BandColMat::apply(const Vec& x, Vec& Ax){ matVec(x,Ax); return false;}
  
real* BandColMat::castToArray(){return M.castToArray();}

const real* BandColMat::castToConstArray() const
{return M.castToConstArray();}

BandColMat::BandColMat(int upperDiagonals,int lowerDiagonals,int Neq):
  LinearOperator(Neq,Neq),
  M(2*lowerDiagonals+upperDiagonals+1,Neq,0.0),
  ku(upperDiagonals),
  kl(lowerDiagonals),
  neq(Neq),
  rowOffSet(kl-1)
{
  pivots = new int[Neq];
}

void BandColMat::newsize(int lowerDiagonals,int upperDiagonals,int Neq)
{
  ku=upperDiagonals;
  kl=lowerDiagonals;
  rowOffSet=kl-1;
  neq=Neq;
  M.newsize(2*kl+ku+1,neq);
  delete [] pivots;
  pivots = new int[neq];
}


ostream& operator<<(ostream& str,BandColMat& A)
{
  for (int i=0;i<A.neq;i++)
    {
      str<<"row "<<i<<": ";
      int bl(max(0,i-A.kl)),bu(min(A.neq-1,i+A.ku));
      for (int j=bl;j<=bu;j++)
        {
          if (fabs(A.M(A.rowOffSet+A.ku+1+i-j,j)) !=0.0)
            {
              str<<j<<" "<<A.M(A.rowOffSet+A.ku+1+i-j,j)<<"  ";
            }
        }
      str<<endl;
    }

  return str;
}
 
void BandColMat::matVec(const Vec& x, Vec& Ax)
{
  Ax=0.0;
  for (int i=0;i<neq;i++)
    {
      int bl(max(0,i-kl)),bu(min(neq-1,i+ku));
      for (int j=bl;j<=bu;j++)
        {
          if (M(rowOffSet+ku+1+i-j,j) !=0)
            Ax(i)+=M(rowOffSet+ku+1+i-j,j)*x(j);
        }
    }
}

}//Daetk
