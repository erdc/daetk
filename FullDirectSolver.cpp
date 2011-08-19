#include "FullDirectSolver.h"

namespace Daetk 
{

FullDirectSolver::FullDirectSolver(Mat& Min): 
  neq(Min.dim(0)),
  x(this,Min.dim(0)),
  b(this,Min.dim(0)),
  scaleFactor(Min.dim(0)),
  M(&Min)
{
  permutationVector=new int[neq];
}

FullDirectSolver::~FullDirectSolver()
{
  delete [] permutationVector;
}

bool FullDirectSolver::prepare()
{
  scaleFactor=0.0;
  real pivot,temp;
  int i,j,k,pivotIndex(0);
  for (i=0;i<neq;i++)
    {
      for(j=0;j<neq;j++)
	{
	temp = ::fabs((*M)(i,j));
	if (temp > scaleFactor(i))
	  scaleFactor(i)=temp;
	}
      if (scaleFactor(i)==0.0) 
	  return true;//matrix is singular
    }
  for (i=0;i<neq;i++)
    {
      pivot=0.0;
      for (j=i;j<neq;j++)
	{
	  temp=(*M)(j,i);
	  if (::fabs(temp/scaleFactor(j)) > ::fabs(pivot))
	    {
	      pivot=temp;
	      pivotIndex=j; 
	    }
	}
      if (i!=pivotIndex)
	{
	  for (k=0;k<neq;k++) 
	    {
	      temp=(*M)(i,k);
	      (*M)(i,k)=(*M)(pivotIndex,k);
	      (*M)(pivotIndex,k)=temp;   
	    }
	}      
      permutationVector[i]=pivotIndex;
      scaleFactor(pivotIndex)=scaleFactor(i);
      for(j=i+1;j<neq;j++)
	{
	  (*M)(j,i)=(*M)(j,i)/(*M)(i,i); 
	  for (k=i+1;k<neq;k++) 
	    {
	      (*M)(j,k)-=(*M)(j,i)*(*M)(i,k);
	    }
	}
    }
  return false;
}

bool FullDirectSolver::solve(const Vec& bIn,Vec& xIn)
{
  int i,j,indexOfFirstNonzeroInb;
  real temp;
  b.attachToTarget(bIn);
  x.attachToTarget(xIn);

  indexOfFirstNonzeroInb=-1; //exploit initial zeroes
  for (i=0;i<neq;i++)         //for rows $i=1,...,n:
    {
      temp=b.v_(permutationVector[i]);
      b.v_(permutationVector[i])=b.v_(i);
      if (indexOfFirstNonzeroInb!=-1)
	{
	  for (j=indexOfFirstNonzeroInb;j<i;j++)
	    {
	      temp-=(*M)(i,j)*x.v_(j);
	    }
	}
      else if (temp!=0)
	indexOfFirstNonzeroInb=i; 
      x.v_(i)=temp;              
    }
  for (i=neq-1;i>=0;i--)
    {
      temp=x.v_(i);
      for (j=i+1;j<neq;j++) 
	{
	  temp-=(*M)(i,j)*x.v_(j);
	}
      x.v_(i)=temp/(*M)(i,i);
    }
  x.restoreToTarget();
  return false;
}

}//Daetk
