#include "WeightedRMSNorm.h"

namespace Daetk 
{

using std::cerr;
using std::endl;
using std::cout;
using std::flush;

const Vec&  WeightedRMSNorm::getScaling(){return scaling;}

real* WeightedRMSNorm::getWeightBegin(){return weight.begin();}
  
real* WeightedRMSNorm::getWeightEnd(){return weight.end();}
  
void WeightedRMSNorm::scale(const Vec& x,Vec& y)
  {
    int ldim=x.ldim();
    for (i=0;i<ldim;i++)
      y[i]=x[i]*scaling[i];
//      y=x;
  }
  
void WeightedRMSNorm::deScale(const Vec& x, Vec& y)
  {
    int ldim=x.ldim();
    for (i=0;i<ldim;i++)
      y[i]=x[i]/scaling[i];
//      y=x;
  }
  
void WeightedRMSNorm::print(){std::cout<<"weight"<<weight<<"atol"<<atol<<"rtol"<<rtol<<std::endl;}

WeightedRMSNorm::WeightedRMSNorm(const int& Neq):
  neq(Neq),
  rneq(real(Neq)),
  oneOver_sqrt_neq(1.0/sqrt(real(Neq))),
  weight(Neq),
  scaling(Neq),
  atol(Neq),
  rtol(Neq),
  tmp(Neq)
{
  Tracer tr("WeightedRMSNorm::WeightedRMSNorm(const int& Neq)");
}
  
WeightedRMSNorm::WeightedRMSNorm(WeightedRMSNorm& WN):
  neq(WN.neq),
  rneq(WN.rneq),
  oneOver_sqrt_neq(WN.oneOver_sqrt_neq),
  weight(WN.weight),
  scaling(WN.weight),
  atol(WN.atol),
  rtol(WN.rtol),
  tmp(WN.tmp)
{
  Tracer tr("WeightedRMSNorm::WeightedRMSNorm(WeightedRMSNorm& WN)");
}
 
WeightedRMSNorm::~WeightedRMSNorm()
{
  Tracer tr("WeightedRMSNorm::~WeightedRMSNorm()");
}

void WeightedRMSNorm::setWeight(const Vec& y)
{
  int ldim=weight.ldim();
  real r;
  for (int i=0;i<ldim;i++)
    {
      weight[i]=atol[i] + fabs(y[i])*rtol[i];
      if (weight[i])
        weight[i]=1.0/weight[i];
      else
        {
          cerr<<"Encountered a zero while calculating weight"<<weight<<endl
              <<"Set atol to a nonzero value if the solution is not"<<endl
              <<"bounded away from zero"<<endl;
          exit(1);
        }
    }
  r = (*this)(y)*100.0*MACHINE_EPSILON; //norm of y
  if (r > 1.0)
    {
      cerr<<"The requested tolerance exceeds machine precision"<<endl;
      cerr<<"Altering weightvector to appropriate precision"<<endl;
      for (int i=0;i<ldim;i++)
        {
          weight[i]=1.0/(r*atol[i] + r*fabs(y[i])*rtol[i]);
        }
    }
//   scaling = oneOver_sqrt_neq*weight;
  scaling = weight;
  scal(oneOver_sqrt_neq,weight);
}

void WeightedRMSNorm::setTolerances(const Vec& atolIn,const Vec& rtolIn)
{
  atol = atolIn;
  rtol = rtolIn;
}

void WeightedRMSNorm::setTolerances(const real& atolR,const real& rtolR)
{
  atol = atolR;
  rtol = rtolR;
}

real WeightedRMSNorm::operator()(const Vec& y)
{
  int ldim=weight.ldim();
  for (int i=0;i<ldim;i++)
    {
      tmp[i] = y[i]*weight[i];
      tmp[i]*=tmp[i];
    }
  return sqrt(tmp.sum()/rneq);
}

  //new weighed L2 norm
  
  WeightedL2Norm::WeightedL2Norm(const int& Neq):
    WeightedRMSNorm(Neq),
    dV(Neq)
  {
    oneOver_sqrt_neq = 1.0; //cek just in case this gest used somewhere
    Tracer tr("WeightedL2Norm::WeightedL2Norm(const int& Neq)");
  }
  void WeightedL2Norm::setTolerances(const real& atolR,const real& rtolR, const Vec& dV_in)
  {
    atol = atolR;
    rtol = rtolR;
    dV = dV_in;//deep copy
    V = dV.sum();
  }
  void WeightedL2Norm::setWeight(const Vec& y)
  {
    int ldim=weight.ldim();
    real r;
    for (int i=0;i<ldim;i++)
      {
        weight[i]=atol[i] + fabs(y[i])*rtol[i];
        if (weight[i])
          weight[i]=1.0/weight[i];
        else
          {
            cerr<<"Encountered a zero while calculating weight"<<weight<<endl
                <<"Set atol to a nonzero value if the solution is not"<<endl
                <<"bounded away from zero"<<endl;
            exit(1);
          }
      }
    r = (*this)(y)*100.0*MACHINE_EPSILON; //norm of y
    if (r > 1.0)
      {
        cerr<<"The requested tolerance exceeds machine precision"<<endl;
        cerr<<"Altering weightvector to appropriate precision"<<endl;
        for (int i=0;i<ldim;i++)
          {
            weight[i]=1.0/(r*atol[i] + r*fabs(y[i])*rtol[i]);
          }
      }
    scaling = weight;
    //don't scale by 1/sqrt(neq) for L2
  }

  real WeightedL2Norm::operator()(const Vec& y)
  {
    int ldim=weight.ldim();
    for (int i=0;i<ldim;i++)
      {
        tmp[i] = y[i]*weight[i];
        tmp[i] *= tmp[i];
        tmp[i] *= dV[i];
      }
    return sqrt(fabs(tmp.sum())/V);
    //return sqrt(fabs(tmp.sum()));
  }
  
}//Daetk
