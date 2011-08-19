#include "DaeDefinition.h"

namespace Daetk 
{

using std::cerr;
using std::endl;

bool DaeDefinition::yPrimeValue(const real& t,const Vec& y,Vec& yp)
{/*not required*/ return true;}

bool DaeDefinition::evaluateDaeJacobian(const real& t,const Vec& y,const Vec& yp, 
                                const real& alpha)
{
  std::cerr<<"You must override evaluateDaeJacobian to use an analytical jacobian"<<std::endl;
  return true;
}

bool DaeDefinition::jacVec(const Vec& x, Vec& Jx)
{
  std::cerr<<"You must override jacVec to use an analytical jacobian with an iterative method"<<std::endl;
  
  return true;
}

void DaeDefinition::initializeFunction(const Vec& y, const Vec& yp, const Vec& f)
{
#ifndef USE_BLAS
  yDaeDef = y;
  ypDaeDef= yp;
  Fcurrent= f;
#else
  copy(y,yDaeDef);
  copy(yp,ypDaeDef);
  copy(f,Fcurrent);
#endif  
  updateF=false;
}

DaeDefinition::~DaeDefinition()
{
  Tracer tr("DaeDefinition:~DaeDefinition()");
}

DaeDefinition::DaeDefinition(int dim,DataCollector* dataIn):
    VectorFunction(dim,dim),
    updateF(true),
    updateJac(true),
    alphaDaeDef(0.0),
    tDaeDef(0.0),
    yDaeDef(dim,12345),
    ypDaeDef(dim,12345),
    Fcurrent(dim,12345),
    yLast(dim,12345),
    ypLast(dim,12345),
    Flast(dim,12345),
    del(dim,12345),
    betaDaeDef(dim,0.0),
    data(dataIn),
    //mwf added
    computeOwnJacobianDelta(false)
{
  Tracer tr("DaeDefinition::DaeDefinition(Data& dataIn)");
}


const Vec& DaeDefinition::argument()
{
  return yDaeDef;
}
  
const Vec& DaeDefinition::value(bool& evalError)
{
  evalError = false;
  if (updateF)
    {
      data->functionEvaluation();
      evalError = residual(tDaeDef,yDaeDef,ypDaeDef,Fcurrent);
      if (!evalError)
	updateF=false;
    }
  return Fcurrent;
}

void DaeDefinition::correctArgument(Vec& correction)
{
  yLast = yDaeDef;
  ypLast = ypDaeDef;
  Flast = Fcurrent;

  updateJac=true;
  updateF=true;
#ifndef USE_BLAS
  yDaeDef-=correction;
#else
  axpy(-1.0,correction,yDaeDef);
#endif
#ifndef USE_BLAS
  ypDaeDef-=alphaDaeDef*correction;
#else
  axpy(-alphaDaeDef,correction,ypDaeDef);
#endif
}

void DaeDefinition::resetFunction()
{
  updateF=true;
  updateJac=true;
}

void DaeDefinition::unCorrect()
{
  updateF = false;
  yDaeDef = yLast;
  ypDaeDef = ypLast;
  Fcurrent = Flast;
}

bool DaeDefinition::evaluateAnalyticalJacobian()
{
  return evaluateDaeJacobian(tDaeDef,yDaeDef,ypDaeDef,alphaDaeDef);
}

bool DaeDefinition::numericalJacVec(const Vec& v, Vec& Jv)
{
  bool evalError=false;
  real delFac=1.0e-5;
  value(evalError);//sets Fcurrent
  if (evalError)
    return evalError;
  else
    {
//       del = -delFac*v;
      del = v;
      scal(-delFac,del);
      correctArgument(del);//sets Flast to Fcurrent -- remebmer correction is -del = delFac*v
      value(evalError);
      while (evalError)
        {
          cerr<<"cutting back on del in numerical jac vec"<<endl;
          unCorrect();//sets Flast to Fcurrent
          delFac*=0.1;
//           del = -delFac*v;
          del = v;
          scal(-delFac,del);
          correctArgument(del);
          value(evalError);//sets Fcurrent
        }
//       Jv = (Fcurrent - Flast)/delFac;
      Jv = Fcurrent;
      axpy(-1.0,Flast,Jv);
      scal(1.0/delFac,Jv);
      unCorrect();
    }
  return evalError;
}

bool DaeDefinition::analyticalJacVec(const Vec& v, Vec& Jv)
{
  bool evalError = false;
  if (updateJac)
    {
      evalError = evaluateAnalyticalJacobian();
      if (evalError)
        return evalError;
      else
        updateJac=false;
    }
  jacVec(v,Jv);
  return evalError;
}
//mwf added to allow other integrators that don't compute delta
void DaeDefinition::computeDeltaForJacobian()
{

  if (computeOwnJacobianDelta)
    {
      deltaVF = 0.0;
      for (int i=0;i<yDaeDef.ldim();i++)
        {
          //mwf this is what it should be 
          deltaVF[i] = SQRT_MACHINE_EPSILON*fabs(yDaeDef[i])+SQRT_MACHINE_EPSILON;
          deltaVF[i]= yDaeDef[i] - (yDaeDef[i] - deltaVF[i]);
        }
    }
}

}//Daetk
