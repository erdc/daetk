#include "GlobalDaeDef.h"

namespace Daetk 
{

const real& GlobalDaeDef::getT0(){return *t0;}
  
const Vec& GlobalDaeDef::getY0(){return *y0;}

const Vec& GlobalDaeDef::getY0prime(){return *y0prime;}

bool GlobalDaeDef::residual(const real& t,const Vec& y,
				   const Vec& yp, Vec& r)
{
  return resGlobal(t,y,yp,r);
}

bool GlobalDaeDef::yPrimeValue(const real& t,const Vec& y,Vec& yp) 
{
  return yPrimeGlobal(t,y,yp);
}

GlobalDaeDef::GlobalDaeDef():
  DaeDefinition(0,0),
  t0(0),
  y0(0),
  y0prime(0),
  resGlobal(0),
  yPrimeGlobal(0),
  indexResGlobal(0)
{
  //default constructor called for vectors
  Tracer tr("GlobalDaeDef::GlobalDaeDef()");
}

GlobalDaeDef::GlobalDaeDef(ResFunction residualIn, const real& t0In, 
                           const Vec& y0In, const Vec& y0primeIn,DataCollector* dataIn): 
  DaeDefinition(y0In.dim(),dataIn),
  t0(&t0In),
  y0(&y0In),
  y0prime(&y0primeIn),
  resGlobal(residualIn),
  yPrimeGlobal(0),
  indexResGlobal(0)
{
  Tracer tr("GlobalDaeDef::GlobalDaeDef(ResFunction residualIn, const real& t0, const Vec& y0, const Vec& y0prime)");
}

GlobalDaeDef::GlobalDaeDef(ResFunction residualIn,
                           YprimeFunction yPrimeValueIn, const real& t0In, 
                           const Vec& y0In,const Vec& y0primeIn,DataCollector* dataIn):
  DaeDefinition(y0In.dim(),dataIn),
  t0(&t0In),
  y0(&y0In),
  y0prime(&y0primeIn),
  resGlobal(residualIn),
  yPrimeGlobal(yPrimeValueIn),
  indexResGlobal(0)
{
  Tracer tr("GlobalDaeDef::GlobalDaeDef(ResFunction residualIn, yPrimeFunction yPrimeValueIn, const real& t0, const Vec& y0,const Vec& y0prime)");
}

GlobalDaeDef::GlobalDaeDef(ResFunction residualIn, 
                           IndexResFunction indexResIn, const real& t0in, 
                           const Vec& y0in,const Vec& y0primein,DataCollector* dataIn):
  DaeDefinition(y0in.dim(),dataIn),
  t0(&t0in),
  y0(&y0in),
  y0prime(&y0primein),
  resGlobal(residualIn),
  indexResGlobal(indexResIn)
{
  Tracer tr("GlobalDaeDef::GlobalDaeDef(ResFunction residualIn, IndexResFunction indexResIn, const real& t0in, const Vec& y0in,const Vec& y0primein)");
}

GlobalDaeDef::~GlobalDaeDef()
{
  Tracer tr("~GlobalDaeDef::GlobalDaeDef()");
}

}//Daetk
