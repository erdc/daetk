#ifndef GLOBALDAEDEF_H
#define GLOBALDAEDEF_H

#include "Definitions.h"
#include "DaeDefinition.h"

namespace Daetk 
{

typedef bool (*ResFunction)(const real& t,const Vec& y,const Vec& yp, Vec& r);

typedef bool (*YprimeFunction)(const real& t,const Vec& y,Vec& yp);

typedef  bool (*IndexResFunction)(const int& node,const real& t, const Vec& y, 
                  const Vec& yp, real& resAtNode);

class GlobalDaeDef : public DaeDefinition 
{
public:
  GlobalDaeDef();

  GlobalDaeDef(ResFunction residualIn, const real& t0, const Vec& y0,
               const Vec& y0prime,DataCollector* dataIn);

  GlobalDaeDef(ResFunction residualIn,YprimeFunction yPrimeValueIn,
               const real& t0, const Vec& y0,const Vec& y0prime,DataCollector* dataIn);

  GlobalDaeDef(ResFunction residualIn, IndexResFunction indexResIn,
               const real& t0, const Vec& y0,const Vec& y0prime,DataCollector* dataIn);

  virtual ~GlobalDaeDef();

  bool residual(const real& t,const Vec& y,const Vec& yp, Vec& r);
  bool yPrimeValue(const real& t,const Vec& y,Vec& yp);
  const real& getT0();
  const Vec& getY0();
  const Vec& getY0prime();
private:
  const real *t0;
  const Vec *y0,*y0prime;
  ResFunction resGlobal;
  YprimeFunction yPrimeGlobal;
  IndexResFunction indexResGlobal;
};

}//Daetk
#endif
