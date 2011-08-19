#ifndef MAPPING_H
#define MAPPING_H

#include "Definitions.h"
#include "ParameterDatabase.h"

//some simple classes to just define (x,y) = F(xHat,yHat)
namespace Daetk
{

class IdentityMap
{
public:
  IdentityMap(ParameterDatabase& pd) {}
  ~IdentityMap() {}
  //map ref to phys
  real Fx(real x,real y)
    { return x; }
  real Fy(real x,real y)
    { return y; }
  //map phys to ref
  real FxInv(real x,real y)
    { return x; }
  real FyInv(real x,real y)
    { return y; }
  //deriv of ref to phys
  real dFxdX(real x,real y)
    { return 1.0; }
  real dFxdY(real x,real y)
    { return 0.0; }
  real dFydX(real x,real y)
    { return 0.0; }
  real dFydY(real x,real y)
    { return 1.0; }

};

//only does 45 deg for now
class RotationMap
{
public:
  RotationMap(ParameterDatabase& pd) {}
  ~RotationMap() {}
  //map ref to phys
  real Fx(real x,real y)
    { return (x-y)/sqrt(2.0);  }
  real Fy(real x,real y)
    { return (x+y)/sqrt(2.0); }
  //map phys to ref
  real FxInv(real x,real y)
    { return (x+y)/sqrt(2.0); }
  real FyInv(real x,real y)
    { return (y-x)/sqrt(2.0); }
  //deriv of ref to phys
  real dFxdX(real x,real y)
    { return 1.0/sqrt(2.0);  }
  real dFxdY(real x,real y)
    { return -1.0/sqrt(2.0); }
  real dFydX(real x,real y)
    { return 1.0/sqrt(2.0); }
  real dFydY(real x,real y)
    { return 1.0/sqrt(2.0);  }

};


//only does 45 deg for now
class PiecewiseLinear1DMap
{
public:
  PiecewiseLinear1DMap(ParameterDatabase& pd): 
    physX(pd.i("nxNodes")*pd.i("nyNodes")*pd.i("nzNodes")),
    refDx(pd.r("xRight") - pd.r("xLeft")/(real(pd.i("nxNodes"))-1.0))
  {
    //initialize to identity
    for (int i=0;i<physX.dim();i++)
      {
        physX(i) = i*refDx;
      }
  }
  ~PiecewiseLinear1DMap() {}
  Vec physX;
  real refDx;
  //map ref to phys
  real Fx(real x,real y)
  { 
    real theta,ir;
    int i;
    theta = std::modf(x/refDx,&ir);
    i=int(ir);
    i=std::min(i,physX.dim()-1);
    if (i== physX.dim()-1)
      return physX[i];
    else
      return theta*(physX[i+1]-physX[i]) + physX[i];
  }
  real Fy(real x,real y)
    { return y; }
  //map phys to ref
  real FxInv(real x,real y)
  {
    int i=0;
    while (x > physX[i+1] && (i+1) < (physX.dim()-1))
      i++;
    return (x - physX[i])/(physX[i+1] - physX[i])*refDx + i*refDx;
  }
  real FyInv(real x,real y)
    { return y; }
  //deriv of ref to phys
  real dFxdX(real x,real y)
  {     
    real theta,ir;
    int i;
    theta = std::modf(x/refDx,&ir);
    i=int(ir);
    i=std::min(i,physX.dim()-2);
    return (physX[i+1] - physX[i])/refDx ;
  }
  real dFxdY(real x,real y)
    { return 0.0; }
  real dFydX(real x,real y)
    { return 0.0; }
  real dFydY(real x,real y)
    { return 1.0;  }

};

class YotovMap
{
public:
  YotovMap(ParameterDatabase& pd):slopeConst(0.1),slopeDenom(M_PI/6.0) {}
  ~YotovMap() {}
  //map ref to phys
  real Fx(real x,real y)
    { return x; }
  real Fy(real x,real y)
    {   
      //this is y+r*sin(pi*x/D) (yotov)
      return y + slopeConst*sin(M_PI*x/slopeDenom);

    }
  //map phys to ref
  real FxInv(real x,real y)
    { return x; }
  real FyInv(real x,real y)
    {   
      //this is y+r*sin(pi*x/D)(yotov)
      return y - slopeConst*sin(M_PI*x/slopeDenom);
    }
  //deriv of ref to phys
  real dFxdX(real x,real y)
    { return 1.0; }
  real dFxdY(real x,real y)
    { return 0.0; }
  real dFydX(real x,real y)
    {
      //this is y+r*sin(pi*x/D) (yotov)
      return slopeConst*M_PI/slopeDenom*cos(M_PI*x/slopeDenom);
    }
  real dFydY(real x,real y)
    { return 1.0; }

  real slopeConst,slopeDenom;
};

class HillSlope1Map
{
public:
  //mwf may 7 was slopeConst(0.1), slopeDenom(25.0e1) 
  //mwf then tried slopeConst(0.5), slopeDenom(25.0e0) 
  HillSlope1Map(ParameterDatabase& pd):slopeConst(0.5),slopeDenom(25.0e0) {}
  ~HillSlope1Map() {}
  //map ref to phys
  real Fx(real x,real y)
    { return x; }
  real Fy(real x,real y)
    {   
      //this is y+r*x*sin(pi*x/D)
      return y + slopeConst*x*sin(M_PI*x/slopeDenom);
    }
  //map phys to ref
  real FxInv(real x,real y)
    { return x; }
  real FyInv(real x,real y)
    {   
      //this is y+r*x*sin(pi*x/D)
      return y - slopeConst*x*sin(M_PI*x/slopeDenom);
    }
  //deriv of ref to phys
  real dFxdX(real x,real y)
    { return 1.0; }
  real dFxdY(real x,real y)
    { return 0.0; }
  real dFydX(real x,real y)
    {
      //this is y+r*x*sin(pi*x/D)
      return slopeConst*sin(M_PI*x/slopeDenom) 
	+ slopeConst*M_PI/slopeDenom*x*cos(M_PI*x/slopeDenom);
    }
  real dFydY(real x,real y)
    { return 1.0; }

  real slopeConst,slopeDenom;
};

class HillSlope2Map
{
public:
  HillSlope2Map(ParameterDatabase& pd):slopeConst(0.1),slopeDenom(25.0e0),
    xCutOff(16.0e0),yOffSet(0.0)
    {
      yOffSet = slopeConst*xCutOff*sin(M_PI*xCutOff/slopeDenom);
    }
  ~HillSlope2Map() {}
  //map ref to phys
  real Fx(real x,real y)
    { return x; }
  real Fy(real x,real y)
    {   
	//this is y+r*x*sin(pi*x/D) for x <=16.0
      real eps(1.0e-10);
      if (x > xCutOff+eps)
	return y + yOffSet;
      else
	return y + slopeConst*x*sin(M_PI*x/slopeDenom);
    }
  //map phys to ref
  real FxInv(real x,real y)
    { return x; }
  real FyInv(real x,real y)
    {   
      //this is y+r*x*sin(pi*x/D) for x <=16.0
      real eps(1.0e-10);
      if (x > xCutOff+eps)
	return y - yOffSet;
      else
	return y - slopeConst*x*sin(M_PI*x/slopeDenom);
    }
  //deriv of ref to phys
  real dFxdX(real x,real y)
    { return 1.0; }
  real dFxdY(real x,real y)
    { return 0.0; }
  real dFydX(real x,real y)
    {
      //this is y+r*x*sin(pi*x/D) for x <=16.0
      real eps(1.0e-10);
      if (x > xCutOff+eps)
	return 0.0;
      else
	return slopeConst*sin(M_PI*x/slopeDenom) 
	  + slopeConst*M_PI/slopeDenom*x*cos(M_PI*x/slopeDenom);
    }
  real dFydY(real x,real y)
    { return 1.0; }

  real slopeConst,slopeDenom,xCutOff,yOffSet;

};

}//Daetk

#endif
