#ifndef DENSITY_H
#define DENSITY_H

#include "Definitions.h"
#include "Vec.h"

namespace Daetk
{
class Density
{
public:
    virtual void setHead(const real& psi)=0;
    virtual void setHead(const Vec& p, Vec& r, Vec& Dr)=0;
    virtual void setHeadD(const Vec& p, Vec& r, Vec& Dr, Vec& DDr)=0;
    virtual void setPressure(const real& p)=0;
    virtual const real& getRho()=0;
    virtual const real& getDrho()=0;
    virtual const real& getDDrho()=0;
};
}//Daetk
#endif
