#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

#include "Definitions.h"
#include "Vec.h"

namespace Daetk 
{
class Preconditioner
{
public:
  virtual ~Preconditioner();
  virtual bool prepare()=0;
  virtual bool apply(const Vec& x, Vec& Mx)=0;
};
}
#endif
