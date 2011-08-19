#ifndef LINEAROPERATOR_H
#define LINEAROPERATOR_H

#include "Definitions.h" 
#include "Vec.h"

namespace Daetk 
{
class LinearOperator
{
public:
  LinearOperator(unsigned int domain, unsigned int range);
  void newsize(unsigned int domain, unsigned int range);
  virtual ~LinearOperator();
  virtual bool apply(const Vec& x, Vec& Ax)=0;
  int dimDomain(){return dimDomain_;}
  int dimRange(){return dimRange_;}
  inline void finalizeRow(int){}
  inline void finalizeBlockRow(int){}
  inline void finalizeAddRow(int){}
  inline void finalizeAddBlockRow(int){}
protected:
  unsigned int dimDomain_, dimRange_;
};
}//Daetk
#endif
