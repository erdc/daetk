#include "LinearOperator.h"

namespace Daetk
{

LinearOperator::LinearOperator(unsigned int domain, unsigned int range):
  dimDomain_(domain),
  dimRange_(range)
{}

void LinearOperator::newsize(unsigned int domain, unsigned int range){dimDomain_=domain; dimRange_=range;}

LinearOperator::~LinearOperator(){}

}//Daetk
