#ifndef ERROR_NORMS_TEST_H
#define ERROR_NORMS_TEST_H

#include "Definitions.h"
#include "WeightedRMSNorm.h"


namespace Daetk
{

class WeightedMaxNorm: public WeightedRMSNorm
{
  /***********************************************************************

    Calculate max_i |y_i|/(\eps_r|yw_i| + \eps_a)


  **********************************************************************/
public:
  WeightedMaxNorm(const int& Neq = 0);
  virtual ~WeightedMaxNorm();

  virtual real operator()(const Vec& y);

};//end Weighted MaxNorm

}//end Daetk
#endif
