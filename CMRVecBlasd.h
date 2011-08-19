#ifndef CMRVECBLASD_H
#define CMRVECBLASD_H

#include "Definitions.h"
#include "CMRVec.h"

namespace Daetk 
{

typedef CMRVec<double> CMRVecd;
///
void rot(CMRVecd& X, CMRVecd& Y, double& C, double& S);
/** Apply plane rotation */


///
void swap(CMRVecd& X, CMRVecd& Y);
/** X <-> Y */

///
void scal(const double& ALPHA, CMRVecd& X);
/** X <- alpha * X */

///
void copy(const CMRVecd& X, CMRVecd& Y);
/** Y <- X */

///
void axpy(const double& ALPHA, const CMRVecd& X, CMRVecd& Y);
/** Y <- alpha * X + Y*/

///
double dot(const CMRVecd& X,const CMRVecd& Y);
/** ddot <- X^T * Y */

///
double nrm2(const CMRVecd& X);
/** dnrm2 <- ||X||_2 */

///
double asum(const CMRVecd& X);
/** dasum <- ||re(X)||_1 + ||im(X)||_1 */

///
int imax(const CMRVecd& X);
/** idamax <- 1st k s.t. X(k) = max { X(i)| X.base() <= i < X.dim() }. */
}//Daetk
#endif


