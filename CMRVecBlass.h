#ifndef CMRVECBLASS_H
#define CMRVECBLASS_H
#include "Definitions.h"
#include "CMRVec.h"
namespace Daetk 
{

typedef CMRVec<float> CMRVecs;
///
void rot(CMRVecs& X, CMRVecs& Y, float& C, float& S);
/** Apply plane rotation */


///
void swap(CMRVecs& X, CMRVecs& Y);
/** X <-> Y */

///
void scal(const float& ALPHA, CMRVecs& X);
/** X <- alpha * X */

///
void copy(const CMRVecs& X, CMRVecs& Y);
/** Y <- X */

///
void axpy(const float& ALPHA, const CMRVecs& X, CMRVecs& Y);
/** Y <- alpha * X + Y*/

///
float dot(const CMRVecs& X,const CMRVecs& Y);
/** ddot <- X^T * Y */

///
float nrm2(const CMRVecs& X);
/** dnrm2 <- ||X||_2 */

///
float asum(const CMRVecs& X);
/** dasum <- ||re(X)||_1 + ||im(X)||_1 */

///
int imax(const CMRVecs& X);
/** idamax <- 1st k s.t. X(k) = max { X(i)| X.base() <= i < X.dim() }. */
}//Daetk
#endif
