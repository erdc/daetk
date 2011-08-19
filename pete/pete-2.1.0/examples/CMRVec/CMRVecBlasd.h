#include "CMRVec.h"

typedef CMRVec<double> CMRVecd;
///
void drot(CMRVecd& X, CMRVecd& Y, double& C, double& S);
/** Apply plane rotation */


///
void dswap(CMRVecd& X, CMRVecd& Y);
/** X <-> Y */

///
void dscal(const double& ALPHA, CMRVecd& X);
/** X <- alpha * X */

///
void dcopy(const CMRVecd& X, CMRVecd& Y);
/** Y <- X */

///
void daxpy(const double& ALPHA, const CMRVecd& X, CMRVecd& Y);
/** Y <- alpha * X + Y*/

///
double ddot(const CMRVecd& X,const CMRVecd& Y);
/** ddot <- X^T * Y */

///
double dnrm2(const CMRVecd& X);
/** dnrm2 <- ||X||_2 */

///
double dasum(const CMRVecd& X);
/** dasum <- ||re(X)||_1 + ||im(X)||_1 */

///
int idamax(const CMRVecd& X);
/** idamax <- 1st k s.t. X(k) = max { X(i)| X.base() <= i < X.dim() }. */




