#ifndef DASPKJAC_H
#define DASPKJAC_H

#include "Definitions.h"
#include "Jacobian.h"
#include "Preconditioner.h"
#include "DaeDefinition.h"
 
namespace Daetk 
{
 
typedef void (*GlobDaspkRes)(const real& t,real* y,real* yp,
                             real& cj,real* delta, int* ires, real* rpar, 
                             int *ipar);
typedef void (*jac_direct_p)( real& t,  real* y,  real* yprime, 
                         real* pd, real& cj, real* rpar, int* ipar);
typedef void (jac_krylov_p)(GlobDaspkRes res,int* ires,int* neq,  real& t,  real* y,  real* yprime,
                   real* weight,  real* residual,real& work, 
                   real& h,  real& CJ, 
                  real* wp, int* iwp,  int& ier, real* rpar,int* ipar);  

class DaspkJacobian
{
public:
  typedef void (*jac_direct_p)( real& t,  real* y,  real* yprime, 
                                real* pd, real& cj, real* rpar, int* ipar);
  typedef void (*jac_krylov_p)(GlobDaspkRes res,int* ires,int* neq,  real& t,  real* y,  real* yprime,
                              real* weight,  real* residual,real& work, 
                              real& h,  real& CJ, 
                              real* wp, int* iwp,  int& ier, real* rpar,int* ipar);  
  DaspkJacobian(){}
  virtual ~DaspkJacobian(){}
  static void jac_direct( real& t,  real* y,  real* yprime, 
                         real* pd, real& cj, real* rpar, int* ipar);
  static void jac_krylov(GlobDaspkRes res,int* ires,int* neq,  real& t,  real* y,  real* yprime,
                   real* weight,  real* residual,real& work, 
                   real& h,  real& CJ, 
                  real* wp, int* iwp,  int& ier, real* rpar,int* ipar);  

  static DaeDefinition* theDae;
  static Preconditioner* thePrec;
  static Jacobian* theJac;


};
}//daetk
#endif
