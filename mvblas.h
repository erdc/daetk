
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                   MV++ Numerical Matrix/Vector C++ Library                */
/*                             MV++ Version 1.5                              */
/*                                                                           */
/*                                  R. Pozo                                  */
/*               National Institute of Standards and Technology              */
/*                                                                           */
/*                                  NOTICE                                   */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that this permission notice appear in all copies and             */
/* supporting documentation.                                                 */
/*                                                                           */
/* Neither the Institution (National Institute of Standards and Technology)  */
/* nor the author makes any representations about the suitability of this    */
/* software for any purpose.  This software is provided ``as is''without     */
/* expressed or implied warranty.                                            */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#ifndef _MV_BLAS1_H_
#define _MV_BLAS1_H_

#include "Definitions.h"
#include "mvv.h"

namespace Daetk 
{
MV_Vector& operator*=(MV_Vector &x, const real &a);
MV_Vector operator*(const real &a, const MV_Vector &x);
MV_Vector operator*(const MV_Vector &x, const real &a);
MV_Vector operator+(const MV_Vector &x, 
    const MV_Vector &y);
MV_Vector operator-(const MV_Vector &x, 
    const MV_Vector &y);
MV_Vector& operator+=(MV_Vector &x, const MV_Vector &y);
MV_Vector& operator-=(MV_Vector &x, const MV_Vector &y);

real dot(const MV_Vector &x, const MV_Vector &y);
real norm(const MV_Vector &x);
}
#endif
// _MV_BLAS1_H_

