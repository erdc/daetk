
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

//
//      mvvind.h        MV_Vector Index class

#ifndef _MV_VEC_INDEX_H_
#define _MV_VEC_INDEX_H_

// A MV_VecIndex is an ordered pair (start,end) denoting a subvector
//  region, similar to a Fortran 90 or Matlab colon notation.  For example, 
//
//  MV_Vector_double A(10), B(20);
//  MV_VecIndex I(2,4);
//
//  A(I) = B(MV_VecIndex(0,2); 
//
//  sets the thrid through fifth elements of A to the first two elements
//  of B.  There is no stride argument, only contiguous regions are allowed.
//
#include "Definitions.h"

namespace Daetk 
{
class MV_VecIndex
{
    private:
            unsigned int start_;
            unsigned int end_;       
            char all_;      // true if this index refers to the complete
                            // vector range.  start_ and end_ are ignored.
    public:
  MV_VecIndex();
  MV_VecIndex(unsigned int i1);
  MV_VecIndex(unsigned int i1, unsigned int i2);
  MV_VecIndex(const MV_VecIndex &s);


  int start() const;
  int end() const;
  int length() const;
  int all() const;
  MV_VecIndex& operator=(const MV_VecIndex& I);
  MV_VecIndex operator+(int i);
  MV_VecIndex& operator+=(int i);
  MV_VecIndex operator-(int i);
  MV_VecIndex& operator-=(int i);

};
}
#endif  
//  _INDEX_H_
