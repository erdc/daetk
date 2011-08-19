
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
//      mv_vector.h       Basic vector class (real precision)
//

#ifndef _MV_VECTOR_H
#define _MV_VECTOR_H    



#include <iostream>       // for formatted printing of matrices
#include "Definitions.h"

#include "mvvind.h"

// this is really used as a sort of global constant. The reason
// for creating its own type is that so it can be overloaded to perform
// a deep or shallow assignement.  (Any variable of type MV_Vector_::ref_type
// has only one possible value: one.)
//   It is included as a seperate file to avoid multiple definitions.

#include "mvvrf.h"

namespace Daetk 
{
class MV_Vector
{                                                                      
    protected:                                                           
           real *p_;
           unsigned int dim_;
           int ref_;  // 0 or 1; does this own its own memory space?
    public:                                                            


        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
    MV_Vector();                             
    MV_Vector(unsigned int);                             
    MV_Vector(unsigned int, const real&);   
                                                     
                                                    
    MV_Vector(real*, unsigned int);      
    MV_Vector(const real*, unsigned int);        
    MV_Vector(const MV_Vector &); 
    
    // reference of an exisiting data structure
    //
    // note that ref() is initalized with i rather than 1.
    // this is so compilers will not generate a warning that i was
    // not used in the construction.  (MV_Vector::ref_type is an enum that
    // can *only* have the value of 1.
    //
  MV_Vector(real* d, unsigned int N, MV_Vector_::ref_type i);

  MV_Vector(const MV_Vector &V, MV_Vector_::ref_type i);

    ~MV_Vector();                              
                                                                       
        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
                                                                       

  inline real&      operator()(unsigned int i);
  inline const  real&       operator()(unsigned int i) const ;

  inline real&      operator[](unsigned int i);
  inline const  real&       operator[](unsigned int i) const ;


    MV_Vector operator()(const MV_VecIndex &I) ;
    MV_Vector operator()(void);
    const MV_Vector operator()(void) const;
    const MV_Vector operator()(const MV_VecIndex &I) const;

  inline unsigned int             size() const;
  inline unsigned int             dim() const;
  inline int                      ref() const;
  inline int                      null() const;
            //
            // Create a new *uninitalized* vector of size N
            MV_Vector & newsize(unsigned int );
                                                                       
        /*::::::::::::::*/                                             
        /*  Assignment  */                                             
        /*::::::::::::::*/                                             
                                                                       
    MV_Vector & operator=(const MV_Vector&);
    MV_Vector & operator=(const real&);


    friend std::ostream& operator<<(std::ostream &s, const MV_Vector &A);

  real* castToArray();
  const real* castToConstArray() const ;
};                                                                     
    
real&      MV_Vector::operator()(unsigned int i)
                  {
#                   ifdef MV_VECTOR_BOUNDS_CHECK
                    assert(i < dim_);
#                   endif
                    return p_[i];
                  }
    const  real&       MV_Vector::operator()(unsigned int i) const 
                  {
#                   ifdef MV_VECTOR_BOUNDS_CHECK
                    assert(i < dim_);
#                   endif
                    return p_[i];
                  }

    real&      MV_Vector::operator[](unsigned int i)
                  {
#                   ifdef MV_VECTOR_BOUNDS_CHECK
                    assert(i < dim_);
#                   endif
                    return p_[i];
                  }
    const  real&       MV_Vector::operator[](unsigned int i) const 
                  {
#                   ifdef MV_VECTOR_BOUNDS_CHECK
                    assert(i < dim_);
#                   endif
                    return p_[i];
                  }

    inline unsigned int             MV_Vector::size() const { return dim_;}
    inline unsigned int             MV_Vector::dim() const { return dim_;}
    inline int                      MV_Vector::ref() const { return  ref_;}
    inline int                      MV_Vector::null() const {return dim_== 0;}
}//Daetk
#endif
