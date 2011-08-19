
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
//      mvmtp.h  : basic templated numerical matrix class, storage
//                  by columns (Fortran oriented.)
//
//
//


#ifndef _MV_MATRIX_H_
#define _MV_MATRIX_H_    

#include "Definitions.h"
#include "LinearOperator.h"
#include "mvv.h"
#include "mvmrf.h"



#include <iostream>       // for formatted printing of matrices


namespace Daetk 
{
class MV_ColMat  : public LinearOperator
{                                                                      
    private:                                                           
           MV_Vector v_;
           int dim0_;   // perferred to using dim_[2]. some compilers
           int dim1_;   // refuse to initalize these in the constructor.
           int lda_;
           int ref_;   // true if this is declared as a reference vector,
                        // i.e. it does not own the memory space, but 
                        // rather it is a view to another vector or array.
           int *pivots;  //added by cek

    public:                                                            
           int* pivotPtr();//added by cek
                                                                 
        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
            MV_ColMat();                             
            MV_ColMat(unsigned int, unsigned int); 

    // some compilers have difficulty with inlined 'for' statements.
    MV_ColMat(unsigned int, unsigned int, const real&);   

    // usual copy by value
    // (can't use default parameter lda=m, because m is not a constant...)
    //
    MV_ColMat(real*, unsigned int m, unsigned int n);
    MV_ColMat(real*, unsigned int m, unsigned int n, unsigned int lda);

    // the "reference" versions
    //
    //
    MV_ColMat(MV_ColMat &A, MV_Matrix_::ref_type i);
    MV_ColMat(real*, unsigned int m, unsigned int n, MV_Matrix_::ref_type i);
    MV_ColMat(real*, unsigned int m, unsigned int n, unsigned int lda,
                MV_Matrix_::ref_type i);

    MV_ColMat(const MV_ColMat&); 
    virtual ~MV_ColMat();                              
                                                                       
        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
    
  bool apply(const Vec& x, Vec& Ax);
                                                               
    inline real&       operator()(unsigned int, unsigned int); 
    inline const real& operator()(unsigned int, unsigned int) const; 
    MV_ColMat operator()(const MV_VecIndex &I, const MV_VecIndex &J) ;
    const MV_ColMat operator()(const MV_VecIndex &I, const MV_VecIndex &J) const;
    unsigned int            dim(int i) const; 
    unsigned int            lda(void) const;
    unsigned int            size(int i) const;
    MV_ColMat&        newsize(unsigned int, unsigned int);
    int ref() const;
                                                                       
        /*::::::::::::::*/                                             
        /*  Assignment  */                                             
        /*::::::::::::::*/                                             
                                                                       
    MV_ColMat & operator=(const MV_ColMat&);
    MV_ColMat & operator=(const real&);


    friend std::ostream& operator<<(std::ostream &s, const MV_ColMat &A);

  //functions for compatability with UNC_SIM
  void setSize(int m, int n);
  void Solve(MV_Vector, MV_Vector);

  //unc_sim functions
  int rowNum();
  int colNum();

  //extensions for wrappers
  real* castToArray();
  const real* castToConstArray() const;
};                                                                     

inline real& MV_ColMat::operator()(unsigned int i, unsigned int j)
{
#ifdef MV_MATRIX_BOUNDS_CHECK
    assert(0<=i && i<dim(0));
    assert(0<=j && j<dim(1));
#endif
    return v_(j*lda_ + i);      // could use indirect addressing
                                // instead...
}

inline const real& MV_ColMat::operator()
                    (unsigned int i, unsigned int j) const
{
#ifdef MV_MATRIX_BOUNDS_CHECK
    assert(0<=i && i<dim(0));
    assert(0<=j && j<dim(1));
#endif
    return v_(j*lda_ + i);
}

inline MV_ColMat::MV_ColMat(real* d, unsigned int m, 
                            unsigned int n, MV_Matrix_::ref_type i ): LinearOperator(m,n),
            v_(d,m*n, MV_Vector_::ref), dim0_(m), dim1_(n), lda_(m), ref_(i) {}

inline MV_ColMat::MV_ColMat( MV_ColMat &A, 
                             MV_Matrix_::ref_type i ):LinearOperator(A),
                v_(&A(0,0), A.dim(0)*A.dim(1), MV_Vector_::ref), 
                dim0_(A.dim(0)), dim1_(A.dim(1)), lda_(A.lda()), ref_(i) {}

inline MV_ColMat::MV_ColMat(real* d, unsigned int m, unsigned int n,
                            unsigned int lda, MV_Matrix_::ref_type i) : LinearOperator(m,n),
            v_(d, lda*n, MV_Vector_::ref), dim0_(m), dim1_(n), lda_(lda),
            ref_(i) {}

}//Daetk
#endif
