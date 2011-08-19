// -*- C++ -*-
// ACL:license
// ----------------------------------------------------------------------
// This software and ancillary information (herein called "SOFTWARE")
// called PETE (Portable Expression Template Engine) is
// made available under the terms described here.  The SOFTWARE has been
// approved for release with associated LA-CC Number LA-CC-99-5.
// 
// Unless otherwise indicated, this SOFTWARE has been authored by an
// employee or employees of the University of California, operator of the
// Los Alamos National Laboratory under Contract No.  W-7405-ENG-36 with
// the U.S. Department of Energy.  The U.S. Government has rights to use,
// reproduce, and distribute this SOFTWARE. The public may copy, distribute,
// prepare derivative works and publicly display this SOFTWARE without 
// charge, provided that this Notice and any statement of authorship are 
// reproduced on all copies.  Neither the Government nor the University 
// makes any warranty, express or implied, or assumes any liability or 
// responsibility for the use of this SOFTWARE.
// 
// If SOFTWARE is modified to produce derivative works, such modified
// SOFTWARE should be clearly marked, so as not to confuse it with the
// version available from LANL.
// 
// For more information about PETE, send e-mail to pete@acl.lanl.gov,
// or visit the PETE web page at http://www.acl.lanl.gov/pete/.
// ----------------------------------------------------------------------
// ACL:license

#ifndef PETE_EXAMPLES_MIXED_EVAL_H
#define PETE_EXAMPLES_MIXED_EVAL_H

//-----------------------------------------------------------------------------
// Includes
//-----------------------------------------------------------------------------

#include <iostream.h>

#include "MixedOperators.h"

//-----------------------------------------------------------------------------
// Eval.h
//
// This file contains several class definitions that are use to evaluate
// expressions containing STL vectors.  The main function defined at the end
// is evaluate(lhs,op,rhs), which allows the syntax:
// vector<int> a,b,c;
// evaluate(a,OpAssign(),b+c);
//
// evaluate() is called by all the global assignment operator functions
// defined in MixedOperators.h
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// We need to specialize CreateLeaf<T> for our classes, so that operators
// know what to stick in the leaves of the expression tree.
//-----------------------------------------------------------------------------

template<class T>
struct CreateLeaf<vector<T> >
{
  typedef vector<T>::const_iterator Leaf_t;
  inline static
  Leaf_t make(const vector<T> &a) { return a.begin(); }
};

template<class T>
struct CreateLeaf<list<T> >
{
  typedef list<T>::const_iterator Leaf_t;
  inline static
  Leaf_t make(const list<T> &a) { return a.begin(); }
};


//-----------------------------------------------------------------------------
// evaluate(lhs,op,rhs) (vector on LHS)
//
// Loop over size and evaluate the expression at each location.
//-----------------------------------------------------------------------------

template<class T,class Op,class RHS>
inline
void evaluate(vector<T> &lhs,const Op &op,const RHS &rhs)
{
  for (int i = 0; i < lhs.size(); ++i)
    {
      // The actual assignment operation is performed here.
      // PETE operator tags all define operator() to perform the operation.
      // (In this case op performs an assignment operation.)
      // forEach is used to compute the rhs value.  DereferenceLeaf gets the
      // values at each node by deferencing an iterator, and the tag 
      // OpCombine tells forEach to use the operator tags in the expression 
      // to combine values together.

      op(lhs[i], forEach(rhs, DereferenceLeaf(), OpCombine()));
      
      // Now, we have to increment the iterators for everything on the RHS.
      // The results don't need to be combined so we use a NullCombiner.
      
      forEach(rhs, IncrementLeaf(), NullCombine());
    }
}


//-----------------------------------------------------------------------------
// evaluate(lhs,op,rhs) (list on LHS)
//
// Loop over size and evaluate the expression at each location.
//-----------------------------------------------------------------------------

template<class T,class Op,class RHS>
inline
void evaluate(list<T> &lhs,const Op &op,const RHS &rhs)
{
  typename list<T>::iterator i = lhs.begin();
  while (i != lhs.end())
    {
      // The actual assignment operation is performed here.
      // PETE operator tags all define operator() to perform the operation.
      // (In this case op performs an assignment operation.)
      // forEach is used to compute the rhs value.  DereferenceLeaf gets the
      // values at each node, and the tag OpCombine tells forEach to use the
      // operator tags in the expression to combine values together.

      op(*i++, forEach(rhs, DereferenceLeaf(), OpCombine()));
      
      // Now, we have to increment the iterators for everything on the RHS.
      
      forEach(rhs, IncrementLeaf(), NullCombine());
    }
}

#endif // PETE_EXAMPLES_MIXED_EVAL_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Eval.h,v $   $Author: ckees $
// $Revision: 1.1 $   $Date: 2001/08/10 15:22:47 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
