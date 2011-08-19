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
//=============================================================================
// Vec3.h
//-----------------------------------------------------------------------------
// DESCRIPTION:
//
// This example serves two purposes. First, it demonstrates the
// construction of an expression-template-aware array class. Second, it 
// demonstrates the flexibility of the ForEach type of template-metaprogram.
//
//=============================================================================

#ifndef PETE_EXAMPLES_VEC3_VEC3_H
#define PETE_EXAMPLES_VEC3_VEC3_H

//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"

#include <iostream>

//-----------------------------------------------------------------------------
//
// CLASS NAME 
//   Vec3
//
// DESCRIPTION
//   A "tiny" three-element expression-template (ET) array class. 
//
//-----------------------------------------------------------------------------

class Vec3
{
public:

  //---------------------------------------------------------------------------
  // Constructors and Destructor
  //---------------------------------------------------------------------------

  Vec3() { d[0] = d[1] = d[2] = 0; }

  Vec3(int i, int j, int k) 
  {
    d[0] = i; d[1] = j; d[2] = k;
  }

  Vec3(const Vec3 &v) 
  {
    d[0] = v.d[0];  d[1] = v.d[1];  d[2] = v.d[2];
  }

  ~Vec3() {}


  //---------------------------------------------------------------------------
  // Vec3 and scalar assigment operators
  //---------------------------------------------------------------------------

  Vec3 &operator=(const Vec3 &v) 
  {
    d[0] = v.d[0];  d[1] = v.d[1];  d[2] = v.d[2];
    return *this;
  }

  Vec3 &operator=(int i) 
  {
    d[0] = d[1] = d[2] = i;
    return *this;
  }

  //---------------------------------------------------------------------------
  // Assignment operator taking expression:
  //---------------------------------------------------------------------------

  template<class RHS>
  Vec3 &operator=(const Expression<RHS> &rhs)
  {
    d[0] = forEach(rhs, EvalLeaf1(0), OpCombine());
    d[1] = forEach(rhs, EvalLeaf1(1), OpCombine());
    d[2] = forEach(rhs, EvalLeaf1(2), OpCombine());

    return *this;
  }

  //---------------------------------------------------------------------------
  // Indexing operators
  //---------------------------------------------------------------------------

  int &operator[](int i)      { return d[i]; }
  int operator[](int i) const { return d[i]; }

  //---------------------------------------------------------------------------
  // Print method used by operator<< free function.
  //---------------------------------------------------------------------------

  void print(ostream &os) const 
  { 
    os << "{" << d[0] << "," << d[1] << "," << d[2] << "}";
  }

private:

  // The underlying complicated data structure

  int d[3];

};


//-----------------------------------------------------------------------------
// We need to specialize CreateLeaf<T> for our class, so that operators
// know what to stick in the leaves of the expression tree.
//-----------------------------------------------------------------------------

template<>
struct CreateLeaf<Vec3>
{
  typedef Reference<Vec3> Leaf_t;
  inline static
  Leaf_t make(const Vec3 &a) { return Leaf_t(a); }
};

//-----------------------------------------------------------------------------
// ostream inserter for Vec3s
//-----------------------------------------------------------------------------

ostream &operator<<(ostream &os, const Vec3 &a)
{
  a.print(os);
  return os;
}

//-----------------------------------------------------------------------------
// Specialization of LeafFunctor class for applying the EvalLeaf1
// tag to a Vec3. The apply method simply returns the array
// evaluated at the point.
//-----------------------------------------------------------------------------

template<>
struct LeafFunctor<Vec3, EvalLeaf1>
{
  typedef int Type_t;
  static Type_t apply(const Vec3 &a, const EvalLeaf1 &f)
    { return a[f.val1()]; }
};

//-----------------------------------------------------------------------------
// Specialization of LeafFunctor class for applying the CountLeaf
// tag to a Vec3. The apply method simply returns 1 for a Vec3 and 0 for
// anything else.
//-----------------------------------------------------------------------------

struct CountLeaf { };

template<>
struct LeafFunctor<Vec3, CountLeaf>
{
  typedef int Type_t;
  static Type_t apply(const Vec3 &, const CountLeaf &)
    { return 1; }
};

template<class T>
struct LeafFunctor<T, CountLeaf>
{
  typedef int Type_t;
  static Type_t apply(const T &a, const CountLeaf &)
    { return 0; }
};

// We put this include at the end because
// the operators can't be defined until after Vec3 and
// CreateLeaf<Vec3> have been defined.
// (Since Vec3 isn't templated the operators aren't just
// templates.)

#include "Vec3Operators.h"

#endif // PETE_EXAMPLES_VEC3_VEC3_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Vec3.h,v $   $Author: ckees $
// $Revision: 1.1 $   $Date: 2001/08/10 15:22:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
