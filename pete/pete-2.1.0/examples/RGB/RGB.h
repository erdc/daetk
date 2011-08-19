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
// This example shows how to use combiners to synthesize a type at compile
// time using ForEach.
//=============================================================================

#ifndef PETE_EXAMPLES_RGB_RGB_H
#define PETE_EXAMPLES_RGB_RGB_H

//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------

#include <iostream.h>

#include "PETE/PETE.h"

//-----------------------------------------------------------------------------
// Tag classes representing colors. Also, a functor for getting a color from
// a leaf and a combiner for combining colors.
//-----------------------------------------------------------------------------

struct Red { };
struct Green { };
struct Blue { };
struct GetColor { };
struct ColorCombine { };

//-----------------------------------------------------------------------------
// A few overloaded functions that print out color names given a type.
//-----------------------------------------------------------------------------

inline void calc(const Red &)
{
  cout << "This expression is red." << endl;
}

inline void calc(const Blue &)
{
  cout << "This expression is blue." << endl;
}

inline void calc(const Green &)
{
  cout << "This expression is green." << endl;
}

//-----------------------------------------------------------------------------
// A function that deduces a color at compile time and calls a special
// function based on the value.
//
//-----------------------------------------------------------------------------

template <class Expr>
void printColor(const Expression<Expr> &expr)
{
  typedef typename ForEach<Expression<Expr>, GetColor, ColorCombine>::Type_t 
    DeducedColor_t;
    
  calc(DeducedColor_t());
}

//-----------------------------------------------------------------------------
// A set of combiners that produce new colors according to some arbitrary
// rules: red & green give blue, red & blue give green, blue and green give 
// red.
//-----------------------------------------------------------------------------

template<class Op>
struct Combine2<Red, Green, Op, ColorCombine>
{
  typedef Blue Type_t;
};

template<class Op>
struct Combine2<Red, Blue, Op, ColorCombine>
{
  typedef Green Type_t;
};

template<class Op>
struct Combine2<Green, Blue, Op, ColorCombine>
{
  typedef Red Type_t;
};

template<class Op>
struct Combine2<Green, Red, Op, ColorCombine>
{
  typedef Blue Type_t;
};

template<class Op>
struct Combine2<Blue, Green, Op, ColorCombine>
{
  typedef Red Type_t;
};

template<class Op>
struct Combine2<Blue, Red, Op, ColorCombine>
{
  typedef Green Type_t;
};

template<class C1, class C2, class Op>
struct Combine2<C1, C2, Op, ColorCombine>
{
  typedef C1 Type_t;
};

//-----------------------------------------------------------------------------
// A class that has a single template parameter that specifies a color.
//-----------------------------------------------------------------------------

template<class ColorTag>
class Operand
{
};

//-----------------------------------------------------------------------------
// We need to specialize CreateLeaf<T> for Operand, so that operators
// know what to stick in the leaves of the expression tree.
//-----------------------------------------------------------------------------

template<class ColorTag>
struct CreateLeaf<Operand<ColorTag> >
{
  typedef Operand<ColorTag> Leaf_t;
  inline static
  const Leaf_t &make(const Operand<ColorTag> &a) { return a; }
};

//-----------------------------------------------------------------------------
// Specialization of LeafFunctor class for applying the getting the "color"
// of an operand.
//-----------------------------------------------------------------------------

template<class Color>
struct LeafFunctor<Operand<Color>, GetColor>
{
  typedef Color Type_t;
};

#include "RGBOperators.h"

#endif // PETE_EXAMPLES_RGB_RGB_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: RGB.h,v $   $Author: ckees $
// $Revision: 1.1 $   $Date: 2001/08/10 15:22:47 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
