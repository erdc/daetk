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
// TArray.h
//-----------------------------------------------------------------------------
// DESCRIPTION:
//
// This example serves two purposes. First, it demonstrates the
// construction of an expression-template-aware array. Second, it 
// demonstrates the flexibility of the ForEach type of template-metaprogram.
//
//=============================================================================

#ifndef PETE_EXAMPLES_TARRAY_TARRAY_H
#define PETE_EXAMPLES_TARRAY_TARRAY_H

//-----------------------------------------------------------------------------
// Include files
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"
#include "ForEachInOrder.h"

#include <iostream.h>

//-----------------------------------------------------------------------------
// Forward declarations
//-----------------------------------------------------------------------------

struct PrintTag;

//-----------------------------------------------------------------------------
//
// CLASS NAME
//   TArray
//
// DESCRIPTION
//   A "tiny" three-element expression-template (ET) array class.
//
//-----------------------------------------------------------------------------

class TArray
{
public:

  //---------------------------------------------------------------------------
  // Constructors and Destructor
  //---------------------------------------------------------------------------

  TArray()
    {
      data[0] = data[1] = data[2] = 0;
    }

  TArray(int i, int j, int k)
    {
      data[0] = i; data[1] = j; data[2] = k;
    }

  TArray(const TArray &model)
    {
      data[0] = model.data[0];
      data[1] = model.data[1];
      data[2] = model.data[2];
    }

  ~TArray() {}


  //---------------------------------------------------------------------------
  // TArray and scalar assigment operators
  //---------------------------------------------------------------------------

  TArray &operator=(const TArray &rhs)
    {
      data[0] = rhs.data[0];
      data[1] = rhs.data[1];
      data[2] = rhs.data[2];

      return *this;
    }

  TArray &operator=(int i)
    {
      data[0] = data[1] = data[2] = i;
      return *this;
    }

  //---------------------------------------------------------------------------
  // Indexing operators
  //---------------------------------------------------------------------------

  int &operator[](int i) 
    { return data[i]; }

  int operator[](int i) const 
    { return data[i]; }

  //---------------------------------------------------------------------------
  // Print method used by operator<< free function.
  //---------------------------------------------------------------------------

  void print(ostream &os) const
    { 
      os << "{" << data[0] << ", " 
                << data[1] << ", " 
	        << data[2] << "}";
    }

  //===========================================================================
  // PETE Requirements
  //===========================================================================

  // Assignment operator taking expression:

  template <class Expr>
  TArray &operator=(const Expr &expr)
  {
    typedef EvalLeaf1   FTag_t;
    typedef OpCombine      CTag_t;
    typedef NullCombine    VTag_t;

    typedef CreateLeaf<Expr>::Leaf_t Expr_t;
    
    typedef ForEachInOrder<Expr_t, FTag_t, VTag_t, CTag_t> ForEach_t;

    const Expr_t &e = CreateLeaf<Expr>::make(expr);

    data[0] = ForEach_t::apply(e, FTag_t(0), VTag_t(), CTag_t());
    data[1] = ForEach_t::apply(e, FTag_t(1), VTag_t(), CTag_t());
    data[2] = ForEach_t::apply(e, FTag_t(2), VTag_t(), CTag_t());

    return *this;
  }

  // printAssign method performs an assignment and the traverses the
  // parse-tree and prints the expression to cout.

  template <class Expr>
  void printAssign(const Expr &expr)
  {
    typedef EvalLeaf1   FTag_t;
    typedef OpCombine      CTag_t;
    typedef NullCombine    VTag_t;
    
    typedef CreateLeaf<Expr>::Leaf_t Expr_t;
    
    typedef ForEachInOrder<Expr_t, FTag_t, VTag_t, CTag_t> ForEach_t;

    const Expr_t &e = CreateLeaf<Expr>::make(expr);

    data[0] = ForEach_t::apply(e, FTag_t(0), VTag_t(), CTag_t());
    data[1] = ForEach_t::apply(e, FTag_t(1), VTag_t(), CTag_t());
    data[2] = ForEach_t::apply(e, FTag_t(2), VTag_t(), CTag_t());

    typedef BinaryNode<OpAssign, TArray, Expr_t> Assign_t;

    typedef ForEachInOrder<Assign_t, PrintTag, PrintTag, NullTag> Print_t;

    OpAssign opAssign;
    Assign_t t(opAssign, *this, e);

    Print_t::apply(t, PrintTag(cout), PrintTag(cout), NullTag());

    cout << endl;

  }

private:

  // The underlying complicated data structure

  int data[3];

};

//-----------------------------------------------------------------------------
// We need to specialize CreateLeaf<T> for our class, so that operators
// know what to stick in the leaves of the expression tree.
// TArray's are containers so we want to store these by reference in the
// expression tree.
//-----------------------------------------------------------------------------

template<>
struct CreateLeaf<TArray>
{
  typedef Reference<TArray> Leaf_t;
  inline static
  Leaf_t make(const TArray &a) { return a; }
};

//-----------------------------------------------------------------------------
// ostream inserter for TArrays
//-----------------------------------------------------------------------------

ostream &operator<<(ostream &os, const TArray &a)
{
  a.print(os);
  return os;
}

//=============================================================================
// Functor definitions required for evaluation loop
//=============================================================================

//
// struct LeafFunctor<TArray, EvalLeaf1
//
// Specialication of LeafFunctor class for applying the EvalLeaf1
// tag to a TArray. The apply method simply returns the array
// evaluated at the point.
//

template<>
struct LeafFunctor<TArray, EvalLeaf1>
{
  typedef int Type_t;
  static Type_t apply(const TArray &a, const EvalLeaf1 &f)
    { return a[f.i1_m]; }
};


//=============================================================================
// Functor definitions required for "print" loop
//=============================================================================

//
// struct PrintTag
//
// Tag struct to carry an ostream ref through the parse tree.
//

struct PrintTag
{
  mutable ostream &os_m;
  PrintTag(ostream &os) : os_m(os) {}
};


//
// struct LeafFunctor<TArray, PrintTag>
//
// Specialization of LeafFunctor class for applying the PrintTag tag to
// a TArray. The apply method simply prints the result to the ostream
// carried by PrintTag.
//

template<>
struct LeafFunctor<TArray, PrintTag>
{
  typedef int Type_t;
  static int apply(const TArray &a, const PrintTag &f)
    { 
      f.os_m << a; 
      return 0;
    }
};

//
// struct LeafFunctor<Scalar<T>, PrintTag>
//
// Specialication of TagFunctor class for applying the PrintTag tag to
// a Scalar<T>. The apply method simply prints the result to the
// ostream carried by PrintTag.
//

template<class T>
struct LeafFunctor<Scalar<T>, PrintTag>
{
  typedef int Type_t;
  static int apply(const Scalar<T> &s, const PrintTag &f)
    { 
      f.os_m << s.value(); 
      return 0;
    }
};


//
// struct ParenPrinter
//
// Utility class to provide start and finish functions that print the
// open paren and closed paren for an expression as the ForEach moves
// down and back up an edge of the parse-tree.
//

template<class Op>
struct ParenPrinter
{
  static void start(Op,PrintTag) {};

  static void finish(Op,PrintTag) {};
};


// Only need to parenthesize assignment

template<>
struct ParenPrinter<OpAdd>
{
  static void start(OpAdd,PrintTag p)
    { p.os_m << "("; }

  static void finish(OpAdd,PrintTag p)
    { p.os_m << ")"; }
};

//
// struct TagVisitor<OpAdd, PrintTag>
//
// Specialication of TagVisitor class for applying the PrintTag
// tag to OpAdd nodes. The visit method simply prints a symbol
// to the ostream carried by PrintTag.
//

template <>
struct TagVisitor<OpAdd, PrintTag> : public ParenPrinter<OpAdd>
{ 
  static void visit(OpAdd op, PrintTag t) 
    {
      t.os_m << " + ";
    }
};

//
// struct TagVisitor<OpMultiply, rintTag>
//
// Specialication of TagVisitor class for applying the PrintTag
// tag to OpMultiply nodes. The visit method simply prints a symbol
// to the ostream carried by PrintTag.
//

template <>
struct TagVisitor<OpMultiply, PrintTag> : public ParenPrinter<OpMultiply>
{ 
  static void visit(OpMultiply op, PrintTag t) 
    {
      t.os_m << " * ";
    }
};

//
// struct TagVisitor<OpAssign, PrintTag>
//
// Specialication of TagVisitor class for applying the PrintTag
// tag to OpAssign nodes. The visit method simply prints a symbol
// to the ostream carried by PrintTag.
//

template <>
struct TagVisitor<OpAssign, PrintTag> : public ParenPrinter<OpAssign>
{ 
  static void visit(OpAssign op, PrintTag t) 
    {
      t.os_m << " = ";
    }
};

#include "TArrayOperators.h"

#endif // PETE_EXAMPLES_TARRAY_TARRAY_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TArray.h,v $   $Author: ckees $
// $Revision: 1.1 $   $Date: 2001/08/10 15:22:48 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
