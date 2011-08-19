#include <iostream.h>
#include <list.h>
#include <vector.h>
#include <math.h>

#include "PETE/PETE.h"
#include "Listing2Operators.h"

//-----------------------------------------------------------------------------
// We need to specialize the CreateLeaf traits class for STL vectors and lists
// so that operators know what to stick in the leaves of the expression tree. 
// In these cases, we store const_iterators at the leaves.
//-----------------------------------------------------------------------------

template<class T, class Allocator>
struct CreateLeaf<vector<T, Allocator> >
{
  typedef typename vector<T, Allocator>::const_iterator Leaf_t;
  inline static Leaf_t make(const vector<T, Allocator> &a) { return a.begin(); }
};

template<class T, class Allocator>
struct CreateLeaf<list<T, Allocator> >
{
  typedef typename list<T, Allocator>::const_iterator Leaf_t;
  inline static Leaf_t make(const list<T, Allocator> &a) { return a.begin(); }
};

//-----------------------------------------------------------------------------
// Loop over vector and evaluate the expression at each location.
//-----------------------------------------------------------------------------

template<class T, class Allocator, class Op, class RHS>
inline void evaluate(vector<T, Allocator> &lhs, const Op &op, const RHS &rhs)
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
// Loop over list and evaluate the expression at each location.
//-----------------------------------------------------------------------------

template<class T, class Allocator, class Op, class RHS>
inline void evaluate(list<T, Allocator> &lhs, const Op &op, const RHS &rhs)
{
  typename list<T, Allocator>::iterator i = lhs.begin();
  while (i != lhs.end())
    {
      // The actual assignment operation is performed here.
      // PETE operator tags all define operator() to perform the operation.
      // (In this case op performs an assignment operation.)
      // forEach is used to compute the rhs value.  DereferenceLeaf gets the
      // values at each node by deferencing an iterator, and the tag 
      // OpCombine tells forEach to use the operator tags in the expression 
      // to combine values together.

      op(*i++, forEach(rhs, DereferenceLeaf(), OpCombine()));
      
      // Now, we have to increment the iterators for everything on the RHS.
      // The results don't need to be combined so we use a NullCombiner.
      
      forEach(rhs, IncrementLeaf(), NullCombine());
    }
}

int main()
{
  int i;
  const int n = 10;
  vector<double> A, C;
  vector<int> B;
  list<int> D;
  list<double> E;

  for (i = 0; i < n; ++i)
    {
      A.push_back(i);
      B.push_back(2*i);
      C.push_back(3*i);
      D.push_back(i);
      E.push_back(0.0);
    }

  A += -B + 2 * C;
  
  assign(B, 2);
  assign(D, A + B * C);
  A += where(D < 30, B, C);

  assign(E, C);
  E += E - 4 / (sin(C) + 1);

  list<int>::const_iterator di = D.begin();
  list<double>::const_iterator ei = E.begin();
  for (i = 0; i < n; ++i)
    { 
      cout << " A[" << i << "] = " << A[i]
	   << " B[" << i << "] = " << B[i]
	   << " C[" << i << "] = " << C[i]
	   << " D[" << i << "] = " << *di++
	   << " E[" << i << "] = " << *ei++
	   << endl;
    }
    
  return 0;
}
