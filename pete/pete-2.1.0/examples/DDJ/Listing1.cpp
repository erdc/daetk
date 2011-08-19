#include <iostream.h>
#include <vector.h>
#include <math.h>

#include "PETE/PETE.h"
#include "Listing1Operators.h"

//-----------------------------------------------------------------------------
// We need to specialize the CreateLeaf traits class for STL vectors so that
// operators know what to stick in the leaves of the expression tree. In this
// case, we store a const_iterator at the leaves.
//-----------------------------------------------------------------------------

template<class T, class Allocator>
struct CreateLeaf<vector<T, Allocator> >
{
  typedef typename vector<T, Allocator>::const_iterator Leaf_t;
  inline static Leaf_t make(const vector<T, Allocator> &a) { return a.begin(); }
};

//-----------------------------------------------------------------------------
// Loop over vector and evaluate the expression at each location.
//-----------------------------------------------------------------------------

template<class T, class Allocator, class AssignOp, class RHS>
inline void evaluate(vector<T, Allocator> &lhs, const AssignOp &op, const RHS &rhs)
{
  typename vector<T, Allocator>::iterator iterLHS = lhs.begin();
  while (iterLHS != lhs.end())
    {
      // The actual assignment operation is performed here.
      // PETE operator tags all define operator() to perform the operation.
      // (In this case op performs an assignment operation.)
      // forEach is used to compute the rhs value.  DereferenceLeaf gets the
      // values at each node by deferencing an iterator, and the tag 
      // OpCombine tells forEach to use the operator tags in the expression 
      // to combine values together.

      op(*iterLHS, forEach(rhs, DereferenceLeaf(), OpCombine()));

      // Increment the LHS iterator.
      
      iterLHS++;
            
      // Now, we have to increment the iterators for everything on the RHS.
      // The results don't need to be combined so we use a NullCombiner.
      
      forEach(rhs, IncrementLeaf(), NullCombine());
    }
}

int main()
{
  int i;
  const int n = 10;
  vector<double> A, C, E(n);
  vector<int> B, D;

  for (i = 0; i < n; ++i)
    {
      A.push_back(i);
      B.push_back(2*i);
      C.push_back(3*i);
      D.push_back(i);
    }

  A += -B + 2 * C;
  
  assign(B, 2);
  assign(D, A + B * C);
  A += where(D < 30, B, C);

  assign(E, C);
  E += E - 4 / (sin(C) + 1);

  for (i = 0; i < n; ++i)
    { 
      cout << " A[" << i << "] = " << A[i]
	   << " B[" << i << "] = " << B[i]
	   << " C[" << i << "] = " << C[i]
	   << " D[" << i << "] = " << D[i]
	   << " E[" << i << "] = " << E[i]
	   << endl;
    }
    
  return 0;
}
