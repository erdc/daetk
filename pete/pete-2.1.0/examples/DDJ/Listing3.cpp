#include <iostream.h>
#include <vector.h>
#include <math.h>

#include "PETE/PETE.h"
#include "Listing1Operators.h"

//-----------------------------------------------------------------------------
// We need to specialize the CreateLeaf traits class for the STL vector
// so that operators know what to stick in the leaves of the expression tree. 
// We're going to use a different evaluation mechanism (operator[]()) so
// store references to the vector operands at the leaves.
//-----------------------------------------------------------------------------

template<class T, class Allocator>
struct CreateLeaf<vector<T, Allocator> >
{
  typedef Reference<vector<T, Allocator> > Leaf_t;
  inline static Leaf_t make(const vector<T, Allocator> &a) { return Leaf_t(a); }
};

//-----------------------------------------------------------------------------
// EvalLeaf1 is used to evaluate the leaves consisting of STL vectors.
// (It's already defined for Scalar values.)
//-----------------------------------------------------------------------------

template<class T, class Allocator>
struct LeafFunctor<vector<T, Allocator>, EvalLeaf1>
{
  typedef T Type_t;
  inline static
  Type_t apply(const vector<T, Allocator>& vec, const EvalLeaf1 &f)
  {
    return vec[f.val1()];
  }
};

//-----------------------------------------------------------------------------
// We need to write a functor that is capable of comparing the size of
// the vector with a stored value. Then, we supply LeafFunctor specializations
// for Scalar<T> and STL vector leaves.
//-----------------------------------------------------------------------------

class SizeLeaf
{
public:

  SizeLeaf(int s) : size_m(s) { }
  SizeLeaf(const SizeLeaf &model) : size_m(model.size_m) { }
  bool operator()(int s) const { return size_m == s; }
  
private:
  
  int size_m;
  
};

template<class T>
struct LeafFunctor<Scalar<T>, SizeLeaf>
{
  typedef bool Type_t;
  inline static
  bool apply(const Scalar<T> &, const SizeLeaf &) 
  {
    // Scalars always conform.
    
    return true;
  }
};

template<class T, class Allocator>
struct LeafFunctor<vector<T, Allocator>, SizeLeaf>
{
  typedef bool Type_t;
  inline static
  bool apply(const vector<T, Allocator> &v, const SizeLeaf &s) 
  {
    return s(v.size());
  }
};

//-----------------------------------------------------------------------------
// Loop over vector and evaluate the expression at each location.
//-----------------------------------------------------------------------------

template<class T, class Allocator, class Op, class RHS>
inline void evaluate(vector<T, Allocator> &lhs, const Op &op, const RHS &rhs)
{
  if (forEach(rhs, SizeLeaf(lhs.size()), AndCombine()))
    {
      // We get here if the vectors on the RHS are the same size as those on
      // the LHS.
      
      for (int i = 0; i < lhs.size(); ++i)
        {
          // The actual assignment operation is performed here.
          // PETE operator tags all define operator() to perform the operation.
          // (In this case op performs an assignment operation.)
          // forEach is used to compute the rhs value.  EvalLeaf1 gets the
          // values at each node using random access, and the tag 
          // OpCombine tells forEach to use the operator tags in the expression 
          // to combine values together.

          op(lhs[i], forEach(rhs, EvalLeaf1(i), OpCombine()));
        }
    }
  else
    {
      cerr << "Error: LHS and RHS don't conform." << endl;
      exit(1);
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
