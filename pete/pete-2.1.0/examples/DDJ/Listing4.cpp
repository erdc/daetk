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
// We need to specialize LeafFunctor to extract iterator tags. By default, 
// it uses iterator_traits to extract the iterator tag. We also need to supply
// a specialization for Scalar<T>, which simply returns 
// random_access_iterator_tag because scalars don't place any constraints on
// iteration.
//-----------------------------------------------------------------------------

struct GetIteratorTagLeaf { };

template<class Iterator>
struct LeafFunctor<Iterator, GetIteratorTagLeaf>
{
  typedef std::iterator_traits<Iterator>::iterator_category Type_t;
};

template<class T>
struct LeafFunctor<Scalar<T>, GetIteratorTagLeaf>
{
  typedef std::random_access_iterator_tag Type_t;
};

#if defined(__MWERKS__)
template<class T>
struct LeafFunctor<const T *, GetIteratorTagLeaf>
{
  typedef std::random_access_iterator_tag Type_t;
  inline static
  Type_t apply(const T *, const GetIteratorTagLeaf &)
    { return Type_t(); }
};
#endif

//-----------------------------------------------------------------------------
// We now need to define some combiners for synthesizing an iterator tag
// from two others. The rule is that we return the tag with the fewest
// capabilities. We're only dealing with vectors and lists, so we need to
// consider only bidirectional and random access tags.
//-----------------------------------------------------------------------------

struct IteratorTagCombine { };

template<class Op>
struct Combine2<std::random_access_iterator_tag, 
                std::random_access_iterator_tag, 
                Op, IteratorTagCombine>
{
  typedef std::random_access_iterator_tag Type_t;
};

template<class Op>
struct Combine2<std::bidirectional_iterator_tag, 
                std::bidirectional_iterator_tag, 
                Op, IteratorTagCombine>
{
  typedef std::bidirectional_iterator_tag Type_t;
};

template<class Op>
struct Combine2<std::random_access_iterator_tag, 
                std::bidirectional_iterator_tag, 
                Op, IteratorTagCombine>
{
  typedef std::bidirectional_iterator_tag Type_t;
};

template<class Op>
struct Combine2<std::bidirectional_iterator_tag, 
                std::random_access_iterator_tag, 
                Op, IteratorTagCombine>
{
  typedef std::bidirectional_iterator_tag Type_t;
};

//-----------------------------------------------------------------------------
// EvalLeaf1 is used to evaluate the leaves consisting of STL vectors.
// (It's already defined for Scalar values.)
//-----------------------------------------------------------------------------

template<class Iterator>
struct LeafFunctor<Iterator, EvalLeaf1>
{
  typedef typename std::iterator_traits<Iterator>::value_type Type_t;
  inline static
  Type_t apply(const Iterator &i, const EvalLeaf1 &f)
  {
    return *(i + f.val1());
  }
};

#if defined(__MWERKS__)
template<class T>
struct LeafFunctor<const T *, EvalLeaf1>
{
  typedef T Type_t;
  inline static
  Type_t apply(const T *i, const EvalLeaf1 &f)
  {
    return *(i + f.val1());
  }
};
#endif

//-----------------------------------------------------------------------------
// Return a particular element in the most efficient way possible. We extract
// the most basic iterator tag from the expression and then dispatch (at
// compile time) to a function that does the work.
//-----------------------------------------------------------------------------

template<class Expr>
inline 
typename ForEach<Expression<Expr>, DereferenceLeaf, OpCombine>::Type_t
elem(const Expression<Expr> &expr, int n)
{
  typedef typename 
    ForEach<Expression<Expr>, GetIteratorTagLeaf, IteratorTagCombine>::Type_t
      IteratorTag_t;
  
  return elem(expr, n, IteratorTag_t());
}

template<class Expr>
inline 
typename ForEach<Expression<Expr>, DereferenceLeaf, OpCombine>::Type_t
elem(const Expression<Expr> &expr, int n, 
  const std::random_access_iterator_tag &)
{
  return forEach(expr, EvalLeaf1(n), OpCombine());
}

template<class Expr>
inline 
typename ForEach<Expression<Expr>, DereferenceLeaf, OpCombine>::Type_t
elem(const Expression<Expr> &expr, int n, 
  const std::bidirectional_iterator_tag &)
{
  for (int i = 0; i < n; ++i)
    forEach(expr, IncrementLeaf(), NullCombine());
    
  return forEach(expr, DereferenceLeaf(), OpCombine());
}

int main()
{
  int n = 10;
  vector<int> a, b, c;
  list<double> d;

  for (int i = 0; i < n;++i)
    {
      a.push_back(i);
      b.push_back(2*i);
      c.push_back(3*i);
      d.push_back(4*i);
    }

  cout << elem(2.0 * a + b * c, 4) << endl;
  cout << elem(2.5 * a + b * c + sin(d), 4) << endl;
}
