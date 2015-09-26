#ifndef CMRVECINDEX_H
#define CMRVECINDEX_H

#include "Definitions.h"

namespace Daetk 
{
///
class CMRVecIndex

/** This class is used to index unit stride subvectors of CMRVec's.  */

{
public:
  ///
   CMRVecIndex();
  /** Create a CMRVecIndex that accesses the entire vector */

  ///
   CMRVecIndex( int i1);
  /** Create a CMRVecIndex that accesses the element i1 */

  ///
   CMRVecIndex( int i1,  int i2);
  /** Create a CMRVecIndex that accesses a subvector starting at i1
      and ending at i2 */

  ///
   CMRVecIndex(const CMRVecIndex &s);
  /** Create a copy of s */

  ///
   int start() const;
  /** Return the first index of subvectors accessed by this CMRVecIndex */
  
  ///
   int end() const;
  /** Return the last index of subvectors accessed by this CMRVecIndex */

  ///
   int length() const;
  /** Return the length of a subvector accessed by this CMRVecIndex */

  ///
   bool all() const;
  /** Return true if the CMRVecIndex accesses an entire vector */

  ///
   CMRVecIndex& operator=(const CMRVecIndex& VI);
  /** Copy I */

  ///
   CMRVecIndex operator+(int i);
  /** Return a CMRVecIndex with start=this->start_+i and
      end=this->end_+i */

  ///
   CMRVecIndex& operator+=(int i);
  /** Add i to start and end */

  ///
   CMRVecIndex operator-(int i);
  /** Return a CMRVecIndex with start= this->start_-i and
      end=this->end_-i */
  
  ///
   CMRVecIndex& operator-=(int i);
  /** Subtract i from start and end */

private:
  int start_;
  int end_;       
  bool all_;
};


}//Daetk
#endif  



