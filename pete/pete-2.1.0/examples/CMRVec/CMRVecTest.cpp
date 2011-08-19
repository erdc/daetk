#include "CMRVec.h"
#include "CMRVecConstIterator.h"
#include "CMRVecIterator.h"
#include "CMRVecOperators.h"
#include "Chronograph.h"
#include "CMRVecBlasd.h"

int main()
{
  CMRVec<double> a(80),b(80),a2(80),b2(80);
  a = 1.397;
  b = 1.397;
  a2 = 1.397;
  b2 = 1.397;
  Chronograph time;
  for (int i=0;i<100000;i++)
    a += 2.0*b;
  time.stop();
  cout<<a<<endl;
  cout<<time.elapsed()<<endl;
  a-=2.0*b;
  a=2.0*a+2.0*b+4.0*b;
  a*=2.1;
  a/=2.1;
  a+=b;
  a-=b;
  time.reset();
  time.start();
  for (int i=0;i<100000;i++)
    daxpy(2.0,b2,a2);
  time.stop();
  time.elapsed();
  cout<<time.elapsed()<<endl;;
  cout<<a2(0)<<endl;
}
