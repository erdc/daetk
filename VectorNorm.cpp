#include "VectorNorm.h"

namespace Daetk
{
  VectorNorm::VectorNorm()
    {
      Tracer tr("VectorNorm::VectorNorm()");
    }
   VectorNorm::~VectorNorm()
    {
      Tracer tr("VectorNorm::~VectorNorm()");
    }
  Norm2::Norm2()
    {
      Tracer tr("Norm2::Norm2()");
    }
  Norm2::~Norm2()
    {
      Tracer tr("Norm2::~Norm2()");
    }
   real Norm2::operator()(const Vec& y){return norm(y);}
  void Norm2::setWeight(const Vec& y){}
  const Vec& Norm2::getWeight() const {static Vec tmp;return tmp;}
   const Vec& Norm2::getScaling(){std::cerr<<"You shouldn't scale when using Norm2"<<std::endl;static Vec r; return r;}
   void Norm2::setTolerances(const Vec& atol, const Vec& rtol){}
   void Norm2::setTolerances(const real& atol,const real& rtol){}
   real* Norm2::getWeightBegin(){return 0;}
   real* Norm2::getWeightEnd(){return 0;}
   void Norm2::scale(const Vec& x,Vec& y)
    {
#ifndef USE_BLAS
      y=x;
#else
      copy(x,y);
#endif
    }
   void Norm2::deScale(const Vec& x, Vec& y)
    {
#ifndef USE_BLAS
      y=x;
#else
      copy(x,y);
#endif
    }

}//Daetk
