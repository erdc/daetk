#include "Poisson2p.h"

namespace Daetk
{
  Poisson2p::Poisson2p(){}
  Poisson2p::~Poisson2p(){}
  void Poisson2p::readParameters(ParameterDatabase& pd)
  {
    Psk2p::readParameters(pd);
  }
  
  void Poisson2p::millerSimilarScaling(Vec& delta){}

}//Daetk
