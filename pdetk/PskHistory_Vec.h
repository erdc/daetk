#ifndef PSKHISTORY_VEC
#define PSKHISTORY_VEC

#include "Definitions.h"
#include "CMRVec.h"
#include "PskHistory.h"

namespace Daetk
{
class PskHistory_Vec : public PskHistory
{
public:
  PskHistory_Vec(){}
  PskHistory_Vec(unsigned int nNodes);  
  ~PskHistory_Vec(){}
  void newsize(int i);
  int operator()(int i);
  int node(int i);
private:
    CMRVec<real> _lastWaterSatVec,
    _lastNaplSatVec,
    _lastAirSatVec,
    _lastTotalSatVec,
    _lastAirWaterHeadVec,
    _lastAirNaplHeadVec,
    _lastNaplWaterHeadVec,
    _waterAppSatVec,
    _totalAppSatVec,    
    _airWaterHeadVec,
    _airNaplHeadVec,
    _naplWaterHeadVec,
    _lastWaterAppSatVec,
    _lastTotalAppSatVec,
    _histMin2pWaterAppSatVec,
    _histMin3pWaterAppSatVec,
    _histMinTotalAppSatVec,
    _maxResNaplEffSatVec,
    _maxResAirByWaterEffSatVec,
    _maxResAirByNaplEffSatVec,
    _maxResAirOnMainNaplImbibitionVec,
    _maxResAirOnMainWaterImbibitionVec,
    _maxResNaplOnMainWaterImbibitionVec,
    _irreducibleWaterSatVec;

  CMRVec< std::vector<ReversalPoint> > _ItoDStack_waterVec,
    _DtoIStack_waterVec,
    _ItoDStack_totalVec,
    _DtoIStack_totalVec;
  
  CMRVec<SaturationPath> _waterPathVec,_totalPathVec;
  CMRVec<SystemState> _stateVec,_lastStateVec;  
};

}//Daetk
#endif
