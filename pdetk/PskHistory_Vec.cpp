#include "PskHistory_Vec.h"
namespace Daetk
{
  void PskHistory_Vec::newsize(int nNodes)
  {
    _lastWaterSatVec.newsize(nNodes);
    _lastNaplSatVec.newsize(nNodes);
    _lastAirSatVec.newsize(nNodes);
    _lastTotalSatVec.newsize(nNodes);
    _lastAirWaterHeadVec.newsize(nNodes);
    _lastAirNaplHeadVec.newsize(nNodes);
    _lastNaplWaterHeadVec.newsize(nNodes);
    _waterAppSatVec.newsize(nNodes);
    _totalAppSatVec.newsize(nNodes);    
    _airWaterHeadVec.newsize(nNodes);
    _airNaplHeadVec.newsize(nNodes);
    _naplWaterHeadVec.newsize(nNodes);
    _lastWaterAppSatVec.newsize(nNodes);
    _lastTotalAppSatVec.newsize(nNodes);
    _histMin2pWaterAppSatVec.newsize(nNodes);
    _histMin3pWaterAppSatVec.newsize(nNodes);
    _histMinTotalAppSatVec.newsize(nNodes);
    _maxResNaplEffSatVec.newsize(nNodes);
    _maxResAirByWaterEffSatVec.newsize(nNodes);
    _maxResAirByNaplEffSatVec.newsize(nNodes);
    _maxResAirOnMainNaplImbibitionVec.newsize(nNodes);
    _maxResAirOnMainWaterImbibitionVec.newsize(nNodes);
    _maxResNaplOnMainWaterImbibitionVec.newsize(nNodes);
    _irreducibleWaterSatVec.newsize(nNodes);
    _ItoDStack_waterVec.newsize(nNodes);
    _DtoIStack_waterVec.newsize(nNodes);
    _ItoDStack_totalVec.newsize(nNodes);
    _DtoIStack_totalVec.newsize(nNodes);
    _waterPathVec.newsize(nNodes);
    _totalPathVec.newsize(nNodes);
    _stateVec.newsize(nNodes);
    _lastStateVec.newsize(nNodes);

    _lastWaterSat = &_lastWaterSatVec(0);
    _lastNaplSat = &_lastNaplSatVec(0);
    _lastAirSat = &_lastAirSatVec(0);
    _lastTotalSat = &_lastTotalSatVec(0);
    _lastAirWaterHead = &_lastAirWaterHeadVec(0);
    _lastAirNaplHead = &_lastAirNaplHeadVec(0);
    _lastNaplWaterHead = &_lastNaplWaterHeadVec(0);
    _airWaterHead = &_airWaterHeadVec(0);
    _airNaplHead = &_airNaplHeadVec(0);
    _naplWaterHead = &_naplWaterHeadVec(0);
    _waterAppSat = &_waterAppSatVec(0);
    _totalAppSat = &_totalAppSatVec(0);
    _lastWaterAppSat = &_lastWaterAppSatVec(0);
    _lastTotalAppSat = &_lastTotalAppSatVec(0);
    _histMin2pWaterAppSat = &_histMin2pWaterAppSatVec(0);
    _histMin3pWaterAppSat = &_histMin3pWaterAppSatVec(0);
    _histMinTotalAppSat = &_histMinTotalAppSatVec(0);
    _maxResNaplEffSat = &_maxResNaplEffSatVec(0);
    _maxResAirByWaterEffSat = &_maxResAirByWaterEffSatVec(0);
    _maxResAirByNaplEffSat = &_maxResAirByNaplEffSatVec(0);
    _maxResAirOnMainNaplImbibition = &_maxResAirOnMainNaplImbibitionVec(0);
    _maxResAirOnMainWaterImbibition = &_maxResAirOnMainWaterImbibitionVec(0);
    _maxResNaplOnMainWaterImbibition = &_maxResNaplOnMainWaterImbibitionVec(0);
    _irreducibleWaterSat = &_irreducibleWaterSatVec(0);
    _ItoDStack_water = &_ItoDStack_waterVec(0);
    _DtoIStack_water = &_DtoIStack_waterVec(0);
    _ItoDStack_total = &_ItoDStack_totalVec(0);
    _DtoIStack_total = &_DtoIStack_totalVec(0);
    _waterPath = &_waterPathVec(0);
    _totalPath = &_totalPathVec(0);
    _state = &_stateVec(0);
    _lastState = &_lastStateVec(0);
  }

  PskHistory_Vec::PskHistory_Vec(unsigned int nNodes):
    _lastWaterSatVec(nNodes),
    _lastNaplSatVec(nNodes),
    _lastAirSatVec(nNodes),
    _lastTotalSatVec(nNodes),
    _lastAirWaterHeadVec(nNodes),
    _lastAirNaplHeadVec(nNodes),
    _lastNaplWaterHeadVec(nNodes),
  _waterAppSatVec(nNodes),
  _totalAppSatVec(nNodes),    
  _airWaterHeadVec(nNodes),
  _airNaplHeadVec(nNodes),
  _naplWaterHeadVec(nNodes),
  _lastWaterAppSatVec(nNodes),
  _lastTotalAppSatVec(nNodes),
  _histMin2pWaterAppSatVec(nNodes),
  _histMin3pWaterAppSatVec(nNodes),
  _histMinTotalAppSatVec(nNodes),
  _maxResNaplEffSatVec(nNodes),
  _maxResAirByWaterEffSatVec(nNodes),
  _maxResAirByNaplEffSatVec(nNodes),
  _maxResAirOnMainNaplImbibitionVec(nNodes),
  _maxResAirOnMainWaterImbibitionVec(nNodes),
  _maxResNaplOnMainWaterImbibitionVec(nNodes),
  _irreducibleWaterSatVec(nNodes),
  _ItoDStack_waterVec(nNodes),
  _DtoIStack_waterVec(nNodes),
  _ItoDStack_totalVec(nNodes),
  _DtoIStack_totalVec(nNodes),
  _waterPathVec(nNodes),
  _totalPathVec(nNodes),
  _stateVec(nNodes),
  _lastStateVec(nNodes)
{  
  _lastWaterSat = &_lastWaterSatVec(0);
  _lastNaplSat = &_lastNaplSatVec(0);
  _lastAirSat = &_lastAirSatVec(0);
  _lastTotalSat = &_lastTotalSatVec(0);
  _lastAirWaterHead = &_lastAirWaterHeadVec(0);
  _lastAirNaplHead = &_lastAirNaplHeadVec(0);
  _lastNaplWaterHead = &_lastNaplWaterHeadVec(0);
  _airWaterHead = &_airWaterHeadVec(0);
  _airNaplHead = &_airNaplHeadVec(0);
  _naplWaterHead = &_naplWaterHeadVec(0);
  _waterAppSat = &_waterAppSatVec(0);
  _totalAppSat = &_totalAppSatVec(0);
  _lastWaterAppSat = &_lastWaterAppSatVec(0);
  _lastTotalAppSat = &_lastTotalAppSatVec(0);
  _histMin2pWaterAppSat = &_histMin2pWaterAppSatVec(0);
  _histMin3pWaterAppSat = &_histMin3pWaterAppSatVec(0);
  _histMinTotalAppSat = &_histMinTotalAppSatVec(0);
  _maxResNaplEffSat = &_maxResNaplEffSatVec(0);
  _maxResAirByWaterEffSat = &_maxResAirByWaterEffSatVec(0);
  _maxResAirByNaplEffSat = &_maxResAirByNaplEffSatVec(0);
  _maxResAirOnMainNaplImbibition = &_maxResAirOnMainNaplImbibitionVec(0);
  _maxResAirOnMainWaterImbibition = &_maxResAirOnMainWaterImbibitionVec(0);
  _maxResNaplOnMainWaterImbibition = &_maxResNaplOnMainWaterImbibitionVec(0);
  _irreducibleWaterSat = &_irreducibleWaterSatVec(0);
  _ItoDStack_water = &_ItoDStack_waterVec(0);
  _DtoIStack_water = &_DtoIStack_waterVec(0);
  _ItoDStack_total = &_ItoDStack_totalVec(0);
  _DtoIStack_total = &_DtoIStack_totalVec(0);
  _waterPath = &_waterPathVec(0);
  _totalPath = &_totalPathVec(0);
  _state = &_stateVec(0);
  _lastState = &_lastStateVec(0);
}

int PskHistory_Vec::operator()(int i)
{
  _lastWaterSat = &_lastWaterSatVec(i);
  _lastNaplSat = &_lastNaplSatVec(i);
  _lastAirSat = &_lastAirSatVec(i);
  _lastTotalSat = &_lastTotalSatVec(i);
  _lastAirWaterHead = &_lastAirWaterHeadVec(i);
  _lastAirNaplHead = &_lastAirNaplHeadVec(i);
  _lastNaplWaterHead = &_lastNaplWaterHeadVec(i);
  _waterAppSat = &_waterAppSatVec(i);
  _totalAppSat = &_totalAppSatVec(i);    
  _airWaterHead = &_airWaterHeadVec(i);
  _airNaplHead = &_airNaplHeadVec(i);
  _naplWaterHead = &_naplWaterHeadVec(i);
  _lastWaterAppSat = &_lastWaterAppSatVec(i);
  _lastTotalAppSat = &_lastTotalAppSatVec(i);
  _histMin2pWaterAppSat = &_histMin2pWaterAppSatVec(i);
  _histMin3pWaterAppSat = &_histMin3pWaterAppSatVec(i);
  _histMinTotalAppSat = &_histMinTotalAppSatVec(i);
  _maxResNaplEffSat = &_maxResNaplEffSatVec(i);
  _maxResAirByWaterEffSat = &_maxResAirByWaterEffSatVec(i);
  _maxResAirByNaplEffSat = &_maxResAirByNaplEffSatVec(i);
  _maxResAirOnMainNaplImbibition = &_maxResAirOnMainNaplImbibitionVec(i);
  _maxResAirOnMainWaterImbibition = &_maxResAirOnMainWaterImbibitionVec(i);
  _maxResNaplOnMainWaterImbibition = &_maxResNaplOnMainWaterImbibitionVec(i);
  _irreducibleWaterSat = &_irreducibleWaterSatVec(i);
  
  _ItoDStack_water = &_ItoDStack_waterVec(i);
  _DtoIStack_water = &_DtoIStack_waterVec(i);
  _ItoDStack_total = &_ItoDStack_totalVec(i);
  _DtoIStack_total = &_DtoIStack_totalVec(i);

  _waterPath = &_waterPathVec(i);
  _totalPath = &_totalPathVec(i);
  _state = &_stateVec(i);
  _lastState = &_lastStateVec(i);  
  return i;
}


int PskHistory_Vec::node(int i)
{
  _lastWaterSat = &_lastWaterSatVec(i);
  _lastNaplSat = &_lastNaplSatVec(i);
  _lastAirSat = &_lastAirSatVec(i);
  _lastTotalSat = &_lastTotalSatVec(i);
  _lastAirWaterHead = &_lastAirWaterHeadVec(i);
  _lastAirNaplHead = &_lastAirNaplHeadVec(i);
  _lastNaplWaterHead = &_lastNaplWaterHeadVec(i);
  _waterAppSat = &_waterAppSatVec(i);
  _totalAppSat = &_totalAppSatVec(i);    
  _airWaterHead = &_airWaterHeadVec(i);
  _airNaplHead = &_airNaplHeadVec(i);
  _naplWaterHead = &_naplWaterHeadVec(i);
  _lastWaterAppSat = &_lastWaterAppSatVec(i);
  _lastTotalAppSat = &_lastTotalAppSatVec(i);
  _histMin2pWaterAppSat = &_histMin2pWaterAppSatVec(i);
  _histMin3pWaterAppSat = &_histMin3pWaterAppSatVec(i);
  _histMinTotalAppSat = &_histMinTotalAppSatVec(i);
  _maxResNaplEffSat = &_maxResNaplEffSatVec(i);
  _maxResAirByWaterEffSat = &_maxResAirByWaterEffSatVec(i);
  _maxResAirByNaplEffSat = &_maxResAirByNaplEffSatVec(i);
  _maxResAirOnMainNaplImbibition = &_maxResAirOnMainNaplImbibitionVec(i);
  _maxResAirOnMainWaterImbibition = &_maxResAirOnMainWaterImbibitionVec(i);
  _maxResNaplOnMainWaterImbibition = &_maxResNaplOnMainWaterImbibitionVec(i);
  _irreducibleWaterSat = &_irreducibleWaterSatVec(i);
  
  _ItoDStack_water = &_ItoDStack_waterVec(i);
  _DtoIStack_water = &_DtoIStack_waterVec(i);
  _ItoDStack_total = &_ItoDStack_totalVec(i);
  _DtoIStack_total = &_DtoIStack_totalVec(i);
  
  _waterPath = &_waterPathVec(i);
  _totalPath = &_totalPathVec(i);
  _state = &_stateVec(i);
  _lastState = &_lastStateVec(i);  
  return i;
}
}//Daetk
