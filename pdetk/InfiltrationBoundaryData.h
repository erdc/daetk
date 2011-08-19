#ifndef INFILTRATION_BOUNDARY_DATA_H
#define INFILTRATION_BOUNDARY_DATA_H

#include "ParameterDatabase.h"
#include "Definitions.h"

//mwf dumb little class for doing infiltration problems
//mwf try to add simple class for administering infiltration
//boundary?
namespace Daetk
{

class InfiltrationBoundaryDataSimple
{
public:
  inline InfiltrationBoundaryDataSimple();
  ~InfiltrationBoundaryDataSimple() {}

  real peak_value,tFreq;
  std::set<int> infiltrationBoundary;

};

inline
InfiltrationBoundaryDataSimple::InfiltrationBoundaryDataSimple():
  peak_value(0.1),tFreq(-1.0),infiltrationBoundary()
{
  
}

class InfiltrationBoundaryData: public InfiltrationBoundaryDataSimple
{
public:
  inline InfiltrationBoundaryData(ParameterDatabase& pd);
  ~InfiltrationBoundaryData() {}

  inline real scaleBoundaryValue(real Xc);
  //uses infilStartLoc and infilStopLoc to set other values
  inline void updateParameters();
  
  //calculate infiltration rate as function of time
  inline real calculateInfiltration(real t);

  int infilStartIndex,infilStopIndex;
    real infilStartLoc,infilStopLoc,
    infilMiddleLoc,infilBoundLength;
  
  real peak_flux;


  int streamNode1,streamNode2;

  //might want this later? at least use for
  //analytical test problems
  real Kxx,Kyy,Kxy;
  //mwf put in for max/min of boundary condition?
  real psiMin,psiMax;
  //mwf put in for 3rd kind boundary condition?
  real Kss;
  //for switching between types of fluxes
  bool useConstantInfiltration;

protected:
  //for time varying infiltration profile
  real sinusoidalInfiltration(real t);
};

inline
InfiltrationBoundaryData::InfiltrationBoundaryData(ParameterDatabase& pd):
  InfiltrationBoundaryDataSimple(),
  infilStartIndex(pd.i("startBC")),infilStopIndex(pd.i("endBC")),
  infilStartLoc(0.0),infilStopLoc(0.0),infilMiddleLoc(0.0),
  infilBoundLength(0.0),peak_flux(pd.r("peakFlux")),
  streamNode1(pd.i("streamNode1")),streamNode2(pd.i("streamNode2")),
  Kxx(pd.r("Kxx")),Kyy(pd.r("Kyy")),Kxy(pd.r("Kxy")),
  psiMin(pd.r("yBack")),psiMax(0.0),Kss(pd.r("Kss")),
  useConstantInfiltration(true)
{
  
}

inline 
real InfiltrationBoundaryData::scaleBoundaryValue(real Xc)
{
#ifndef LINEAR_SCALE_BOUNDARY_VALUES
  return 1.0;
#else
  return fabs(fabs(Xc-infilMiddleLoc)
	      -0.5*infilBoundLength);
#endif
}

inline 
void InfiltrationBoundaryData::updateParameters()
{
  infilMiddleLoc  = 0.5*(infilStartLoc+infilStopLoc);
  infilBoundLength= infilStopLoc-infilStartLoc;
}

inline 
real InfiltrationBoundaryData::calculateInfiltration(real t)
{
  if (useConstantInfiltration)
    return peak_flux;
  else
    return sinusoidalInfiltration(t);
}
  
inline
real InfiltrationBoundaryData::sinusoidalInfiltration(real t)
{
  real flux,hour,minute,time;
  int hour_o_day;
  time = t*24.0;
  hour = int(time);
  minute =(time - hour);//not really a minute
  hour_o_day = int(hour) % 24;
  //  	  cout<<" in InfilBndData "<<endl;
  //  	  cout<<hour<<":"<<minute<<'\t'<<hour_o_day<<endl;
  if (hour_o_day >= 0 && hour_o_day < 1)
      flux = peak_flux * sin(minute*M_PI);
  else
    flux = 0.0;
  
  return flux;
}


}//Daetk
#endif
