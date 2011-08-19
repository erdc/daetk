// chrono.cpp - Chronograph class member
//              function.
// Greg Messer

#include "Chronograph.h"

#define ULONG unsigned long

double Chronograph::elapsedHMS(double &Hours,
                               double &Mins,
                               double &Secs)
{
    double dHour   = 60.0 * 60.0;
    double dMinute = 60.0;
    double dElapsed;
    double dWork;
    double dHours;
    double dMins;
    double dSecs;

    dElapsed = elapsed();

    if(dElapsed < dMinute)
    {
        dHours = 0.0;
        dMins = 0.0;
        dSecs = dElapsed;
    }
    else if(dElapsed < dHour)
    {
        dWork = dElapsed;

        dHours = 0.0;

        dMins = (ULONG) (dElapsed / dMinute);
        dWork -= dMins * dMinute;

        dSecs = dWork;
    }
    else
    {
        dWork = dElapsed;

        dHours = (ULONG) (dElapsed / dHour);
        dWork -= dHours * dHour;

        dMins = (ULONG) (dWork / dMinute);
        dWork -= dMins * dMinute;

        dSecs = dWork;
    }

    Hours = dHours;
    Mins  = dMins;
    Secs  = dSecs;

    return dElapsed;
}

// end of chrono.cpp
