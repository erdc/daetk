/* Translated by the Edison Design Group C++/C front end (version 2.35) */
/* Tue Nov  4 19:41:38 1997 */
int __KAI_KCC_3_2;
#line 1 "Chronograph.C"
#line 45 "Chronograph.h"
struct Chronograph;
#line 53 "/usr/include/sys/_inttypes.h"
typedef unsigned uint32_t;
#line 214 "/usr/include/sys/types.h"
typedef uint32_t clock_t;
#line 45 "Chronograph.h"
struct Chronograph {




clock_t _start;



clock_t _stopped;



clock_t _lastlap;


double _elapsed;};
#line 306 "/usr/include/sys/time.h"
extern clock_t clock(void);
#line 9 "Chronograph.C"
extern double elapsedHMS__11ChronographFRdN21(struct Chronograph *const, double *, double *, double *); double elapsedHMS__11ChronographFRdN21( struct Chronograph *const this,  double *__2264_40_Hours, 
double *__2265_40_Mins, 
double *__2266_40_Secs)
{ auto double __T1074465212; auto double __T1074465268; auto double __T1074465324; auto double __T1074465380; auto double __T1074522136; auto double __T1074522216; auto clock_t __T1074522320; auto double __T1074522424; auto double __T1074522504; auto double __T1074522656; auto double __T1074522784;
#line 235 "Chronograph.h"
{

auto clock_t __T1074455952; __T1074522320 = ((this->_stopped));

if (__T1074522320)
__T1074455952 = __T1074522320;

else  {
__T1074455952 = (clock());

{if (0U == __T1074455952)
__T1074455952 = 1U;} } __T1074522424 = (((double)(__T1074455952 - ((this->_start)))) / (1.e+6));


(this->_elapsed) = __T1074522424;

__T1074465380 = __T1074522424; } __T1074522136 = ((this->_elapsed));
#line 23 "Chronograph.C"
if (__T1074522136 < (6.e+1))
{
__T1074465324 = (0.);
__T1074465268 = (0.);
__T1074465212 = __T1074522136;
}
else  if (__T1074522136 < (3.6e+3))
{


__T1074465324 = (0.); __T1074522504 = ((double)((unsigned long)(((this->_elapsed)) / (6.e+1))));

__T1074465268 = __T1074522504;


__T1074465212 = (__T1074522136 - ((6.e+1) * __T1074522504));
}

else  { __T1074522784 = ((double)((unsigned long)(((this->_elapsed)) / (3.6e+3))));


__T1074465324 = __T1074522784; __T1074522656 = (((this->_elapsed)) - ((3.6e+3) * ((double)((unsigned long)(((this->_elapsed)) / (3.6e+3)))))); __T1074522216 = ((double)((unsigned long)((((this->_elapsed)) - ((3.6e+3) * ((double)((unsigned long)(((this->_elapsed)) / (3.6e+3)))))) / (6.e+1))));


__T1074465268 = __T1074522216;


__T1074465212 = (__T1074522656 - ((6.e+1) * __T1074522216));
}

(*__2264_40_Hours) = __T1074465324;
(*__2265_40_Mins) = __T1074465268;
(*__2266_40_Secs) = __T1074465212;

return __T1074465380;
}
