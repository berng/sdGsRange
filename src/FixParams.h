#ifndef __FIX_PARAMS_H__
#define __FIX_PARAMS_H__
#include "zenith_angle.h"
extern double H_MAX_EXTERN;
extern double DH_EXTERN;

typedef struct {
 float range;
 float P;
 float t;
 float elv_normal;
 float elv_low;
 float elv_high;
 long rec_id;
 float foF2;
 float hmin;
 float hmax;
 float f0;
 double zenith_cos;
 int beam;
 double noise;
} rec_type;

rec_type FixParams(double Xlat0,double Xlong0, float DayNo, rec_type* history,long idx, double foF2max,double foF2min,double betta,double alpha,double dt,long history_len);

#endif
