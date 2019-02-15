#ifndef __ZENITH_ANGLE_H__
#define __ZENITH_ANGLE_H__

#include <stdio.h>
#include <math.h>

#ifndef PI
#define  PI 3.1415926
#endif


double get_zenith_angle(float LAT,float DAYNO,float LT);
float get_cos_zenith(float LAT,float LON,float UT,float DAYNO);
float model_foF2(float lat,float lon,float ut,float dayno,float foF2max,float foF2min,float betta,float alpha);
float model_normalized(float lat,float lon,float ut,float dayno,float betta,float alpha);
#endif
