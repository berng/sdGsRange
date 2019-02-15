#include "zenith_angle.h"

#include <stdio.h>
#include <math.h>



double get_zenith_angle(float LAT,float DAYNO,float LT)
{
 float lat=LAT;
// float DAYNO;
 float h;
 float Fi;
 h=LT;
 Fi=lat/180.*PI;
 float decl;
 decl=-0.40915*cos(2.*PI*(DAYNO+8.)/365.);
 double  h_angle;
 h_angle=(h-12.)/12.*PI;
 double  cos_zen_angle;
 cos_zen_angle=sin(Fi)*sin(decl)+cos(Fi)*cos(decl)*cos(h_angle);
// fprintf(stderr,"cos: (h=%lf) %lf",h,cos_zen_angle);
  return cos_zen_angle;
}

float get_cos_zenith(float LAT,float LON,float UT,float DAYNO)
{
 float LT_max;
 LT_max=UT+LON/15;
 while(LT_max>24)
  { LT_max-=24; }
 while(LT_max<0)
  { LT_max+=24; }
// fprintf(stderr,"cos: (UT=%lf)",UT);
 return (float)get_zenith_angle(LAT,DAYNO,LT_max);
}

/*
float model_foF2(float lat,float lon,float ut,float dayno,float foF2max,float foF2min,float betta,float alpha)
 {
  double zen,zen1,zen2;

  zen=get_cos_zenith(lat,lon,ut,dayno);
  float ut1;
  ut1=-lon/15;
  while(ut1>24)
   { ut1-=24; }
  while(ut1<0)
   { ut1+=24; }
  zen1=get_cos_zenith(lat,lon,ut1,dayno);

  ut1+=12;
  while(ut1>24)
   { ut1-=24; }

  zen2=get_cos_zenith(lat,lon,ut1,dayno);
  
  return (foF2max-foF2min)*(
                            (atan(betta*(zen-alpha))-atan(betta*(zen1-alpha)))
                             /(atan(betta*(zen2-alpha))-atan(betta*(zen1-alpha)))
			   )
		    +foF2min;
 }
*/

float model_foF2(float lat,float lon,float ut,float dayno,float foF2max,float foF2min,float betta,float alpha)
 {
  float mod_norm;
  mod_norm=model_normalized(lat,lon,ut,dayno,betta,alpha);
  return (foF2max-foF2min)*mod_norm+foF2min;
 }


float model_normalized(float lat,float lon,float ut,float dayno,float betta,float alpha)
 {
  double zen,zen1,zen2;

  zen=get_cos_zenith(lat,lon,ut,dayno);
  float ut1;
  ut1=-lon/15;
  while(ut1>24)
   { ut1-=24; }
  while(ut1<0)
   { ut1+=24; }
  zen1=get_cos_zenith(lat,lon,ut1,dayno);

  ut1+=12;
  while(ut1>24)
   { ut1-=24; }

  zen2=get_cos_zenith(lat,lon,ut1,dayno);
  
  return (atan(betta*(zen-alpha))-atan(betta*(zen1-alpha)))/(atan(betta*(zen2-alpha))-atan(betta*(zen1-alpha)));
 }


/*
main()
{
 float dayno;
 float lat=52;
 float lon=102;
 float ut;

 for(dayno=1;dayno<366;dayno+=1)
  for(ut=0;ut<24;ut+=1)
   {
    printf("%f %f %f\n",dayno+ut/24,get_cos_zenith(lat,lon,ut,dayno),model_foF2(lat,lon,ut,dayno,10.0,3.0,10.0,-0.3));
   }
}
*/
