/// claculate spherical propagation trajectories
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "spherical.h"


#define PI (3.1415926)
#define Re (6371.0)

typedef struct {
 double x,y,z;
}vec;

vec spheric2decart(double r,double lambda,double theta)
{
 vec v;
 v.x=r*sin(theta)*cos(lambda);
 v.y=r*sin(theta)*sin(lambda);
 v.z=r*cos(theta);
 return v;
}


void decart2spheric(vec v,
		    double *r,double *lambda,double *theta
		    )
{
 *r=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
 *theta=atan2(v.z,sqrt(v.x*v.x+v.y*v.y));
// if(v.x>0.0 && v.y>=0.0)
  *lambda=atan2(v.y,v.x);
// else
//  if(v.x<=0.0)
//   *lambda=atan2(v.y,v.x)+PI;
//  else
//   *lambda=atan2(v.y,v.x)+2.*PI;
}


vec rotateOvrX(vec src, double phi)
{
 vec res;
 res.x=1*src.x+ (0*src.y) + (0*src.z);
/*
 res.y=0*src.x+ (cos(phi)*src.y) + (-sin(phi)*src.z);
 res.z=0*src.x+ (sin(phi)*src.y) + (cos(phi)*src.z);
*/
 res.y=0*src.x+ (cos(phi)*src.y) + (sin(phi)*src.z);
 res.z=0*src.x+ (-sin(phi)*src.y) + (cos(phi)*src.z);
 return res;
 
}

vec rotateOvrY(vec src, double phi)
{
 vec res;
 res.x=cos(phi)*src.x + 0*src.y+ (sin(phi)*src.z);
 res.y=0*src.x + 1.0*src.y + (0*src.z);
 res.z=-sin(phi)*src.x + 0*src.y+ (cos(phi)*src.z);
 return res; 
}

vec rotateOvrZ(vec src, double phi)
{
 vec res;
/*
 res.x=cos(phi)*src.x+ (-sin(phi)*src.y)+ 0*src.z;
 res.y=sin(phi)*src.x+ (cos(phi)*src.y) + 0*src.z;
*/
 res.x=cos(phi)*src.x+ (sin(phi)*src.y)+ 0*src.z;
 res.y=-sin(phi)*src.x+ (cos(phi)*src.y) + 0*src.z;
 res.z=0*src.x+ (0*src.y) + 1.0*src.z;
 return res;
}


double Spherical_R_az(double LatSrc,double LonSrc,double R,double az, double* LatDest,double* LonDest,double* k1,double* az1)
 {
    vec res;
    double az0;
    double lambda,theta,r;
    az0=LonSrc-180.0;
/*
    if(LatSrc<0)
     az0+=180.;
*/

    res=spheric2decart(1.0,PI*(az0-az)/180.0,PI*R/(PI*Re));
    res=rotateOvrZ(res,PI*(-(-LonSrc)/180.0)); //Longitude
    res=rotateOvrY(res,PI*((90.0-LatSrc)/180.0)); //Latitude //checked
    res=rotateOvrZ(res,PI*(-(LonSrc)/180.0)); //Longitude
    decart2spheric(res,&r,&lambda,&theta);
    *LatDest=180.*theta/PI;
    *LonDest=180.*lambda/PI;
    *k1=0.0;
    *az1=0.0;
 }
/*
main()
{
 vec res;
 double r,lambda,theta;
//#define LAT0 (56.5)
double LAT0=-76.5;
double LON0=58.5;

 double az0;//-180.0;//-LAT0;

//az0=-90;
//az0=0;
double az=25.0;
double R=1000.0;
double lat,lon;
double k1,az1;

//for(LON0=0;LON0<360;LON0+=20)
for(az=-180.;az<180.0;az+=10.)
{
 for(R=0.0;R<5000.0;R+=300.0)
 {
//  az=-180;
  Spherical_R_az(LAT0,LON0,R,az,&lat,&lon,&k1,&az1);
  printf("%lf %lf\n",lon,lat);
 }
 printf("\n");
}

}
*/