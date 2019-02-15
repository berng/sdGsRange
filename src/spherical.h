/// corect version, works well at every latitudes and ranges
#ifndef __SPHERICAL_H__
#define __SPHERICAL_H__
/// \bug do not calculate k1 and az1
/// \bug not useful for magnetic aspect calculations
double Spherical_R_az(double LatSrc,double LonSrc,double R,double az, double* LatDest,double* LonDest,double* k1,double* az1);
#endif
