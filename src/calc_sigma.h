#ifndef __CALC_SIGMA_H__
#define __CALC_SIGMA_H__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

typedef struct{
double dr;
double r;
double dr_w;
double Ampl;
long idx;
}dr_type;

int order_dr(const void* a,const void* b);
int order_drw(const void* a,const void* b);
int order_Ampl(const void* a,const void* b);
int order_Ampl_inv(const void* a,const void* b);
double calc_sigma(dr_type* dr,long len,dr_type* dr_sve,int sigma_type);  //sigma_type==-1: nonweighted else: weighted exp(/200) or exp(/20)
#endif
