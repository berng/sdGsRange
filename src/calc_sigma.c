#include "calc_sigma.h"

int order_Ampl(const void* a,const void* b)
 {
  double aa,bb;
  aa=((dr_type*)a)->Ampl;
  bb=((dr_type*)b)->Ampl;
  if(aa>bb)
   return -1;
  else
   if(aa<bb)
    return 1;
   else
    return 0;
 }

int order_Ampl_inv(const void* a,const void* b)
 {
  double aa,bb;
  aa=((dr_type*)b)->Ampl;
  bb=((dr_type*)a)->Ampl;
  if(aa>bb)
   return -1;
  else
   if(aa<bb)
    return 1;
   else
    return 0;
 }

int order_dr(const void* a,const void* b)
 {
  double aa,bb;
  aa=((dr_type*)a)->dr;
  bb=((dr_type*)b)->dr;
  if(aa>bb)
   return 1;
  else
   if(aa<bb)
    return -1;
   else
    return 0;
 }

int order_drw(const void* a,const void* b)
 {
  double aa,bb;
  aa=((dr_type*)a)->dr_w;
  bb=((dr_type*)b)->dr_w;

  aa=(aa>0)?aa:-0.3*aa;
  bb=(bb>0)?bb:-0.3*bb;

  if(aa>bb)
   return 1;
  if(aa<bb)
   return -1;
  return 0;

 }


double calc_sigma(dr_type* drt,long len,dr_type* drt_sve,int sigma_type)
{
 long i;
 double sigma=0;
 double n=0;
//berng: not necessary: qsort((void*)drt,len,sizeof(dr_type),order_drw);
 for(i=0;i<len;i++)
  {
   drt_sve[i]=drt[i];
   drt_sve[i].dr=0;
   drt_sve[i].r=0;
   drt_sve[i].dr_w=0;
   drt_sve[i].Ampl=0;
   drt_sve[i].idx=0;
  }
//fprintf(stderr,"calc sigma(%ld)\n",len);
// n=1;
// for(sigma=0.0,i=0;i<len*3/4;i++)
// for(sigma=0.0,i=0;i<len*15/16/*/4*/;i++)
// for(sigma=0.0,i=0;i<len*9/10/*/4*/;i++)
 for(sigma=0.0,i=0;i<len/*4*/;i++)
 {
//  fprintf(stderr,"[%lf]...",drt[i].Ampl);
  if(drt[i].Ampl>0.0 /* && fabs(drt[i].dr_w)<500*/)
  {
   drt_sve[i]=drt[i];
#ifdef DEBUG
   fprintf(stderr,"drw: %lf\n",drt[i].dr_w);
#endif
   if(sigma_type==-1)
   {
    sigma+=drt[i].dr_w*drt[i].dr_w;
    //sigma+=1./(fabs(drt[i].dr_w)+100.0);
    n+=1.0;
   }
   else
   {
   if(drt[i].dr_w<0)
    {
     //sigma+=exp(-fabs(drt[i].dr_w)/500.0);
//     sigma+=1./(fabs(drt[i].dr_w)+45);
#ifdef _SIGMA_MAX
     sigma+=exp(-fabs(drt[i].dr_w)/200.0);
#else
     sigma+=fabs(drt[i].dr_w*drt[i].dr_w);
#endif
//     sigma+=exp(-fabs(drt[i].dr_w)/500.0);
     n+=1.0;
    }
   else
    {
//     sigma+=exp(-fabs(drt[i].dr_w)/200.0);
#ifdef _SIGMA_MAX
     sigma+=exp(-fabs(drt[i].dr_w)/20.0);
#else
     sigma+=fabs(drt[i].dr_w*drt[i].dr_w);
#endif
//    exp(-fabs(drt[i].dr_w)/500.0);
//     sigma+=exp(-fabs(drt[i].dr_w)/500.0);
//    sigma+=exp(-fabs(drt[i].dr_w)/50.0);
//    sigma+=exp(-fabs(drt[i].dr_w)/100.0);
//     sigma+=0;
//    sigma+=exp(-fabs(drt[i].dr_w)/300.0);
//   if(fabs(drt[i].dr_w)>300)
     n+=1.0;
    }
  }
//   n=1.0;
  }
  else
   n+=1;
 }
//fprintf(stderr,"n for sigma(%lf), sigma: %lf\n",n,sigma);
 if(n>0)
  {
#ifdef _SIGMA_MAX
//   sigma/=(double)n;
#endif
//   fprintf(stderr,"sig: %le\n",sigma);
   sigma=sqrt(sigma);
   return sigma;
  }
 else
  if(sigma_type==-1)
   return 1e100;
  else
   return 0.0;
//  return 1e100;
}

/*
#define MAX_ARR_LEN 100000
main()
{
 double* src_x;
 double* src_y;
 double* src_det;
 double* dr;
 double r;
 long len;
 long i;
 double k;
 double old_sigma,opt_k;
 double sigma;
 src_x=calloc(sizeof(double),MAX_ARR_LEN);
 src_y=calloc(sizeof(double),MAX_ARR_LEN);
 src_det=calloc(sizeof(double),MAX_ARR_LEN);
 dr=calloc(sizeof(double),MAX_ARR_LEN);

 for(i=0;!feof(stdin) && i<100000;i++)
  {
   fscanf(stdin,"%lf%lf%lf",src_x+i,src_y+i,src_det+i);
  }
 len=i-1;
 for(i=0;i<len;i++)
  dr[i]=0;
 for(old_sigma=1e50,k=0;k<3;k+=0.1)
  {
   for(i=0;i<len;i++)
    {
     dr[i]=fabs(src_x[i]*k-src_y[i]);
    }
   sigma=calc_sigma(dr,len);
   printf("sigma[%lf]: %lf\n",k,sigma);

   if(sigma<old_sigma)
    {
     old_sigma=sigma;
     opt_k=k;
     printf("OPT sigma[%lf]: %lf\n",k,sigma);
    }
  }
// sigma=calc_sigma(src,len);
 printf("sigma: %lf\n",sigma);
}
*/

