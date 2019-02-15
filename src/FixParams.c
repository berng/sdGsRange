#include "FixParams.h"

rec_type FixParams(double Xlat0, double Xlong0, float DayNo, rec_type* history,long idx, double foF2max,double foF2min,double betta,double alpha,double dt,long history_len)
 {
   long i;
   static double max_hmax=0.;
   static double min_hmax=0.;
   static double max_foF2=0.;
   static double min_foF2=0.;
   float mod_norm;
   rec_type res;
   res=history[idx];
   mod_norm=model_normalized(Xlat0,Xlong0,history[idx].t+dt,DayNo,betta,alpha);   
   res.foF2=(foF2max-foF2min)*mod_norm+foF2min;
/*
   res.hmax=350.0;//(1-mod_norm*320./350.);
   res.hmin=res.hmax-100;
*/
/*
fprintf(stderr,"h:%lf dh:%lf\n",H_MAX_EXTERN,DH_EXTERN);
fprintf(stdout,"h:%lf dh:%lf\n",H_MAX_EXTERN,DH_EXTERN);
fflush(stderr);
fflush(stdout);
*/
   res.hmax=H_MAX_EXTERN;//(1-mod_norm*320./350.);
   res.hmin=res.hmax-DH_EXTERN;
   if(res.hmax>200)
    res.hmin=res.hmax-100;
   else
    res.hmin=100.0;
   return res;
 }


