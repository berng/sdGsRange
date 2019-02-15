#include "get_fof2_iri.h"
#define LN(x) log(x)
void get_muf2_iri(int yyyy,
		  int mm, 
		  int dd,
		  int hh,
		  int mi,
		  double src_lat,
		   double src_long,
		   double dest_lat,
		   double dest_long,
		   double* resfof2,
		   double* hmax,
		   double* hmin)
{


 int n=9;
 long offset=20;

 int i,j;

 char date_str[100];
 double time_ut;
 double r;
 double mof_exp,mof_moh,mof_gnss;
 
 long history_len;
 long k;
 long HistoryLen;
 long date;
 double x1,x2;



 iri_conf IRIconf;
 IRIconf.NeComputed=_FTRUE;  	 	//Ne computed		 Ne not computed                     t
 IRIconf.TeTiComputed=_FTRUE;    		//Te, Ti computed        Te, Ti not computed                 t	   
 IRIconf.NeNiComputed=_FTRUE;	 	//Ne & Ni computed       Ni not computed                     t
 IRIconf.MagneticFieldByTable=_FTRUE;	//B0 - Table option      B0 - Gulyaeva (1987)                t
 IRIconf.foF2_CCIR=_FFALSE;		//foF2 - CCIR            foF2 - URSI                     false
 IRIconf.Ni_old=_FFALSE;    		//Ni - DS-78 & DY-85     Ni - DS-95 & TTS-03             false
 IRIconf.Ne_tops=_FTRUE;    		//Ne - Tops: f10.7<188   f10.7 unlimited                     t            
 IRIconf.foF2_model=_FTRUE;    		//foF2 from model        foF2 or NmF2 - user input           t
 IRIconf.hmF2_model=_FTRUE;    		//hmF2 from model        hmF2 or M3000F2 - user input        t
 IRIconf.TeStandard=_FTRUE;    		//Te - Standard          Te - Using Te/Ne correlation        t
 IRIconf.NeStandard=_FTRUE;    		//Ne - Standard Profile  Ne - Lay-function formalism         t
 IRIconf.MessagesToUnit6=_FFALSE; 		//Messages to unit 6     no messages                         t
 IRIconf.foF1_model=_FTRUE;    		//foF1 from model        foF1 or NmF1 - user input           t
 IRIconf.hmF1_model=_FTRUE;    		//hmF1 from model        hmF1 - user input (only Lay version)t
 IRIconf.foE_model=_FTRUE;    		//foE  from model        foE or NmE - user input             t
 IRIconf.hmE_model=_FTRUE;    		//hmE  from model        hmE - user input                    t
 IRIconf.Rz12_from_file=_FTRUE;  		//Rz12 from file         Rz12 - user input                   t
 IRIconf.IGRF_dip=_FTRUE;  		//IGRF dip, magbr, modip old FIELDG using POGO68/10 for 1973 t
 IRIconf.F1_probability_model=_FTRUE;	//F1 probability model   critical solar zenith angle (old)   t
 IRIconf.standard_F1=_FTRUE;		//standard F1            standard F1 plus L condition        t
 IRIconf.IonDriftComputed=_FFALSE;		//ion drift computed     ion drift not computed          false
 IRIconf.IonDensitiesInPercent=_FTRUE;	//ion densities in %     ion densities in m-3                t
 IRIconf.Te_tops=_FFALSE;			//Te_tops (Aeros,ISIS)   Te_topside (Intercosmos)        false
 IRIconf.Dregion_IRI95=_FTRUE;		//D-region: IRI-95       Special: 3 D-region models          t
 IRIconf.F107D_from_AP_DAT=_FTRUE;	//F107D from AP.DAT      F107D user input (oarr(41))         t
 IRIconf.foF2_storm_model=_FTRUE;		//foF2 storm model       no storm updating                   t
 IRIconf.IG12_from_file=_FTRUE;		//IG12 from file         IG12 - user input		     t
 IRIconf.spreadF_probability=_FFALSE;	//spread-F probability 	 not computed                    false
 IRIconf.IRI01_topside=_FFALSE;		//IRI01-topside          new options as def. by JF(30)   false
 IRIconf.IRI01_topside_corr=_FFALSE;	//IRI01-topside corr.    NeQuick topside model   	 false 

 oarr_type OARR;
 outf_type* OUTF;
 OUTF=(outf_type*)calloc(100,sizeof(outf_type));
 int JMAG;
 float ALATI;
 float ALONG;

 int IYYYY=yyyy;
 int MMDD=mm*100+dd;

 float DHOUR=(float)hh+(float)mi/60.;
 float HEIBEG=80.;
 float HEIEND=400.;
 float HEISTP=10.;

 double OUTF2[20];
 outf_type OUTF3;

 ALATI=60;

 IYYYY=yyyy;
 MMDD=mm*100+dd;


 HEIBEG=250.;//Altitude;

// initialize_();
double Altitude=250;
 HEIBEG=Altitude;
 HEIEND=Altitude;
 HEISTP=1.;


long k2;
 

mof_gnss=mof_exp;
mof_moh=mof_exp;

  double foF2;


// ALONG=ALONG_SRC;
// ALATI=ALATI_SRC;

// double L,az,k_az,aoc;
// int KT=10;
// int NHA,N,M;
// float Lr,dLr;
//


  ALATI=src_lat;
  ALONG=src_long;
double foF2s[100];
double heights[100];
long M,Mmax;

for(i=0;i<100;i++)
 foF2s[i]=heights[i]=0.;
foF2=0.;
*hmax=0.;
Mmax=0;

  for(M=0,HEIBEG=40;HEIBEG<=500;HEIBEG+=10,M++)
  {
   MMDD=mm*100+dd;
   IYYYY=yyyy-22;
//   DHOUR=((float)hh+(float)mi/60.)+25.;

//   DHOUR=((float)hh+(float)mi/60.);
   DHOUR=((float)hh+(float)mi/60.)+(double)src_long/15.0;
   while(DHOUR>=24.)
    {
     DHOUR-=24.;
//     MMDD=mm*100+dd+1;
    }
   while(DHOUR<0)
    {
     DHOUR+=24.;
//     MMDD=mm*100+dd-1;
    }
    
//   DHOUR=((float)hh+(float)mi/60.);

   iri_fast_(&IRIconf,&JMAG,&ALATI,&ALONG,&IYYYY,&MMDD,&DHOUR,&HEIBEG,(void*)(&OUTF2),(&OARR));
   OUTF3=((outf_type*)(&OUTF2))[0];

   foF2s[M]=(double)OUTF3.Ne_m3;
   heights[M]=HEIBEG;

   if(foF2<(double)OUTF3.Ne_m3)
    {
     foF2=(double)OUTF3.Ne_m3;
     *hmax=HEIBEG;
     Mmax=M;
    }
  }

int N=3;
double SRC_M[Len1][Len1];
double SRC_M1[Len1];
double SRC_M2[Len1];
double SRC_B[Len1];
double res;
double dh;

  for(i=0;i<N;i++)
  {
   for(j=0;j<N;j++)
    SRC_M[i][j]=0.;
   SRC_M1[i]=SRC_M2[i]=SRC_B[i]=0.;
  }

 for(M=0;M<=Mmax;M++)
  {
   if(heights[M]<=*hmax
     &&
      foF2s[M]>foF2*0.25
     )
    {
//     printf("found: %ld\n",M);
     dh=fabs(*hmax-heights[M])/100.;
     for(i=0;i<N;i++)
     {
      for(j=0;j<N;j++)
        SRC_M[i][j]+=pow((double)dh,(double)(i+j));
      SRC_B[i]+=pow(dh,(double)i)*foF2s[M];
     }
    }
  }

//     for(i=0;i<N;i++)
//     {
//      for(j=0;j<N;j++)
//       printf("%lf ",SRC_M[i][j]);
//      printf("|%lf\n",SRC_B[i]); 
//     }
//
 QRdecomp(SRC_M,SRC_M1,SRC_M2,N);
 QRsolve(SRC_M,SRC_M1,SRC_M2,SRC_B,N);
 for(M=0;M<=Mmax;M++)
  {
   dh=fabs(*hmax-heights[M])/100.;
   for(res=0.,i=0;i<N;i++)
    res+=SRC_B[i]*pow(dh,(double)i);
   if(res>0) 
    {
     *hmin=heights[M];
     break;
    }
    
//    fprintf(stderr,"%lf %lf %lf\n",heights[M],foF2s[M],res);
  }
//  exit(1);
//fprintf(stderr,"here we are5.... \n");

//berng: fixed thickness
// *hmax/=1.5;
 *hmax=300;
 *hmin=100;
// *hmin=*hmax-100;

 if(!(*hmin>110))
  *hmin=110;
 *resfof2=(double)(sqrt(foF2/1.24e10));
//  fprintf(stderr,"here we are: %d %d %f (%f,%f).. %lf %f ..\n",IYYYY,MMDD,DHOUR,ALATI,ALONG,*resfof2,HEIBEG);
  return;
}



/*
int main(int argn,char* argv[])
{

 initialize_();
int yyyy=2014;
int mm=6;
int dd=15;
int hh;
int mi=0;
double foF2,hmax,hmin;


//for(mm=1;mm<=12;mm++)
// for(dd=1;dd<=28;dd++)
mm=6;
dd=1;
double Ch;
double Cf;

Ch=1.1;
Cf=0.5;
  for(hh=0;hh<24;hh++)
   for(mi=0;mi<60;mi+=15)
{
 float muf;

 get_muf2_iri(yyyy,mm,dd,hh,mi,(float)ALATI_SRC,
		  (float)ALONG_SRC,
		  (float)ALATI_DEST,
		  (float)ALONG_DEST,
		  &foF2,
		  &hmax,
		  &hmin
		  );

 printf("%f %f %f %f %f %f %f %f\n",((float)mm*31.+(float)dd)+(float)hh/24.+(float)mi/(24.*60.),hmin,hmax,foF2,Distance_Range(hmin,hmax,foF2,10.),Range(hmin*Ch,hmax*Ch,foF2*Cf,10.),approx_Distance_Range(foF2),approx_Range(foF2));
}




}


*/
