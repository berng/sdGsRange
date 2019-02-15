#include "main.h"
#include "process_noise_and_range.h"
#include "calc_sigma.h"
#include "genetic_search.h"
#include "claster_analysis.h"

#include "spherical.h"
#define Re (6371.)
#define Pi 3.1415926
#define RAD2DEG (180./Pi)
#define DEG2RAD (Pi/180.)


/*
double H_MAX_EXTERN=350;
double DH_EXTERN=100;
*/
double H_MAX_EXTERN=450;
double DH_EXTERN=100;


double P_MAX=1.0;
long azimuth=0;
rec_type* history=NULL;
rec_type* new_history=NULL;
rec_type* history_sort=NULL;
long history_len=0;
long current_rec_id=0;

model_type model[10000];
long model_len=0;


int compare_records(const void* a, const void* b)
 {
  rec_type* a_rec;
  rec_type* b_rec;
  a_rec=(rec_type*)a;
  b_rec=(rec_type*)b;
  if(a_rec->range>b_rec->range)
   return -1;
  else
   if(a_rec->range<b_rec->range)
     return 1;
   else
    return 0;
 }


void fill_history_array(struct RadarParm prm,struct FitData fit,long select_beam)
{
//  double avrg2,avrg,avrg_n;
  long i;
//  long j,k,l;
//  double ampl=1;
//  long range=-1;
//  double xlat,xlong;
  current_rec_id++;
  if(history==NULL)
   {
    history=calloc(MAX_HISTORY_LEN,sizeof(rec_type));
    new_history=calloc(MAX_HISTORY_LEN,sizeof(rec_type));
    history_sort=calloc(MAX_HISTORY_LEN,sizeof(rec_type));
    if(history==NULL || new_history==NULL || history_sort==NULL)
      fprintf(stderr,"error allocating history, exiting...\n"),
      exit(1);
    history_len=0;
   }
long pos;
  for(i=0;i<prm.nrang && history_len<MAX_HISTORY_LEN;i++)
  {

   double r1;
   r1=(double)prm.frang+(double)i*(double)prm.rsep;
   if(fit.rng[i].p_l>3 && fit.rng[i].gsct==1 && (((long)prm.bmnum==select_beam && select_beam>-1) || select_beam<0 )/*&& (long)prm.channel==2*/  && r1<3500)
   {
#ifdef DEBUG
    fprintf(stderr,"fill. %ld (%ld) %p bm: %ld sbm:%ld\n",i,prm.nrang,fit.rng,prm.bmnum,select_beam);    
#endif    
    history[history_len].range=r1;
    if(fit.elv!=NULL)
    {
     history[history_len].elv_normal=fit.elv[i].normal;
     history[history_len].elv_low=fit.elv[i].low;
     history[history_len].elv_high=fit.elv[i].high;
    }
    else
     {
     history[history_len].elv_normal=-1;
     history[history_len].elv_low=-1;
     history[history_len].elv_high=-1;
     }
//    history[history_len].P=fit.rng[i].p_l;  //std
    history[history_len].P=fit.rng[i].p_l+2.*log(r1/200.0);   //with R^2
    history[history_len].rec_id=current_rec_id;
    history[history_len].f0=(double)prm.tfreq/1000.;
    history[history_len].t=(double)prm.time.hr+(double)prm.time.mt/60.;
    history[history_len].beam=(int)prm.bmnum;
    history[history_len].noise=(int)prm.noise.search;
    pos=prm.time.hr*4+prm.time.mt/15;
//    double offs;
//    offs=(double)(prm.time.mt%15);

//    fprintf(stderr,"fill2. %ld (%ld).",i,prm.nrang);

    if(pos<model_len)
    {
     if(pos<model_len-1)
      {
//       fprintf(stderr,"fill3. %ld (%ld).",i,prm.nrang);
       history[history_len].zenith_cos=model[pos].zenith_cos;//+(model[pos+1].foF2-model[pos].foF2)*(double)offs/15.;
       history[history_len].foF2=model[pos].foF2;//+(model[pos+1].foF2-model[pos].foF2)*(double)offs/15.;
       history[history_len].hmin=model[pos].hmin;//+(model[pos+1].hmin-model[pos].hmin)*(double)offs/15.;
       history[history_len].hmax=model[pos].hmax;//+(model[pos+1].hmax-model[pos].hmax)*(double)offs/15.;
      }
     else
      {
//       fprintf(stderr,"fill4. %ld (%ld).",i,prm.nrang);
       history[history_len].zenith_cos=model[pos].zenith_cos;//+(model[0].foF2-model[pos].foF2)*(double)offs/15.;
       history[history_len].foF2=model[pos].foF2;//+(model[0].foF2-model[pos].foF2)*(double)offs/15.;
       history[history_len].hmin=model[pos].hmin;//+(model[0].hmin-model[pos].hmin)*(double)offs/15.;
       history[history_len].hmax=model[pos].hmax;//+(model[0].hmax-model[pos].hmax)*(double)offs/15.;
      }
      
    }
//     fprintf(stderr,"%lf %lf %lf\n",history[history_len].t,history[history_len].foF2,model[pos].foF2);
    history_len++;
   }
  }
// fprintf(stderr,"history len: %ld\n",history_len);
// exit(1);
 return;
}




int main(int argc,char *argv[]) {
  FILE *fp;
  struct RadarParm prm;
  struct RadarParm *PRM;
  struct FitData fit;
  struct FitData *FIT;
  char* radarname;
double Xlat0=(double)56.5,Xlong0=(double)58.5;
float DayNo=0;
long SELECT_BEAM=0;


  
  if(argc<5)
   {
      fprintf(stderr,"usage %s fname radarlat radarlong beam(-1 or actual beam) radarname\n",argv[0]);
      exit(1);
   }
  fp=fopen(argv[1],"r");
  Xlat0=atof(argv[2]);
  Xlong0=atof(argv[3]);
  SELECT_BEAM=-1;
  if(argc>=5)
   SELECT_BEAM=atol(argv[4]);
  radarname=argv[5];
  H_MAX_EXTERN=atof(argv[6]);
  DH_EXTERN=atof(argv[7]);
  
#ifdef DEBUG
fprintf(stderr,"select beam %ld %ld\n",(long)SELECT_BEAM,(long)argc);
#endif

  if (fp==NULL) 
   {
    fprintf(stderr,"File %s not found.\n",argv[1]);
    exit(-1);
   }

  PRM=RadarParmMake();
  FIT=FitMake();
  prm=*PRM;
  fit=*FIT;
  int old_h=0;
  int yyyy,mm,dd,hh,mi;
  long i;

  int init=0;
  while(FitFread(fp,&prm,&fit) !=-1) 
    {
     if(init==0)
      {
      yyyy=prm.time.yr;
      mm=prm.time.mo;
      dd=prm.time.dy;
      DayNo=mm*30+dd;
      
      dd = 10;
//      dd++;
#ifdef DEBUG
      fprintf(stderr,"prm: %d-%d-%d %d\n",prm.time.yr,prm.time.mo,prm.time.dy,prm.time.hr);
#endif
      init++;
      model_len=0;
      for(hh=0;hh<24;hh+=1)
       for(mi=0;mi<60;mi+=15)
        {
         get_muf2_iri(yyyy,mm,dd,hh,mi,(double)Xlat0,(double)Xlong0,(double)Xlat0,(double)Xlong0,&(model[model_len].foF2),&(model[model_len].hmax),&(model[model_len].hmin));
//berng: fix it         
//         model[model_len].hmax/=1.5;
//         model[model_len].hmin=model[model_len].hmax-100.0;
// /berng

         model[model_len].zenith_cos=(double)get_cos_zenith((float)Xlat0,(float)Xlong0,(float)(hh+mi/60.),(float)(dd+mm*30));
         model[model_len].t=hh+mi/60.; //unnencessary, unused
         model_len++;//=i+1;
//         fprintf(stderr,"%d %d %lf %lf %lf %lf\n",hh,mi,(model[model_len-1].foF2),(model[model_len-1].hmax),(model[model_len-1].hmin),(model[model_len-1].zenith_cos));
        }
//       exit(1);
//      dd++;

      }
     if(prm.time.hr<old_h)
      break;
     else
      old_h=prm.time.hr;
     fill_history_array(prm,fit,SELECT_BEAM);
    }
    
 
double max_f=0.0,min_f=1e100,ref=0.5;
for(i=0;i<model_len;i++)
 {
  max_f=(max_f>(double)model[i].foF2 && (double)model[i].foF2>0)?max_f:(double)model[i].foF2;
  min_f=(min_f<(double)model[i].foF2 && (double)model[i].foF2>0)?min_f:(double)model[i].foF2;
#ifdef DEBUG
  fprintf(stderr,"%lf (%lf %lf)=>%lf\n",ref,min_f,max_f,model[i].foF2);
#endif
 }
ref=min_f/max_f;

#ifdef DEBUG
fprintf(stderr,"ref: %lf (%lf %lf)\n",ref,min_f,max_f);
#endif

ref*=1.2;
ref=(ref>0.7)?0.7:ref;

long jj,pos_max=-1;
FILE* stream_track;
stream_track=fopen("stream.track_1","wt");

  for(jj=0,i=1;i<history_len;i++)
   {
     double max_r,max_a;
      double elv_normal;
      double elv_low,elv_high;
//      int max_found=0;
      long j;
      if(history[i-1].t<history[i].t)
      {
       max_r=-1;//history[i].range;
       max_a=-1;//history[i].P;
       elv_normal=0;
       elv_low=elv_high=0;
       pos_max=-1;
       for(j=i;(history[j+1].range>history[j].range)&&(history[j].range<3500)&&j<history_len-1;j++)
        if(history[j].P>max_a)
        {
         max_r=history[j].range;
         max_a=history[j].P;
         elv_normal=history[j].elv_normal;
         elv_low=history[j].elv_low;
         elv_high=history[j].elv_high;
         pos_max=j;
        }
       if(pos_max>=0)
        {
         new_history[jj++]=history[pos_max];
	 fprintf(stream_track,"%lf %lf %ld %lf\n",new_history[jj-1].t,new_history[jj-1].range,new_history[jj-1].beam,new_history[jj-1].f0);

//         fprintf(stderr,"%ld max %ld: %lf %lf %lf\n",jj,i,max_r,max_a,history[pos_max].t);
        }
      }
    }
   history_len=jj;
fclose(stream_track);


//  double* dr;
  dr_type* drt;
  dr_type* drt_save;
  dr_type* drt_opt;
  double R;

//  dr=calloc(1000000,sizeof(double));
  drt=calloc(1000000,sizeof(dr_type));
  drt_save=calloc(1000000,sizeof(dr_type));
  drt_opt=calloc(1000000,sizeof(dr_type));
//  double sigma;
  double old_sigma=1e100;
  long j;


/*
  double opt_Cf,opt_Cpow,opt_Ch;
  double opt_dCf;
  double Cf,Cpow,dCf,Ch;
*/
  double dt=0;
  double opt_dt=0;

  double foF2max;
  double foF2min;
  double opt_foF2max;
  double opt_foF2min;
  double betta;
  double alpha;
  double opt_betta;
  double opt_alpha;
  

rec_type fixed_params;
#ifdef _SIGMA_MAX
old_sigma=0.;
#else
old_sigma=1e100;
#endif

//int init2=0;
dt=0;
alpha=0;
//foF2max=8.0;
foF2max=max_f;
foF2min=min_f;
betta=10.0;
//long calc_drt_len;

ref=0.9;
#define ANIMAL_NUMBER 100
#define BEST_ANIMAL_NUMBER 10

//find approximation by genetic algorithm
animal_type animals[ANIMAL_NUMBER];
long itter;
CreateInititalGeneration(animals,ANIMAL_NUMBER);
#define SIGMA_LEVEL_0 0
#define SIGMA_LEVEL_1 -1


#define GA_ITTERATIONS 200
#define GA_ITTERATIONS_1 400
for(itter=0;itter<GA_ITTERATIONS;itter++)
 {
  CalcSigma(Xlat0,Xlong0,DayNo,animals,ANIMAL_NUMBER,new_history,history_len,drt,drt_save,SIGMA_LEVEL_0);
  OrderAnimals(animals,ANIMAL_NUMBER,SIGMA_LEVEL_0);
#ifdef DEBUG
  fprintf(stderr,"itter: %ld, sigma: %lf\n",itter,animals[0].sigma);
#endif

  GenerateNewAnimals(animals,ANIMAL_NUMBER,BEST_ANIMAL_NUMBER);
  CalcSigma(Xlat0,Xlong0,DayNo,animals,ANIMAL_NUMBER,new_history,history_len,drt,drt_save,SIGMA_LEVEL_0);
  OrderAnimals(animals,ANIMAL_NUMBER,SIGMA_LEVEL_0);
//  PrintAnimals(animals,BEST_ANIMAL_NUMBER);
 }

//get best approximation by genetic algorithm
OrderAnimals(animals,ANIMAL_NUMBER,SIGMA_LEVEL_0);
CalcSigma(Xlat0,Xlong0,DayNo,animals,0,new_history,history_len,drt,drt_save,SIGMA_LEVEL_0);
long l;
for(l=0;l<animals[0].calc_drt_len;l++)
 drt_opt[l]=drt_save[l];
old_sigma=animals[0].sigma;
opt_foF2max=animals[0].foF2max;
opt_foF2min=animals[0].foF2min;
opt_betta=animals[0].betta;
opt_alpha=animals[0].alpha;
opt_dt=animals[0].dt;
//#ifdef DEBUG
fprintf(stderr,"opt: %lf foF2max:%lf foF2min:%lf betta:%lf alpha:%lf dt:%lf\n",old_sigma,opt_foF2max,opt_foF2min,opt_betta,opt_alpha,opt_dt);
//#endif
//exit(1);


//get best approximation curve
profile_type* profile;
profile_type* profile_src;
//double Ampl;
//double h_to_correct;
//double opt_h2corr;
//double opt_Ampl;
profile=calloc(sizeof(profile_type),history_len);
profile_src=calloc(sizeof(profile_type),history_len);
/*
fprintf(stderr,"calc FixParams\n");
fprintf(stdout,"calc FixParams %ld\n",history_len);
*/

for(i=0;i<history_len;i++)
 {
//  fprintf(stdout,"calcing FixParams %ld\n",i);
  fixed_params=FixParams(Xlat0,Xlong0,DayNo, new_history,i,opt_foF2max,opt_foF2min,opt_betta,opt_alpha,opt_dt,history_len);
  profile_src[i].foF2=fixed_params.foF2;
  profile[i].foF2=profile_src[i].foF2;
  profile_src[i].t=fixed_params.t;
  profile[i].t=profile_src[i].t;
  profile_src[i].hmax=fixed_params.hmax;
  profile[i].hmax=profile_src[i].hmax;
  profile_src[i].hmin=fixed_params.hmin;
  profile[i].hmin=profile_src[i].hmin;  
 }
/*
fprintf(stderr,"calced FixParams\n");
fprintf(stdout,"calced FixParams\n");
*/
/*
free(drt);
free(drt_save);
*/
// ====================== CLASTER ANALSIS AND FILTERING ========================================================

 long k;
 long arr_len=0;
 long pair_arr_len=0;

//read
 double min_y,min_x;
 double max_y,max_x;
 min_y=min_x=1e100;
 max_y=max_x=-1e100;

 double x_width;
 double y_width;
 


// float xx,yy;
//fill initial array
 point_float* point_arr_float;
 point* point_arr;
 pair* pair_arr;
 long median_distance=0;
 long* median_arr;
 
 long median_distance_pos=-1;
 



 point_arr_float=calloc(history_len,sizeof(*point_arr_float));
 if(!point_arr_float)
  fprintf(stderr,"no space to allocate point_arr_float, exit"),exit(1);

 for(k=i=0;i<history_len;i++) 
  {
   point_arr_float[k].r.x=(float)new_history[i].t;
   point_arr_float[k].r.y=(float)new_history[i].range;
   point_arr_float[k].beam=new_history[i].beam;
   point_arr_float[k].noise=new_history[i].noise;
   
   min_x=(point_arr_float[k].r.x<min_x)?point_arr_float[k].r.x:min_x;
   max_x=(point_arr_float[k].r.x>max_x)?point_arr_float[k].r.x:max_x;
   min_y=(point_arr_float[k].r.y<min_y)?point_arr_float[k].r.y:min_y;
   max_y=(point_arr_float[k].r.y>max_y)?point_arr_float[k].r.y:max_y;
   point_arr_float[k].color=0;
   k++;
  }
// arr_len=k-1;
 long full_arr_len;
 full_arr_len=k-1;

//calculate scales(cell size) for matrix
 x_width=(max_x-min_x)/100.0;
 y_width=(max_y-min_y)/100.0;

//fill matrix - convert float coordinates to integer ones
 point_arr=calloc(full_arr_len,sizeof(*point_arr));
 arr_len=GenerateMatrix(point_arr_float,full_arr_len,min_x,x_width,min_y,y_width,point_arr);

#ifdef DEBUG
 fprintf(stderr,"xw:%f yw:%f\n",x_width,y_width);
#endif

#ifdef DEBUG
 fprintf(stderr,"arr len %ld\n",arr_len);
#endif
if(arr_len==0)
 return 0;

//fill pairs array
 pair_arr=calloc((long)history_len*(long)history_len,sizeof(pair));
 if(!pair_arr)
  fprintf(stderr,"no space to allocate pair_arr, exit"),exit(1);
 pair_arr_len=FillPairsArray(point_arr,arr_len,pair_arr);

//calc min len for each element to its neighbour
 median_arr=calloc(arr_len,sizeof(*median_arr));
 median_distance=GetMedianDistance(point_arr,arr_len,median_arr);

//berng fix;
median_distance*=2;

//find compact pairs with distance between points smaller than median distance
 median_distance_pos=GetCompactPairs(pair_arr,pair_arr_len,median_distance);


//Coloring graphs that are close enough to model
 int max_color;
 long found_1=0;
    
    found_1=-1;
    for(max_color=-1,i=0;i<median_distance_pos;i++)
     {
      if(pair_arr[i].color>max_color)
       max_color=pair_arr[i].color;
      if(pair_arr[i].color==0 
        &&
         found_1<0 
	)
       {
        for(j=0;j<animals[0].calc_drt_len;j++)
        {
         k=drt_opt[j].idx;
         if(
           (
             point_arr_float[k].point_idx==pair_arr[i].idx1
           ||
            point_arr_float[k].point_idx==pair_arr[i].idx2
           )
           && new_history[k].P>P_MAX 
           )
          {
           R=Range(profile_src[k].hmin,profile_src[k].hmax,profile_src[k].foF2,new_history[k].f0);
           if(fabs(R-new_history[k].range)<y_width)
    	    {
#ifdef DEBUG
	     fprintf(stderr,"found: %lf (%lf %lf)\n",R,new_history[k].t,new_history[k].range);
#endif
             found_1=i;
             break;
            }
          }
        }
       }
     if(found_1>-1)
      {
#ifdef DEBUG
      fprintf(stderr,"coloring by %ld \n",max_color+1);
#endif
      count_graph(pair_arr,median_distance_pos,median_distance,max_color+1,found_1);
      max_color+=1;
      found_1=-1;
//   break;
     }
    }
#ifdef DEBUG
 fprintf(stderr,"count graph finished\n");
#endif
// exit(1);
 
#ifdef DEBUG
 fprintf(stderr,"median distance: %ld\n",median_distance);
 fprintf(stderr,"median distance pos: %ld\n",median_distance_pos);
#endif

//coloring source array
ColorPointArr(pair_arr,median_distance_pos,point_arr);

//count length of each colored graph
long* colors;
colors=calloc(median_distance_pos,sizeof(*colors));
max_color=CountColors(point_arr,arr_len,median_distance_pos,colors);

//print colored points
fprintf(stderr,"opt: %lf foF2max:%lf foF2min:%lf betta:%lf alpha:%lf dt:%lf\n",old_sigma,opt_foF2max,opt_foF2min,opt_betta,opt_alpha,opt_dt);
fprintf(stderr,"max color: %ld\n",max_color);
double elev;
long id_next_level;
rec_type* next_level_history;
next_level_history=calloc(MAX_HISTORY_LEN+1,sizeof(rec_type));
//FILE* stream_track;
stream_track=fopen("stream.track_2","wt");

 for(id_next_level=0,j=1;(j<=max_color);j++)
  if(colors[j]>1)
   {
    for(i=0;i<arr_len;i++)
     if(j==point_arr[i].color)
      {
       for(k=0;k<full_arr_len;k++)
        if(i==point_arr_float[k].point_idx)
          {
           float foF2,D;
//           calc_Distance((double)point_arr_float[k].r.y,(double)(new_history[k].f0),&D,&foF2);
           fixed_params=FixParams(Xlat0,Xlong0,DayNo,new_history,k,opt_foF2max,opt_foF2min,opt_betta,opt_alpha,opt_dt,history_len);
           R=Range(fixed_params.hmin,fixed_params.hmax,fixed_params.foF2,fixed_params.f0);
           calc_Distance((double)R,(double)(new_history[k].f0),&D,&foF2);
           if(R<D) 
            {
//             D=R;
             elev=0;
            }
           else
            elev=acos(D/R)*180./3.1415;

//           fprintf(stdout,"%lf %lf %lf\n",point_arr_float[k].r.x,point_arr_float[k].r.y,R);
//#ifdef DEBUG
	   fprintf(stream_track,"%lf %lf %lf %ld %lf %lf %lf %lf \t%lf\t%lf %lf\n",new_history[k].t,new_history[k].range,R,new_history[k].beam,new_history[k].f0,new_history[k].elv_normal,new_history[k].elv_low,new_history[k].elv_high,elev,D,R);
//#endif
	   next_level_history[id_next_level++]=new_history[k];
           }
      }
   }
fclose(stream_track);
#ifdef DEBUG
exit(1);
#endif

//fprintf(stderr,"History Len Level 2: %ld\n",id_next_level);

//find approximation by genetic algorithm - NEXT APPROXIMATION
history_len=id_next_level-1;
//animal_type animals[ANIMAL_NUMBER];
//long itter;
//berng:reply best itterations: 
CreateInititalGeneration(animals,ANIMAL_NUMBER);
//for(i=0;i<ANIMAL_NUMBER;i++)
// animals[i].sigma=-1;
for(old_sigma=1e100,itter=0;itter<GA_ITTERATIONS_1;itter++)
//for(old_sigma=0,itter=0;itter<GA_ITTERATIONS_1;itter++)
 {
  CalcSigma(Xlat0,Xlong0,DayNo,animals,ANIMAL_NUMBER,next_level_history,history_len,drt,drt_save,SIGMA_LEVEL_1);
//exit(1);

  OrderAnimals(animals,ANIMAL_NUMBER,SIGMA_LEVEL_1);
#ifdef DEBUG
  fprintf(stderr,"itter: %ld, sigma: %lf\n",itter,animals[0].sigma);
#endif

  GenerateNewAnimals(animals,ANIMAL_NUMBER,BEST_ANIMAL_NUMBER);
  CalcSigma(Xlat0,Xlong0,DayNo,animals,ANIMAL_NUMBER,next_level_history,history_len,drt,drt_save,SIGMA_LEVEL_1);
  OrderAnimals(animals,ANIMAL_NUMBER,SIGMA_LEVEL_1);
  if(old_sigma>animals[0].sigma+0.1)
//  if(old_sigma+0.1<animals[0].sigma)
   {
    itter=0;
   }
  
//  PrintAnimals(animals,BEST_ANIMAL_NUMBER);
  

//  if(old_sigma>animals[0].sigma+0.01)
  if(old_sigma<animals[0].sigma-0.01)
   {
    fprintf(stderr,"strange error, missed order check!!\n\b");
//    PrintAnimals(animals,BEST_ANIMAL_NUMBER);
    break;
   }
  if(animals[0].sigma>0)
   old_sigma=animals[0].sigma;
//  PrintAnimals(animals,BEST_ANIMAL_NUMBER);
 }

//get best approximation by genetic algorithm
OrderAnimals(animals,ANIMAL_NUMBER,SIGMA_LEVEL_1);
CalcSigma(Xlat0,Xlong0,DayNo,animals,0,next_level_history,history_len,drt,drt_save,SIGMA_LEVEL_1);
//long l;
for(l=0;l<animals[0].calc_drt_len;l++)
 drt_opt[l]=drt_save[l];
old_sigma=animals[0].sigma;
opt_foF2max=animals[0].foF2max;
opt_foF2min=animals[0].foF2min;
opt_betta=animals[0].betta;
opt_alpha=animals[0].alpha;
opt_dt=animals[0].dt;
//#ifdef DEBUG
//#endif
//exit(1);



double sigma_model,sigma_model_n;
stream_track=fopen("stream.track_3","wt");
//for(sigma_model_n=sigma_model=0.,k=0;k<full_arr_len;k++)
for(sigma_model_n=sigma_model=0.,k=0;k<history_len;k++)

 {
   float foF2,D;
   calc_Distance((double)point_arr_float[k].r.y,(double)(next_level_history[k].f0),&D,&foF2);
   fixed_params=FixParams(Xlat0,Xlong0,DayNo,next_level_history,k,opt_foF2max,opt_foF2min,opt_betta,opt_alpha,opt_dt,history_len);
   R=Range(fixed_params.hmin,fixed_params.hmax,fixed_params.foF2,fixed_params.f0);
   sigma_model+=(next_level_history[k].range-R)*(next_level_history[k].range-R);
   sigma_model_n+=1.0;
//   fprintf(stdout,"%lf %lf %lf %ld\n",next_level_history[k].t,next_level_history[k].range,R,next_level_history[k].beam);
   fprintf(stream_track,"%lf %lf %lf %ld %lf\n",next_level_history[k].t,next_level_history[k].range,R,next_level_history[k].beam,next_level_history[k].f0);
 }
fclose(stream_track);

if(sigma_model_n>1)
 sigma_model/=(sigma_model_n-1);
sigma_model=sqrt(sigma_model);

FILE* stream_solution;
stream_solution=fopen("solution.out","wt");
fprintf(stderr,"NEXT opt: %lf foF2max:%lf foF2min:%lf betta:%lf alpha:%lf dt:%lf\n",old_sigma,opt_foF2max,opt_foF2min,opt_betta,opt_alpha,opt_dt);
fprintf(stream_solution,"opt: %lf foF2max:%lf foF2min:%lf betta:%lf alpha:%lf dt:%lf sigma_model_km:%lf\n",old_sigma,opt_foF2max,opt_foF2min,opt_betta,opt_alpha,opt_dt,sigma_model);
fclose(stream_solution);


// exit(1);
#ifdef DEBUG
fprintf(stderr,"rewind...\n");
#endif

  rewind(fp);
  rec_type src_data;
#ifdef DEBUG
fprintf(stderr,"recalc...\n");
#endif

  while(FitFread(fp,&prm,&fit) !=-1) 
   {
#ifdef DEBUG
fprintf(stderr,"reading...\n");
#endif


   if(!(((long)prm.bmnum==SELECT_BEAM && SELECT_BEAM>-1) || SELECT_BEAM<0 ))
    continue;


//    if(prm.channel<=1)
//    {
    src_data.f0=(double)prm.tfreq/1000.;
    src_data.t=(double)prm.time.hr+(double)prm.time.mt/60.+(double)prm.time.sc/3600.;
    src_data.beam=(int)prm.bmnum;
    src_data.noise=(int)prm.noise.search;
    fixed_params=FixParams(Xlat0,Xlong0,DayNo,&src_data,0,opt_foF2max,opt_foF2min,opt_betta,opt_alpha,opt_dt,1);
    R=Range(fixed_params.hmin,fixed_params.hmax,fixed_params.foF2,fixed_params.f0);
    float foF2,D;
    double az,CN,CE;
    double k_az,aoc;
    az=(double)prm.bmazm;

//    float foF2_2,D_2;
    double CN_2,CE_2;
#ifdef DEBUG
fprintf(stderr,"calcing...\n");
#endif
    
    calc_Distance((double)R,(double)(src_data.f0),&D,&foF2);
//    az=(double)(src_data.beam-2)*3.;
#ifdef DEBUG
fprintf(stderr,"sphe...\n");
#endif

    Spherical_R_az(Xlat0,Xlong0,D*5./6.,az,&CN,&CE,&k_az,&aoc);
#ifdef DEBUG
fprintf(stderr,"sphe2...\n");
#endif
//    calc_Distance((double)R,(double)(src_data.f0),&D_2,&foF2_2);
//    az=(double)(src_data.beam-2)*3.;
    Spherical_R_az(Xlat0,Xlong0,D*1./6.,az,&CN_2,&CE_2,&k_az,&aoc);

#ifdef DEBUG
fprintf(stderr,"cidx %s...\n",radarname);
#endif
    unsigned long idx;
    unsigned long idx2;
    idx=(radarname[0]-(int)'a');
    idx=(idx<<5)+(radarname[1]-(int)'a');
    idx=(idx<<5)+(radarname[2]-(int)'a');
//    idx=(idx<<5)+src_data.beam;
    idx=(idx<<5);
    idx=(idx<<2);  //    idx=(idx<<2)+prm.channel;
    idx2=idx+1;
    idx2=(idx2<<4)+(long)(prm.tfreq/1000.-8);
    idx=(idx<<4)+(long)(prm.tfreq/1000.-8);

#ifdef DEBUG
fprintf(stderr,"indexing...\n");
#endif    
    
//    +src_data.beam+prm.channel+(long)(prm.tfreq/1000.)
    for(i=0;i<prm.nrang;i++)
    {
//     if(fit.rng[i].p_l>0 && fit.rng[i].gsct>0)
     elev=acos(D/R)*180./3.1415;
/* far range D-layer */
     if(fit.rng[i].p_l>0 && fit.rng[i].gsct>0)
      fprintf(stdout,"%lf %d %lf\t%lf %lf %lf\t%lf %lf %lf\t%lf %lf %lf\t%lf %lf\t%lf %ld %s %ld\t%lx\n",
            src_data.t,
            src_data.beam,
            src_data.noise,
            foF2,D,elev,
            ((src_data.noise>0)?(10.*log(src_data.noise)/log(10.)*sin(elev*3.1415/180.)):0.0),
//            10.*log(src_data.noise)/log(10.),
//            R*5./6.,
            R,
            (double)prm.frang+(double)i*(double)prm.rsep,
            	CN,
		CE,
		az,
		Xlat0,
		Xlong0,
		(double)prm.bmazm,
		(long)prm.channel,
		radarname,
		(long)(prm.tfreq/2000.)*2,
		idx
            );
    }
/* near range D-layer */
/*
      elev=acos(D/R)*180./3.1415;
      fprintf(stdout,"%lf %d %lf\t%lf %lf %lf\t%lf %lf %lf\t%lf %lf %lf\t%lf %lf\t%lf %ld %s %ld\t%lx\n",
            src_data.t,
            src_data.beam,
            src_data.noise,
            foF2,D,elev,
            10.*log(src_data.noise)/log(10.)*sin(elev*3.1415/180.),
//            10.*log(src_data.noise)/log(10.),
            R*1./6.,
            (double)prm.frang+(double)i*(double)prm.rsep,
            	CN_2,
		CE_2,
		az,
		Xlat0,
		Xlong0,
		(double)prm.bmazm,
		(long)prm.channel,
		radarname,
		(long)(prm.tfreq/2000.)*2,
		idx2
            );
*/

//    }
   }

  fclose(fp);


}

