#include "genetic_search.h"

void ConvertGenesToVars(animal_type* animal)
 {
  animal->foF2max=1.0+(double)((animal->genes)&0xFF)/8.0;
//  animal->foF2min=1.0+(animal->foF2max-1.0)*(double)(((animal->genes>>8)&0xF))/16.0;
  animal->foF2min=animal->foF2max/16.0+(animal->foF2max*8.0/16.0)*(double)(((animal->genes>>8)&0xF))/16.0;
//  animal->betta=1.0+(double)((animal->genes>>12)&0xF)*2.0;  
  animal->betta=1.0+(double)((animal->genes>>12)&0xF)/4.0;  
  animal->alpha=((double)((animal->genes>>16)&0xF)-8.0)/8.0;  
//  animal->dt=(double)((animal->genes>>20)&0xF)/2.0-4.0;
  animal->dt=-(double)((animal->genes>>20)&0xF)*3.0/16.0;
//  animal->dt=0;
  animal->sigma=-1.0;
 }

void CreateInititalGeneration(animal_type* animals, long animals_number)
 {
  long i;
  srand(time(NULL));
  for(i=0;i<animals_number;i++)
   {
    animals[i].genes=RANDOM_GENE;
//    animals[i].sigma=-1;
    ConvertGenesToVars(animals+i);
   }
 }

void  CalcSigma(double Xlat0,double Xlong0,float DayNo,animal_type* animals,long animal_num,rec_type* new_history,long history_len,dr_type* drt,dr_type* drt_save,int sigma_type)
{
 long i,j,k;
 static long init2=0;
 long calc_drt_len=0;
 double R;
 rec_type fixed_params;
 P_MAX=0;
#ifdef DEBUG
 fprintf(stderr,"calcing sigma %ld %ld\n",animal_num,history_len);
#endif

 for(init2=0,P_MAX=0.,k=0;k<animal_num;k++)
 {
  if(animals[k].sigma<0 || k==0)
  {
#ifdef DEBUG
   fprintf(stderr,"check animal %ld...\n",k);
#endif

   for(j=i=0;i<history_len;i++/*,j++*/)
    {
#ifdef DEBUG
     if(k==0)
      fprintf(stderr,"fixing params\n");
#endif

     fixed_params=FixParams(Xlat0,Xlong0,DayNo,new_history,i,animals[k].foF2max,animals[k].foF2min,animals[k].betta,animals[k].alpha,animals[k].dt,history_len);
#ifdef DEBUG
     if(k==0)
      fprintf(stderr,"fixing range\n");
#endif
     R=Range(fixed_params.hmin,fixed_params.hmax,fixed_params.foF2,fixed_params.f0);
#ifdef DEBUG
     if(k==0)
      fprintf(stderr,"saving...\n");
#endif

     drt[j].dr=fabs(R-new_history[i].range);
     drt[j].r=R;
     drt[j].dr_w=R-new_history[i].range;
#ifdef DEBUG
     if(k==0)
      fprintf(stderr,"drw[%ld]: %lf %lf %lf\n",j,(double)R,(double)new_history[i].range,R-new_history[i].range);
#endif
     drt[j].Ampl=(new_history[i].P>P_MAX)?new_history[i].P:0.0;
     drt[j].idx=i;
     if(drt[j].Ampl>0.0)
      {
//       fprintf(stderr,"found: %ld\n",j);
       j++;
      }
    }
   if(init2==0)
    {
#ifdef DEBUG
   fprintf(stderr,"drt filled, qsroting...\n");
#endif

     long l;
//if(sigma_type==-1)
//     qsort((void*)drt,j,sizeof(dr_type),order_Ampl_inv);
//else
     qsort((void*)drt,j,sizeof(dr_type),order_Ampl);
#ifdef DEBUG
   fprintf(stderr,"drt filled, qsroted\n");
#endif

     calc_drt_len=j;
/*
     if(sigma_type==-1)
      P_MAX=drt[1].Ampl;
     else
      P_MAX=drt[calc_drt_len-1].Ampl;
*/   P_MAX=0;
//     fprintf(stderr,"PMAX: %lf cdl: %ld\n",P_MAX,calc_drt_len);
     init2=1;
    }
//   fprintf(stderr,"drt filled, qsroted\n");
#ifdef DEBUG
   fprintf(stderr,"calcing sigma %ld...\n",k);
#endif

   animals[k].sigma=calc_sigma(drt,calc_drt_len,drt_save,sigma_type);
   animals[k].calc_drt_len=calc_drt_len;
#ifdef DEBUG
   fprintf(stderr,"calced sigma: %lf for (%lf,%lf,%lf,%lf,%lf) %ld,%ld\n",animals[k].sigma,animals[k].foF2max,animals[k].foF2min,animals[k].betta,animals[k].alpha,animals[k].dt,calc_drt_len,history_len);
#endif
  }
 }
}

int order_Animal(const void* a,const void* b)
 {
  double aa,bb;
  aa=((animal_type*)a)->sigma;
  bb=((animal_type*)b)->sigma;
  if(aa>bb)
   return -1;
  else
   if(aa<bb)
    return 1;
   else
    return 0;
 }

int order_Animal_inv(const void* a,const void* b)
 {
  double aa,bb;
  aa=((animal_type*)b)->sigma;
  bb=((animal_type*)a)->sigma;
  if(aa>bb)
   return -1;
  else
   if(aa<bb)
    return 1;
   else
    return 0;
 }


void OrderAnimals(animal_type* animals,long animal_num,int sigma_type)
 {

  if(sigma_type==-1)
//   qsort((void*)animals,animal_num,sizeof(*animals),order_Animal_inv);
   qsort((void*)animals,animal_num,sizeof(animal_type),order_Animal_inv);
//   qsort((void*)animals,animal_num,sizeof(*animals),order_Animal);
  else
//   qsort((void*)animals,animal_num,sizeof(*animals),order_Animal);
   qsort((void*)animals,animal_num,sizeof(animal_type),order_Animal);

 }

long GetChild(long genes1,long genes2)
 {
  long ret;
  long mix;
  mix=RANDOM_GENE;
  ret = (mix&genes1) | ((~mix)&genes2);  
  long modify;
  modify=(rand()%5);
  mix=1<<(rand()%24);
  long i;
  for(i=0;i<modify;i++)
   mix |= 1<<(rand()%24);
/*   
  mix |= 1<<(rand()%24);
  mix |= 1<<(rand()%24);
  mix |= 1<<(rand()%24);
*/
  ret ^= mix;
  return ret;
 }

void  GenerateNewAnimals(animal_type* animals,long animal_num,long best_animal_num)
{
 long i,parent1,parent2;
/*
 for(i=best_animal_num+1;i<animal_num;i++)
  {
   parent1=rand()%(best_animal_num-1);
   parent2=parent1+1+(rand()%(best_animal_num*2-parent1));
   animals[i].genes=GetChild(animals[parent1].genes,animals[parent2].genes);
//   fprintf(stderr,"mix %ld(%ld) %ld(%ld) res %ld(%ld)\n",parent1,animals[parent1].genes,parent2,animals[parent2].genes,i,animals[i].genes);
   ConvertGenesToVars(animals+i);
  }
*/
// for(i=animal_num-1;i>best_animal_num;i--)
 for(i=animal_num-1;i>1;i--)
  {
   parent1=rand()%(best_animal_num-1);
   parent2=rand()%(best_animal_num-1);
   animals[i].genes=GetChild(animals[parent1].genes,animals[parent2].genes);
//   fprintf(stderr,"mix %ld(%ld) %ld(%ld) res %ld(%ld)\n",parent1,animals[parent1].genes,parent2,animals[parent2].genes,i,animals[i].genes);
   ConvertGenesToVars(animals+i);
  }

}

void  PrintAnimals(animal_type* animals,long animal_num)
 {
  long i;
  for(i=0;i<animal_num;i++)
   {
    printf("%ld: %ld %lf\t",i,animals[i].genes,animals[i].sigma);
//    ConvertGenesToVars(animals+i);
//    printf("%lf\t",i,animals[i].genes,animals[i].sigma);
   }
  printf("\n");
  fflush(stdout);
 }
