#include "claster_analysis.h"

int cmp(const void* a,const void* b)
 {
  int ai,bi;
  ai=*((int*)a);
  bi=*((int*)b);
  if(ai>bi)
   return -1;
  else
   if(ai<bi)
    return 1;
   else
    return 0;
 }

int cmp_long(const void* a,const void* b)
 {
  long ai,bi;
  ai=*((long*)a);
  bi=*((long*)b);
  if(ai>bi)
   return -1;
  else
   if(ai<bi)
    return 1;
   else
    return 0;
 }

int cmp_median(const void* a,const void* b)
 {
  long ai,bi;
  ai=*((long*)a);
  bi=*((long*)b);
  if(ai>bi)
   return 1;
  else
   if(ai<bi)
    return -1;
   else
    return 0;
 }

int cmp_distance(const void* a,const void* b)
 {
  pair* ai;
  pair* bi;
  ai=(pair*)a;
  bi=(pair*)b;
  if(ai->distance>bi->distance)
   return 1;
  else
   if(ai->distance<bi->distance)
    return -1;
   else
    return 0;
 }


long calc_distance(vec r1,vec r2)
 {
  long dist;
  dist=(long)sqrt((double)(r1.x-r2.x)*(double)(r1.x-r2.x)+(double)(r1.y-r2.y)*(double)(r1.y-r2.y));
  return dist;
 }

void count_graph(pair* pair_arr,long median_distance_pos,long median_distance,int color,long pos)
 {
  long i,found;
  long cur_col;
  if(pair_arr[pos].color!=0)
   return;
//  fprintf(stderr,"count %ld with %ld\n",pos,color);
  pair_arr[pos].color=color;
  found=0;
  for(i=0;i<median_distance_pos;i++)
   {
    if(
     (
      pair_arr[i].idx1==pair_arr[pos].idx1
     ||
      pair_arr[i].idx2==pair_arr[pos].idx1
     || 
      pair_arr[i].idx1==pair_arr[pos].idx2
     ||
      pair_arr[i].idx2==pair_arr[pos].idx2
     )
     &&
      pair_arr[i].color==0
    )
     {
//       fprintf(stderr,"recursive call %ld with %ld\n",i,color);
 
      count_graph(pair_arr,median_distance_pos,median_distance,color,i);
      found=1;
     }
   }

 }



long GenerateMatrix(point_float* point_arr_float,long arr_len,double min_x,double x_width, double min_y,double y_width, point* point_arr)
{
 long i,j,l;
 long x_long,y_long;
 for(l=i=0;i<arr_len;i++)
 {
  x_long=(point_arr_float[i].r.x-min_x)/x_width;
  y_long=(point_arr_float[i].r.y-min_y)/y_width;
  for(j=0;j<l;j++)
   if(point_arr[j].r.x==x_long && point_arr[j].r.y==y_long)
    {
     point_arr_float[i].point_idx=j;
     break;
    }
  if(j==l)
   {
    point_arr[j].r.x=x_long;
    point_arr[j].r.y=y_long;
    point_arr[j].color=0;
    point_arr_float[i].point_idx=j;
    l++;
   }
 }
// arr_len=l;
 return l;
}

long FillPairsArray(point* point_arr, long arr_len, pair* pair_arr)
{
 long k,l,i;
 for(l=0,k=0;k<arr_len;k++)
  for(i=k+1;i<arr_len;i++)
   {
    pair_arr[l].idx1=k;
    pair_arr[l].idx2=i;
    pair_arr[l].color=0;
    pair_arr[l].distance=calc_distance(point_arr[k].r,point_arr[i].r);
    l++;
   }
 return l;
}

long GetMedianDistance(point* point_arr,long arr_len,long* median_arr)
{
 
 long min_len;
 long i,k;

 for(k=0;k<arr_len;k++)
 {
  min_len=100000;
  for(i=0;i<arr_len;i++)
   {
    if(calc_distance(point_arr[k].r,point_arr[i].r)<min_len 
	&& calc_distance(point_arr[k].r,point_arr[i].r)>0)
      min_len=calc_distance(point_arr[k].r,point_arr[i].r);
   }
  median_arr[k]=min_len;
 }

//Calculate median distance between points
 qsort(median_arr,arr_len,sizeof(*median_arr),cmp_median);
 long md;
 if(arr_len%2==1)
   md=median_arr[(arr_len-1)/2];
 else
   md=(median_arr[arr_len/2]+median_arr[arr_len/2-1])/2;
 return md*3;

}

long GetCompactPairs(pair* pair_arr,long pair_arr_len,long distance)
{
 long i;
 qsort(pair_arr,pair_arr_len,sizeof(*pair_arr),cmp_distance);
 for(i=0;i<pair_arr_len;i++)
  {
   if(pair_arr[i].distance>distance)
    {
     return i;
     break;
    }
  }
 return pair_arr_len-1;
}

void ColorPointArr(pair* pair_arr,long distance,point* point_arr)
 {
  long i;
  for(i=0;i<distance;i++)
   {
    point_arr[pair_arr[i].idx1].color=(pair_arr[i].color>0)?pair_arr[i].color:point_arr[pair_arr[i].idx1].color;
    point_arr[pair_arr[i].idx2].color=(pair_arr[i].color>0)?pair_arr[i].color:point_arr[pair_arr[i].idx2].color;
   }
 }

long CountColors(point* point_arr,long arr_len,long distance,long* colors)
{
 long max_color;
 long i;
 
 for(i=0;i<distance;i++)
  colors[i]=0;
 max_color=0;
 for(i=0;i<arr_len;i++)
   {
    colors[point_arr[i].color]++;
    max_color=(point_arr[i].color>max_color)?point_arr[i].color:max_color;
   }
 return max_color;
}

