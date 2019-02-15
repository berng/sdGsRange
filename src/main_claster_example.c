#include "claster_analysis.h"
long median_val=1;
int* range_arr;
long range_arr_len=0;
double x_width;
double y_width;
point* point_arr;
point_float* point_arr_float;


main(int argn,char* argv[])
{
long i,j,j2,k,l;
int col;
long arr_len=0;
pair* pair_arr;
long pair_arr_len=0;
long median_distance=0;
long* median_arr;
long median_distance_pos=-1;
vec r;

range_arr=calloc(LEN,sizeof(*range_arr));
point_arr=calloc(LEN,sizeof(*point_arr));
point_arr_float=calloc(LEN,sizeof(*point_arr_float));
median_arr=calloc(LEN,sizeof(*median_arr));
pair_arr=calloc((long)LEN*(long)LEN,sizeof(pair));
if(!pair_arr)
 {
  fprintf(stderr,"no space for pair_arr\n");
  exit(1);
 }


FILE* fin;
fin=fopen(argv[1],"rt");

//read
double min_y,min_x;
double max_y,max_x;

 min_y=min_x=1e100;
 max_y=max_x=-1e100;
 float xx,yy;
//fill initial array
 for(k=i=0;!feof(fin);i++) 
  {
   fscanf(fin,"%f%f",&xx,&yy);
   point_arr_float[k].r.x=(float)xx;
   point_arr_float[k].r.y=(float)yy;
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

 fclose(fin);

//calculate scales(cell size) for matrix
 x_width=(max_x-min_x)/100.0;
 y_width=(max_y-min_y)/100.0;
/*
 x_width=0.05;
 y_width=15;
*/

//fill matrix - convert float coordinates to integer ones
 arr_len=GenerateMatrix(point_arr_float,full_arr_len,min_x,x_width,min_y,y_width,point_arr);

 fprintf(stderr,"xw:%f yw:%f\n",x_width,y_width);
 fprintf(stderr,"arr len %ld\n",arr_len);

//fill pairs array
 pair_arr_len=FillPairsArray(point_arr,arr_len,pair_arr);
//calc min len for each element to its neighbour
 median_distance=GetMedianDistance(point_arr,arr_len,median_arr);
//find compact pairs with distance between points smaller than median distance
 median_distance_pos=GetCompactPairs(pair_arr,pair_arr_len,median_distance);

 fprintf(stderr,"median distance: %ld\n",median_distance);
 fprintf(stderr,"median distance pos: %ld\n",median_distance_pos);

//Coloring graph
int max_color;
long found_1=0;
for(;found_1>-1;)
{
 found_1=-1;
 for(max_color=-1,i=0;i<median_distance_pos;i++)
  {
   if(pair_arr[i].color>max_color)
    max_color=pair_arr[i].color;
   if(pair_arr[i].color==0 && found_1<0)
    found_1=i;
  }
 if(found_1>-1)
  {
   fprintf(stderr,"coloring by %ld \n",max_color+1);
   count_graph(pair_arr,median_distance_pos,median_distance,max_color+1,found_1);
   max_color+=1;
//   break;
  }
}
fprintf(stderr,"count graph finished\n");

//exit(1);
//coloring source array
ColorPointArr(pair_arr,median_distance_pos,point_arr);

//count length of each colored graph

long* colors;
colors=calloc(median_distance_pos,sizeof(*colors));
max_color=CountColors(point_arr,arr_len,median_distance_pos,colors);


//qsort(colors,median_distance_pos,sizeof(*colors),cmp_long);
/*
for(j=1;(j<=max_color);j++)
 fprintf(stderr,"%ld\n",colors[j]);
exit(1);
*/
//print colored points
fprintf(stderr,"max color: %ld\n",max_color);
 for(j=1;(j<=max_color);j++)
  if(colors[j]>1)
   {
    for(i=0;i<arr_len;i++)
     if(j==point_arr[i].color)
      {
       for(k=0;k<full_arr_len;k++)
        if(i==point_arr_float[k].point_idx)
          fprintf(stdout,"%f %f %ld %ld\n",point_arr_float[k].r.x,point_arr_float[k].r.y,point_arr[i].color,colors[j]);
      }
//      fprintf(stdout,"%d %d %ld\n",(int)point_arr[i].r.x,(int)point_arr[i].r.y,(long)point_arr[i].color);
    fprintf(stdout,"\n");
   }

 exit(1);


}

