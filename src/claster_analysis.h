#ifndef __CLASTER_ANALYSIS_H__
#define __CLASTER_ANALYSIS_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define LEN 10000

typedef struct {
 float x;
 float y;
}vec_float;

typedef struct {
 long x;
 long y;
}vec;

typedef struct{
 vec r;
 int color;
}point;

typedef struct{
 vec_float r;
 int color;
 long point_idx;
 int beam;
 double noise;
}point_float;

typedef struct{
 int idx1,idx2;
 long distance;
 int color;
}pair;

int cmp(const void* a,const void* b);
int cmp_long(const void* a,const void* b);
int cmp_median(const void* a,const void* b);
int cmp_distance(const void* a,const void* b);
long calc_distance(vec r1,vec r2);

void count_graph(pair* pair_arr,long median_distance_pos,long median_distance,int color,long pos);
long GenerateMatrix(point_float* point_arr_float,long arr_len,double min_x,double x_width, double min_y,double y_width, point* point_arr);
long FillPairsArray(point* point_arr, long arr_len, pair* pair_arr);
long GetMedianDistance(point* point_arr,long arr_len,long* median_arr);
long GetCompactPairs(pair* pair_arr,long pair_arr_len,long distance);
void ColorPointArr(pair* pair_arr,long distance,point* point_arr);
long CountColors(point* point_arr,long arr_len,long distance,long* colors);

#endif
