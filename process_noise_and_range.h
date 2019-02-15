#ifndef __PROCESS_NOISE_AND_RANGE_H__
#define __PROCESS_NOISE_AND_RANGE_H__
#include "calc_Distance.h"
#include "zenith_angle.h"
#include "get_fof2_iri.h"
#include "FixParams.h"

typedef struct {
 double foF2;
 double zenith_cos;
 double hmin;
 double hmax;
 double t;
} model_type;

#define MAX_HISTORY_LEN 1000000
#define MAX_HISTORY_RECORD_NUM 10

typedef struct{
float foF2;
float hmax;
float hmin;
float t;
} profile_type;

#include "genetic_search.h"


#endif
