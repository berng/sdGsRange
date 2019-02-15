#ifndef __GET_FOF2_IRI_H__
#define __GET_FOF2_IRI_H__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//#define KH_I
//#define M_I
//#define N_I
#include "main.h"

#include "iri2007.h"
#include "spherical.h"
#include "invmag.h"
//#include "libmuf.h"
#include "Solution.h"

#define _FTRUE  1
#define _FFALSE 0



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
		  double* hmin
		  
		  );


#endif
