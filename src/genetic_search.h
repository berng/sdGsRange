#ifndef __GENETIC_SEARCH_H__
#define __GENETIC_SEARCH_H__
#include <stdlib.h>
#include <time.h>
#include "FixParams.h"
#include "calc_sigma.h"
#include "calc_Distance.h"
extern double P_MAX;

typedef struct {
 long int genes;
 double dt;
 double foF2max;
 double foF2min;
 double betta;
 double alpha;
 double sigma;
 long calc_drt_len;
}animal_type;

#define RANDOM_GENE ((long)rand()%(0xFFFFFF))

void ConvertGenesToVars(animal_type* animal);
void CreateInititalGeneration(animal_type* animals, long animals_number);
void CalcSigma(double Xlat0,double Xlong0,float DayNo,animal_type* animals,long animal_num,rec_type* new_history,long history_len,dr_type* drt,dr_type* drt_save,int sigma_type);
int  order_Animal(const void* a,const void* b);
int  order_Animal_inv(const void* a,const void* b);
void OrderAnimals(animal_type* animals,long animal_num,int sigma_type);
long GetChild(long genes1,long genes2);
void GenerateNewAnimals(animal_type* animals,long animal_num,long best_animal_num);
void PrintAnimals(animal_type* animals,long animal_num);
#endif
