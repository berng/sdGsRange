#ifndef __SOLUTION
#define __SOLUTION
#include <math.h>
#include <stdio.h>
#ifndef Len1
#define Len1 100
#endif

double __max(double M[Len1][Len1],int k,int n);
double __modul(double M[Len1][Len1],int k,int n);
double __Multiplicate(double M[Len1][Len1],int k,int n,int j);
int QRdecomp(double M[Len1][Len1],double M1[Len1],double M2[Len1],int n);
double __mul2(double M[Len1][Len1],double b[Len1],int n,int i);
void Rsolve(double M[Len1][Len1],double M2[Len1],double b[Len1],int n);
double __mul1(double M[Len1][Len1],double b[Len1],int n,int j);
void QRsolve(double M[Len1][Len1],double M1[Len1],double M2[Len1],double b[Len1],int n);
#endif
