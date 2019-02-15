#include "Solution.h"

double __max(double M[Len1][Len1],int k,int n)
{
 double maxim=0.0;
 double wrk;
  int i;
 for(maxim=0.0e-108,i=k;i<n;i++)
  maxim=((wrk=fabs(M[i][k]))>maxim)?wrk:maxim;
 return maxim;
}

double __modul(double M[Len1][Len1],int k,int n)
{
 double rez=0.0,wrk;
 int i;
 for(rez=0.0,i=k;i<n;i++)
  {
    wrk=M[i][k];
    rez+=wrk*wrk; 
  }
  return sqrt(rez); 
}
double __Multiplicate(double M[Len1][Len1],int k,int n,int j)
{
 double rez=0.0;
 int i;
 for(i=k;i<n;i++)
  rez+=M[i][k]*M[i][j];
 return rez;
}

int QRdecomp(double M[Len1][Len1],double M1[Len1],double M2[Len1],int n)
{
 int sing=0;
 int i,j,k;
 double nu,wrk,sigma,tou;

 for(k=0;k<n-1;k++)
   {
     nu=__max(M,k,n);
     if(nu<=0.0) /* singularity */
      {
        M1[k]=0.0;
        M2[k]=0.0;
        printf("singularity");
        sing=-1;
      }	
     else
      {
        for(i=k;i<n;i++)
         M[i][k]/=nu;
        wrk=__modul(M,k,n);
        sigma=(M[k][k]>0.0e-108)?wrk:-wrk;
        M[k][k]+=sigma;
        M1[k]=sigma*M[k][k];
        M2[k]=-nu*sigma;
      }	
    for(j=k+1;j<n;j++)
     {
       tou=__Multiplicate(M,k,n,j)/M1[k];
       for(i=k;i<n;i++)
        M[i][j]-=tou*M[i][k];
     }
   }


 if(M[n-1][n-1]<=0.0e-108)
  sing=-1;
 M2[n-1]=M[n-1][n-1];
return sing;
}

double __mul2(double M[Len1][Len1],double b[Len1],int n,int i)
{
 int j;
 double rez=0.0;
 for(rez=0.0,j=i+1;j<n;j++)
  rez+=M[i][j]*b[j];
 return rez;
}

void Rsolve(double M[Len1][Len1],double M2[Len1],double b[Len1],int n)
{
  int i;
  b[n-1]=b[n-1]/M2[n-1];
  for(i=n-2;i>=0;i--)
   b[i]=(b[i]-__mul2(M,b,n,i))/M2[i];
 return;
}

double __mul1(double M[Len1][Len1],double b[Len1],int n,int j)
{
 int i;
 double rez=0.0;
 for(rez=0.0,i=j;i<n;i++)
  rez+=M[i][j]*b[i];
 return rez;
}

void QRsolve(double M[Len1][Len1],double M1[Len1],double M2[Len1],double b[Len1],int n)
{
 int j,i;
 double tou;
 for(j=0;j<n-1;j++)
  {
    tou=__mul1(M,b,n,j)/M1[j];
    for(i=j;i<n;i++)
      b[i]-=tou*M[i][j];
  }
 Rsolve(M,M2,b,n);
 return;
}

/*
main()
{
 double M[Len1][Len1];
 double M1[Len1];
 double M2[Len1];
 double b[Len1];
 int n=7;
 int i,j;
 for(i=0;i<n;i++)
   for(j=0;j<n;j++)
     M[i][j]=(double)(360360.0/(double)(i+j+1));
  for(i=0;i<n;i++) 
    b[i]=0.0;
  b[0]=1.0;
 QRdecomp(M,M1,M2,n);
 QRsolve(M,M1,M2,b,n);
 for(i=0;i<n;i++)
  printf("%d %le\n",i,b[i]);
}
*/
