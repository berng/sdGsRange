#include "calc_Distance.h"
/*
 Решение согласно статье 
В.А. Березовский, И.Д. Золотарев, А.А. Васенина, Ю.К. Свешников, ОПРЕДЕЛЕНИЕ ПАРАМЕТРОВ КВ РАДИОЛИНИИ ПО РЕЗУЛЬТАТАМ ВОЗВРАТНО-НАКЛОННОГО ЗОНДИРОВАНИЯ ИОНОСФЕРЫ//Вестн. Ом. ун-та. 2011. No 2. С. 98–102.
 в плоском случае и для параболического слоя
*/
// Прямая задача - Коээффициент для трансформации группового пуи в дальность по земле
double Distance_Range(double hmin,double hmax,double foF2,double f0)
 {
  double res;
  double ksi=hmin/hmax;
  res=ksi*foF2*foF2/*/(f0*f0)*/;
/*
  if(res<1)
   res=sqrt(1.-res);
  else
   res=-1;
*/
  return res;
 }

// Прямая задача - зависимость групповго пути от параметров параболического слоя
double Range(double hmin,double hmax,double foF2,double f0)
 {
  double res;
  double ksi=hmin/hmax;
#ifdef DEBUG
  fprintf(stdout,"ksi: %lf h: [%lf %lf]\n",ksi,hmin,hmax);
#endif
  res=(2.*hmax*sqrt(ksi)+(hmax-hmin)*LN((1+sqrt(ksi))/(1-sqrt(ksi))))/foF2/**f0*/;
#ifdef DEBUG
  fprintf(stdout,"Range: %lf\n",res*f0);
#endif
  return res*f0;
 }


// Прямая задача - для радара EKB средняя зависимость группового пути от foF2 (расчитана по модели IRI 2008,2010,2012гг).
float approx_Range(double foF2)
 {
//  return 648.411/(foF2+(-0.20275)+(-0.00596609)*foF2*foF2);
  return 648.411/(foF2+(-0.20275));
 }


// Прямая задача - для радара EKB средняя зависимость дальности от группового пути (расчитана по модели IRI 2008,2010,2012гг).

float approx_Distance_Range(double foF2)
 {
    return 0.579129*foF2*foF2;
 }


// Обратная задача - для радара EKB расчитывает foF2 по групповму пути и частоте заондирования (расчитана по модели IRI 2008,2010,2012гг).

float inv_RangeToFoF2(double R,double f0)
 {
  return (f0*648.411)/R + 0.20275;
 }


// Обратная задача - для радара EKB расчитывает дальность по групповму пути и частоте заондирования (расчитана по модели IRI 2008,2010,2012гг).
double calc_Distance(double R,double f0,float* D,float* foF2_ptr)
 {
  double foF2,d;

  foF2=inv_RangeToFoF2(R,f0);
  *foF2_ptr=foF2;
//printf("foF2:%lf\n",foF2);
  d=approx_Distance_Range(foF2)/(f0*f0);
//printf("d:%lf R:%lf\n",d,R);
  if(d<1)
   d=R*sqrt(1-d);
  else
   d=-1;
//printf("res d:%lf\n",d);
  *D=d;
  return d;
 }

