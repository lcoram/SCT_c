#include <stdio.h>
#include "sct_smart_boxes.h"

int main() {
   int n = 1;
   double x[10] = {1,2,3,4,5,6,7,8,9,10};
   double y[10] = {1,2,3,4,5,6,7,8,9,10};
   double z[10] = {1,1,3,1,1,1,1,10,100,1000};
   double t[10] = {1,1,11,4,5,10,100,100,10,10};
   int nmax = 100;
   int nmin = 10;
   int nminprof = 10;
   double dzmin = 100;
   double dhmin = 10000;
   double dz = 30;
   double t2pos[10] = {1,1,1,1,1,1,1,1,1,1};
   double t2neg[10] = {1,1,1,1,1,1,1,1,1,1};
   double eps2[10] = {1,1,1,1,1,1,1,1,1,1};
   int flags[10];
   double sct[10];
   double rep[10];

   for(n = 0; n < 10; n++) {
      sct_smart_boxes(&n, x, y, z, t, &nmax, &nmin, &nminprof, &dzmin, &dhmin, &dz, t2pos, t2neg, eps2, flags, sct, rep);
      printf("n=%d\n", n);
      for(int i = 0; i < n; i++) {
         printf("%d %f %f\n", flags[i], sct[i], rep[i]);
      }
      printf("\n");
   }

   return 0;
}
