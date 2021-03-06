#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_blas.h>
#include "sct_smart_boxes.h"

int main(int argc, const char* argv[])
{
  // Testing functions from main, but eventually have to call it from R code
  FILE *fp;
  const char *filename = "input_for_example.txt";
  if(argc == 2)
     filename = argv[1];
  fp = fopen(filename, "r");
  if(fp == NULL) {
      printf("could not open file: %s \n", filename);
      return 1;
  }
  printf("opened file: %s \n", filename);

  // arrays z (elevation), x,y easting and northing coords
  double *z; // array do not know size of yet
  double *x;
  double *y;
  double *t;
  double *t2pos;
  double *t2neg;
  double *eps2;
  z = malloc(sizeof(double) * 50000);
  x = malloc(sizeof(double) * 50000);
  y = malloc(sizeof(double) * 50000);
  t = malloc(sizeof(double) * 50000);
  t2pos = malloc(sizeof(double) * 50000);
  t2neg = malloc(sizeof(double) * 50000);
  eps2 = malloc(sizeof(double) * 50000);
  int n = 0;
  if (!z || !x || !y || !t) { /* If data == 0 after the call to malloc, allocation failed for some reason */
    printf("Error allocating memory \n");
    return 1;
  }

  int lenLine=150;
  char line[lenLine];
  char delim[] = ";";

  fgets(line,lenLine,fp); // Read header line
  while(fgets(line,lenLine,fp)) { // loop through lines in file
    //itot;pridtot;lon;lat;xtot;ytot;ztot;ttot;laftot;
    // printf("line: %s \n", line);
    char *ptr = strtok(line, delim);
    int j = 0; // index inside line
    bool mybox = false;
    while(ptr != NULL) { // break up line
      // if itot == 266 then in oslo box
      if(j==0 && (strcmp(ptr, "266")==0 || strcmp(ptr, " 266")==0 )) {
        mybox = true;
      }
      if(mybox) {
        if(j == 4) { //easting
          x[n] = atof(ptr);
        }
        if(j == 5) { //northing
          y[n] = atof(ptr);
        }
        if(j == 6) { // elevation
          // keep the elevation values for oslo stations
          z[n] = atof(ptr);
        }
        if(j == 7) { // temperature
          t[n] = atof(ptr);
          n++; // increment size of these arrays
        }
      }
      ptr = strtok(NULL, delim);
      j++;
    }
  }
  printf("Number of stations in box: %i \n", n);
  printf("Mean x: %f y: %f z: %f \n", mean(x,n), mean(y,n), mean(z,n));
  printf("Mean t: %f\n", mean(t,n));

  fclose(fp);

  // do some stuff with the box
  int maxNumStationsInBox = 1000;
  int minNumStationsInBox = 100;
  int nminprof = 100;
  double dzmin = 30;
  double dhmin = 10000;
  double dz = 200;

  // initial variables for VP (passed into function)
  // variables for SCT
  for(int i=0; i<n; i++) {
    t2pos[i] = 4;
    t2neg[i] = 8;
    eps2[i] = 0.5;
  }

  // allocate memory for the indices
  int *flags = malloc(sizeof(int) * n);
  double *rep = malloc(sizeof(double) * n);
  double *sct = malloc(sizeof(double) * n);
  int *boxids = malloc(sizeof(int) * n);

  sct_smart_boxes(&n, x, y, z, t, &maxNumStationsInBox, &minNumStationsInBox, &nminprof, &dzmin, &dhmin, &dz, t2pos, t2neg, eps2, flags, sct, rep, boxids);

  FILE *out1;
  out1 = fopen("output.txt", "w");
  // save the original box for plotting
  for(int j=0; j<n; j++) {
    char str_temp[10];
    char str[100];
    // 0
    sprintf(str_temp,"%f",x[j]);
    strcpy(str,str_temp);
    strcat(str,";");
    // 1
    sprintf(str_temp,"%f",y[j]);
    strcat(str,str_temp);
    strcat(str,";");
    // 2
    sprintf(str_temp,"%f",z[j]);
    strcat(str,str_temp);
    strcat(str,";");
    // 3
    sprintf(str_temp,"%f",t[j]);
    strcat(str,str_temp);
    strcat(str,";");
    // 4
    sprintf(str_temp,"%d",flags[j]);
    strcat(str,str_temp);
    strcat(str,";");
    // 5
    sprintf(str_temp,"%f",rep[j]);
    strcat(str,str_temp);
    strcat(str,";");
    // 6
    sprintf(str_temp,"%f",sct[j]);
    strcat(str,str_temp);
    strcat(str,";");
    // 7
    sprintf(str_temp,"%d",boxids[j]);
    strcat(str,str_temp);
    strcat(str,";");
    fputs(str,out1);
    fputs("\n",out1);
  }
  fclose(out1);

  free(x);
  free(y);
  free(z);
  free(t);
  free(t2pos);
  free(t2neg);

  return 0;
}
