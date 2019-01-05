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

// export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
// compile: gcc -L/user/local/lib SCT_rewrite.c -lm -lgsl -lgslcblas

/*
The command line I've used to run titan is:
$>./titan.R test/privateObs_20180326T08Z.txt ~/data/out.txt --input.files test/kdvh_Norway_wmo_2018032608.txt,test/kdvh_Norway_nowmo_nomet_2018032608.txt,test/kdvh_Norway_nowmo_met_2018032608.txt,test/smhi_Sweden_2018032608.txt -c --varname.elev elev,z,z,z,z --varname.value TA,TA,TA,TA,TA --prid 9,1,2,6,4 -v

which means that I'm using netatmo (privateObs_20180326T08Z.txt, prid set to 9); kdvh_Norway_wmo_ prid set to 1; kdvh_Norway_nowmo_nomet_ prid set to 2; ...
All the files I've used are available on github

header of input_for_example.txt
itot;pridtot;lon;lat;xtot;ytot;ztot;ttot;laftot;
266;9;10.698517;59.520877;-242857;-378966;38;-0.2;1;
the column names are identical to the titan (global) variables.
itot= identifier for the box (Oslo fjord has the box  label = 266)
pridtot= observation provider identifier (9,1,2,6,4)
lon/lat= longitude latitude
xtot/ytot= easting and northing coordinates in km (I've transformed lat/lon into x/y kilometric coordinates so to make easier to compute distances)
ztot= elevation
ttot= observation (temperature)
laftot= land area fraction (the user can decide whether to use it in the interpolation)
*/

// Cristian's functions
void spatial_consistency_test(double *t2, int *box, double *boxCentre, int *numStationsInBox,
                             double *x, double *y, double *z, double *t, double *vp); // this is the interface to R (then need pointers)?

int vertical_profile_optimizer(gsl_vector *input, double **data);
double vertical_profile_optimizer_function(const gsl_vector *v, void *data); // GSL format

void vertical_profile(int nz, double *z,
    double t0, double gamma, double a, double h0, double h1i, double *t_out);


// helpers
double compute_quantile(double quantile, double *array, int sizeArray);
double mean(const double *array, int sizeArray);
void print_vector(double *vector, int size);
void print_gsl_vector(gsl_vector *vector, int size);
void print_matrix(double **matrix, int rows, int columns);
void print_gsl_matrix(gsl_matrix *matrix, int rows, int columns);
void print_sub_gsl_matrix(gsl_matrix *matrix, int start, int stop);

gsl_matrix* inverse_matrix(const gsl_matrix *matrix);

int main()
{
  bool testing = true;
  if(testing == false) {
    // Testing functions from main, but eventually have to call it from R code
    FILE *fp;
    char *filename = "/home/louiseo/Documents/SCT/myrepo/input_for_example.txt";
    fp = fopen(filename, "r");
    if(fp == NULL) {
        printf("could not open file: %s \n", filename);
        return 1;
    }
    printf("opened file: %s \n", filename);

    // write just the oslo box stuff to a file for testing
    //FILE *out;
    //out = fopen("33.txt", "w");

    // arrays z (elevation), x,y easting and northing coords
    double *z; // array do not know size of yet
    double *x;
    double *y;
    double *t;
    z = malloc(sizeof(double) * 50000);
    x = malloc(sizeof(double) * 50000);
    y = malloc(sizeof(double) * 50000);
    t = malloc(sizeof(double) * 50000);
    int n = 0;
    if (!z || !x || !y || !t) { /* If data == 0 after the call to malloc, allocation failed for some reason */
      printf("Error allocating memory \n");
      return 1;
    }
    //memset(z, 0, sizeof(double)*50000);

    int lenLine=150;
    char line[lenLine];
    char delim[] = ";";

    while(fgets(line,lenLine,fp)) { // loop through lines in file
      //itot;pridtot;lon;lat;xtot;ytot;ztot;ttot;laftot;
      //printf("line: %s \n", line);
      char *ptr = strtok(line, delim);
      int j = 0; // index inside line
      bool mybox = false;
      while(ptr != NULL) { // break up line
        // if itot == 266 then in oslo box
        if(j==0 && (strcmp(ptr, "289")==0 || strcmp(ptr, " 289")==0 )) {
          mybox = true;
        }
        if(mybox) {
          if(j == 1) { //easting
            x[n] = atof(ptr);
          }
          if(j == 2) { //northing
            y[n] = atof(ptr);
          }
          if(j == 3) { // elevation
            // keep the elevation values for oslo stations
            z[n] = atof(ptr);
          }
          if(j == 4) { // temperature
            t[n] = atof(ptr);
            n++; // increment size of these arrays
          }
          /*
          fputs(ptr,out);
          if(j!=9) {
            fputs(";",out);
          } */
          //printf("%s\n", ptr);
        }
        ptr = strtok(NULL, delim);
        j++;
      }
    }
    printf("Number of stations in box: %i \n", n);
    printf("Mean x: %f y: %f z: %f \n", mean(x,n), mean(y,n), mean(z,n));
    printf("Mean t: %f\n", mean(t,n));

    fclose(fp);
    //fclose(out);

    // initial variables for VP (first guess, adjusted during optimization)
    double gamma = -0.0065;
    double a = 5;

    double meanT = mean(t,n);

    // make copies of z for the quantile computation
    double * z_temp1;
    double * z_temp2;
    z_temp1 = malloc(sizeof(double) * n);
    z_temp2 = malloc(sizeof(double) * n);
    // allocate for output
    double *t_out = malloc(sizeof(double) * n);
    for(int i=0; i<n; i++) {
      z_temp1[i] = z[i];
      z_temp2[i] = z[i];
      t_out[i] = -999;
    }

    double exact_p10 = compute_quantile(0.10, z_temp1, n);
    double exact_p90 = compute_quantile(0.90, z_temp2, n);
    free(z_temp1);
    free(z_temp2);

    // data (params) that needs to be passed into vp
    double nd = (double) n; // cast + resize to double
    // data (double *n, double *z, double *t, double *t_out)
    double * data[4] = {&nd, z, t, t_out};

    // vector (double t0, double gamma, double a, double h0, double h1i)
    // Starting point for optimization
    gsl_vector *input = gsl_vector_alloc(5);
    gsl_vector_set(input,0,meanT);
    gsl_vector_set(input,1,gamma);
    gsl_vector_set(input,2,a);
    gsl_vector_set(input,3,exact_p10);
    gsl_vector_set(input,4,exact_p90);
    printf ("Input vector set = t0: %.4f gamma: %.4f a: %.4f h0: %.4f h1i: %.4f\n",
            meanT, gamma, a, exact_p10, exact_p90);

    int status = vertical_profile_optimizer(input, data);
    printf("status optimizer: %d\n", status);
    printf ("t0: %.4f gamma: %.4f a: %.4f h0: %.4f h1i: %.4f\n",
              gsl_vector_get(input, 0),
              gsl_vector_get(input, 1),
              gsl_vector_get(input, 2),
              gsl_vector_get(input, 3),
              gsl_vector_get(input, 4));

    vertical_profile(n, z, gsl_vector_get(input, 0), gsl_vector_get(input, 1),
      gsl_vector_get(input, 2), gsl_vector_get(input, 3), gsl_vector_get(input, 4), t_out);
    // now have temperature profile (t_out)
    gsl_vector_free(input);

    for(int i=0; i<n; i++) {
      assert(t_out[i] !=-999);
    }
    printf("first values in z (after vp) %f %f %f %f %f %f\n", z[0], z[1], z[2], z[3], z[4], z[5]);
    printf("first values in t_out (after vp) %f %f %f %f %f %f\n", t_out[0], t_out[1], t_out[2], t_out[3], t_out[4], t_out[5]);

    // variables for SCT
    double t2 = 16; // input by user into SCT function (TITAN seems to use 16? Cristian said 25)

    // 266;-247429.070909252;-365421.526660934;3466;
    // 289;-16140.3479247747;-468700.431167886;90;
    //int box = 266;
    int box = 289;
    double boxCentre[2];
    boxCentre[0] = -16140.3479247747;
    boxCentre[1] = -468700.431167886;
    //boxCentre[0] = -247429.070909252;
    //boxCentre[1] = -365421.526660934;
    int numStationsInBox = n;
    printf("num stations: %d\n", n);

    // void spatial_consistency_test(int *t2, int *box, double *boxCentre, int *numStationsInBox,
                                  //double *x, double *y, double *z, double *vp)
    clock_t start = clock(), diff;
    spatial_consistency_test(&t2, &box, boxCentre, &numStationsInBox, x, y, z, t, t_out);
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("SCT end\n");
    printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);

    free(z);
    free(x);
    free(y);
    free(t);
    free(t_out);
    return 0;
  }
  else {
    FILE *fp;
    char *filename = "/home/louiseo/Documents/SCT/myrepo/testdata.txt";
    fp = fopen(filename, "r");
    if(fp == NULL) {
        printf("could not open file: %s \n", filename);
        return 1;
    }
    printf("opened file: %s \n", filename);

    int n = 10;
    double *z;
    double *x;
    double *y;
    double *t;
    double *vp;
    z = malloc(sizeof(double) * n);
    x = malloc(sizeof(double) * n);
    y = malloc(sizeof(double) * n);
    t = malloc(sizeof(double) * n);
    vp = malloc(sizeof(double) * n);
    if (!z || !x || !y || !t || !vp) { /* If data == 0 after the call to malloc, allocation failed for some reason */
      printf("Error allocating memory \n");
      return 1;
    }
    printf("allocated memory \n");

    int lenLine=150;
    char line[lenLine];
    char delim[] = ";";

    bool first = true;
    int size = 0;
    while(fgets(line,lenLine,fp)) { // loop through lines in file
      //x;y;z;t;vp
      //printf("line: %s \n", line);
      char *ptr = strtok(line, delim);
      int j = 0; // index inside line
      while(ptr != NULL && !first) { // break up line
        //printf("j: %d \n", j);
        if(j == 0) { //easting
          x[size] = atof(ptr);
        }
        if(j == 1) { //northing
          y[size] = atof(ptr);
        }
        if(j == 2) { // elevation
          z[size] = atof(ptr);
        }
        if(j == 3) { // temperature
          t[size] = atof(ptr);
        }
        if(j == 4) { // vertical profile
          vp[size] = atof(ptr);
          size++;
        }
        ptr = strtok(NULL, delim);
        j++;
      }
      first = false;
    }
    printf("Number of stations in box: %i \n", size);
    printf("Mean x: %f y: %f z: %f \n", mean(x,n), mean(y,n), mean(z,n));
    printf("Mean t: %f\n", mean(t,n));

    fclose(fp);

    // variables for SCT
    double t2 = 3; // input by user into SCT function (TITAN seems to use 16? Cristian said 25)
    int box = 1;
    double boxCentre[2];
    boxCentre[0] = 1;
    boxCentre[1] = 1;

    // void spatial_consistency_test(int *t2, int *box, double *boxCentre, int *numStationsInBox,
                                  //double *x, double *y, double *z, double *vp)
    clock_t start = clock(), diff;
    spatial_consistency_test(&t2, &box, boxCentre, &n, x, y, z, t, vp);
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("SCT end\n");
    printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);

    free(z);
    free(x);
    free(y);
    free(t);
    free(vp);
    return 0;
  }

}

//----------------------------------------------------------------------------//
int vertical_profile_optimizer(gsl_vector *input, double **data)
{
  // optimize inputs for VP (using Nelder-Mead Simplex algorithm)
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss;
  gsl_multimin_function vp_optim;

  int iter = 0;
  int status;
  double size;

  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (5);
  gsl_vector_set_all (ss, 1.0);

  /* Initialize method and iterate */
  vp_optim.n = 5;
  vp_optim.f = vertical_profile_optimizer_function;
  vp_optim.params = data;

  s = gsl_multimin_fminimizer_alloc (T, 5);
  gsl_multimin_fminimizer_set (s, &vp_optim, input, ss);
  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status)
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-2);

      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }
      //printf ("iter= %5d f()= %10.3e size= %.3f\n", iter, s->fval, size);
      /*
      printf ("t0: %.4f gamma: %.4f a: %.4f h0: %.4f h1i: %.4f\n",
              gsl_vector_get (s->x, 0),
              gsl_vector_get (s->x, 1),
              gsl_vector_get (s->x, 2),
              gsl_vector_get (s->x, 3),
              gsl_vector_get (s->x, 4));
              */
    }
  while (status == GSL_CONTINUE && iter < 100);

  gsl_vector_set(input,0,gsl_vector_get(s->x, 0));
  gsl_vector_set(input,1,gsl_vector_get(s->x, 1));
  gsl_vector_set(input,2,gsl_vector_get(s->x, 2));
  gsl_vector_set(input,3,gsl_vector_get(s->x, 3));
  gsl_vector_set(input,4,gsl_vector_get(s->x, 4));

  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
  return status;
}

//----------------------------------------------------------------------------//
/*
#+ cost function used for optimization of tvertprof parameter
tvertprof2opt<-function(par) {
te<-tvertprof(z=zopt,t0=par[1],gamma=par[2],a=par[3],h0=par[4],h1i=par[5])
return(log((mean((te-topt)**2))**0.5))
}
*/
// vector (double t0, double gamma, double a, double h0, double h1i)
// data (int n, double *z, double *t, double *t_out)
double vertical_profile_optimizer_function(const gsl_vector *v, void *data)
{
  // my input dat:a double * data[4] = {&nd, z, t, t_out};
  double **p = (double **)data;
  int n = (int) *p[0]; // is of type double but should be an int
  double *z = p[1];
  double *t = p[2];
  double *t_out = p[3];

  // the parameters to mess with
  double t0 = gsl_vector_get(v,0);
  double gamma = gsl_vector_get(v,1);
  double a = gsl_vector_get(v,2);
  double h0 = gsl_vector_get(v,3);
  double h1i = gsl_vector_get(v,4);

  //printf("t0: %f gamma: %f a: %f h0: %f h1i: %f\n", t0, gamma, a, h0, h1i);

  // give everything to vp to compute t_out
  vertical_profile(n, z, t0, gamma, a, h0, h1i, t_out);
  // RMS
  double total = 0;
  for(int i=0; i<n; i++) {
    total += pow((t_out[i]-t[i]),2);
    //printf("%f", t_out[i]);
  }
  //printf("\n");
  double value = log(pow((total / n),0.5));

  //printf("first values in t_out %f %f %f\n", t_out[0], t_out[1], t_out[2]);
  //printf("first values in t %f %f %f\n", t[0], t[1], t[2]);
  //printf("first values in z %f %f %f\n", z[0], z[1], z[2]);
  //printf("optimizer value: %f\n", value);
  return value;
}


/*
#+ vertical profile of temperature (Frei, 2014)
tvertprof<-function(z,t0,gamma,a,h0,h1i) {
# ref:
# Frei, C. (2014). Interpolation of temperature in a mountainous region
#  using nonlinear profiles and non‐Euclidean distances.
#  International Journal of Climatology, 34(5), 1585-1605.
# input
#  z= array. elevations [m amsl]
#  t0= numeric. temperature at z=0 [K or degC]
#  gamma=numeric. temperature lapse rate [K/m]
#  a= numeric. inversion (spatial) length
#  h0= numeric. z where inversion starts [m]
#  h1i= numeric. h0+h1i is z where inversion stops [m]
#       (Frei uses h1 directly, I use an increment to h0 so to avoid ending
#        up with h1<=h0 during the optimization)
# Output
#  t= array. temperature [K or degC]
*/
void vertical_profile(int nz, double *z,
    double t0, double gamma, double a, double h0, double h1i, double *t_out)
{
  /*
  t<-z
  t[]<-NA
  h1<-h0+abs(h1i)
  z.le.h0<-which(z<=h0)
  z.ge.h1<-which(z>=h1)
  z.in<-which(z>h0 & z<h1)
  if (length(z.le.h0)>0)
   t[z.le.h0]<-t0-gamma*z[z.le.h0]-a
  if (length(z.ge.h1)>0)
   t[z.ge.h1]<-t0-gamma*z[z.ge.h1]
  if (length(z.in)>0)
   t[z.in]<-t0-gamma*z[z.in]-a/2*(1+cos(pi*(z[z.in]-h0)/(h1-h0)))
  return(t)
  */
  double h1 = h0 + fabs(h1i); // h1<-h0+abs(h1i)
  //printf("h1: %f\n", h1);
  // loop over the array of elevations (z)
  // t_out is an empty array of length nz (e.g. the same length as z)
  for(int i=0; i<nz; i++) {
    // define some bools
    bool z_le_h0 = z[i] <= h0; // z.le.h0<-which(z<=h0)
    bool z_ge_h1 = z[i] >= h1; // z.ge.h1<-which(z>=h1)
    bool z_in = (z[i]>h0 && z[i]<h1); // z.in<-which(z>h0 & z<h1)
    if(z_le_h0) {
      // t[z.le.h0]<-t0-gamma*z[z.le.h0]-a
      t_out[i] = t0-gamma*z[i]-a;
      //printf("calling vp 1: %f %f\n", z[i], t_out[i]);
    }
    if(z_ge_h1) {
      // t[z.ge.h1]<-t0-gamma*z[z.ge.h1]
      t_out[i] = t0-gamma*z[i];
      //printf("calling vp 2: %f %f\n", z[i], t_out[i]);
    }
    if(z_in) {
      // t[z.in]<-t0-gamma*z[z.in]-a/2*(1+cos(pi*(z[z.in]-h0)/(h1-h0)))
      t_out[i] = t0-gamma*z[i]-a/2*(1+cos(M_PI*(z[i]-h0)/(h1-h0)));
      //printf("calling vp 3: %f %f %f %f %f %f %f %f\n", t0, gamma, z[i], a, h0, h1, h1i, t_out[i]);
    }
  }
}

/*
# ref:
#  Lussana, C., Uboldi, F., & Salvati, M. R. (2010). A spatial consistency
#   test for surface observations from mesoscale meteorological networks.
#   Quarterly Journal of the Royal Meteorological Society, 136(649), 1075-1088.
# input
#  ixynp= vector(4). 1=box identifier on the grid;
#                   2/3=easting/northing coord (center of the box);
#                   4=number of stations within the box
#  NOTE: stations are associated to a box before running this function
#  nmin= numeric. minimum number of stations to fit a vertical profile
#  dzmin= numeric. minimum elevation range to fit a vertical profile [m]
#  Dhmin= numeric. minimum value for OI horizontal decorellation length [m]
#  Dz= numeric. OI vertical decorellation length [m]
#  Dz.bg= numeric. OI vertical decorellation length
#         for the construction of the background (RH,RR,SD) [m]
#  eps2.bg= numeric. OI ratio between obs_err_variance/backg_err_variance
#  eps2= numeric. OI ratio between obs_err_variance/backg_err_variance
#  T2=numeric. SCT threshold. (obs-pred)^2/(varObs+varPred)^2 > T2, suspect!
#  T2pos=numeric. SCT threshold. (obs-pred)>=0 AND
#                                (obs-pred)^2/(varObs+varPred)^2 > T2, suspect!
#  T2neg=numeric. SCT threshold. (obs-pred) <0 AND
#                                (obs-pred)^2/(varObs+varPred)^2 > T2, suspect!
#  NOTE: if T2pos and T2neg are specified then they are used, otherwise use T2
#  sus.code=numeric. identifier code for suspect observation
# output
#  number of rejected stations. (special cases: (i) NA if the function has not
#   been applied; (ii) -1 if just one station in the domain
#
*/
void spatial_consistency_test(double *t2, int *box, double *boxCentre, int *numStationsInBox,
                              double *x, double *y, double *z, double *t, double *vp)
                              //int *nmin, int *dzmin, int *dhmin, int *dz,
                              //int *dz_bg, double *eps2_bg, double *eps2,
                              //double *T2)
{
  int n = numStationsInBox[0];
  printf("SCT - number stations %i in box %i \n", n, box[0]);
  /*
  # distance matrices
  disth<-(outer(xtot[j],xtot[j],FUN="-")**2.+ outer(ytot[j],ytot[j],FUN="-")**2.)**0.5
  distz<-abs(outer(ztot[j],ztot[j],FUN="-"))
  */
  double** disth = malloc(sizeof(double*)*n);
  double** distz = malloc(sizeof(double*)*n);
  double *Dh = malloc(sizeof(double)*n);

  print_vector(x,n);
  print_vector(y,n);
  print_vector(z,n);
  print_vector(t,n);
  print_vector(vp,n);

  // no need to select j since already only have those for a particular box
  // outer product of the matrices
  for(int i=0; i<n; i++) {
    disth[i] = malloc(sizeof(double)*n);
    distz[i] = malloc(sizeof(double)*n);
    double *Dh_vector = malloc(sizeof(double)*(n-1)); // need to remove one since not considering the diagonal
    for(int j=0; j<n; j++) {
      //printf("i %i j %i \n", i, j);
      disth[i][j] = pow((pow((x[i]-x[j]),2)+pow((y[i]-y[j]),2)),0.5);
      distz[i][j] = abs(z[i]-z[j]);
      if(i != j) { // do not want to consider the diagonal
        if(i < j) {
          Dh_vector[j-1] = disth[i][j];
        }
        else if(i > j) {
          Dh_vector[j] = disth[i][j];
        }
      }
    }
    Dh[i] = compute_quantile(0.10, Dh_vector, n-1);
    free(Dh_vector);
  }
  print_matrix(disth,n,n);
  print_matrix(distz,n,n);
  /*
  # set to optimal Dh for the SCT
  # either Dhmin or the average 10-percentile of the set of distances between a
  # station and all the others
  Dh<-max(Dhmin,
  mean(apply(cbind(1:nrow(disth),disth),MARGIN=1,FUN=function(x){
     as.numeric(quantile(x[which(!((1:length(x))%in%c(1,(x[1]+1))))],probs=0.1))})))
  */
  double Dh_mean = mean(Dh,n);
  printf("Dh: %f\n", Dh_mean);
  // TODO: what number should this be? Dhmin
  if(Dh_mean < 10000) {
    Dh_mean = 10000;
  }
  printf("Dh_mean: %f\n", Dh_mean);

  // background error correlation matrix
  int Dz = 200; // good value (use this default)
  gsl_matrix *S, *Sinv;
  S = gsl_matrix_alloc(n,n);
  Sinv = gsl_matrix_alloc(n,n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      //  S = exp(-0.5*(disth/Dh)**2.-0.5*(distz/Dz)**2.)
      double value = exp(-.5*pow((disth[i][j]/Dh_mean),2)-.5*pow((distz[i][j]/Dz),2));
      if(i == j) { // weight the diagonal?? (0.5 default)
        // TODO: make this added value user provided
        value = value + 0.5;
      }
      gsl_matrix_set(S,i,j,value);
      //gsl_matrix_set(Sinv,i,j,value); // not yet inverted, but need a copy of S
    }
  }
  printf("created the S matrix - size1 %lu size2 %lu \n", S->size1, S->size2);
  print_gsl_matrix(S,n,n);

  gsl_vector *stationFlags;
  stationFlags = gsl_vector_alloc(n);
  // d<-topt-tb
  gsl_vector *d;
  d = gsl_vector_alloc(n);
  for(int i=0; i<n; i++) {
    // initialize with 0 (flag no stations)
    gsl_vector_set(stationFlags,i,0);
    gsl_vector_set(d,i,(t[i]-vp[i])); // difference between actual temp and temperature from vertical profile
  }

  bool first = true;
  int current_n = n;
  int throwOut = 0;
  // loop for SCT
  // note that current_n should be used inside this loop!!!
  while(1) {
    if(first) {
      // if first time then invert matrix
      //gsl_matrix_memcpy (Sinv, S);
      clock_t start = clock(), diff;
      Sinv = inverse_matrix(S);
      // if use this, then make Sinv a copy of S
      //gsl_linalg_cholesky_decomp1(Sinv);
      //gsl_linalg_cholesky_invert(Sinv); // in place
      diff = clock() - start;
      int msec = diff * 1000 / CLOCKS_PER_SEC;
      printf("Time taken to invert matrix %d seconds %d milliseconds \n", msec/1000, msec%1000);

      // UN-weight the diagonal of S
      for(int i=0; i<current_n; i++) {
        // TODO: make this added value user provided
        double value = gsl_matrix_get(S,i,i) - 0.5;
        gsl_matrix_set(S,i,i,value);
      }
      printf("S first\n");
      print_gsl_matrix(S, current_n, current_n); //(int rows, int columns, gsl_matrix *matrix)
      printf("Sinv first\n");
      print_gsl_matrix(Sinv, current_n, current_n);

      // no longer the first iteration
      first = false;
    }
    else { // not first time
      if (throwOut > 0) { // change to "else if ( >1 ) if implement the shortcut
        /*
        S<-S[-indx,-indx]
        eps2.vec<-eps2.vec[-indx]
        diag(S)<-diag(S)+eps2.vec
        SRinv<-chol2inv(chol(S))
        # from S+R go back to S
        diag(S)<-diag(S)-eps2.vec
        */
        int newSize = current_n - throwOut;
        gsl_vector *sf_temp = gsl_vector_alloc(newSize);
        gsl_vector *d_temp = gsl_vector_alloc(newSize);
        gsl_matrix *s_temp = gsl_matrix_alloc(newSize, newSize);

        printf("d before: ");
        for(int vec=0; vec<current_n; vec++) {
          printf(" %f", gsl_vector_get(d, vec));
        }
        printf("\n");

        int counter_i = 0;
        for(int i=0; i<current_n; i++) {
          int sf = gsl_vector_get(stationFlags,i);
          if(sf == 1) {
            printf("Removing column - counter_i: %i, i: %i \n", counter_i, i);
          }
          else if(sf == 0){ // add all rows and columns that we want to keep
            // update stationFlags
            gsl_vector_set(sf_temp,counter_i,0);
            // update d
            gsl_vector_set(d_temp,counter_i,gsl_vector_get(d,i));
            int counter_j = 0;
            for(int j=0; j<current_n; j++) {
              int sfj = gsl_vector_get(stationFlags,j);
              if(sfj == 0) {
                // update S
                gsl_matrix_set(s_temp, counter_i, counter_j, gsl_matrix_get(S, i, j));
                counter_j++;
              }
              else if (sfj == 1){
                //printf("Removing row - counter_j: %i, j: %i \n", counter_j, j);
              }
            }
            assert(counter_j == newSize);
            counter_i++;
          }
        }
        assert(counter_i == newSize);
        current_n = newSize;
        gsl_vector_free(stationFlags);
        stationFlags = sf_temp;
        gsl_vector_free(d);
        d = d_temp;
        gsl_matrix_free(S);
        S = s_temp;
        printf("d after: ");
        for(int vec=0; vec<current_n; vec++) {
          printf(" %f", gsl_vector_get(d, vec));
        }
        printf("\n");
        assert(stationFlags->size == current_n);
        assert(d->size == current_n);
        assert(S->size1 == current_n);
        assert(S->size2 == current_n);

        // weight the diagonal again
        for(int i=0; i<current_n; i++) {
          // TODO: make this added value user provided
          double value = gsl_matrix_get(S,i,i) + 0.5;
          gsl_matrix_set(S,i,i,value);
        }
        gsl_matrix_free(Sinv); // free the old Sinv
        Sinv = gsl_matrix_alloc(current_n,current_n);
        // invert the matrix again
        //gsl_matrix_memcpy (Sinv, S);
        clock_t start = clock(), diff;
        Sinv = inverse_matrix(S);
        // if use this, then make Sinv a copy of S
        //gsl_linalg_cholesky_decomp1(Sinv);
        //gsl_linalg_cholesky_invert(Sinv); // in place
        diff = clock() - start;
        int msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("Time taken to invert matrix %d seconds %d milliseconds \n", msec/1000, msec%1000);

        // UN-weight the diagonal of S
        for(int i=0; i<current_n; i++) {
          // TODO: make this added value user provided
          double value = gsl_matrix_get(S,i,i) - 0.5;
          gsl_matrix_set(S,i,i,value);
        }

      }
      //printf("S\n");
      //print_gsl_matrix(S, current_n, current_n); //(int rows, int columns, gsl_matrix *matrix)
      //printf("Sinv\n");
      //print_gsl_matrix(Sinv, current_n, current_n);

      // implent to remove only one station at once (faster?)
      /*
      else { //throwout == 1
        # Update inverse matrix (Uboldi et al 2008, Appendix AND erratum!)
        aux<-SRinv
        SRinv<-aux[-indx,-indx]-(tcrossprod(aux[indx,-indx],aux[-indx,indx]))*Zinv[indx]
        S<-S[-indx,-indx]
        eps2.vec<-eps2.vec[-indx]
        rm(aux)
        aux = gsl_matrix_alloc(current_n,current_n);
        gsl_matrix_transpose_memcpy(aux, Sinv); // aux is now the transpose of Sinv
        gsl_matrix_free(aux);
      }
      */
    }
    printf("Current n (end of matrix and vector updates) %i \n", current_n);
    printf("d size %lu \n", d->size);
    printf("S size1 %lu size2 %lu \n", S->size1, S->size2);
    assert(stationFlags->size == current_n);
    assert(d->size == current_n);
    assert(S->size1 == current_n);
    assert(S->size2 == current_n);
    assert(Sinv->size1 == current_n);
    assert(Sinv->size2 == current_n);
    assert(current_n > 0); // Should hopefully not throw out all the stations...

    gsl_vector *Zinv, *Sinv_d, *ares_temp, *ares;
    Zinv = gsl_vector_alloc(current_n);
    Sinv_d = gsl_vector_alloc(current_n);
    ares_temp = gsl_vector_alloc(current_n);
    ares = gsl_vector_alloc(current_n);

    // (SRinv.d<-crossprod(SRinv,d[sel]))
    // this function does not appear to work properly!!!
    //gsl_blas_dgemv(CblasNoTrans, 1, Sinv, d, 1, Sinv_d); // crossprod Sinv & d to create Sinv_d
    //gsl_blas_dgemv(CblasNoTrans, 1, S, Sinv_d, 1, ares_temp); // crossprod S and Sinv_d to create ares_temp
    for(int i=0; i<current_n; i++) {
      double acc = 0;
      for(int j=0; j<current_n; j++) {
        acc += gsl_matrix_get(Sinv,i,j)*gsl_vector_get(d,j);
      }
      gsl_vector_set(Sinv_d, i, acc);
    }
    for(int i=0; i<current_n; i++) {
      double acc = 0;
      for(int j=0; j<current_n; j++) {
        acc += gsl_matrix_get(S,i,j)*gsl_vector_get(Sinv_d,j);
      }
      gsl_vector_set(ares_temp, i, acc);
    }

    for(int i=0; i<current_n; i++) {
      gsl_vector_set(Zinv,i,(1/gsl_matrix_get(Sinv,i,i))); //Zinv<-1/diag(SRinv)
      gsl_vector_set(ares,i,(gsl_vector_get(ares_temp,i)-gsl_vector_get(d,i))); // ares<-crossprod(S,SRinv.d)-d[sel]
    }
    gsl_vector_free(ares_temp);
    //printf("Zinv: ");
    //print_gsl_vector(Zinv,current_n);
    //printf("Sinv_d: ");
    //print_gsl_vector(Sinv_d,current_n);

    // cvres<--Zinv*SRinv.d
    gsl_vector *cvres;
    cvres = gsl_vector_alloc(current_n);
    for(int i=0; i<current_n; i++) {
      gsl_vector_set(cvres,i,-1*gsl_vector_get(Zinv,i));
    }
    gsl_vector_mul(cvres,Sinv_d); // multiplies -Zinv(initial cvres) by Sinv_d (result stored in cvres)

    // sig2o<-mean(d[sel]*(-ares))
    gsl_vector *sig2o_temp, *negAres_temp;
    sig2o_temp = gsl_vector_alloc(current_n);
    negAres_temp = gsl_vector_alloc(current_n);
    for(int i=0; i<current_n; i++) {
      gsl_vector_set(negAres_temp,i,-1*gsl_vector_get(ares,i));
    }
    gsl_vector_memcpy(sig2o_temp,d); // copies d into sig2o_temp
    gsl_vector_mul(sig2o_temp,negAres_temp); // multiplies d by -ares
    double sig2o = 0;
    for(int i=0; i<current_n; i++) {
      sig2o = sig2o + gsl_vector_get(sig2o_temp,i);
    }
    //printf("d: ");
    //print_gsl_vector(d,current_n);
    //printf("neg ares: ");
    //print_gsl_vector(negAres_temp,current_n);
    sig2o = sig2o/current_n;
    printf("sig2o: %f\n", sig2o);
    gsl_vector_free(sig2o_temp);
    gsl_vector_free(negAres_temp);

    assert(sig2o > 0); // really should never have negative sig2o
    //if (sig2o<0.01) sig2o<-0.01       # safe threshold
    if(sig2o < 0.01) {
      sig2o = 0.01;
      printf("using backup sig2o: %f\n", sig2o);
    }

    // pog[sel]<-(ares*cvres)/sig2o
    gsl_vector *pog, *pog_temp;
    pog = gsl_vector_alloc(current_n);
    pog_temp = gsl_vector_alloc(current_n);
    gsl_vector_memcpy(pog_temp,ares); // copies ares into pog_temp
    gsl_vector_mul(pog_temp,cvres); // multiplies ares by cvres
    printf("pog: ");
    for(int i=0; i<current_n; i++) {
      gsl_vector_set(pog,i,(gsl_vector_get(pog_temp,i)/sig2o));
      printf(" %f", gsl_vector_get(pog,i));
    }
    printf("\n");
    gsl_vector_free(pog_temp);

    // figure out if we should flag a station
    throwOut = 0; // reset this number
    for(int i=0; i<current_n; i++) {
      // haven't already thrown it out
      int sf = gsl_vector_get(stationFlags,i);
      if(sf != 1) {
        // does it fail the test
        if(gsl_vector_get(pog,i) > t2[0]) {
          //printf("throw out this piece of data: %f\n", gsl_vector_get(pog,i));
          throwOut = throwOut + 1;
          gsl_vector_set(stationFlags,i,1);
        }
      }
    }
    printf("throw out: %i \n",throwOut);

    // FREE ALL THE MEMORY !!! (but not S or d)
    gsl_vector_free(Zinv);
    gsl_vector_free(Sinv_d);
    gsl_vector_free(ares);
    gsl_vector_free(cvres);
    gsl_vector_free(pog);

    if(throwOut == 0) { // no stations to remove, done SCT iterations
      gsl_vector_free(d);
      gsl_matrix_free(S);
      gsl_matrix_free(Sinv);
      break;
    }
  }
  // end of while SCT loop
}

//----------------------------------------------------------------------------//
// HELPER FUNCTIONS
double compute_quantile(double quantile, double *array, int sizeArray)
{
  // compute quantiles
  gsl_rstat_quantile_workspace *work_q = gsl_rstat_quantile_alloc(quantile);
  double exact_q;
  double val_q;
  size_t i;
  /* add data to quantile accumulators; also store data for exact comparisons */
  for (i = 0; i < sizeArray; ++i) {
      gsl_rstat_quantile_add(array[i], work_q);
  }
  /* exact values*/
  gsl_sort(array, 1, sizeArray);
  exact_q = gsl_stats_quantile_from_sorted_data(array, 1, sizeArray, quantile);
  /* estimated values */
  val_q = gsl_rstat_quantile_get(work_q);
  //printf ("%.5f quartile: exact = %.5f, estimated = %.5f, error = %.6e\n",
  //        quantile, exact_q, val_q, (val_q - exact_q) / exact_q);
  gsl_rstat_quantile_free(work_q);
  return exact_q;
}

double mean(const double *array, int sizeArray)
{
  double sum = 0;
  for(int i=0; i<sizeArray; i++) {
    sum = sum + array[i];
  }
  double mean = sum/sizeArray;
  //printf("Mean: %f\n", mean);
  return mean;
}

void print_vector(double *vector, int size) {
  for (int s=0; s<size; s++)
  {
    printf("%.2f ", vector[s]);
  }
  printf("\n");
}

void print_gsl_vector(gsl_vector *vector, int size) {
  for (int s=0; s<size; s++)
  {
    printf("%.2f ", gsl_vector_get(vector,s));
  }
  printf("\n");
}

void print_matrix(double **matrix, int rows, int columns) {
  for (int r=0; r<rows; r++)
  {
      for(int c=0; c<columns; c++)
          {
           printf("%.2f ", matrix[r][c]);
          }
      printf("\n");
   }
}

void print_gsl_matrix(gsl_matrix *matrix, int rows, int columns) {
  for (int r=0; r<rows; r++)
  {
      for(int c=0; c<columns; c++)
          {
           printf("%.2f ", gsl_matrix_get(matrix,r,c));
          }
      printf("\n");
   }
}

void print_sub_gsl_matrix(gsl_matrix *matrix, int start, int stop) {
  for (int r=start; r<=stop; r++)
  {
      for(int c=start; c<=stop; c++)
          {
           printf("%.2f ", gsl_matrix_get(matrix,r,c));
          }
      printf("\n");
   }
}

gsl_matrix* inverse_matrix(const gsl_matrix *matrix) {
  int s;
  gsl_matrix *return_matrix = gsl_matrix_alloc(matrix->size1,matrix->size2);
  gsl_matrix *temp_matrix = gsl_matrix_alloc(matrix->size1,matrix->size2);
  gsl_matrix_memcpy(temp_matrix, matrix);
  gsl_permutation * p = gsl_permutation_alloc (matrix->size1);

  gsl_linalg_LU_decomp (temp_matrix, p, &s);
  gsl_linalg_LU_invert (temp_matrix, p, return_matrix);
  gsl_matrix_free(temp_matrix);
  gsl_permutation_free (p);
  return return_matrix;
}
