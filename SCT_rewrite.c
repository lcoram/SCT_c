#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_blas.h>

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
void spatial_consistency_test(int *t2, int *box, double *boxCentre, int *numStationsInBox,
                              double *x, double *y, double *z, double *t, double *vp); // this is the interface to R (then need pointers)?

int vertical_profile_optimizer(gsl_vector *input, double **data);
double vertical_profile_optimizer_function(const gsl_vector *v, void *data); // GSL format

void vertical_profile(int nz, double *z,
    double t0, double gamma, double a, double h0, double h1i, double *t_out);


// helpers
gsl_matrix* invert_matrix(gsl_matrix *matrix, int size);
double compute_quantile(double quantile, double *array, int sizeArray);
double mean(double *array, int sizeArray);

int main()
{
  // Testing functions from main, but eventually have to call it from R code
  FILE *fp;
  char *filename = "/home/louiseo/Documents/SCT/myrepo/input_for_example.txt";
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
    bool oslo = false;
    while(ptr != NULL) { // break up line
      // if itot == 266 then in oslo box
      if(j==0 && (strcmp(ptr, "266")==0 || strcmp(ptr, " 266")==0 )) {
        oslo = true;
      }
      if(oslo) {
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
        //printf("%s\n", ptr);
      }
      ptr = strtok(NULL, delim);
      j++;
    }
  }
  printf("Number of stations in Oslo box: %i \n", n);
  printf("Mean x: %f y: %f z: %f \n", mean(x,n), mean(y,n), mean(z,n));
  printf("Mean t: %f\n", mean(t,n));

  fclose(fp);

  // initial variables for VP (first guess, adjusted during optimization)
  double gamma = 0.0065; // default?
  double a = 5;

  double meanT = mean(t,n);

  double exact_p10 = compute_quantile(0.10, z, n);
  double exact_p90 = compute_quantile(0.90, z, n);

  // allocate for output
  double *t_out = malloc(sizeof(double) * n);

  /* data (params) that needs to be passed into vp */
  double nd = (double) n; // cast + resize to double
  // data (double *n, double *z, double *t, double *t_out)
  double * data[4] = {&nd, z, t, t_out};

  // vector (double t0, double gamma, double a, double h0, double h1i)
  /* Starting point for optimization */
  gsl_vector *input = gsl_vector_alloc(5);
  gsl_vector_set(input,0,meanT);
  gsl_vector_set(input,1,gamma);
  gsl_vector_set(input,2,a);
  gsl_vector_set(input,3,exact_p10);
  gsl_vector_set(input,4,exact_p90);
  printf ("Input vector set = t0: %.4f gamma: %.4f a: %.4f h0: %.4f h1i: %.4f\n",
          meanT, gamma, a, exact_p10, exact_p90);

  int status = vertical_profile_optimizer(input, data);
  printf("status: %d\n", status);
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

  // variables for SCT
  int t2 = 25; // input by user into SCT function
  
    // 266;-247429.070909252;-365421.526660934;3466;
  int box = 266;
  double boxCentre[2];
  boxCentre[0] = -247429.070909252;
  boxCentre[1] = -365421.526660934;
  int numStationsInBox = n;

  // void spatial_consistency_test(int *t2, int *box, double *boxCentre, int *numStationsInBox,
                                //double *x, double *y, double *z, double *vp)
  clock_t start = clock(), diff;
  spatial_consistency_test(&t2, &box, boxCentre, &numStationsInBox, z, x, y, t, t_out);
  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("SCT end\n");
  printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);

  free(t);
  return 0;
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
  // my input data
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

  // give everything to vp to compute t_out
  vertical_profile(n, z, t0, gamma, h0, h0, h1i, t_out);
  // RMS
  double total = 0;
  for(int i=0; i<n; i++) {
    total += pow((t[i]-t_out[i]),2);
  }
  double value = total / n;
  // Do not need the log?
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
  double h1 = h0 + fabs(h1i); // h1<-h0+abs(h1i)
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
      //printf("calling vp 1: %f \n", t_out[i]);
    }
    if(z_ge_h1) {
      // t[z.ge.h1]<-t0-gamma*z[z.ge.h1]
      t_out[i] = t0-gamma*z[i];
      //printf("calling vp 2: %f \n", t_out[i]);
    }
    if(z_in) {
      // t[z.in]<-t0-gamma*z[z.in]-a/2*(1+cos(pi*(z[z.in]-h0)/(h1-h0)))
      t_out[i] = t0-gamma*z[i]-a/2*(1+cos(M_PI*(z[i]-h0)/(h1-h0)));
      //printf("calling vp 3: %f \n", t_out[i]);
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
void spatial_consistency_test(int *t2, int *box, double *boxCentre, int *numStationsInBox,
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
  // no need to select j since already only have those for a particular box
  // outer product of the matrices
  for(int i=0; i<n; i++) {
    disth[i] = malloc(sizeof(double)*n);
    distz[i] = malloc(sizeof(double)*n);
    double *Dh_vector = malloc(sizeof(double)*n);
    for(int j=0; j<n; j++) {
      //printf("i %i j %i \n", i, j);
      disth[i][j] = pow((pow((x[i]-x[j]),2)+pow((y[i]-y[j]),2)),0.5);
      distz[i][j] = abs(z[i]-z[j]);
      Dh_vector[j] = disth[i][j];
    }
    Dh[i] = compute_quantile(0.10, Dh_vector, n);
  }
  /*
  # set to optimal Dh for the SCT
  # either Dhmin or the average 10-percentile of the set of distances between a
  # station and all the others
  Dh<-max(Dhmin,
  mean(apply(cbind(1:nrow(disth),disth),MARGIN=1,FUN=function(x){
     as.numeric(quantile(x[which(!((1:length(x))%in%c(1,(x[1]+1))))],probs=0.1))})))
  */
  int Dz = 200; // good value (use this default)
  double Dh_mean = mean(Dh,n);
  printf("Dh: %f\n", Dh_mean);

  // background error correlation matrix
  gsl_matrix *S, *Sinv;
  S = gsl_matrix_alloc(n,n);
  Sinv = gsl_matrix_alloc(n,n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      //  S = exp(-0.5*(disth/Dh)**2.-0.5*(distz/Dz)**2.)
      double value = exp(-.5*pow((disth[i][j]/Dh_mean),2)-.5*pow((distz[i][j]/Dz),2));
      if(i == j) { // weight the diagonal?? (0.5 default)
        //value = value *0.5;
      }
      gsl_matrix_set(S,i,j,value);
      gsl_matrix_set(Sinv,i,j,value); // not yet inverted, but need a copy of S
    }
  }
  printf("created the S matrix\n");


  bool first = true;
  int stationFlags[n];
  //  d<-topt-tb
  gsl_vector *d;
  d = gsl_vector_alloc(n);
  for(int i=0; i<n; i++) {
    // initialize with 0 (flag no stations)
    stationFlags[i] = 0;
    gsl_vector_set(d,i,(t[i]-vp[i])); // difference between actual temp and temperature from vertical profile
  }

  // now loop for SCT
  while(1) {

    if(first) {
      // if first time then invert matrix
      //gsl_matrix* invert_matrix(gsl_matrix *matrix, int size)
      clock_t start = clock(), diff;
      //gsl_matrix *S_gsl_inv = invert_matrix(S_gsl,n);
      gsl_linalg_cholesky_invert(Sinv); // in place
      diff = clock() - start;
      int msec = diff * 1000 / CLOCKS_PER_SEC;
      printf("Time taken to invert matrix %d seconds %d milliseconds \n", msec/1000, msec%1000);
      first = false;
    }
    else { // not first time
      // Update inverse matrix (Uboldi et al 2008, Appendix AND erratum!)
      /*
      aux<-SRinv
      SRinv<-aux[-indx,-indx]-(tcrossprod(aux[indx,-indx],aux[-indx,indx]))*Zinv[indx]
      S<-S[-indx,-indx]
      eps2.vec<-eps2.vec[-indx]
      rm(aux)
      */
      // take out the stations that we have decided to remove
    }

    gsl_vector *Zinv, *Sinv_d, *ares_temp, *ares;
    Zinv = gsl_vector_alloc(n);
    Sinv_d = gsl_vector_alloc(n);
    ares_temp = gsl_vector_alloc(n);
    ares = gsl_vector_alloc(n);
    // (SRinv.d<-crossprod(SRinv,d[sel]))
    gsl_blas_dgemv(CblasNoTrans, 1, Sinv, d, 1, Sinv_d); // crossprod Sinv & d to create Sinv_d
    gsl_blas_dgemv(CblasNoTrans, 1, S, Sinv_d, 1, ares_temp); // crossprod S and Sinv_d to create ares_temp
    for(int i=0; i<n; i++) {
      gsl_vector_set(Zinv,i,(1/gsl_matrix_get(Sinv,i,i))); //Zinv<-1/diag(SRinv)
      gsl_vector_set(ares,i,(gsl_vector_get(ares_temp,i)-gsl_vector_get(d,i))); // ares<-crossprod(S,SRinv.d)-d[sel]
      if(abs(gsl_vector_get(ares,i)) > 1) {
        printf("abs(ares)>1: %i, %f\n", i, gsl_vector_get(ares,i));
      }
    }
    gsl_vector_free(ares_temp);

    // cvres<--Zinv*SRinv.d
    gsl_vector *cvres;
    cvres = gsl_vector_alloc(n);
    for(int i=0; i<n; i++) {
      gsl_vector_set(cvres,i,-1*gsl_vector_get(Zinv,i));
    }
    gsl_vector_mul(cvres,Sinv_d); // multiplies -Zinv(initial cvres) by Sinv_d (result stored in cvres)

    // sig2o<-mean(d[sel]*(-ares))
    gsl_vector *sig2o_temp, *negAres_temp;
    sig2o_temp = gsl_vector_alloc(n);
    negAres_temp = gsl_vector_alloc(n);
    for(int i=0; i<n; i++) {
      gsl_vector_set(negAres_temp,i,-1*gsl_vector_get(ares,i));
    }
    gsl_vector_memcpy(sig2o_temp,d); // copies d into sig2o_temp
    gsl_vector_mul(sig2o_temp,negAres_temp); // multiplies d by -ares
    double sig2o = 0;
    for(int i=0; i<n; i++) {
      sig2o = sig2o + gsl_vector_get(sig2o_temp,i);
    }
    sig2o = sig2o/n;
    printf("sig2o: %f\n", sig2o);
    gsl_vector_free(sig2o_temp);
    gsl_vector_free(negAres_temp);

    // pog[sel]<-(ares*cvres)/sig2o
    gsl_vector *pog, *pog_temp;
    pog = gsl_vector_alloc(n);
    pog_temp = gsl_vector_alloc(n);
    gsl_vector_memcpy(pog_temp,ares); // copies ares into pog_temp
    gsl_vector_mul(pog_temp,cvres); // multiplies ares by cvres
    //printf("pog: ");
    for(int i=0; i<n; i++) {
      gsl_vector_set(pog,i,(gsl_vector_get(pog_temp,i)/sig2o));
      //printf("%f", gsl_vector_get(pog,i));
    }
    //printf("\n");
    gsl_vector_free(pog_temp);

    // figure out if we should flag a station
    int throwOut = 0;
    for(int i=0; i<n; i++) {
      // does it fail the test
      if(gsl_vector_get(pog,i) > t2[0]) {
        //printf("throw out this piece of data: %f\n", gsl_vector_get(pog,i));
        throwOut = throwOut + 1;
      }
    }
    printf("throw out: %i \n",throwOut);

  // FREE ALL THE MEMORY !!!
  gsl_matrix_free(S);
  gsl_matrix_free(Sinv);
  gsl_vector_free(d);
  gsl_vector_free(Zinv);
  gsl_vector_free(Sinv_d);
  gsl_vector_free(ares);
  gsl_vector_free(cvres);
  gsl_vector_free(pog);

  break; // have not currently implemented more than one loop...
  }
  // end of while SCT loop

}

//----------------------------------------------------------------------------//
// HELPER FUNCTIONS
gsl_matrix* invert_matrix(gsl_matrix *matrix, int size)
{
  gsl_permutation *p = gsl_permutation_alloc(size);
  gsl_matrix *matrix_inv = gsl_matrix_alloc(size,size);
  gsl_linalg_LU_invert(matrix,p,matrix_inv);
  gsl_permutation_free(p);
  return matrix_inv;
}

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

double mean(double *array, int sizeArray)
{
  double sum = 0;
  for(int i=0; i<sizeArray; i++) {
    sum = sum + array[i];
  }
  double mean = sum/sizeArray;
  //printf("Mean: %f\n", mean);
  return mean;
}
