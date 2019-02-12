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
// compile: gcc -L/user/local/lib SCT_wrapper.c -lm -lgsl -lgslcblas


/*
 * Runs the SCT by splitting into boxes first
 *
 * Arguments:
 *    n: Number of stations
 *    x: Array of x-position (in meters in some projection)
 *    y: Array of y-position (in meters in some projection)
 *    z: Array of station altitudes (in meters)
 *    t: Array of temperatures (in any units)
 *    nmax: Split a box into two if it contains more than this number of stations
 *    nmin: Don't allow boxes to have fewer than this  a box into two if it contains more than this number of stations
 *
 * Outputs:
 *    flags: Output array of flags, one for each station. 0 means passed the SCT, 1 means fails.
 * */
void sct_wrapper(int *n, double *x, double *y, double *z, double *t, int *nmax, int *nmin, int *nminprof, double *gam, double *as, double *t2pos, double *t2neg, int *flags, double *corep, double *pog);

/*
 * Structure to contain station information for a box
 * */
struct box {
  int  n;
  double *x;
  double *y;
  double *z;
  double *t;
  int *i;
};

void spatial_consistency_test(struct box *currentBox, int *nminprof, double *gam, double *as, double *t2pos, double *t2neg, int *flags, double *corep, double *pog);

int vertical_profile_optimizer(gsl_vector *input, struct box *currentBox, int nminprof, double *vp);
// optimizer functions
double basic_vertical_profile_optimizer_function(const gsl_vector *v, void *data);
double vertical_profile_optimizer_function(const gsl_vector *v, void *data); // GSL format
// vp profile calculations
void basic_vertical_profile(int nz, double *z, double t0, double gamma, double a, double h0, double h1i, double *t_out);
void vertical_profile(int nz, double *z, double t0, double gamma, double a, double h0, double h1i, double *t_out);

// box division controller (other functions used by this controller not forward declared)
struct box * control_box_division(int maxNumStationsInBox, int minNumStationsInBox, struct box inputBox);

// Helper functions
double compute_quantile(double quantile, double *array, int sizeArray);
double mean(const double *array, int sizeArray);
double max(const double *array, int sizeArray);
double min(const double *array, int sizeArray);
void print_vector(double *vector, int size);
void print_gsl_vector(gsl_vector *vector, int size);
void print_matrix(double **matrix, int rows, int columns);
void print_gsl_matrix(gsl_matrix *matrix, int rows, int columns);
void print_sub_gsl_matrix(gsl_matrix *matrix, int start, int stop);

gsl_matrix* inverse_matrix(const gsl_matrix *matrix);


//----------------------------------------------------------------------------//
void sct_wrapper(int *n, double *x, double *y, double *z, double *t, int *nmax, int *nmin, int *nminprof, double *gam, double *as, double *t2pos, double *t2neg, int *flags, double *corep, double *pog) {

  // put the input box in the struct
  struct box inputBox;
  inputBox.n = n[0];
  inputBox.x = x;
  inputBox.y = y;
  inputBox.z = z;
  inputBox.t = t;
  inputBox.i = malloc(sizeof(int) * n[0]);
  for(int i = 0; i < n[0]; i++)
     inputBox.i[i] = i;

  // split the box if needed
  struct box * nAndBoxes = control_box_division(nmax[0], nmin[0], inputBox);

  int nB = nAndBoxes[0].n;
  printf("SCT wrapper - number of boxes after division: %i \n", nB);

  // loop over the boxes to call SCT
  for(int i=1; i<nB+1; i++) {
    int box_n = nAndBoxes[i].n;
    double * box_x = nAndBoxes[i].x;
    double * box_y = nAndBoxes[i].y;
    double * box_z = nAndBoxes[i].z;
    double * box_t = nAndBoxes[i].t;
    int * box_i = nAndBoxes[i].i;
    printf("box: %i \n", (i-1));

    // Run the SCT on the current box
    printf("num stations before SCT: %d \n", box_n);
    clock_t start = clock(), diff;
    int* local_flags = malloc(sizeof(int) * box_n);
    double* local_t2pos = malloc(sizeof(double) * box_n);
    double* local_t2neg = malloc(sizeof(double) * box_n);
    double* local_corep = malloc(sizeof(double) * box_n);
    double* local_pog = malloc(sizeof(double) * box_n);
    for(int r = 0; r < box_n; r++) {
       if(!(box_i[r] < n[0]))
          printf("%d %d\n", box_i[r], n[0]);
       assert(box_i[r] < n[0]);
       local_t2pos[r] = t2pos[box_i[r]];
       local_t2neg[r] = t2neg[box_i[r]];
       local_corep[r] = corep[box_i[r]];
       local_pog[r] = pog[box_i[r]];
    }
    spatial_consistency_test(&nAndBoxes[i], nminprof, gam, as, local_t2pos, local_t2neg, local_flags, local_corep, local_pog);

    for(int r = 0; r < box_n; r++) {
       assert(box_i[r] < n[0]);
       flags[box_i[r]] = local_flags[r];
       corep[box_i[r]] = local_corep[r];
       pog[box_i[r]] = local_pog[r];
    }
    free(local_t2pos);
    free(local_t2neg);
    free(local_flags);
    free(local_corep);
    free(local_pog);

    // Merge flags back into global array
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("SCT end - Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
    printf("num stations after SCT: %d \n", box_n);

  } // end of looping over boxes
  free(inputBox.i);
  return;
}


//----------------------------------------------------------------------------//
int vertical_profile_optimizer(gsl_vector *input, struct box *currentBox, int nminprof, double *vp)
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

  // data (params) that needs to be passed into vp
  double nd = (double) currentBox[0].n; // cast + resize to double
  double *z = currentBox[0].z;
  double *t = currentBox[0].t;
  double * data[4] = {&nd, z, t, vp};

  /* Initialize method and iterate */
  vp_optim.n = 5;
  if(currentBox[0].n < nminprof) {
    vp_optim.f = basic_vertical_profile_optimizer_function;
  }
  else {
    vp_optim.f = vertical_profile_optimizer_function;
  }
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

  // now actually calculate the vertical profile using these optimized variables
  if(currentBox[0].n < nminprof) { // basic vp
    basic_vertical_profile(currentBox[0].n, currentBox[0].z, gsl_vector_get(input, 0), gsl_vector_get(input, 1),
            gsl_vector_get(input, 2), gsl_vector_get(input, 3), gsl_vector_get(input, 4), vp);
  }
  else { // more complicated vp
    vertical_profile(currentBox[0].n, currentBox[0].z, gsl_vector_get(input, 0), gsl_vector_get(input, 1),
            gsl_vector_get(input, 2), gsl_vector_get(input, 3), gsl_vector_get(input, 4), vp);
  }

  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return status;
}

//----------------------------------------------------------------------------//
/*
#+ cost function used for optimization of tvertprof parameter
tvertprofbasic2opt<-function(par) {
  te<-tvertprof_basic(z=zopt,t0=par[1],gamma=argv$gamma.standard)
  return(log((mean((te-topt)**2))**0.5))
}
*/
double basic_vertical_profile_optimizer_function(const gsl_vector *v, void *data)
{
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
  basic_vertical_profile(n, z, t0, gamma, a, h0, h1i, t_out);
  // RMS
  double total = 0;
  for(int i=0; i<n; i++) {
    total += pow((t_out[i]-t[i]),2);
  }
  double value = log(pow((total / n),0.5));

  return value;
}

/*
#+ vertical profile of temperature (linear)
tvertprof_basic<-function(z,t0,gamma) {
# input
#  z= array. elevations [m amsl]
#  t0= numeric. temperature at z=0 [K or degC]
#  gamma=numeric. temperature lapse rate [K/m]
# Output
#  t= array. temperature [K or degC]
#------------------------------------------------------------------------------
  return(t0+gamma*z)
}
*/
void basic_vertical_profile(int nz, double *z,
    double t0, double gamma, double a, double h0, double h1i, double *t_out)
{
  for(int i=0; i<nz; i++) {
    t_out[i] = t0 + gamma*z[i];
  }
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
void spatial_consistency_test(struct box *currentBox, int *nminprof, double *gam, double *as, double *t2pos, double *t2neg, int *flags, double* corep_out, double* pog_out)
{
  // break out the box for simplicity
  int n = currentBox[0].n;
  double *x = currentBox[0].x;
  double *y = currentBox[0].y;
  double *z = currentBox[0].z;
  double *t = currentBox[0].t;
  printf("SCT - number stations %i \n", n);

  // fill with 0 to keep track of stations that are flagged
  for(int i=0; i<n; i++) {
    flags[i] = 0;
  }

  /*
  Stuff for VP
  */
  double gamma = gam[0];
  double a = as[0];
  double meanT = mean(t,n);
  // make copies of z for the quantile computation
  double * z_temp1;
  double * z_temp2;
  z_temp1 = malloc(sizeof(double) * n);
  z_temp2 = malloc(sizeof(double) * n);
  // allocate for output
  double *vp = malloc(sizeof(double) * n);
  for(int i=0; i<n; i++) {
    z_temp1[i] = z[i];
    z_temp2[i] = z[i];
    vp[i] = -999;
  }
  double exact_p10 = compute_quantile(0.10, z_temp1, n);
  double exact_p90 = compute_quantile(0.90, z_temp2, n);
  free(z_temp1);
  free(z_temp2);
  // vector (double t0, double gamma, double a, double h0, double h1i)
  // Starting point for optimization
  gsl_vector *vp_input = gsl_vector_alloc(5);
  gsl_vector_set(vp_input,0,meanT);
  gsl_vector_set(vp_input,1,gamma);
  gsl_vector_set(vp_input,2,a);
  gsl_vector_set(vp_input,3,exact_p10);
  gsl_vector_set(vp_input,4,exact_p90);
  printf ("VP input vector set = t0: %.4f gamma: %.4f a: %.4f h0: %.4f h1i: %.4f\n",
          meanT, gamma, a, exact_p10, exact_p90);
  // calculate VP
  // int vertical_profile_optimizer(gsl_vector *input, struct box *currentBox, int nminprof, double *vp);
  int status = vertical_profile_optimizer(vp_input, currentBox, nminprof[0], vp);
  printf("status optimizer: %d\n", status);
  printf ("t0: %.4f gamma: %.4f a: %.4f h0: %.4f h1i: %.4f\n",
          gsl_vector_get(vp_input, 0),
          gsl_vector_get(vp_input, 1),
          gsl_vector_get(vp_input, 2),
          gsl_vector_get(vp_input, 3),
          gsl_vector_get(vp_input, 4));
  // now have temperature profile (vp)
  for(int i=0; i<n; i++) {
    assert(vp[i] !=-999);
  }
  // done calculating VP for the first time
  int sizeWhenProfileCalculated = currentBox[0].n;
  // END of calculating vp

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
  double **S_global = malloc(sizeof(double*)*n);
  for(int i=0; i<n; i++) {
    S_global[i] = malloc(sizeof(double)*n);
  }
  Sinv = gsl_matrix_alloc(n,n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      //  S = exp(-0.5*(disth/Dh)**2.-0.5*(distz/Dz)**2.)
      double value = exp(-.5*pow((disth[i][j]/Dh_mean),2)-.5*pow((distz[i][j]/Dz),2));
      S_global[i][j] = value;
      if(i == j) { // weight the diagonal?? (0.5 default)
        // TODO: make this added value user provided
        value = value + 0.5;
      }
      gsl_matrix_set(S,i,j,value);
      //gsl_matrix_set(Sinv,i,j,value); // not yet inverted, but need a copy of S
    }
  }
  printf("created the S matrix - size1 %lu size2 %lu \n", S->size1, S->size2);
  //print_gsl_matrix(S,n,n);

  // d<-topt-tb
  gsl_vector *d;
  d = gsl_vector_alloc(n);
  double *d_global = malloc(sizeof(double)*n);
  for(int i=0; i<n; i++) {
    gsl_vector_set(d,i,(t[i]-vp[i])); // difference between actual temp and temperature from vertical profile
    d_global[i] = (t[i]-vp[i]);
  }

  /* ---------------------------------------------------
  Beginning of real SCT looping
  ------------------------------------------------------*/
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
      //printf("S first\n");
      //print_gsl_matrix(S, current_n, current_n); //(int rows, int columns, gsl_matrix *matrix)
      //printf("Sinv first\n");
      //print_gsl_matrix(Sinv, current_n, current_n);

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
        /* --------------------------------------------------
        have 3 sizes:
        n = original
        current_n = current size of vectors / matrices
        newSize = size that will be current n after throwing things out
        ---------------------------------------------------- */
        int newSize = current_n - throwOut;
        gsl_vector *sf_temp = gsl_vector_alloc(newSize);
        gsl_vector *d_temp = gsl_vector_alloc(newSize);
        gsl_matrix *s_temp = gsl_matrix_alloc(newSize, newSize);

        //printf("d before: ");
        //for(int vec=0; vec<current_n; vec++) {
          //printf(" %f", gsl_vector_get(d, vec));
        //}
        //printf("\n");
        // loop over the original length
        //for(int original=0; original<n; original++) {}

        int li = 0;
        for(int i=0; i<n; i++) { // this needs to loop over current_n
          if(flags[i] == 2) { // stations flagged for removal
            printf("Removing column - li: %i, i: %i \n", li, i);
            flags[i] = 1; // newly removed station now set to 1
          }
          else if(flags[i] == 0){ // all rows and columns that we want to keep
            // update d
            gsl_vector_set(d_temp,li, d_global[i]);
            int lj = 0;
            for(int j=0; j<n; j++) { // this needs to loop over current_n
              if(flags[j] == 0) {
                // update S
                gsl_matrix_set(s_temp, li, lj, S_global[i][j]);
                lj++;
              }
              else if (flags[j] == 1){
                //printf("Removing row - lj: %i, j: %i \n", lj, j);
              }
            }
            assert(lj == newSize);
            li++;
          }
        }
        assert(li == newSize);
        current_n = newSize;
        gsl_vector_free(d);
        d = d_temp;
        gsl_matrix_free(S);
        S = s_temp;
        //printf("d after: ");
        //for(int vec=0; vec<current_n; vec++) {
          //printf(" %f", gsl_vector_get(d, vec));
        //}
        //printf("\n");
        assert(d->size == current_n);
        assert(S->size1 == current_n);
        assert(S->size2 == current_n);

        // we now have an empty box!!!
        if(current_n == 0) {
          // FREE ALL THE MEMORY !!!
          gsl_vector_free(d);
          gsl_matrix_free(S);
          gsl_matrix_free(Sinv);
          break; // exit sct loop
        }
        // chech size has not changed too much, if it has then recompute VP
        float percentageSizeChange = (float)(sizeWhenProfileCalculated - current_n) / sizeWhenProfileCalculated;
        if(percentageSizeChange > 0.1) {
          printf("size of box has changed significantly (recalculate vp): %f \n", percentageSizeChange);
          // gsl_vector *vp_input, double *vp
          int status = vertical_profile_optimizer(vp_input, currentBox, nminprof[0], vp);
          printf("status optimizer: %d\n", status);
          sizeWhenProfileCalculated = current_n;
        }

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
    } // end else
    printf("Current n (end of matrix and vector updates) %i \n", current_n);
    printf("d size %lu \n", d->size);
    printf("S size1 %lu size2 %lu \n", S->size1, S->size2);
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
    int li = 0;
    for(int i=0; i<n; i++) {
      double acc = 0;
      int lj = 0;
      if(flags[i] == 0) {
        for(int j=0; j<n; j++) {
          if(flags[j] == 0) {
            acc += gsl_matrix_get(Sinv,li,lj)*gsl_vector_get(d,lj);
            lj++;
          }
        }
        gsl_vector_set(Sinv_d, li, acc);
        li++;
      }
    }
    li = 0;
    for(int i=0; i<n; i++) {
      double acc = 0;
      int lj = 0;
      if(flags[i] == 0) {
        for(int j=0; j<n; j++) {
          if(flags[j] == 0) {
            acc += gsl_matrix_get(S,li,lj)*gsl_vector_get(Sinv_d,lj);
            lj++;
          }
        }
        gsl_vector_set(ares_temp, li, acc);
        li++;
      }
    }
    li = 0;
    for(int i=0; i<n; i++) {
      if(flags[i] == 0) {
        gsl_vector_set(Zinv,li,(1/gsl_matrix_get(Sinv,li,li))); //Zinv<-1/diag(SRinv)
        gsl_vector_set(ares,li,(gsl_vector_get(ares_temp,li)-gsl_vector_get(d,li))); // ares<-crossprod(S,SRinv.d)-d[sel]
        li++;
      }
    }
    gsl_vector_free(ares_temp);
    //printf("Zinv: ");
    //print_gsl_vector(Zinv,current_n);
    //printf("Sinv_d: ");
    //print_gsl_vector(Sinv_d,current_n);

    // cvres<--Zinv*SRinv.d
    gsl_vector *cvres;
    cvres = gsl_vector_alloc(current_n);
    li = 0;
    for(int i=0; i<current_n; i++) {
      if(flags[i] == 0) {
        gsl_vector_set(cvres,li,-1*gsl_vector_get(Zinv,li));
        li++;
      }
    }
    gsl_vector_mul(cvres,Sinv_d); // multiplies -Zinv(initial cvres) by Sinv_d (result stored in cvres)

    // sig2o<-mean(d[sel]*(-ares))
    gsl_vector *sig2o_temp, *negAres_temp;
    sig2o_temp = gsl_vector_alloc(current_n);
    negAres_temp = gsl_vector_alloc(current_n);
    li = 0;
    for(int i=0; i<n; i++) {
      if(flags[i] == 0) {
        gsl_vector_set(negAres_temp,li,-1*gsl_vector_get(ares,li));
        li++;
      }
    }
    gsl_vector_memcpy(sig2o_temp,d); // copies d into sig2o_temp
    gsl_vector_mul(sig2o_temp,negAres_temp); // multiplies d by -ares
    double sig2o = 0;
    li = 0;
    for(int i=0; i<n; i++) {
      if(flags[i] == 0) {
        sig2o = sig2o + gsl_vector_get(sig2o_temp,li);
        li++;
      }
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
    //printf("pog: ");

    for(int li=0; li<current_n; li++) {
      gsl_vector_set(pog,li,(gsl_vector_get(pog_temp,li)/sig2o));
      //printf(" %f", gsl_vector_get(pog,i));
    }
    //printf("\n");
    gsl_vector_free(pog_temp);

    // figure out if we should flag a station
    throwOut = 0; // reset this number
    li = 0;
    for(int i=0; i<n; i++) {
      if(flags[i] == 0) {
         assert(sig2o > 0);
         corep_out[i] = gsl_vector_get(d, li)*gsl_vector_get(ares, li) * -1 /sig2o;
         pog_out[i] = gsl_vector_get(pog, li);
        // does it fail the test
        if((gsl_vector_get(cvres,li) < 0 && gsl_vector_get(pog,li) > t2pos[i]) ||
          (gsl_vector_get(cvres,li) >= 0 && gsl_vector_get(pog,li) > t2neg[i])) {
          printf("throw out this piece of data: %f cvres=%f pog=%f corep=%f\n", t[i], gsl_vector_get(cvres, li), gsl_vector_get(pog,li), corep_out[i]);
          throwOut = throwOut + 1;
          flags[i] = 2; // temporarily set to 2 so we know its a newly flagged station
        }
        li++;
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
      free(vp);
      break;
    }
  } // end of while SCT loop

}

struct box merge_boxes(struct box box1, struct box box2) {

  struct box mergedBox;
  mergedBox.n = (box1.n+box2.n); // set initial n
  mergedBox.x = malloc(sizeof(double) * mergedBox.n);
  mergedBox.y = malloc(sizeof(double) * mergedBox.n);
  mergedBox.z = malloc(sizeof(double) * mergedBox.n);
  mergedBox.t = malloc(sizeof(double) * mergedBox.n);
  mergedBox.i = malloc(sizeof(double) * mergedBox.n);

  for(int i=0; i<box1.n; i++) {
    mergedBox.x[i] = box1.x[i];
    mergedBox.y[i] = box1.y[i];
    mergedBox.z[i] = box1.z[i];
    mergedBox.t[i] = box1.t[i];
    mergedBox.i[i] = box1.i[i];
  }
  for(int j=box1.n; j<mergedBox.n; j++) {
    mergedBox.x[j] = box2.x[j-box1.n];
    mergedBox.y[j] = box2.y[j-box1.n];
    mergedBox.z[j] = box2.z[j-box1.n];
    mergedBox.t[j] = box2.t[j-box1.n];
    mergedBox.i[j] = box2.i[j-box1.n];
  }

  //printf("Merged box size: %i \n", mergedBox.n);
  return mergedBox;
}

// recursive function that will keep splitting boxes until they are the right size
// returns a box at the "leaf"
void split_box(int maxNumStationsInBox, int minNumStationsInBox, struct box inputBox, int * finalNumBoxes, struct box ** finalBoxes) {

  struct box * boxes;
  // allocate memory, currently for 4 boxes
  boxes = malloc(sizeof(struct box) * 4);
  for(int i=0; i<4; i++) {
    //boxes[i].n = malloc(sizeof(int));
    boxes[i].n = 0; // set initial n
    boxes[i].x = malloc(sizeof(double) * inputBox.n);
    boxes[i].y = malloc(sizeof(double) * inputBox.n);
    boxes[i].z = malloc(sizeof(double) * inputBox.n);
    boxes[i].t = malloc(sizeof(double) * inputBox.n);
    boxes[i].i = malloc(sizeof(double) * inputBox.n);
  }
  double maxX = max(inputBox.x,inputBox.n);
  double maxY = max(inputBox.y,inputBox.n);
  double minX = min(inputBox.x,inputBox.n);
  double minY = min(inputBox.y,inputBox.n);
  // halfway between min and max
  double halfwayX = minX + abs(abs(maxX)-abs(minX))/2;
  double halfwayY = minY + abs(abs(maxY)-abs(minY))/2;
  //printf("halfway x: %f y: %f \n", halfwayX, halfwayY);

  // new boxes
  for(int i=0; i<inputBox.n; i++) {
    // (0,0)
    if(inputBox.x[i] < halfwayX && inputBox.y[i] < halfwayY) {
      boxes[0].x[boxes[0].n] = inputBox.x[i];
      boxes[0].y[boxes[0].n] = inputBox.y[i];
      boxes[0].z[boxes[0].n] = inputBox.z[i];
      boxes[0].t[boxes[0].n] = inputBox.t[i];
      boxes[0].i[boxes[0].n] = inputBox.i[i];
      boxes[0].n++;
    }
    // (0,1)
    if(inputBox.x[i] >= halfwayX && inputBox.y[i] < halfwayY) {
      boxes[1].x[boxes[1].n] = inputBox.x[i];
      boxes[1].y[boxes[1].n] = inputBox.y[i];
      boxes[1].z[boxes[1].n] = inputBox.z[i];
      boxes[1].t[boxes[1].n] = inputBox.t[i];
      boxes[1].i[boxes[1].n] = inputBox.i[i];
      boxes[1].n++;
    }
    // (1,0)
    if(inputBox.x[i] < halfwayX && inputBox.y[i] >= halfwayY) {
      boxes[2].x[boxes[2].n] = inputBox.x[i];
      boxes[2].y[boxes[2].n] = inputBox.y[i];
      boxes[2].z[boxes[2].n] = inputBox.z[i];
      boxes[2].t[boxes[2].n] = inputBox.t[i];
      boxes[2].i[boxes[2].n] = inputBox.i[i];
      boxes[2].n++;
    }
    // (1,1)
    if(inputBox.x[i] >= halfwayX && inputBox.y[i] >= halfwayY) {
      boxes[3].x[boxes[3].n] = inputBox.x[i];
      boxes[3].y[boxes[3].n] = inputBox.y[i];
      boxes[3].z[boxes[3].n] = inputBox.z[i];
      boxes[3].t[boxes[3].n] = inputBox.t[i];
      boxes[3].i[boxes[3].n] = inputBox.i[i];
      boxes[3].n++;
    }
  }
  printf("4 way split - 0: %i 1: %i 2: %i 3: %i \n", boxes[0].n, boxes[1].n, boxes[2].n, boxes[3].n);
  int numBoxes = 4;

  // what kind of aspect ratio do the boxes have
  // is the same for all boxes currently...
  //for(int i=0; i<numBoxes; i++) {
    double maX = max(boxes[0].x,boxes[0].n);
    double maY = max(boxes[0].y,boxes[0].n);
    double miX = min(boxes[0].x,boxes[0].n);
    double miY = min(boxes[0].y,boxes[0].n);
    double diffX = abs(maX - miX);
    double diffY = abs(maY - miY);
    //printf("diff: x %f y %f \n", diffX, diffY);
  //}

  // first check which way it makes more sense to merge
  // ok to merge the way that keeps the boxes squarer?
  if(diffY < diffX) {
    // wide boxes, so preferable to merge 0,2 + 1,3 (vertically)
    int n1 = boxes[0].n + boxes[2].n;
    int n2 = boxes[1].n + boxes[3].n;
    if(n1 > minNumStationsInBox && n2 > minNumStationsInBox) {
      // merge vertically
      printf("best to merge vertically \n");
      boxes[0] = merge_boxes(boxes[0],boxes[2]);
      boxes[1] = merge_boxes(boxes[1],boxes[3]);
    }
    else {
      // merge horizontally
      printf("had to merge horizontally \n");
      boxes[0] = merge_boxes(boxes[0],boxes[1]);
      boxes[1] = merge_boxes(boxes[2],boxes[3]);
    }
  }
  else { // diffY > diffX
    // tall boxes, so preferable to merge 0,1 + 2,3 (horizontally)
    int n1 = boxes[0].n + boxes[1].n;
    int n2 = boxes[2].n + boxes[3].n;
    if(n1 > minNumStationsInBox && n2 > minNumStationsInBox) {
      // merge horizontally
      printf("best to merge horizontally \n");
      boxes[0] = merge_boxes(boxes[0],boxes[1]);
      boxes[1] = merge_boxes(boxes[2],boxes[3]);
    }
    else {
      // merge vertically
      printf("had to merge vertically \n");
      boxes[0] = merge_boxes(boxes[0],boxes[2]);
      boxes[1] = merge_boxes(boxes[1],boxes[3]);
    }
  }
  // don't need these boxes anymore
  for(int i=2; i<numBoxes; i++) {
    printf("free box: %i \n", i);
    free(boxes[i].x);
    free(boxes[i].y);
    free(boxes[i].z);
    free(boxes[i].t);
    free(boxes[i].i);
  }
  printf("2 way split - 0: %i 1: %i \n", boxes[0].n, boxes[1].n);
  numBoxes = 2;

  // loop over the boxes
  for(int i=0; i<numBoxes; i++) {
    int n_temp = boxes[i].n;
    //printf("still too big or return? %i i: %i \n", n_temp, i);
    if(n_temp > maxNumStationsInBox) {
      printf("box still too big %i (being further recursively split)\n", n_temp);
      // split the box further
      split_box(maxNumStationsInBox, minNumStationsInBox, boxes[i], finalNumBoxes, finalBoxes);
    }
    else {
      printf("box size: %i, being returned \n", n_temp);
      int current_n = *finalNumBoxes;
      //printf("current total number of boxes: %i \n", current_n);
      // add to list of boxes (finalBoxes)
      ((*finalBoxes)[current_n]).n = boxes[i].n;
      ((*finalBoxes)[current_n]).x = boxes[i].x;
      ((*finalBoxes)[current_n]).y = boxes[i].y;
      ((*finalBoxes)[current_n]).z = boxes[i].z;
      ((*finalBoxes)[current_n]).t = boxes[i].t;
      ((*finalBoxes)[current_n]).i = boxes[i].i;
      // increment the number of boxes
      (*finalNumBoxes)++;
    }
  }
}

struct box * control_box_division(int maxNumStationsInBox, int minNumStationsInBox, struct box inputBox) {

  // check is isn't already smaller than the max
  if(inputBox.n < maxNumStationsInBox) {
    struct box * nAndBoxes = malloc(sizeof(struct box) * 2);
    // just return
    nAndBoxes[0].n = 1;
    nAndBoxes[1].n = inputBox.n;
    nAndBoxes[1].x = inputBox.x;
    nAndBoxes[1].y = inputBox.y;
    nAndBoxes[1].z = inputBox.z;
    nAndBoxes[1].t = inputBox.t;
    nAndBoxes[1].i = inputBox.i;

    return nAndBoxes;
  }

  struct box * totalBoxes;
  int maxNumBoxes = floor(inputBox.n/minNumStationsInBox);
  printf("allocating memory for potential max number of boxes: %i \n",  maxNumBoxes);
  totalBoxes = malloc(sizeof(struct box) * maxNumBoxes);
  int totalNumBoxes = 0;

  // pass the outputs in by reference (so the recursive function can add to them)
  split_box(maxNumStationsInBox, minNumStationsInBox, inputBox, &totalNumBoxes, &totalBoxes);

  printf("total number of boxes: %i \n", totalNumBoxes);
  struct box * nAndBoxes = malloc(sizeof(struct box) * (totalNumBoxes+1));
  // (TODO: package up a pointer to also return the number of boxes)
  // easiest just to use first box to encode how many boxes there are
  nAndBoxes[0].n = totalNumBoxes;
  for(int b=1; b<(totalNumBoxes+1); b++) {
    nAndBoxes[b].n = totalBoxes[b-1].n;
    nAndBoxes[b].x = totalBoxes[b-1].x;
    nAndBoxes[b].y = totalBoxes[b-1].y;
    nAndBoxes[b].z = totalBoxes[b-1].z;
    nAndBoxes[b].t = totalBoxes[b-1].t;
    nAndBoxes[b].i = totalBoxes[b-1].i;
  }
  free(totalBoxes);
  // return a list of boxes
  return nAndBoxes;
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
  return mean;
}

double max(const double *array, int sizeArray)
{
  double max = array[0];
  for(int i=0; i<sizeArray; i++) {
    if(array[i] > max) {
        max = array[i];
    }
  }
  return max;
}

double min(const double *array, int sizeArray)
{
  double min = array[0];
  for(int i=0; i<sizeArray; i++) {
    if(array[i] < min) {
        min = array[i];
    }
  }
  return min;
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
