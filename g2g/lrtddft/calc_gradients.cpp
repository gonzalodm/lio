#include <iostream>
#include <stdio.h>
#include <string.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"
#include "../libxc/libxcproxy.h"

#include "calc_VXC.h"
#include "calc_FXC.h"

using namespace G2G;
using namespace std;

void calc_gradients(double* dens, double* diff, double* trad, double sigma)
{

// LIBXC INITIALIZATION
   const int nspin = XC_POLARIZED;
   const int functionalExchange = fortran_vars.ex_functional_id; //101;
   const int functionalCorrelation = fortran_vars.ec_functional_id; // 130;
   double pd = dens[0];
   LibxcProxy<double,3> libxcProxy(functionalExchange, functionalCorrelation, nspin);

// OUTPUTS FOR LIBXC
   double* vrho        = (double*)malloc(4*sizeof(double));
   double* vsigma      = (double*)malloc(5*sizeof(double));
   double* v2rho2      = (double*)malloc(5*sizeof(double));
   double* v2rhosigma  = (double*)malloc(8*sizeof(double));
   double* v2sigma2    = (double*)malloc(8*sizeof(double));
   double* v3rho3      = (double*)malloc(6*sizeof(double));
   double* v3rho2sigma = (double*)malloc(11*sizeof(double));
   double* v3rhosigma2 = (double*)malloc(14*sizeof(double));
   double* v3sigma3    = (double*)malloc(12*sizeof(double));

/*
// DEBUG
   // punto 1
   dens[0] = 0.20136241320638337;
   dens[1] = -1.4257645526048659E-003;
   dens[2] = 0.0000000000000000;
   dens[3] = -7.6834757112023463E-002;

   diff[0] = -4.6926890412794395E-002;
   diff[1] = 2.6949332892050394E-004;
   diff[2] = 0.0000000000000000;
   diff[3] = 4.0465227355638206E-002;
   pd = dens[0] * 2.0f;
   sigma = 5.9056127050231052E-003 * 4.0f;
   trad[0] = -4.8246479671649509E-002;
   trad[1] = 2.7974624993282616E-004;
   trad[2] = 0.0000000000000000;
   trad[3] = -2.8474564875021414E-002;
// END DEBUG
*/


// LIBXC CALCULATE DERIVATIVES
   libxcProxy.coefZv(&pd,&sigma,
                     vrho,vsigma,
                     v2rho2,v2rhosigma,v2sigma2,
                     v3rho3,v3rho2sigma,v3rhosigma2,v3sigma3);
   
   double* DDUM = (double*)malloc(4*sizeof(double));
   double* PDUM = (double*)malloc(4*sizeof(double));
   double* VDUM = (double*)malloc(4*sizeof(double));
   memset(DDUM,0.0f,4*sizeof(double));
   memset(PDUM,0.0f,4*sizeof(double));
   memset(VDUM,0.0f,4*sizeof(double));

   calc_VXC(dens,diff,
            vrho,vsigma,
            v2rho2,v2rhosigma,v2sigma2,
            v3rho3,v3rho2sigma,v3rhosigma2,v3sigma3,
            DDUM,PDUM);

   calc_FXC(dens,trad,
            vrho,vsigma,
            v2rho2,v2rhosigma,v2sigma2,
            v3rho3,v3rho2sigma,v3rhosigma2,v3sigma3,
            DDUM,VDUM);
  
   // COPY RESULTS
   for(int i=0; i<4; i++) {
      dens[i] = DDUM[i];
      diff[i] = PDUM[i];
      trad[i] = VDUM[i];
   }

   free(DDUM); DDUM = NULL;
   free(PDUM); PDUM = NULL;
   free(VDUM); VDUM = NULL;
}
