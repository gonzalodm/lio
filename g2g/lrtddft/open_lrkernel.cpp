#include <iostream>
#include <xc.h>
#include <omp.h>

#include <stdio.h>
#include <string.h> 

#include "../common.h"
#include "../init.h"
#include "../partition.h"

#include "obtain.h"
#include "eri.h"
#include "centros.h"
#include "../libxc/libxcproxy.h"

#define DENSMIN 1e-5

using namespace G2G;
extern Partition partition;

//######################################################################
//######################################################################
extern "C" void g2g_open_calculate2e_(double* Ta,double* Tb,double* Cbas,
                         int& numvec,double* Fa,double* Fb,int& int2elec)
{
   int M = fortran_vars.m;
   int M2 = M*M;
   int M3 = M2*M;
   int s_func = fortran_vars.s_funcs;
   int p_func = fortran_vars.p_funcs;
   int d_func = fortran_vars.d_funcs;
   double* aContr = &fortran_vars.a_values(0,0);
   double* pos = &fortran_vars.atom_positions_pointer(0,0);
   double timeI, timeF;
   uint* ncont = &fortran_vars.contractions(0);
   uint* nuc = &fortran_vars.nucleii(0);

   timeI = timeF = 0.0f;
   if (int2elec == 0 ) {
     timeI = omp_get_wtime();
     eri(M,fortran_vars.atoms,ncont,Cbas,aContr,pos,nuc,
         s_func,p_func,d_func);
     timeF = omp_get_wtime();
     printf(" ERI SUBROUTINE %f\n",timeF-timeI);
   }
   Obtain fock(numvec,M,fortran_vars.dim);
   fock.calculate(Ta,Tb,fortran_vars.Kmat,Fa,Fb);

   fflush(stdout); // NOT BUFFERED
}
//######################################################################
//######################################################################
extern "C" void g2g_open_calculatedft_(double* Ta,double* Tb,double* Ca,double* Cb,
                                       double* Fa,double* Fb, int& NCOa,int& NCOb)
{
   partition.open_solve_lr(Ta,Tb,Ca,Cb,Fa,Fb,NCOa,NCOb);
   fflush(stdout); // NOT BUFFERED
}

//######################################################################
//######################################################################

namespace G2G {

void Partition::open_solve_lr(double* Ta,double* Tb,double* Ca,double* Cb,
                              double* Fa,double* Fb,int& NCOa,int& NCOb)
{

#pragma omp parallel for schedule(static)
    for(uint i=0;i<work.size();i++) {
      for(uint j=0;j<work[i].size();j++) {
         int ind = work[i][j];
         if(ind >= cubes.size()) {
           spheres[ind - cubes.size()]->open_solve_closed_lr(Ta,Tb,Ca,Cb,Fa,Fb,NCOa,NCOb);
         } else {
           cubes[ind]->open_solve_closed_lr(Ta,Tb,Ca,Cb,Fa,Fb,NCOa,NCOb);
         }
      }
   }
   fflush(stdout);
}
//######################################################################
//######################################################################

//######################################################################
//######################################################################

template<class scalar_type> void PointGroupCPU<scalar_type>::
               open_solve_closed_lr(double* Ta,double* Tb,double* Ca,double* Cb,
                                    double* Fa,double* Fb,int& NCOa,int& NCOb)
{
//=================================================
// INPUTS
   // Ta = Trial Vectors of alpha
   // Tb = Trial Vectors of alpha
   // Ca = SCF Coeficients of alpha
   // Cb = SCF Coeficients of beta
   // NCOa = number of occupied orbitals alpha
   // NCOb = number of occupied orbitals beta

// OUTPUTS
   // Fa = Dft Fock alpha
   // Fb = Dft Fock beta
//=================================================

   const uint group_m = this->total_functions();
   const int npoints = this->points.size();
   bool lda = false;
   bool compute_forces = true;
   compute_functions(compute_forces,!lda);

   // Densities alpha and beta
   HostMatrix<scalar_type> rmm_input_a(group_m, group_m);
   HostMatrix<scalar_type> rmm_input_b(group_m, group_m);

   // Obtain reduced density matrix for this group
   get_rmm_input(rmm_input_a, rmm_input_b);

   // Obtain indexes of basis functions for this group
   int* numeros = new int[group_m];
   get_coef_input(rmm_input_a,numeros); // rmm_input_a is not used

   // Obtain reduced Coefficients and Vectors for this group
   HostMatrix<scalar_type> Cred_a(NCOa,group_m);
   HostMatrix<scalar_type> Cred_b(NCOb,group_m);
   HostMatrix<scalar_type> Tred_a(NCOa,group_m);
   HostMatrix<scalar_type> Tred_b(NCOb,group_m);
   int M = fortran_vars.m;

/*
   for(int ii=0; ii<group_m; ii++) {
     int col = numeros[ii];
     // For alpha Coeff
     for(int jj=0; jj<NCOa; jj++) {
       Tred_a(jj,ii) = Ta[col*M+jj];
       Cred_a(jj,ii) = Ca[jj*M+col];
     }
     // For beta Coeff
     for(int jj=0; jj<NCOb; jj++) {
       Tred_b(jj,ii) = Tb[col*M+jj];
       Cred_b(jj,ii) = Cb[jj*M+col];
     }
   }
*/

   for(int ii=0; ii<group_m; ii++) {
     int col = numeros[ii];
     // For alpha Coeff
     for(int jj=0; jj<NCOa; jj++) {
       Tred_a(jj,ii) = Ta[col*M+jj];
       Cred_a(jj,ii) = Ca[jj*M+col];

       Tred_b(jj,ii) = Tb[col*M+jj];
       Cred_b(jj,ii) = Cb[jj*M+col];
     }
     // For beta Coeff
     for(int jj=NCOa; jj<NCOb; jj++) {
       Tred_b(jj,ii) = Tb[col*M+jj];
       Cred_b(jj,ii) = Cb[jj*M+col];
     }
   }

   double* grdMO_a = (double*)malloc(4*NCOa*sizeof(double));
   double* grdMO_b = (double*)malloc(4*NCOb*sizeof(double));

// Reduced Fock alpha and beta for this group
   double* smallFock_a  = (double*)malloc(group_m*group_m*sizeof(double));
   memset(smallFock_a,0.0f,group_m*group_m*sizeof(double));
   double* smallFock_b  = (double*)malloc(group_m*group_m*sizeof(double));
   memset(smallFock_b,0.0f,group_m*group_m*sizeof(double));

   double* rho = (double*)malloc(8*sizeof(double));
   double* tra = (double*)malloc(8*sizeof(double));
   double* lrCoef = (double*)malloc(16*sizeof(double));
   double* precond_1 = (double*)malloc(group_m*sizeof(double));
   double* precond_2 = (double*)malloc(group_m*sizeof(double));

//LIBXC INITIALIZATION
   const int nspin = XC_POLARIZED;
   const int functionalExchange = fortran_vars.ex_functional_id; //101;
   const int functionalCorrelation = fortran_vars.ec_functional_id; // 130;
   LibxcProxy<scalar_type,3> libxcProxy(functionalExchange, functionalCorrelation, nspin);

   // Loop over grid points of this group
   for (int point = 0; point < npoints; point++) {
     scalar_type pd_a, pd_b;
     scalar_type tdx_a, tdy_a, tdz_a, tdx_b, tdy_b, tdz_b;  // dx, dy, dz
     pd_a = tdx_a = tdy_a = tdz_a = 0.0f;
     pd_b = tdx_b = tdy_b = tdz_b = 0.0f;

     // Evaluated basis functions, derivatives and gradient for the point
     const scalar_type* fv = function_values.row(point);
     const scalar_type* gxv = gX.row(point);    // dx
     const scalar_type* gyv = gY.row(point);    // dy
     const scalar_type* gzv = gZ.row(point);    // dz

     // Loop over basis functions
     for (int i = 0; i < group_m; i++) {
       scalar_type w_a, w_b;
       scalar_type w3x_a, w3y_a, w3z_a, w3x_b, w3y_b, w3z_b;

       w_a = w3x_a = w3y_a = w3z_a = 0.0f;
       w_b = w3x_b = w3y_b = w3z_b = 0.0f;

       const scalar_type* rmm_a = rmm_input_a.row(i);
       const scalar_type* rmm_b = rmm_input_b.row(i);

       for (int j = 0; j <= i; j++) {
         const scalar_type rmm_aj = rmm_a[j];
         const scalar_type rmm_bj = rmm_b[j];
         // alpha
         w_a += fv[j] * rmm_aj;
         w3x_a += gxv[j] * rmm_aj;
         w3y_a += gyv[j] * rmm_aj;
         w3z_a += gzv[j] * rmm_aj;
         // beta
         w_b += fv[j] * rmm_bj;
         w3x_b += gxv[j] * rmm_bj;
         w3y_b += gyv[j] * rmm_bj;
         w3z_b += gzv[j] * rmm_bj;
       }
       const scalar_type Fi = fv[i];
       const scalar_type gx = gxv[i], gy = gyv[i], gz = gzv[i];
       // Alpha results
       pd_a += Fi * w_a;
       tdx_a += gx * w_a + w3x_a * Fi;
       tdy_a += gy * w_a + w3y_a * Fi;
       tdz_a += gz * w_a + w3z_a * Fi;
       // Beta results
       pd_b += Fi * w_b;
       tdx_b += gx * w_b + w3x_b * Fi;
       tdy_b += gy * w_b + w3y_b * Fi;
       tdz_b += gz * w_b + w3z_b * Fi;
     }

     // Transform IA --> AO
     memset(grdMO_a,0.0f,NCOa*4*sizeof(double));
     memset(grdMO_b,0.0f,NCOb*4*sizeof(double));
/*
     for(int ii=0; ii<group_m; ii++) {
       // ALPHA
       for(int jj=0; jj<NCOa; jj++) {
          grdMO_a[jj*4]   += fv[ii]  * Cred_a(jj,ii);
          grdMO_a[jj*4+1] += gxv[ii] * Cred_a(jj,ii);
          grdMO_a[jj*4+2] += gyv[ii] * Cred_a(jj,ii);
          grdMO_a[jj*4+3] += gzv[ii] * Cred_a(jj,ii);
       }
       // BETA
       for(int jj=0; jj<NCOb; jj++) {
          grdMO_b[jj*4]   += fv[ii]  * Cred_b(jj,ii);
          grdMO_b[jj*4+1] += gxv[ii] * Cred_b(jj,ii);
          grdMO_b[jj*4+2] += gyv[ii] * Cred_b(jj,ii);
          grdMO_b[jj*4+3] += gzv[ii] * Cred_b(jj,ii);
       }
     }
*/

     for(int ii=0; ii<group_m; ii++) {
       // ALPHA
       for(int jj=0; jj<NCOa; jj++) {
          grdMO_a[jj*4]   += fv[ii]  * Cred_a(jj,ii);
          grdMO_a[jj*4+1] += gxv[ii] * Cred_a(jj,ii);
          grdMO_a[jj*4+2] += gyv[ii] * Cred_a(jj,ii);
          grdMO_a[jj*4+3] += gzv[ii] * Cred_a(jj,ii);

          grdMO_b[jj*4]   += fv[ii]  * Cred_b(jj,ii);
          grdMO_b[jj*4+1] += gxv[ii] * Cred_b(jj,ii);
          grdMO_b[jj*4+2] += gyv[ii] * Cred_b(jj,ii);
          grdMO_b[jj*4+3] += gzv[ii] * Cred_b(jj,ii);
       }
       // BETA
       for(int jj=NCOa; jj<NCOb; jj++) {
          grdMO_b[jj*4]   += fv[ii]  * Cred_b(jj,ii);
          grdMO_b[jj*4+1] += gxv[ii] * Cred_b(jj,ii);
          grdMO_b[jj*4+2] += gyv[ii] * Cred_b(jj,ii);
          grdMO_b[jj*4+3] += gzv[ii] * Cred_b(jj,ii);
       }
     }


     double D0_a, DX_a, DY_a, DZ_a; D0_a = DX_a = DY_a = DZ_a = 0.0f;
     double D0_b, DX_b, DY_b, DZ_b; D0_b = DX_b = DY_b = DZ_b = 0.0f;
     double red_a, redx_a, redy_a, redz_a; red_a = redx_a = redy_a = redz_a = 0.0f;
     double red_b, redx_b, redy_b, redz_b; red_b = redx_b = redy_b = redz_b = 0.0f;
/*
     for(int ii=0; ii<group_m; ii++) {
        double Fi = fv[ii];
        double gxi = gxv[ii];
        double gyi = gyv[ii];
        double gzi = gzv[ii];
        // ALPHA
        D0_a = DX_a = DY_a = DZ_a = 0.0f;
        for(int jj=0; jj<NCOa; jj++) {
           D0_a += Tred_a(jj,ii) * grdMO_a[jj*4];
           DX_a += Tred_a(jj,ii) * grdMO_a[jj*4+1];
           DY_a += Tred_a(jj,ii) * grdMO_a[jj*4+2];
           DZ_a += Tred_a(jj,ii) * grdMO_a[jj*4+3];
        }
        // BETA
        D0_b = DX_b = DY_b = DZ_b = 0.0f;
        for(int jj=0; jj<NCOb; jj++) {
           D0_b += Tred_b(jj,ii) * grdMO_b[jj*4];
           DX_b += Tred_b(jj,ii) * grdMO_b[jj*4+1];
           DY_b += Tred_b(jj,ii) * grdMO_b[jj*4+2];
           DZ_b += Tred_b(jj,ii) * grdMO_b[jj*4+3];
        }
        // ALPHA
        red_a  += D0_a * Fi;
        redx_a += D0_a * gxi + DX_a * Fi;
        redy_a += D0_a * gyi + DY_a * Fi;
        redz_a += D0_a * gzi + DZ_a * Fi;
        // BETA
        red_b  += D0_b * Fi;
        redx_b += D0_b * gxi + DX_b * Fi;
        redy_b += D0_b * gyi + DY_b * Fi;
        redz_b += D0_b * gzi + DZ_b * Fi;
     }
*/

     for(int ii=0; ii<group_m; ii++) {
        double Fi = fv[ii];
        double gxi = gxv[ii];
        double gyi = gyv[ii];
        double gzi = gzv[ii];
        // ALPHA
        D0_a = DX_a = DY_a = DZ_a = 0.0f;
        D0_b = DX_b = DY_b = DZ_b = 0.0f;
        for(int jj=0; jj<NCOa; jj++) {
           D0_a += Tred_a(jj,ii) * grdMO_a[jj*4];
           DX_a += Tred_a(jj,ii) * grdMO_a[jj*4+1];
           DY_a += Tred_a(jj,ii) * grdMO_a[jj*4+2];
           DZ_a += Tred_a(jj,ii) * grdMO_a[jj*4+3];

           D0_b += Tred_b(jj,ii) * grdMO_b[jj*4];
           DX_b += Tred_b(jj,ii) * grdMO_b[jj*4+1];
           DY_b += Tred_b(jj,ii) * grdMO_b[jj*4+2];
           DZ_b += Tred_b(jj,ii) * grdMO_b[jj*4+3];
        }
        // BETA
        for(int jj=NCOa; jj<NCOb; jj++) {
           D0_b += Tred_b(jj,ii) * grdMO_b[jj*4];
           DX_b += Tred_b(jj,ii) * grdMO_b[jj*4+1];
           DY_b += Tred_b(jj,ii) * grdMO_b[jj*4+2];
           DZ_b += Tred_b(jj,ii) * grdMO_b[jj*4+3];
        }
        // ALPHA
        red_a  += D0_a * Fi;
        redx_a += D0_a * gxi + DX_a * Fi;
        redy_a += D0_a * gyi + DY_a * Fi;
        redz_a += D0_a * gzi + DZ_a * Fi;
        // BETA
        red_b  += D0_b * Fi;
        redx_b += D0_b * gxi + DX_b * Fi;
        redy_b += D0_b * gyi + DY_b * Fi;
        redz_b += D0_b * gzi + DZ_b * Fi;
     }


     // COPY RESULTS TO POINTER
     rho[0] = pd_a; rho[1] = tdx_a; rho[2] = tdy_a; rho[3] = tdz_a;
     rho[4] = pd_b; rho[5] = tdx_b; rho[6] = tdy_b; rho[7] = tdz_b;

     tra[0] = red_a; tra[1] = redx_a; tra[2] = redy_a; tra[3] = redz_a;
     tra[4] = red_b; tra[5] = redx_b; tra[6] = redy_b; tra[7] = redz_b;

     // Derivatives Calculate: Obtain Coefficients of integrals
     libxcProxy.coefLR(rho,tra,lrCoef);

     const scalar_type wp = this->points[point].weight;
     double temp1, temp2, temp3, temp4, precii_1, precii_2, result;

     // ALPHA 
     for(int ii=0; ii<group_m; ii++) {
        // alpha-alpha
        temp1 = lrCoef[0] * 0.5f * fv[ii];
        temp2 = lrCoef[1] * (rho[1]*gxv[ii]+rho[2]*gyv[ii]+rho[3]*gzv[ii]);
        temp3 = lrCoef[2] * (rho[5]*gxv[ii]+rho[6]*gyv[ii]+rho[7]*gzv[ii]);
        temp4 = lrCoef[3] * (tra[1]*gxv[ii]+tra[2]*gyv[ii]+tra[3]*gzv[ii]);
        precond_1[ii] = (temp1+temp2+temp3+temp4) * wp;
        // alpha-beta
        temp1 = lrCoef[12] * 0.5f * fv[ii];
        temp2 = lrCoef[13] * (rho[1]*gxv[ii]+rho[2]*gyv[ii]+rho[3]*gzv[ii]);
        temp3 = lrCoef[14] * (rho[5]*gxv[ii]+rho[6]*gyv[ii]+rho[7]*gzv[ii]);
        temp4 = lrCoef[15] * (tra[5]*gxv[ii]+tra[6]*gyv[ii]+tra[7]*gzv[ii]);
        precond_2[ii] = (temp1+temp2+temp3+temp4) * wp;
        precii_1 = precond_1[ii];
        precii_2 = precond_2[ii];

        for(int jj=0; jj<=ii; jj++) {
           result = precond_1[jj] * fv[ii] + precii_1 * fv[jj];
           result += precond_2[jj] * fv[ii] + precii_2 * fv[jj];
           smallFock_a[ii*group_m+jj] += result;
        }
     }

     // BETA
     for(int ii=0; ii<group_m; ii++) {
        // beta-beta
        temp1 = lrCoef[8] * 0.5f * fv[ii];
        temp2 = lrCoef[9] * (rho[5]*gxv[ii]+rho[6]*gyv[ii]+rho[7]*gzv[ii]);
        temp3 = lrCoef[10] * (rho[1]*gxv[ii]+rho[2]*gyv[ii]+rho[3]*gzv[ii]);
        temp4 = lrCoef[11] * (tra[5]*gxv[ii]+tra[6]*gyv[ii]+tra[7]*gzv[ii]);
        precond_1[ii] = (temp1+temp2+temp3+temp4) * wp;
        // beta-alpha
        temp1 = lrCoef[4] * 0.5f * fv[ii];
        temp2 = lrCoef[5] * (rho[5]*gxv[ii]+rho[6]*gyv[ii]+rho[7]*gzv[ii]);
        temp3 = lrCoef[6] * (rho[1]*gxv[ii]+rho[2]*gyv[ii]+rho[3]*gzv[ii]);
        temp4 = lrCoef[7] * (tra[1]*gxv[ii]+tra[2]*gyv[ii]+tra[3]*gzv[ii]);
        precond_2[ii] = (temp1+temp2+temp3+temp4) * wp;
        precii_1 = precond_1[ii];
        precii_2 = precond_2[ii];
        for(int jj=0; jj<=ii; jj++) {
           result = precond_1[jj] * fv[ii] + precii_1 * fv[jj];
           result += precond_2[jj] * fv[ii] + precii_2 * fv[jj];
           smallFock_b[ii*group_m+jj] += result;
        }
     }

   } // END POINTS

   // OBTAIN FOCK GLOBAL
   for(int ii=0; ii<group_m; ii++) {
     int row = numeros[ii];
#pragma omp critical
     for(int jj=0; jj<=ii; jj++) {
       int col = numeros[jj];
       // alpha
       Fa[row*M+col] += smallFock_a[ii*group_m+jj];
       Fa[col*M+row] = Fa[row*M+col];
       // beta
       Fb[row*M+col] += smallFock_b[ii*group_m+jj];
       Fb[col*M+row] = Fb[row*M+col];
     }
   }

   // Free Memory
   free(grdMO_a);       grdMO_a = NULL;
   free(grdMO_b);       grdMO_b = NULL;
   free(smallFock_a);   smallFock_a = NULL;
   free(smallFock_b);   smallFock_b = NULL;
   free(precond_1);     precond_1 = NULL;
   free(precond_2);     precond_2 = NULL;
   free(lrCoef);        lrCoef = NULL;
   delete[] numeros;    numeros = NULL;
}
//######################################################################
//######################################################################
#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#endif
}
