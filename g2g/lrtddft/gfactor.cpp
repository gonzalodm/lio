#include <iostream>
#include <xc.h>
#include <omp.h>

#include <stdio.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"
#include "../libxc/libxcproxy.h"

#include "eri.h"
#include "couplingform.h"
#include "calc_coef.h"

#define DENSMIN 1e-5

using namespace G2G;
extern Partition partition;

//######################################################################
//######################################################################
extern "C" void g2g_calculateg_(double* Tmat,double* F,int& DER)
{
   printf("ESTOY EN CALCULATE GXC, DER=%d\n",DER);
   partition.solve_Gxc(Tmat,F,DER);
   fflush(stdout); // NOT BUFFERED
}
//######################################################################
//######################################################################

namespace G2G {

void Partition::solve_Gxc(double* Tmat,double* F,int& DER)
{

   double timeI, timeF;
   timeI = omp_get_wtime();
#pragma omp parallel for schedule(static)
   for(uint i=0;i<work.size();i++) {
      for(uint j=0;j<work[i].size();j++) {
         int ind = work[i][j];
         if(ind >= cubes.size()) {
           spheres[ind-cubes.size()]->solve_3rd_der(Tmat,F,DER);
         } else {
           cubes[ind]->solve_3rd_der(Tmat,F,DER);
         }
      }
   }
   timeF = omp_get_wtime();
   //printf("SOLVE_CLOSED_LR SUBROUTINE %f\n",timeF-timeI);
   fflush(stdout);
}
//######################################################################
//######################################################################

//######################################################################
//######################################################################

template<class scalar_type> void PointGroupCPU<scalar_type>::
               solve_3rd_der(double* Tmat,double* F,int& DER)
{
   const uint group_m = this->total_functions();
   const int npoints = this->points.size();
   bool lda = false;
   bool compute_forces = true;
   compute_functions(compute_forces,!lda);
   HostMatrix<scalar_type> rmm_input(group_m,group_m);
   int* numeros = new int[group_m];
   int M = fortran_vars.m;
   get_rmm_input(rmm_input);
   get_coef_input(rmm_input,numeros);

   double* zcoef = (double*) malloc(3*sizeof(double));
   double* precond = (double*) malloc(group_m*sizeof(double));

   int row, col;

// FORMAMOS LA TRANSITION DENSITY REDUCIDA
   //cout << "group_m: " << group_m << endl;
   HostMatrix<double> tred(group_m,group_m);
   for(int i=0; i<group_m; i++) {
     row = numeros[i];
     tred(i,i) = Tmat[row*M+row];
     //cout << i << " " <<Tmat[row*M+row] << endl;
     for(int j=0; j<i; j++) {
       col = numeros[j];
       tred(i,j) = Tmat[row*M+col] + Tmat[col*M+row];
       //cout << i << " " << j << " " <<Tmat[row*M+col] << " " << Tmat[col*M+row] << endl;
     }
   }

   for(int point=0;point<npoints;point++) {
      scalar_type pd, pdx, pdy, pdz; pd = pdx = pdy = pdz = 0.0f;
      scalar_type red, redx, redy, redz; red = redx = redy = redz = 0.0f;
      const scalar_type* fv = function_values.row(point);
      const scalar_type* gxv = gX.row(point);
      const scalar_type* gyv = gY.row(point);
      const scalar_type* gzv = gZ.row(point);
      for(int i=0;i<group_m;i++) {
         double w3xc, w3yc, w3zc, w; w3xc = w3yc = w3zc = w = 0.0f;
         double z3xc, z3yc, z3zc, z; z3xc = z3yc = z3zc = z = 0.0f;
         const scalar_type* rm = rmm_input.row(i);
         for(int j=0;j<=i;j++) {
            const double rmj = rm[j];
            // Density
            w += fv[j] * rmj;
            w3xc += gxv[j] * rmj;
            w3yc += gyv[j] * rmj;
            w3zc += gzv[j] * rmj;
            // Transition Density
            z += fv[j] * tred(i,j);
            z3xc += gxv[j] * tred(i,j);
            z3yc += gyv[j] * tred(i,j);
            z3zc += gzv[j] * tred(i,j);
         }
         const double Fi = fv[i];
         const double gx = gxv[i], gy = gyv[i], gz = gzv[i];
         // Density
         pd += Fi * w;
         pdx += gx * w + w3xc * Fi;
         pdy += gy * w + w3yc * Fi;
         pdz += gz * w + w3zc * Fi;
         // Transition Density
         red += Fi * z;
         redx += gx * z + z3xc * Fi;
         redy += gy * z + z3yc * Fi;
         redz += gz * z + z3zc * Fi;
      }
      double sigma = pdx * pdx + pdy * pdy + pdz * pdz;
/*
//    DEBUG 
      pd = 8.6792344676238335E-006 * 2.0f;
      sigma = 6.1329471537828528E-010 * 4.0f;
      pdx = -1.2195835377375016E-005;
      pdy = 1.6363583240860429E-005; pdz = -1.4028166606712547E-005;
      red = -1.8661364453244286E-005;
      redx = 3.3089162612377751E-006;
      redy = 1.8056171389542832E-005; redz = 5.4889385076820334E-005;
*/
      pdx *= 0.5f; pdy *= 0.5f; pdz *= 0.5f;

      if (DER == 2 ) {
          // INIT LIBXC
          double cruz = (redx * pdx + redy * pdy + redz * pdz);
          const int nspin = XC_POLARIZED;
          const int functionalExchange = fortran_vars.ex_functional_id; //101;
          const int functionalCorrelation = fortran_vars.ec_functional_id; // 130;
          LibxcProxy<scalar_type,3> libxcProxy(functionalExchange, functionalCorrelation, nspin);
          // RUN LIBXC
          libxcProxy.coefLR(&pd,&sigma,red,cruz,zcoef);
      } else {
          coef_calculator(pd,sigma,pdx,pdy,pdz,
                          red,redx,redy,redz,zcoef);
      }

      const scalar_type wp = this->points[point].weight;
      double term1, term2, term3, term4, precondii, result;
      for(int i=0; i<group_m; i++) {
        row = numeros[i];
        term1 = zcoef[0] * 0.5f * fv[i] + zcoef[1] * pdx * gxv[i];
        term2 = zcoef[1] * pdy * gyv[i] + zcoef[1] * pdz * gzv[i];
        term3 = zcoef[2] * redx * gxv[i] + zcoef[2] * redy * gyv[i];
        term4 = zcoef[2] * redz * gzv[i];
        precond[i] = (term1 + term2 + term3 + term4) * wp;
        precondii = precond[i];
        for(int j=0; j<=i; j++) {
           col = numeros[j];
           result = fv[i] * precond[j] + precondii * fv[j];
           F[row*M+col] += result;
        }
      }

   }  // END points loop

   // Free Memory
   delete[] numeros;
   free(precond);
   free(zcoef);
}
//######################################################################
//######################################################################

//######################################################################
//######################################################################
/*
template<class scalar_type>
void PointGroupCPU<scalar_type>::get_coef_input(HostMatrix<scalar_type>&
               rmm_input,int* vecnum) const
{
   const int indexes = this->rmm_bigs.size();
   std::vector<int> row;
   for(int i = 0; i < indexes; i++) {
      int bi = this->rmm_bigs[i];
      for(int l=1; l<=fortran_vars.m; l++) {
         for(int k=1; k<=(l-1); k++) {
            int idx = l + (2*fortran_vars.m-k)*(k-1)/2;
            if(bi+1 == idx) {
              row.push_back(l-1);
              row.push_back(k-1);
            }
         }
         int idx = l+(2*fortran_vars.m-l)*(l-1)/2;
         if(bi+1 == idx) {
           row.push_back(l-1);
         }
      }
   }
// delete of repeated indexes
   int numeros[row.size()];
   int r = 0, i, k;
   for(i=0,k=0;i<row.size(),k<row.size();i++,k++) {
      numeros[i]=row[r];r++;
      for(int j =i-1; j>=0;j--) {
        if(numeros[i] == numeros[j])
        {
           i--;break;
        }
     }
  }
  for(int wx=0; wx<i; wx++)
        vecnum[wx] = numeros[wx];
}
*/
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
