#include <iostream>
#include <xc.h>
#include <omp.h>

#include <stdio.h>
#include <string.h> 

#include "../common.h"
#include "../init.h"
#include "../partition.h"

#include "eri.h"
#include "couplingform.h"
#include "../libxc/libxcproxy.h"

#define DENSMIN 1e-5

using namespace G2G;
extern Partition partition;

//######################################################################
//######################################################################
extern "C" void g2g_calculate2e_(double* Tmat,double* K2eAO,double* Cbas,
                                 int& numvec,double* F,int& int2elec)
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

   timeI = timeF = 0.0;
   if (int2elec == 0 ) {
     timeI = omp_get_wtime();
     eri(K2eAO,M,fortran_vars.atoms,ncont,Cbas,aContr,pos,nuc,
         s_func,p_func,d_func);
     timeF = omp_get_wtime();
     printf("POR PASAS UNA VEZ\n");
     //printf("ERI SUBROUTINE %f\n",timeF-timeI);
   }

   double valor, Dens;
   valor = 0.0f;
   for(int ivec=0; ivec<numvec; ivec++) {
     for(int u=0; u<M; u++) {
       for(int v=0; v<=u; v++) {
         for(int k=0; k<M; k++) {
           for(int l=0;l<k; l++) {
              Dens = (Tmat[ivec*M2+k*M+l] + Tmat[ivec*M2+l*M+k])*2.0f; // Tmat es una matriz NO SIMETRICA
              valor += Dens * K2eAO[u*M3+v*M2+k*M+l];
           }
           valor += Tmat[ivec*M2+k*M+k]*2.0f*K2eAO[u*M3+v*M2+k*M+k];
         }
         F[ivec*M2+u*M+v] = valor;
         F[ivec*M2+v*M+u] = valor;
         valor = 0.0f;
       }
     }
   }

   fflush(stdout); // NOT BUFFERED
}
//######################################################################
//######################################################################
extern "C" void g2g_calculatedft_(double* Tmat,double* C,double* Fv, int& NCO)
{
   partition.solve_lr(Tmat,C,Fv,NCO);
   fflush(stdout); // NOT BUFFERED
}

//######################################################################
//######################################################################

namespace G2G {

void Partition::solve_lr(double* T,double* C,double* F,int& NCO)
{

   double timeI, timeF;
   timeI = omp_get_wtime();
#pragma omp parallel for schedule(static)
   for(uint i=0;i<work.size();i++) {
      for(uint j=0;j<work[i].size();j++) {
         int ind = work[i][j];
         if(ind >= cubes.size()) {
           spheres[ind-cubes.size()]->solve_closed_lr(T,C,F,NCO);
         } else {
           cubes[ind]->solve_closed_lr(T,C,F,NCO);
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
               solve_closed_lr(double* T,double* C,double* F,int& NCO)
{
// aqui T no es la transicion density, sino es B * C**T
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

// Formamos B * C**T para este grupo
   HostMatrix<scalar_type> tred(NCO,group_m);
   HostMatrix<scalar_type> Cred(NCO,group_m);

   double* grdMO  = (double*)malloc(4*NCO*sizeof(double));
   double* smallFock  = (double*)malloc(group_m*group_m*sizeof(double));
   memset(smallFock,0.0f,group_m*group_m*sizeof(double));

// Formamos matrices reducidas para el grupo
   for(int i=0; i<NCO; i++) {
     for(int j=0; j<group_m; j++) {
       int jj = numeros[j];
       tred(i,j) = T[i*M+jj]; // con esto utilizo el triangulo inferior de tred
       Cred(i,j) = C[i*M+jj];
     }
   } 

//LIBXC INITIALIZATION
   const int nspin = XC_POLARIZED;
   const int functionalExchange = fortran_vars.ex_functional_id; //101;
   const int functionalCorrelation = fortran_vars.ec_functional_id; // 130;
   double* lrCoef = new double[3];
   double* precond = new double[group_m];
   LibxcProxy<scalar_type,3> libxcProxy(functionalExchange, functionalCorrelation, nspin);
   
   for(int point=0;point<npoints;point++) {
      scalar_type pd, fxc, tdx, tdy, tdz; pd = fxc = tdx = tdy = tdz = 0.0f;
      const scalar_type* fv = function_values.row(point);
      const scalar_type* gxv = gX.row(point);
      const scalar_type* gyv = gY.row(point);
      const scalar_type* gzv = gZ.row(point);
      for(int i=0;i<group_m;i++) {
         double w3xc, w3yc, w3zc, w; w3xc = w3yc = w3zc = w = 0.0f;
         const scalar_type* rm = rmm_input.row(i);
         for(int j=0;j<=i;j++) {
            const double rmj = rm[j];
            w += fv[j] * rmj;
            w3xc += gxv[j] * rmj;
            w3yc += gyv[j] * rmj;
            w3zc += gzv[j] * rmj;
         }
         const double Fi = fv[i];
         const double gx = gxv[i], gy = gyv[i], gz = gzv[i];
         pd += Fi * w;
         tdx += gx * w + w3xc * Fi;
         tdy += gy * w + w3yc * Fi;
         tdz += gz * w + w3zc * Fi;
      }

      memset(grdMO,0.0f,NCO*4*sizeof(double));
      for(int i=0; i<group_m; i++) {
        int row = numeros[i];
        for(int j=0; j<NCO; j++) {
           grdMO[j*4] += fv[i] * Cred(j,i);
           grdMO[j*4+1] += gxv[i] * Cred(j,i);
           grdMO[j*4+2] += gyv[i] * Cred(j,i);
           grdMO[j*4+3] += gzv[i] * Cred(j,i);
        }
      }

      double D0, DX, DY, DZ; D0 = DX = DY = DZ =0.0f;
      double red, redx, redy, redz; red = redx = redy = redz = 0.0f;
      for(int i=0; i<group_m; i++) {
         double Fi = fv[i];
         double gxi = gxv[i];
         double gyi = gyv[i];
         double gzi = gzv[i];
         D0=DX=DY=DZ=0.0f;
         for(int j=0; j<NCO; j++) {
            D0 += tred(j,i) * grdMO[j*4];
            DX += tred(j,i) * grdMO[j*4+1];
            DY += tred(j,i) * grdMO[j*4+2];
            DZ += tred(j,i) * grdMO[j*4+3];
         }
         red += D0 * Fi;
         redx += D0 * gxi + DX * Fi;
         redy += D0 * gyi + DY * Fi;
         redz += D0 * gzi + DZ * Fi;
      }

      double sigma = tdx * tdx + tdy * tdy + tdz * tdz;
      double cruz = redx * tdx + redy * tdy + redz * tdz;
      cruz *= 0.50f;tdx *= 0.5f;tdy *= 0.5f;tdz *= 0.5f;

// LIBXC RUN FOR LINEAR RESPONSE
      libxcProxy.coefLR(&pd,&sigma,red,cruz,lrCoef);
      const scalar_type wp = this->points[point].weight;
      double term1, term2, term3, term4, precondii, result;

      for(int i=0; i<group_m; i++) {
        term1 = lrCoef[0] * 0.5f * fv[i] + lrCoef[1] * tdx * gxv[i];
        term2 = lrCoef[1] * tdy * gyv[i] + lrCoef[1] * tdz * gzv[i];
        term3 = lrCoef[2] * redx * gxv[i] + lrCoef[2] * redy * gyv[i];
        term4 = lrCoef[2] * redz * gzv[i];
        precond[i] = (term1 + term2 + term3 + term4) * wp;
        precondii = precond[i];
        for(int j=0; j<=i; j++) {
           result = fv[i] * precond[j] + precondii * fv[j];
           smallFock[i*group_m+j] += result;
        }
      }

   }  // END points loop

   // ARMAMOS LA FOCK GLOBAL CON LA DEL GRUPO DE PUNTOS
   int row, col;
   for(int i=0; i<group_m; i++) {
     row = numeros[i];
     for(int j=0; j<=i; j++) {
       col = numeros[j];
       F[row*M+col] += smallFock[i*group_m+j];
     }
   }

   // Free Memory
   free(grdMO); grdMO = NULL;
   free(smallFock); smallFock = NULL;
   delete[] numeros; numeros = NULL;
   delete[] precond; precond = NULL;
   delete[] lrCoef; lrCoef = NULL;
}
//######################################################################
//######################################################################

//######################################################################
//######################################################################
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
