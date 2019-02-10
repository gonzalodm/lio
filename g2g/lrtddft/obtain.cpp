#include <iostream>
#include <omp.h>
#include <stdio.h>

#include "obtain.h"

Obtain::Obtain(int vec, int size)
{
   total_vec = vec;
   M  = size;
   M2 = size * size;
   M3 = size * size * size;
   timerI = timerF = 0.0f;
}

Obtain::~Obtain()
{
   printf(" TIME FOCK CALCULATION %f\n",timerF-timerI);
}

void Obtain::calculate(double* Tmat, double* Kmat, double* Fmat)
{
   method_A(Tmat,Kmat,Fmat);
}

void Obtain::method_A(double* T, double* K, double* F)
{
// T is not a symmetric matrix

   double valor, Dens;
   valor = Dens = 0.0f;
   
   timerI = omp_get_wtime();
#pragma omp parallel for private(Dens,valor)
   for(int ivec=0; ivec<total_vec; ivec++) {
     for(int u=0; u<M; u++) {
       for(int v=0; v<=u; v++) {
         for(int k=0; k<M; k++) {
           for(int l=0;l<k; l++) {
              Dens = ( T[ivec*M2+k*M+l] + T[ivec*M2+l*M+k] )*2.0f;
              valor += Dens * K[u*M3+v*M2+k*M+l];
           }
           valor += T[ivec*M2+k*M+k] * 2.0f * K[u*M3+v*M2+k*M+k];
         }
         F[ivec*M2+u*M+v] = valor;
         F[ivec*M2+v*M+u] = valor;
         valor = 0.0f;
       }
     }
   }
   timerF = omp_get_wtime();
}
