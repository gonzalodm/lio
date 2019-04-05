#include <iostream>
#include <omp.h>
#include <stdio.h>

#include "obtain.h"
#include "centros.h"

// borrar luego
using namespace std;
#include <cstdlib>

Obtain::Obtain(int vec, int size, int dimension)
{
   total_vec = vec;
   int_total = dimension;
   M  = size;
   M2 = size * size;
   M3 = size * size * size;
   timerI = timerF = 0.0f;
}

Obtain::~Obtain()
{
   printf(" TIME FOCK CALCULATION %f\n",timerF-timerI);
}

void Obtain::calculate(double* Tmat, FourCenter* Kmat, double* Fmat)
{
   method_A(Tmat,Kmat,Fmat);
}

void Obtain::method_A(double* T, FourCenter* K, double* F)
{
// T is not a symmetric matrix
   double Dens = 0.0f;
   int p1, p2, p3, p4;

   timerI = omp_get_wtime();
#pragma omp parallel for private(p1,p2,p3,p4,Dens)
   for(int ivec=0; ivec<total_vec; ivec++) {
      for(int i=0; i<int_total; i++) {
        p1 = ivec*M2+K[i].k*M+K[i].l;
        p2 = ivec*M2+K[i].l*M+K[i].k;
        p3 = ivec*M2+K[i].u*M+K[i].v;
        p4 = ivec*M2+K[i].v*M+K[i].u;

        Dens = ( T[p1] + T[p2] ) * 2.0f;
        F[p3] += K[i].result * Dens;
        F[p4] += K[i].result * Dens;

        Dens = ( T[p3] + T[p4] ) * 2.0f;
        F[p1] += K[i].result * Dens;
        F[p2] += K[i].result * Dens;
     }
   }
   timerF = omp_get_wtime();
}

// THIS IS FOR OPEN SHELL
void Obtain::calculate(double* TA, double* TB, FourCenter* K, double* FA, double* FB)
{
   double Dens; Dens = 0.0f;
   double value; value = 0.0f;
   int p1, p2, p3, p4;

   timerI = omp_get_wtime();
#pragma omp parallel for private(p1,p2,p3,p4,Dens,value)
   for(int ivec=0; ivec<total_vec; ivec++) {
      for(int i=0; i<int_total; i++) {
        p1 = ivec*M2+K[i].k*M+K[i].l;
        p2 = ivec*M2+K[i].l*M+K[i].k;
        p3 = ivec*M2+K[i].u*M+K[i].v;
        p4 = ivec*M2+K[i].v*M+K[i].u;

        Dens = TA[p1] + TA[p2] + TB[p1] + TB[p2];

        value = K[i].result * Dens;
        FA[p3] += value;
        FA[p4] += value;
        
        Dens = TA[p3] + TA[p4] + TB[p3] + TB[p4];
        value = K[i].result * Dens;
        FA[p1] += value;
        FA[p2] += value;

        Dens = TA[p1] + TA[p2] + TB[p1] + TB[p2];
        value = K[i].result * Dens;
        FB[p3] += value;
        FB[p4] += value;
        
        Dens = TA[p3] + TA[p4] + TB[p3] + TB[p4];
        value = K[i].result * Dens;
        FB[p1] += value;
        FB[p2] += value;
     }
   }
   timerF = omp_get_wtime();
}
