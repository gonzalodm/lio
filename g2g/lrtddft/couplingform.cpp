#include <iostream>
#include <cstdlib>
#include <math.h>
#include <omp.h>

#define DENSMIN 1e-10

using namespace std;

void ObtainFock_A(double* T, double* K, double* Fock, int vec, int M)
{
// In general, T is not a symmetry matrix
   int M2 = M * M;
   int M3 = M2 * M;
   double valor, Dens;
   valor = 0.0f;

#pragma omp parallel for private(Dens,valor)
   for(int ivec=0; ivec<vec; ivec++) {
     for(int u=0; u<M; u++) {
       for(int v=0; v<=u; v++) {
         for(int k=0; k<M; k++) {
           for(int l=0;l<k; l++) {
              Dens = ( T[ivec*M2+k*M+l] + T[ivec*M2+l*M+k] )*2.0f;
              valor += Dens * K[u*M3+v*M2+k*M+l];
           }
           valor += T[ivec*M2+k*M+k] * 2.0f * K[u*M3+v*M2+k*M+k];
         }
         Fock[ivec*M2+u*M+v] = valor;
         Fock[ivec*M2+v*M+u] = valor;
         valor = 0.0f;
       }
     }
   }
}
