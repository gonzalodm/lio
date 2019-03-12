#include <iostream>
#include <omp.h>

#include "../common.h"
#include "../init.h"
#include "../partition.h"
#include "calc_gradients.h"

#include <stdio.h>
#include <string.h>


using namespace G2G;

extern Partition partition;

extern "C" void g2g_calcgrdexc_(double* P,double* V, 
                                double* C,double* F, int& NCO)
{
   cout << "en forces c++ " << endl;
   partition.solveForcesExc(P,V,C,F,NCO);
}

namespace G2G {

void Partition::solveForcesExc(double*P,double*V,
                               double*C,double*F,int& NCO)
{
#pragma omp parallel for schedule(static)
    for(uint i=0;i<work.size();i++) {
      for(uint j=0;j<work[i].size();j++) {
         int ind = work[i][j];
         if(ind >= cubes.size()) {
           spheres[ind - cubes.size()]->solve_for_exc(P,V,C,F,NCO);
         } else {
           cubes[ind]->solve_for_exc(P,V,C,F,NCO);
         }
      }
   }
   fflush(stdout);
}

template<class scalar_type> void PointGroupCPU<scalar_type>::
               solve_for_exc(double*P,double*V,
                             double*C,double*F,int& NCO)
{
   const uint group_m = this->total_functions();
   const int npoints = this->points.size();
   bool lda = false;
   bool compute_forces = true;
   compute_functions(compute_forces,!lda);
   int* numeros = new int[group_m];
   int M = fortran_vars.m;
   int natom = fortran_vars.atoms;
   int local_atoms = this->total_nucleii();
 
// Reduced Matrix for this group
   HostMatrix<scalar_type> rmm_input(group_m,group_m);
   HostMatrix<scalar_type> Pred(group_m,group_m);
   HostMatrix<scalar_type> Vred(group_m,group_m);
   HostMatrix<scalar_type> Cred(NCO,group_m);

   double* dens = (double*) malloc(4*sizeof(double));
   double* trad = (double*) malloc(4*sizeof(double));
   double* diff = (double*) malloc(4*sizeof(double));

// Force for this group
   int dim = 3*local_atoms;
   double* sForce = (double*) malloc(dim*sizeof(double));
   memset(sForce,0.0f,dim*sizeof(double));

   get_rmm_input(rmm_input);
   get_coef_input(rmm_input,numeros);

   int row, col;
   for(int i=0; i<group_m; i++) {
     row = numeros[i];
     Pred(i,i) = P[row*M+row];
     Vred(i,i) = V[row*M+row];
     for(int j=0; j<i; j++) {
       col = numeros[j];
       Pred(i,j) = P[row*M+col] + P[col*M+row];
       Pred(j,i) = Pred(i,j);
       Vred(i,j) = V[row*M+col] + V[col*M+row];
       Vred(j,i) = Vred(i,j);
     }
     for(int k=0; k<NCO; k++) {
        Cred(k,i) = C[k*M+row];
     }
   }
   delete[] numeros; numeros = NULL;

/*
   for(int i=0;i<NCO;i++)
     for(int j=0; j<group_m; j++)
        cout << i << " " << j << " " << Cred(i,j) << endl;
*/

   // Gradients for nuclei
   HostMatrix<scalar_type> ddx, ddy, ddz;
   ddx.resize(local_atoms, 1);
   ddy.resize(local_atoms, 1);
   ddz.resize(local_atoms, 1);
   
   for(int point=0;point<npoints;point++) {
    double pd, pdx, pdy, pdz; pd = pdx = pdy = pdz = 0.0f;
    double pp, ppx, ppy, ppz; pp = ppx = ppy = ppz = 0.0f;
    double pt, ptx, pty, ptz; pt = ptx = pty = ptz = 0.0f;
    // functions and gradients values
    const scalar_type* fv = function_values.row(point);
    const scalar_type* gfx = gX.row(point);
    const scalar_type* gfy = gY.row(point);
    const scalar_type* gfz = gZ.row(point);
    // Hessianos de las funciones bases
    const scalar_type* hxx = hPX.row(point); //XX
    const scalar_type* hyy = hPY.row(point); //YY
    const scalar_type* hzz = hPZ.row(point); //ZZ
    const scalar_type* hxy = hIX.row(point); //XY
    const scalar_type* hxz = hIY.row(point); //XZ
    const scalar_type* hyz = hIZ.row(point); //YZ

    for(int i=0;i<group_m;i++) {
       double w3xc, w3yc, w3zc, w; w3xc = w3yc = w3zc = w = 0.0f;
       double z3xc, z3yc, z3zc, z; z3xc = z3yc = z3zc = z = 0.0f;
       double q3xc, q3yc, q3zc, q; q3xc = q3yc = q3zc = q = 0.0f;
       const scalar_type* rm = rmm_input.row(i);
       for(int j=0;j<=i;j++) {
          const scalar_type rmj = rm[j];
          // Ground State Density
          w += fv[j] * rmj;
          w3xc += gfx[j] * rmj;
          w3yc += gfy[j] * rmj;
          w3zc += gfz[j] * rmj;
          // Difference Relaxed Excited State Density
          z += fv[j] * Pred(i,j);
          z3xc += gfx[j] * Pred(i,j);
          z3yc += gfy[j] * Pred(i,j);
          z3zc += gfz[j] * Pred(i,j);
          // Transition Density
          q += fv[j] * Vred(i,j);
          q3xc += gfx[j] * Vred(i,j);
          q3yc += gfy[j] * Vred(i,j);
          q3zc += gfz[j] * Vred(i,j);
       }
       const double Fi = fv[i];
       const double gx = gfx[i], gy = gfy[i], gz = gfz[i];
       // Ground State Density
       pd += Fi * w;
       pdx += gx * w + w3xc * Fi;
       pdy += gy * w + w3yc * Fi;
       pdz += gz * w + w3zc * Fi;
       // Difference Relaxed Excited State Density
       pp += Fi * z;
       ppx += gx * z + z3xc * Fi;
       ppy += gy * z + z3yc * Fi;
       ppz += gz * z + z3zc * Fi;
       // Transition Density
       pt += Fi * q;
       ptx += gx * q + q3xc * Fi;
       pty += gy * q + q3yc * Fi;
       ptz += gz * q + q3zc * Fi;
    }
    double sigma = (pdx * pdx) + (pdy * pdy) + (pdz * pdz);
    pdx *= 0.5f; pdy *= 0.5f; pdz *= 0.5f;
    dens[0] = pd; dens[1] = pdx; dens[2] = pdy; dens[3] = pdz;
    diff[0] = pp; diff[1] = ppx; diff[2] = ppy; diff[3] = ppz;
    trad[0] = pt; trad[1] = ptx; trad[2] = pty; trad[3] = ptz;

    // esto sobreescribe dens, diff, trad
    calc_gradients(dens,diff,trad,sigma);

    // FORCES CALCULATE
    double DJII, PJII, VJII;
    double temp[4];
    // Gradients for the nuclei
    ddx.zero();
    ddy.zero();
    ddz.zero();

    const scalar_type wp = this->points[point].weight;

    for (int i = 0, ii = 0; i < this->total_functions_simple(); i++) {// cantidad de bases, s+p+d (sin tener en cuenta las 3p y 5d
       uint nuc = this->func2local_nuc(ii);// da a que nucleo LOCAL pertenece la base ii
       uint inc_i = this->small_function_type(i);// da cuantas funciones tiene ese tipo: s->1, p->3, d->5
       double grdx, grdy, grdz; grdx = grdy = grdz = 0.0f;

       for (uint k = 0; k < inc_i; k++, ii++) {
          for (uint j = 0; j < group_m; j++) {
             double factor = (ii == j ? 2.0f : 1.0f);
             DJII = rmm_input(j,ii) * factor * 0.5f;
             PJII = Pred(j,ii) * factor;
             VJII = Vred(j,ii) * factor;
             //cout << j << " " << ii << " " << PJII << endl;
             temp[0] = 2.0f * (dens[0]*DJII + diff[0]*PJII + trad[0]*VJII);
             temp[1] = 2.0f * (dens[1]*DJII + diff[1]*PJII + trad[1]*VJII);
             temp[2] = 2.0f * (dens[2]*DJII + diff[2]*PJII + trad[2]*VJII);
             temp[3] = 2.0f * (dens[3]*DJII + diff[3]*PJII + trad[3]*VJII);

             grdx += temp[0]*gfx[ii]*fv[j];
             grdx += temp[1]*( hxx[ii]*fv[j] + gfx[ii]*gfx[j] );
             grdx += temp[2]*( hxy[ii]*fv[j] + gfx[ii]*gfy[j] );
             grdx += temp[3]*( hxz[ii]*fv[j] + gfx[ii]*gfz[j] );

             grdy += temp[0]*gfy[ii]*fv[j];
             grdy += temp[1]*( hxy[ii]*fv[j] + gfy[ii]*gfx[j] );
             grdy += temp[2]*( hyy[ii]*fv[j] + gfy[ii]*gfy[j] );
             grdy += temp[3]*( hyz[ii]*fv[j] + gfy[ii]*gfz[j] );

             grdz += temp[0]*gfz[ii]*fv[j];
             grdz += temp[1]*( hxz[ii]*fv[j] + gfz[ii]*gfx[j] );
             grdz += temp[2]*( hyz[ii]*fv[j] + gfz[ii]*gfy[j] );
             grdz += temp[3]*( hzz[ii]*fv[j] + gfz[ii]*gfz[j] );
          }
       }
       ddx(nuc) += grdx;
       ddy(nuc) += grdy;
       ddz(nuc) += grdz;
    }
    for (int i = 0; i < local_atoms; i++) {
       sForce[i*3] -= ddx(i) * wp;
       sForce[i*3+1] -= ddy(i) * wp;
       sForce[i*3+2] -= ddz(i) * wp;
    }
   } // END POINTS
  
   // Acumulate forces for this group
   for (int i = 0; i < local_atoms; i++) {
      uint global_atom = this->local2global_nuc[i]; // da el indice del atomo LOCAL al GLOBAL
      F[global_atom*3] += sForce[i*3];
      F[global_atom*3 + 1] += sForce[i*3+1];
      F[global_atom*3 + 2] += sForce[i*3+2];
      //cout << "at x y z " << global_atom << " " << F[global_atom*3] << " " << F[global_atom*3+1] << " " << F[global_atom*3+2] <<endl;
   }

   // Free Memory
   free(dens); dens = NULL;
   free(diff); diff = NULL;
   free(trad); trad = NULL;
   free(sForce); sForce = NULL;
}

#if FULL_DOUBLE
template class PointGroup<double>;
template class PointGroupCPU<double>;
#else
template class PointGroup<float>;
template class PointGroupCPU<float>;
#endif
}


