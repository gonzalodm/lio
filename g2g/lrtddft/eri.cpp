#include <iostream>
#include <omp.h>
#include <libint2.hpp>
#include "../init.h"
#include "centros.h"

using libint2::Atom;
using libint2::BasisSet;
using libint2::Shell;
using libint2::Engine;
using libint2::Operator;
using namespace std;
using namespace G2G;

typedef unsigned int uint;

std::vector<Atom> libint_geom(double*,int);
std::vector<libint2::Shell> make_basis(const std::vector<Atom>&,double*,double*,uint*,uint*,
                            int,int,int,int);

size_t max_nprim(const std::vector<libint2::Shell>&);

int max_l(const std::vector<libint2::Shell>&);

std::vector<int> map_shell(const std::vector<libint2::Shell>&);

void eri(int M,uint natoms,uint*ncont,
         double*cbas,double*a,double*r,uint*nuc,
         int sfunc,int pfunc, int dfunc)
{
  int M3, M2;
  M3 = M*M*M; M2 = M*M;
  libint2::initialize();
  Shell::do_enforce_unit_normalization(false);
  std::vector<Atom> atoms = libint_geom(r,natoms);
  auto obs = make_basis(atoms,a,cbas,ncont,nuc,sfunc,pfunc,dfunc,M);

// Check basis normalization with libint
//TODO: ver como hacer para que libint no normalize
/*
  cout << "BASIS SET LIBINT" << endl; // print SET BASIS
  std::copy(begin(obs), end(obs),
         std::ostream_iterator<Shell>(std::cout, "\n"));
*/

  Engine eri_engine(Operator::coulomb, max_nprim(obs), max_l(obs), 0);
  const auto& buf = eri_engine.results();
  std::vector<int> shell2bf = map_shell(obs);

  // Temporal arrays
  int ind_tot = 0;
  std::vector<unsigned int> indu, indv, indk, indl;
  std::vector<double> intvalue4;

  // We use FULL SYMMETRY
  for(int s1=0; s1<obs.size(); ++s1) {
    int bf1_first = shell2bf[s1]; // first basis function in this shell
    int n1 = obs[s1].size();   // number of basis function in this shell

    for(int s2=0; s2<=s1; ++s2) {
      int bf2_first = shell2bf[s2];
      int n2 = obs[s2].size();

      for(int s3=0; s3<=s1; ++s3) {
        int bf3_first = shell2bf[s3];
        int n3 = obs[s3].size();

        int s4_lim = (s1 == s3) ? s2 : s3;
        for(int s4=0; s4<=s4_lim; ++s4) {
          int bf4_first = shell2bf[s4];
          int n4 = obs[s4].size();

          // degeneration factor
          double s12_deg = (s1 == s2) ? 2.0f : 1.0f;
          double s34_deg = (s3 == s4) ? 2.0f : 1.0f;
          double s12_34_deg = (s1 == s3) ? ( (s2 == s4) ? 2.0f : 1.0f ) : 1.0f;
          double s1234_deg = s12_deg * s34_deg * s12_34_deg;

          eri_engine.compute(obs[s1], obs[s2], obs[s3], obs[s4]);
          const auto* buf_1234 = buf[0];
          if (buf_1234 == nullptr) continue;
          for(int f1=0, f1234=0; f1<n1; ++f1){
            const int bf1 = f1 + bf1_first;
            for(int f2=0; f2<n2; ++f2){
              const int bf2 = f2 + bf2_first;
              for(int f3=0; f3<n3; ++f3){
                const int bf3 = f3 + bf3_first;
                for(int f4=0; f4<n4; ++f4, ++f1234){
                  const int bf4 = f4 + bf4_first;
                  double value = buf_1234[f1234];
                  
                  ind_tot += 1;
                  value = value / s1234_deg;
                  intvalue4.push_back(value);
                  indu.push_back(bf1);
                  indv.push_back(bf2);
                  indk.push_back(bf3);
                  indl.push_back(bf4);

                }
              }
            }
          }

        }
      }
    }
  }

  libint2::finalize();
  // Allocate memory
  fortran_vars.dim = ind_tot;
  fortran_vars.Kmat = NULL;
  fortran_vars.Kmat = (FourCenter*) malloc(ind_tot*sizeof(FourCenter));
  if ( fortran_vars.Kmat == NULL ) {
     printf("Can't allocate memory of Kmat in eri\n");
     exit(-1);
  }
  
  for(int i=0; i<ind_tot; i++) {
      fortran_vars.Kmat[i].u = indu[i];
      fortran_vars.Kmat[i].v = indv[i];
      fortran_vars.Kmat[i].k = indk[i];
      fortran_vars.Kmat[i].l = indl[i];
      fortran_vars.Kmat[i].result = intvalue4[i];
  }

  // Free Memory
  std::vector<double>().swap(intvalue4);
  std::vector<unsigned int>().swap(indu);
  std::vector<unsigned int>().swap(indv);
  std::vector<unsigned int>().swap(indk);
  std::vector<unsigned int>().swap(indl);

}

std::vector<Atom> libint_geom(double* r,int natoms)
{
   std::vector<Atom> atoms(natoms);

   for(int i=0; i<natoms; i++) {
     atoms[i].x = r[i];
     atoms[i].y = r[i+natoms];
     atoms[i].z = r[i+natoms*2];
   }
   return atoms;
}

std::vector<libint2::Shell> make_basis(const std::vector<Atom>& atoms,
     double*a,double*c,uint*ncont,uint*nuc,int s_func, int p_func, int d_func,int M)
{
   std::vector<Shell> obs;
   int from = 0;
   int to = s_func;

   for(int i=from; i<to; i++) {  // for s functions
     int centro = nuc[i]-1;
     int tam = ncont[i];
     std::vector<double> exp(tam);
     std::vector<double> coef(tam);
     for(int cont=0; cont<tam; cont++) {
        exp[cont] = a[i+M*cont];
        coef[cont] = c[i+M*cont];
     }
     obs.push_back(
        {
          exp,
          {
            {0, false, coef}
          },
          {{atoms[centro].x,atoms[centro].y,atoms[centro].z}}
        }
     );
   }

   from = s_func;
   to = s_func+p_func*3;
   for(int i=from; i<to; i=i+3) { // for p functions
      int centro = nuc[i]-1;   
      int tam = ncont[i];
      std::vector<double> exp(tam);
      std::vector<double> coef(tam);
      for(int cont=0; cont<tam; cont++) {
         exp[cont] = a[i+M*cont];
         coef[cont] = c[i+M*cont];
      }
      obs.push_back(
         {
           exp,
           {
             {1, false, coef}
           },
           {{atoms[centro].x,atoms[centro].y,atoms[centro].z}}
         }
      );
   }

   from = s_func+p_func*3;
   to = M;
   for(int i=from; i<to; i=i+6) { // for d functions
      int centro = nuc[i]-1;
      int tam = ncont[i];
      std::vector<double> exp(tam);
      std::vector<double> coef(tam);
      for(int cont=0; cont<tam; cont++) {
         exp[cont] = a[i+M*cont];
         coef[cont] = c[i+M*cont];
      }
      obs.push_back(
         {
           exp,
           {
             {2, false, coef}
           },
           {{atoms[centro].x,atoms[centro].y,atoms[centro].z}}
         }
      );
    }

   return obs;
}

size_t max_nprim(const std::vector<libint2::Shell>& obs) {
  size_t n = 0;
  for (auto shell: obs)
    n = std::max(shell.nprim(), n);
  return n;
}

int max_l(const std::vector<libint2::Shell>& obs) {
  int l = 0;
  for (auto shell: obs)
    for (auto c: shell.contr)
      l = std::max(c.l, l);
  return l;
}

std::vector<int> map_shell(const std::vector<libint2::Shell>& obs) {
  std::vector<int> result;
  result.reserve(obs.size());

  int n = 0;
  for (auto shell: obs) {
    result.push_back(n);
    n += shell.size();
  }
  return result;
}
