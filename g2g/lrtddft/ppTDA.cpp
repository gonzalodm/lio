#include <iostream>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <libint2.hpp>

#include "../common.h"
#include "../init.h"

using namespace G2G;

using libint2::Atom;
using libint2::BasisSet;
using libint2::Shell;
using libint2::Engine;
using libint2::Operator;

typedef unsigned int uint;

std::vector<Atom> libint_geom_ppTDA(double*,int);
std::vector<libint2::Shell> make_basis_ppTDA(const std::vector<Atom>&,double*,double*,uint*,uint*,
                            int,int,int,int);

size_t max_nprim_ppTDA(const std::vector<libint2::Shell>&);

int max_l_ppTDA(const std::vector<libint2::Shell>&);

std::vector<int> map_shell_ppTDA(const std::vector<libint2::Shell>&);

void doIntegrals(int&,uint,uint*, double*,double*,double*,uint*,int,int,int,int&,double*, double*);

extern "C" void g2g_calcpptda_(double* T,double* Cbas,int& numvec,double* F)
{
   int M = fortran_vars.m;
   int s_func = fortran_vars.s_funcs;
   int p_func = fortran_vars.p_funcs;
   int d_func = fortran_vars.d_funcs;
   double* aContr = &fortran_vars.a_values(0,0);
   double* pos = &fortran_vars.atom_positions_pointer(0,0);
   uint* ncont = &fortran_vars.contractions(0);
   uint* nuc = &fortran_vars.nucleii(0);

   int timeI = omp_get_wtime();
   doIntegrals(M,fortran_vars.atoms,ncont,Cbas,aContr,pos,nuc,
       s_func,p_func,d_func,numvec,T,F);
   int timeF = omp_get_wtime();
}

void doIntegrals(int& M,uint natoms,uint*ncont, double*cbas,
     double*a,double*r,uint*nuc,int sfunc,int pfunc, int dfunc,
     int& totvec,double* T, double* F) // outputs
{
  int M4, M3, M2; 
  M2 = M*M;
  M3 = M2*M;
  M4 = M3*M;

  libint2::initialize();
  Shell::do_enforce_unit_normalization(false);
  std::vector<Atom> atoms = libint_geom_ppTDA(r,natoms);
  auto obs = make_basis_ppTDA(atoms,a,cbas,ncont,nuc,sfunc,pfunc,dfunc,M);

  double* Kmat = (double*)malloc(M4*sizeof(double));

// Check basis normalization with libint
//TODO: ver como hacer para que libint no normalize
  std::cout << "BASIS SET LIBINT" << std::endl; // print SET BASIS
  std::copy(begin(obs), end(obs),
         std::ostream_iterator<Shell>(std::cout, "\n"));

  Engine eri_engine(Operator::coulomb, max_nprim_ppTDA(obs), max_l_ppTDA(obs), 0);
  const auto& buf = eri_engine.results();
  std::vector<int> shell2bf = map_shell_ppTDA(obs);

  // We use FULL SYMMETRY
  for(int s1=0; s1<obs.size(); ++s1) {
    int bf1_first = shell2bf[s1]; // first basis function in this shell
    int n1 = obs[s1].size();   // number of basis function in this shell

    for(int s2=0; s2<obs.size(); ++s2) {
      int bf2_first = shell2bf[s2];
      int n2 = obs[s2].size();

      for(int s3=0; s3<obs.size(); ++s3) {
        int bf3_first = shell2bf[s3];
        int n3 = obs[s3].size();

        for(int s4=0; s4<obs.size(); ++s4) {
          int bf4_first = shell2bf[s4];
          int n4 = obs[s4].size();

          // this solve (ac|bd)
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
                  const int ptr = bf1*M3+bf2*M2+bf3*M+bf4;
                  Kmat[ptr] = value;
                  //std::cout << bf1 << " " << bf2 << " " << bf3;
                  //std::cout << " " << bf4 << "    " << Kmat[ptr] << std::endl;
                }
              }
            }
          }

        }
      }
    }
  }

  double temp = 0.0f;
 
  for(int ivec=0; ivec<totvec; ivec++) {
    for(int u=0; u<M; u++) {
      for(int l=0; l<M; l++) {
        const int pF = ivec*M2+u*M+l;
        for(int v=0; v<M; v++) {
          for(int s=0; s<M; s++) {
             const int pT1 = ivec*M2+v*M+s;
             const int pT2 = ivec*M2+s*M+v;
             const int pK1 = u*M3+v*M2+l*M+s;
             const int pK2 = u*M3+s*M2+l*M+v;
             temp += ( Kmat[pK1] - Kmat[pK2] )*( T[pT1]+T[pT2] )*2.0f;
             //std::cout << Kmat[pK1] << "  " << Kmat[pK2] << std::endl;
             //std::cout << u << " " << l << " " << v << " " << s << std::endl;
             //std::cout << T[pT1] + T[pT2] << std::endl;
             //std::cout << Kmat[pK1] - Kmat[pK2] << std::endl;
          }
        }
        F[pF] += temp;
        temp = 0.0f;
      }
    }
  }

  libint2::finalize();

// Free Memory
  std::vector<Atom>().swap(atoms);
  std::vector<int>().swap(shell2bf);
  std::vector<Shell>().swap(obs);
  free(Kmat); Kmat=NULL;
}

std::vector<Atom> libint_geom_ppTDA(double* r,int natoms)
{
   std::vector<Atom> atoms(natoms);

   for(int i=0; i<natoms; i++) {
     atoms[i].x = r[i];
     atoms[i].y = r[i+natoms];
     atoms[i].z = r[i+natoms*2];
   }
   return atoms;
}

std::vector<libint2::Shell> make_basis_ppTDA(const std::vector<Atom>& atoms,
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
     std::vector<double>().swap(exp);
     std::vector<double>().swap(coef);
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
      std::vector<double>().swap(exp);
      std::vector<double>().swap(coef);
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
      std::vector<double>().swap(exp);
      std::vector<double>().swap(coef);
    }

   return obs;
}

size_t max_nprim_ppTDA(const std::vector<libint2::Shell>& obs) {
  size_t n = 0;
  for (auto shell: obs)
    n = std::max(shell.nprim(), n);
  return n;
}

int max_l_ppTDA(const std::vector<libint2::Shell>& obs) {
  int l = 0;
  for (auto shell: obs)
    for (auto c: shell.contr)
      l = std::max(c.l, l);
  return l;
}

std::vector<int> map_shell_ppTDA(const std::vector<libint2::Shell>& obs) {
  std::vector<int> result;
  result.reserve(obs.size());

  int n = 0;
  for (auto shell: obs) {
    result.push_back(n);
    n += shell.size();
  }
  return result;
}
