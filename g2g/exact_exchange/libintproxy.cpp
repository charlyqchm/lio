// TODO: hacer que cuando calcule gradientes no tenga que volver a calcular
//       obs atoms shell2bf y shell2atom. Hay que guardarlos en un Struct o namespace
#include <iostream>
#include "libintproxy.h"

using namespace G2G;
using namespace std;

namespace libint2 {
  int nthreads;
}

LIBINTproxy::LIBINTproxy(int M,uint natoms,uint*ncont,
               double*cbas,double*a,double*r,uint*nuc,
               int sfunc,int pfunc,int dfunc,int DER)
{
// LIBINT initialization
  cout << " LIBINT initialization" << endl;
  libint2::initialize();
  Shell::do_enforce_unit_normalization(false);

// Put coordinates in a object
  atoms = libint_geom(r,natoms);

// Basis in LIBINT format
  obs = make_basis(atoms,a,cbas,ncont,nuc,sfunc,pfunc,dfunc,M);

// First basis in a shell and atom centre in the shell
  shell2bf = map_shell(); // this return shell2atom too

// Counter Integrals only if no derivative is required
  if(DER == 0)
    counter_integrals();
}

LIBINTproxy::~LIBINTproxy()
{
  // LIBINT DEINITIALIZATION AND FREE MEMORIES
  libint2::finalize();
  vector<int>().swap(shell2bf);
  vector<Shell>().swap(obs); 
  vector<Atom>().swap(atoms);
  vector<int>().swap(shell2atom); 
}

void LIBINTproxy::PrintBasis()
{
// Check basis normalization with libint
  std::cout << "BASIS SET LIBINT" << std::endl; // print SET BASIS
  std::copy(begin(obs), end(obs),
         std::ostream_iterator<Shell>(std::cout, "\n"));
}

void LIBINTproxy::counter_integrals()
{
  Engine eri_engine(Operator::coulomb, max_nprim(), max_l(), 0);
  const auto& buf = eri_engine.results();

  int ind_tot = 0;

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

                  ind_tot += 1;

                }
              }
            }
          }

        }
      }
    }
  }

  // Allocate memory of Integrals
  double mem = sizeof(ind_tot) * 1.25e-10;
  printf(" Memory of integrals: %f Gb\n",mem);
  fortran_vars.dim = ind_tot;
  fortran_vars.Kmat = NULL;
  fortran_vars.Kmat = (FourCenter*) malloc(ind_tot*sizeof(FourCenter));
  if ( fortran_vars.Kmat == NULL ) {
     printf(" Can't allocate memory of Kmat in eri\n");
     exit(-1);
  }
}

void LIBINTproxy::calculate()
{
  Engine eri_engine(Operator::coulomb, max_nprim(), max_l(), 0);
  const auto& buf = eri_engine.results();

  int ii = 0;

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

                  // Save indexes and values in a Object
                  value = value / s1234_deg;
                  fortran_vars.Kmat[ii].u = bf1;
                  fortran_vars.Kmat[ii].v = bf2;
                  fortran_vars.Kmat[ii].k = bf3;
                  fortran_vars.Kmat[ii].l = bf4;
                  fortran_vars.Kmat[ii].result = value;
                  ii += 1;
                }
              }
            }
          }

        }
      }
    }
  }

}

vector<Matrix_E> LIBINTproxy::compute_deriv(const Matrix_E& D, int M)
{
  using libint2::nthreads;

#pragma omp parallel
  nthreads = omp_get_num_threads();

  // nderiv = 3 * N_atom, 3 = x y z
  int nderiv = libint2::num_geometrical_derivatives(atoms.size(),1);
  vector<Matrix_E> W(nthreads*nderiv,Matrix_E::Zero(M,M));
  double precision = numeric_limits<double>::epsilon();

  // Set precision values
  vector<Engine> engines(nthreads);
  engines[0] = Engine(Operator::coulomb, max_nprim(), max_l(), 1);
  engines[0].set_precision(precision);
  for(int i=1; i<nthreads; i++)
    engines[i] = engines[0];

  // Run LIBINT in parallel
  auto lambda = [&](int thread_id) {
     auto& engine = engines[thread_id];
     const auto& buf = engine.results();
     int shell_atoms[4];

     // loop over shells
     for(int s1=0, s1234=0; s1<obs.size(); ++s1) {
       int bf1_first = shell2bf[s1];
       int n1 = obs[s1].size();
       shell_atoms[0] = shell2atom[s1];

       for(int s2=0; s2<=s1; ++s2) {
         int bf2_first = shell2bf[s2];
         int n2 = obs[s2].size();
         shell_atoms[1] = shell2atom[s2];

         for(int s3=0; s3<=s1; ++s3) {
           int bf3_first = shell2bf[s3];
           int n3 = obs[s3].size();
           shell_atoms[2] = shell2atom[s3];
           
           int s4_max = (s1 == s3) ? s2 : s3;
           for(int s4=0; s4<=s4_max; ++s4) {
              if( (s1234++) % nthreads != thread_id) continue;
             
              int bf4_first = shell2bf[s4];
              int n4 = obs[s4].size();
              shell_atoms[3] = shell2atom[s4];

              int s12_deg = (s1 == s2) ? 1.0 : 2.0;
              int s34_deg = (s3 == s4) ? 1.0 : 2.0;
              int s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
              int s1234_deg = s12_deg * s34_deg * s12_34_deg;

              auto add_shellset_to_dest = [&] (int op, int idx, int coord1, int coord2) {
                 auto& g = W[op];
                 auto shset = buf[idx];
                 const int weight = s1234_deg;

                 for (int f1=0,f1234=0; f1<n1; ++f1) {
                   const int bf1 = f1 + bf1_first;
                   for (int f2=0; f2<n2; ++f2) {
                     const int bf2 = f2 + bf2_first;
                     for (int f3=0; f3<n3; ++f3) {
                       const int bf3 = f3 + bf3_first;
                       for (int f4=0; f4<n4; ++f4,++f1234) {
                         const int bf4 = f4 + bf4_first;

                         const double value = shset[f1234];
                         const double wvalue = value * weight;

                         // COULOMB DER
                         //g(bf1, bf2) += D(bf3, bf4) * wvalue;
                         //g(bf3, bf4) += D(bf1, bf2) * wvalue;
                         // EXCHANGE DER
                         g(bf1, bf3) -= 0.25f * D(bf2, bf4) * wvalue;
                         g(bf2, bf4) -= 0.25f * D(bf1, bf3) * wvalue;
                         g(bf1, bf4) -= 0.25f * D(bf2, bf3) * wvalue;
                         g(bf2, bf3) -= 0.25f * D(bf1, bf4) * wvalue;
                       }
                     }
                   }
                 }
              }; // END add_shellset_to_dest

              engine.compute2<Operator::coulomb, BraKet::xx_xx,1>
                             (obs[s1],obs[s2],obs[s3],obs[s4]);
              if (buf[0] == nullptr)
                  continue; // if all integrals screened out, skip to next quartet

              for(int d=0; d<12; ++d) {
                 const int a = d / 3;
                 const int xyz = d % 3;
            
                 int coord = shell_atoms[a] * 3 + xyz;
                 auto& g = W[thread_id * nderiv + coord];
                 int coord1 = 0, coord2 = 0;
  
                 add_shellset_to_dest(thread_id*nderiv+coord,d,coord1,coord2);
              } // END for d

           } // end s4
         } // end s3
       } // end s2
     } //end s1
  
  }; // END OF LAMBDA
  
  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for(int t=1; t<nthreads; ++t) {
    for(int d=0; d<nderiv; ++d) {
       W[d] += W[t * nderiv + d];
    }
  }

  vector<Matrix_E> WW(nderiv);
  for(int d=0; d<nderiv; ++d) {
     WW[d] = 0.5f * ( W[d] + W[d].transpose() );
  }
  vector<Matrix_E>().swap(W);
  vector<Engine>().swap(engines);
  return WW;
}

void LIBINTproxy::derivative(double* Dens, double* For)
{
  int M = fortran_vars.m;
// Pass Pointer to Eigen Matrix
  Matrix_E DE = Matrix_E::Zero(M,M);
  for(int ii=0;ii<M;ii++) {
    for(int jj=0;jj<=ii;jj++) {
      DE(ii,jj) = Dens[ii*M+jj];
      DE(jj,ii) = DE(ii,jj);
    }
  }

// Calculate Matrix G of derivatives
  vector<Matrix_E> G = compute_deriv(DE,M);
  int dim_atom = atoms.size();
  double force;

  for(int atom=0,ii=0; atom<dim_atom; atom++) {
    for(int xyz=0; xyz<3; xyz++,ii++) {
       force = G[ii].cwiseProduct(DE).sum();
       For[xyz*dim_atom+atom] += 0.25f*force;
    }
  }
  for(int atom=0;atom<dim_atom;atom++)
     for(int xyz=0;xyz<3;xyz++)
       cout << atom << " " << xyz << " " << For[xyz*dim_atom+atom] << endl;

  DE.resize(0,0);
  vector<Matrix_E>().swap(G);
}

vector<Matrix_E> LIBINTproxy::compute_deriv(const Matrix_E& D,const Matrix_E& P,
                                            const Matrix_E& T,int M)
{
  cout << "NEW deriv" << endl;
  using libint2::nthreads;
#pragma omp parallel
  nthreads = omp_get_num_threads();

  // nderiv = 3 * N_atom, 3 = x y z
  int nderiv = libint2::num_geometrical_derivatives(atoms.size(),1);
  vector<Matrix_E> W(nthreads*nderiv,Matrix_E::Zero(M,M));
  //vector<Matrix_E> W(nthreads*nderiv);
  double precision = numeric_limits<double>::epsilon();

  // Set precision values
  vector<Engine> engines(nthreads);
  engines[0] = Engine(Operator::coulomb, max_nprim(), max_l(), 1);
  engines[0].set_precision(precision);
  for(int i=1; i<nthreads; i++)
    engines[i] = engines[0];

  // Run LIBINT in parallel
  auto lambda = [&](int thread_id) {
     auto& engine = engines[thread_id];
     const auto& buf = engine.results();
     int shell_atoms[4];

     // loop over shells
     for(int s1=0, s1234=0; s1<obs.size(); ++s1) {
       int bf1_first = shell2bf[s1];
       int n1 = obs[s1].size();
       shell_atoms[0] = shell2atom[s1];

       for(int s2=0; s2<obs.size(); ++s2) {
         int bf2_first = shell2bf[s2];
         int n2 = obs[s2].size();
         shell_atoms[1] = shell2atom[s2];

         for(int s3=0; s3<obs.size(); ++s3) {
           int bf3_first = shell2bf[s3];
           int n3 = obs[s3].size();
           shell_atoms[2] = shell2atom[s3];

           int s4_max = (s1 == s3) ? s2 : s3;
           for(int s4=0; s4<obs.size(); ++s4) {
              if( (s1234++) % nthreads != thread_id) continue;

              int bf4_first = shell2bf[s4];
              int n4 = obs[s4].size();
              shell_atoms[3] = shell2atom[s4];

/*
              int s12_deg = (s1 == s2) ? 1.0 : 2.0;
              int s34_deg = (s3 == s4) ? 1.0 : 2.0;
              int s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
              int s1234_deg = s12_deg * s34_deg * s12_34_deg;
*/

              int s1234_deg = 1;
              auto add_shellset_to_dest = [&] (int op, int idx, int coord1, int coord2) {
                 auto& g = W[op];
                 auto shset = buf[idx];
                 const int weight = s1234_deg;

                 for (int f1=0,f1234=0; f1<n1; ++f1) {
                   const int bf1 = f1 + bf1_first;
                   for (int f2=0; f2<n2; ++f2) {
                     const int bf2 = f2 + bf2_first;
                     for (int f3=0; f3<n3; ++f3) {
                       const int bf3 = f3 + bf3_first;
                       for (int f4=0; f4<n4; ++f4,++f1234) {
                         const int bf4 = f4 + bf4_first;

                         const double value = shset[f1234];
                         const double wvalue = value * weight;
                         double matrix = D(bf1,bf3) * D(bf2, bf4);
                         matrix += 0.5f * P(bf1,bf4) * D(bf3,bf2);
                         matrix += 0.5f * P(bf3,bf1) * D(bf4,bf2);
                         matrix += 2.0f * T(bf3,bf1) * T(bf4,bf2);
                         //cout << bf1 << "  " << bf3 << "  " << T(bf1,bf3) << endl;

                         // COULOMB DER
                         //g(bf1, bf2) += D(bf3, bf4) * wvalue;
                         //g(bf3, bf4) += D(bf1, bf2) * wvalue;
                         // EXCHANGE DER
                         //g(bf1, bf3) -= D(bf2, bf4) * wvalue;
                         g(0,0) -= matrix * wvalue;
                         //g(bf2, bf4) -= 0.25f * D(bf1, bf3) * wvalue;
                         //g(bf1, bf4) -= 0.25f * D(bf2, bf3) * wvalue;
                         //g(bf2, bf3) -= 0.25f * D(bf1, bf4) * wvalue;
                       }
                     }
                   }
                 }
              }; // END add_shellset_to_dest

              engine.compute2<Operator::coulomb, BraKet::xx_xx,1>
                             (obs[s1],obs[s2],obs[s3],obs[s4]);
              if (buf[0] == nullptr)
                  continue; // if all integrals screened out, skip to next quartet

              for(int d=0; d<12; ++d) {
                 const int a = d / 3;
                 const int xyz = d % 3;

                 int coord = shell_atoms[a] * 3 + xyz;
                 auto& g = W[thread_id * nderiv + coord];
                 int coord1 = 0, coord2 = 0;

                 add_shellset_to_dest(thread_id*nderiv+coord,d,coord1,coord2);
              } // END for d

           } // end s4
         } // end s3
       } // end s2
     } //end s1

  }; // END OF LAMBDA

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for(int t=1; t<nthreads; ++t) {
    for(int d=0; d<nderiv; ++d) {
       W[d] += W[t * nderiv + d];
    }
  }

  vector<Matrix_E> WW(nderiv);
  for(int d=0; d<nderiv; ++d) {
     WW[d] = 0.5f * ( W[d] + W[d].transpose() );
  }
  vector<Matrix_E>().swap(W);
  vector<Engine>().swap(engines);
  return WW;
}

void LIBINTproxy::derivative(double* D,double* Diff,
                             double* X, double* For)
{
  int M = fortran_vars.m;
// Pass Pointer to Eigen Matrix
  // Ground state Density: DG
  // Difference relaxed density: P
  // Transition density: T
  Matrix_E DG = Matrix_E::Zero(M,M);
  Matrix_E P = Matrix_E::Zero(M,M);
  Matrix_E T = Matrix_E::Zero(M,M);

  for(int ii=0;ii<M;ii++) {
    for(int jj=0;jj<=ii;jj++) {
      // Ground state Density
      DG(ii,jj) = D[ii*M+jj];
      DG(jj,ii) = DG(ii,jj);
      // Difference relaxed density
      //P(ii,jj) = Diff[ii*M+jj] + Diff[jj*M+ii];
      //P(jj,ii) = P(ii,jj);
      P(ii,jj) = Diff[ii*M+jj];
      P(jj,ii) = Diff[jj*M+ii];
      // Transition density
      //T(ii,jj) = X[ii*M+jj] + X[jj*M+ii];
      //T(jj,ii) = T(ii,jj);
      T(ii,jj) = X[ii*M+jj];
      T(jj,ii) = X[jj*M+ii];
    }
  }

// Calculate Matrix G of derivatives
  vector<Matrix_E> G = compute_deriv(DG,P,T,M);
  int dim_atom = atoms.size();
  double force;

  for(int atom=0,ii=0; atom<dim_atom; atom++) {
    for(int xyz=0; xyz<3; xyz++,ii++) {
       //force = G[ii].cwiseProduct(DG).sum();
       force = G[ii](0,0);
       //force += 2.0f * G[ii].cwiseProduct(P).sum();
       //force += 2.0f * GT[ii].cwiseProduct(T).sum();
       For[xyz*dim_atom+atom] += 0.25f * force;
    }
  }

  for(int atom=0;atom<dim_atom;atom++)
     for(int xyz=0;xyz<3;xyz++)
       cout << atom << " " << xyz << " " << For[xyz*dim_atom+atom] << endl;

  DG.resize(0,0);
  P.resize(0,0);
  T.resize(0,0);
  vector<Matrix_E>().swap(G);
  //vector<Matrix_E>().swap(GT);
}

vector<Atom> LIBINTproxy::libint_geom(double* r,int natoms)
{
   vector<Atom> atoms(natoms);

   for(int i=0; i<natoms; i++) {
     atoms[i].x = r[i];
     atoms[i].y = r[i+natoms];
     atoms[i].z = r[i+natoms*2];
   }
   return atoms;
}

vector<libint2::Shell> LIBINTproxy::make_basis(
     const vector<Atom>& atoms,double*a,double*c,
     uint*ncont,uint*nuc,int s_func,int p_func,int d_func,int M)
{
   vector<Shell> obs_in;
   int from = 0;
   int to = s_func;

   for(int i=from; i<to; i++) {  // for s functions
     int centro = nuc[i]-1;
     int tam = ncont[i];
     vector<double> exp(tam);
     vector<double> coef(tam);
     for(int cont=0; cont<tam; cont++) {
        exp[cont] = a[i+M*cont];
        coef[cont] = c[i+M*cont];
     }
     obs_in.push_back(
        {
          exp,
          {
            {0, false, coef}
          },
          {{atoms[centro].x,atoms[centro].y,atoms[centro].z}}
        }
     );
     shell2atom.push_back(centro);
     vector<double>().swap(exp);
     vector<double>().swap(coef);
   }

   from = s_func;
   to = s_func+p_func*3;
   for(int i=from; i<to; i=i+3) { // for p functions
      int centro = nuc[i]-1;   
      int tam = ncont[i];
      vector<double> exp(tam);
      vector<double> coef(tam);
      for(int cont=0; cont<tam; cont++) {
         exp[cont] = a[i+M*cont];
         coef[cont] = c[i+M*cont];
      }
      obs_in.push_back(
         {
           exp,
           {
             {1, false, coef}
           },
           {{atoms[centro].x,atoms[centro].y,atoms[centro].z}}
         }
      );
      shell2atom.push_back(centro);
      vector<double>().swap(exp);
      vector<double>().swap(coef);
   }

   from = s_func+p_func*3;
   to = M;
   for(int i=from; i<to; i=i+6) { // for d functions
      int centro = nuc[i]-1;
      int tam = ncont[i];
      vector<double> exp(tam);
      vector<double> coef(tam);
      for(int cont=0; cont<tam; cont++) {
         exp[cont] = a[i+M*cont];
         coef[cont] = c[i+M*cont];
      }
      obs_in.push_back(
         {
           exp,
           {
             {2, false, coef}
           },
           {{atoms[centro].x,atoms[centro].y,atoms[centro].z}}
         }
      );
      shell2atom.push_back(centro);
      vector<double>().swap(exp);
      vector<double>().swap(coef);
    }

   return obs_in;
}

size_t LIBINTproxy::max_nprim() {
  size_t n = 0;
  for (auto shell: obs)
    n = std::max(shell.nprim(), n);
  return n;
}

int LIBINTproxy::max_l() {
  int l = 0;
  for (auto shell: obs)
    for (auto c: shell.contr)
      l = std::max(c.l, l);
  return l;
}

vector<int> LIBINTproxy::map_shell() {
  vector<int> result;
  result.reserve(obs.size());

  int n = 0;
  for (auto shell: obs) {
    result.push_back(n);
    n += shell.size();
  }
  return result;
}
