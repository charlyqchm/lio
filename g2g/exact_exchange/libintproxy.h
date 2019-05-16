#ifndef LIBINTPROXY_H
#define LIBINTPROXY_H

#include <libint2.hpp>
#include "centros.h"
#include "../init.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix_E;

// namespace LIBINT
using libint2::Atom;
using libint2::BasisSet;
using libint2::Shell;
using libint2::Engine;
using libint2::Operator;
using libint2::BraKet;

// namespace STD
using std::vector;

typedef unsigned int uint;

namespace libint2 {
extern int nthreads;

/// fires off \c nthreads instances of lambda in parallel
template <typename Lambda>
void parallel_do(Lambda& lambda) {
#ifdef _OPENMP
#pragma omp parallel
  { 
    auto thread_id = omp_get_thread_num();
    lambda(thread_id);
  }
#endif
}
}

class LIBINTproxy 
{
private:
       // Variables
       vector<libint2::Shell> obs; // Object Basis
       vector<int> shell2bf; // first basis function
       vector<Atom> atoms;
       vector<int> shell2atom; // atom centre of shell
  

       // Functions
       vector<Atom> libint_geom(double*,int);

       vector<int> map_shell();

       vector<libint2::Shell> make_basis(const 
            vector<Atom>&,double*,double*,uint*,
            uint*,int,int,int,int);

       size_t max_nprim();
    
       int max_l();

       void counter_integrals();

       vector<Matrix_E> compute_deriv(const Matrix_E& D,int M);
       vector<Matrix_E> compute_deriv(const Matrix_E& D,const Matrix_E& P,
                                      const Matrix_E& T,int M);

public:
       LIBINTproxy(int,uint,uint*,double*,double*,
                   double*,uint*,int,int,int,int); // Constructor
       ~LIBINTproxy(); // Destructor

       void PrintBasis();

       void calculate(); // Save four center integrals

       void derivative(double* Dens, double* For);

       void derivative(double* D, double* Diff, 
                       double* X, double* For);
};

#endif
