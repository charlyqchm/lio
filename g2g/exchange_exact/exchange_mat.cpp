#include <iostream>
#include <stdio.h>

#include "../init.h"
#include "../common.h"
#include <stdlib.h>

using namespace G2G;

//##############################################################################
//##############################################################################
extern "C" void g2g_exchange_mat_(double* rho_mat, double* K_mat, double* K_int)

{
   int M = fortran_vars.m;
   int M2, M3;
   int ii, jj, kk, ll;
   M3 = M*M*M; M2=M*M;

   for (ii=0; ii<M; ++ii) {
   for (jj=0; jj<M; ++jj) {
   for (kk=0; kk<M; ++kk) {
   for (ll=0; ll<M; ++ll) {
      K_mat[ii*M+jj] += 0.5*rho_mat[kk*M+ll]*K_int[ii*M3+kk*M2+ll*M+jj];
   }
   }
   }
   }
}
