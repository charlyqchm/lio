#include <iostream>
#include <stdio.h>

#include "../init.h"
#include "../common.h"
//#include "../partition.h"

#include "../lrtddft/eri.h"
#include <stdlib.h>

using namespace G2G;
//extern Partition partition;

//##############################################################################
//##############################################################################
extern "C" void g2g_ex_exact_(double* rho_mat, double* Cbas ,double* K_mat,double* K_int)

{
   int M = fortran_vars.m;
   int nco = fortran_vars.nco;
   int s_func = fortran_vars.s_funcs;
   int p_func = fortran_vars.p_funcs;
   int d_func = fortran_vars.d_funcs;
   double* aContr = &fortran_vars.a_values(0,0);
   double* pos = &fortran_vars.atom_positions_pointer(0,0);
   uint* ncont = &fortran_vars.contractions(0);
   uint* nuc = &fortran_vars.nucleii(0);
//   double* K_int = (double*) malloc(M*M*M*M*sizeof(double)) ;
   int M2, M3;
   int ii, jj, kk, ll;
   M3 = M*M*M; M2=M*M;

//carlos: K_int is calculated, this matrix stored the 4 center integrals, and
//        for thar reason it has 4 dimensions
//   K_int[M*M*M*M] = {0.0};
   // for (ii=0; ii<M; ++ii) {
   // for (jj=0; jj<M; ++jj) {
   // for (kk=0; kk<M; ++kk) {
   // for (ll=0; ll<M; ++ll) {
   //    printf("%d %d %d %d %f \n", ii, kk, ll,jj, K_int[ii*M3+kk*M2+ll*M+jj] );
   // }
   // }
   // }
   // }
   printf("%s %d %d %f", "HOLA", M, fortran_vars.atoms, Cbas[0]);
   eri(K_int,M,fortran_vars.atoms,ncont,Cbas,aContr,pos,nuc,
       s_func,p_func,d_func);
   // printf("%s \n","--------------------------------");
   // for (ii=0; ii<M; ++ii) {
   // for (jj=0; jj<M; ++jj) {
   // for (kk=0; kk<M; ++kk) {
   // for (ll=0; ll<M; ++ll) {
   //     printf("%d %d %d %d %f \n", ii, kk, ll,jj, K_int[ii*M3+kk*M2+ll*M+jj] );
   // }
   // }
   // }
   // }



//carlos: now we try to suit this results in the exchange matrix:
   // for (ii=0; ii<M; ++ii) {
   // for (jj=0; jj<M; ++jj) {
   //    printf("%d %d %f \n", ii, jj, rho_mat[ii*M+jj]);
   // }
   // }


   for (ii=0; ii<M; ++ii) {
   for (jj=0; jj<M; ++jj) {
   for (kk=0; kk<M; ++kk) {
   for (ll=0; ll<M; ++ll) {
      K_mat[ii*M+jj] += 0.5*rho_mat[kk*M+ll]*K_int[ii*M3+kk*M2+ll*M+jj];
   }
   }
   }
   }

   // for (ii=0; ii<M; ++ii) {
   // for (jj=0; jj<M; ++jj) {
   // for (kk=0; kk<M; ++kk) {
   // for (ll=0; ll<M; ++ll) {
   //    printf("%d %d %d %d %f %f %f %f %f %f \n", ii, kk, ll,jj, K_mat[ii*M+jj], K_mat[jj*M+ii],rho_mat[kk*M+ll],rho_mat[ll*M+kk], K_int[ii*M3+kk*M2+ll*M+jj], K_int[jj*M3+kk*M2+ll*M+ii] );
   // }
   // }
   // }
   // }
   // printf("%s %f %f","hola", K_int[0*M3+0*M2+0*M+2], K_int[2*M3+0*M2+0*M+0]);
}
