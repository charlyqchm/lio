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
extern "C" void g2g_ex_exact_(double* rho_mat, double* Cbas ,double* Fk_mat,
                              int& numvec, int& int2elec)
//carlos: esto esta en el orden para trabajar, sin embargo algunas variables de
//        obtein no estan definidas. Volver a checkear eso.
{
   int M = fortran_vars.m;
   int M2, M3;
   M3 = M*M*M; M2=M*M;
   int s_func = fortran_vars.s_funcs;
   int p_func = fortran_vars.p_funcs;
   int d_func = fortran_vars.d_funcs;
   double* aContr = &fortran_vars.a_values(0,0);
   double* pos = &fortran_vars.atom_positions_pointer(0,0);
   uint* ncont = &fortran_vars.contractions(0);
   uint* nuc = &fortran_vars.nucleii(0);
   double timeI, timeF;

   timeI = timeF = 0.0f;
   if (int2elec == 0 ) {
     timeI = omp_get_wtime();
     // INITIALIZATION LIBINT
     LIBINTproxy libintproxy(M,fortran_vars.atoms,ncont,
                 Cbas,aContr,pos,nuc,s_func,p_func,d_func,0);
     // LIBINT CALCULATE INTEGRALS
     libintproxy.calculate();
     timeF = omp_get_wtime();
     printf(" LIBINT CALCULATION %f\n",timeF-timeI);
   }
   Obtain fock(numvec,M,fortran_vars.dim,0,fortran_vars.a_PBE0);
   fock.calculate(rho_mat,fortran_vars.Kmat,Fk_mat);

   fflush(stdout); // NOT BUFFERED

}
