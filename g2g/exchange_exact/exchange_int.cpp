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
extern "C" void g2g_ex_int_(double* Cbas ,double* K_int)

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
   int M2, M3;
   int ii, jj, kk, ll;
   M3 = M*M*M; M2=M*M;

   eri(K_int,M,fortran_vars.atoms,ncont,Cbas,aContr,pos,nuc,
       s_func,p_func,d_func);

}
