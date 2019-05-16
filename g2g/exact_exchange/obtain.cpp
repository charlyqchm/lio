// TODO: checkear la cuenta de intercambio exacto en Excited State
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

#include "obtain.h"
#include "centros.h"


Obtain::Obtain(int vec, int size, int dimension, int MET, double a_PBE0)
{
   total_vec = vec;
   int_total = dimension;
   M  = size;
   M2 = size * size;
   M3 = size * size * size;
   timerI = timerF = 0.0f;
   method = MET;
   cfactor = 1.0f - a_PBE0 ;
}

Obtain::~Obtain()
{
   if (method == 1)
      printf(" TIME FOCK METHOD: %f\n",timerF-timerI);
}

void Obtain::calculate(double* Tmat, FourCenter* Kmat, double* Fmat)
{
   if ( method == 0 ) {
      method_B(Tmat,Kmat,Fmat); // Exact Exchange - PBE0 - SCF
   } else {
      method_A(Tmat,Kmat,Fmat); // Exact Exchange - PBE0 - Excited State
   }
}

void Obtain::method_A(double* T, FourCenter* K, double* F)
{
// T is not a symmetric matrix
   double Dens = 0.0f;
   double alpha = 0.0f;
   int p1, p2, p3, p4;

   alpha = cfactor;
   timerI = omp_get_wtime();
#pragma omp parallel for private(p1,p2,p3,p4,Dens)
   for(int ivec=0; ivec<total_vec; ivec++) {
      for(int i=0; i<int_total; i++) {
        // FOR EXCHANGE EXACT
        p1 = ivec*M2+K[i].u*M+K[i].k;
        p2 = ivec*M2+K[i].k*M+K[i].u;
        p3 = ivec*M2+K[i].v*M+K[i].l;
        p4 = ivec*M2+K[i].l*M+K[i].v;

        Dens = alpha * T[p3];
        F[p1] -= K[i].result * Dens;
        Dens = alpha * T[p4];
        F[p2] -= K[i].result * Dens;
        Dens = alpha * T[p1];
        F[p3] -= K[i].result * Dens;
        Dens = alpha * T[p2];
        F[p4] -= K[i].result * Dens;

        p1 = ivec*M2+K[i].u*M+K[i].l;
        p2 = ivec*M2+K[i].l*M+K[i].u;
        p3 = ivec*M2+K[i].v*M+K[i].k;
        p4 = ivec*M2+K[i].k*M+K[i].v;

        Dens = alpha * T[p3];
        F[p1] -= K[i].result * Dens;
        Dens = alpha * T[p4];
        F[p2] -= K[i].result * Dens;
        Dens = alpha * T[p1];
        F[p3] -= K[i].result * Dens;
        Dens = alpha * T[p2];
        F[p4] -= K[i].result * Dens;

        // FOR TYPE COULOMB
        p1 = ivec*M2+K[i].k*M+K[i].l;
        p2 = ivec*M2+K[i].l*M+K[i].k;
        p3 = ivec*M2+K[i].u*M+K[i].v;
        p4 = ivec*M2+K[i].v*M+K[i].u;

        Dens = ( T[p1] + T[p2] ) * 2.0f;
        F[p3] += K[i].result * Dens;
        F[p4] += K[i].result * Dens;

        Dens = ( T[p3] + T[p4] ) * 2.0f;
        F[p1] += K[i].result * Dens;
        F[p2] += K[i].result * Dens;
     }
   }
   timerF = omp_get_wtime();
}

void Obtain::method_B(double* T, FourCenter* K, double* F)
{
   double Dens = 0.0f;
   double alpha = cfactor * 2.0f;
   int p1, p2, p3, p4;

   timerI = omp_get_wtime();
  for(int i=0; i<int_total; i++) {
     p1 = K[i].u*M+K[i].k;
     p2 = K[i].k*M+K[i].u;
     p3 = K[i].v*M+K[i].l;
     p4 = K[i].l*M+K[i].v;

     Dens = alpha * T[p3];
     F[p1] += K[i].result * Dens;
     Dens = alpha * T[p4];
     F[p2] += K[i].result * Dens;
     Dens = alpha * T[p1];
     F[p3] += K[i].result * Dens;
     Dens = alpha * T[p2];
     F[p4] += K[i].result * Dens;

     p1 = K[i].u*M+K[i].l;
     p2 = K[i].l*M+K[i].u;
     p3 = K[i].v*M+K[i].k;
     p4 = K[i].k*M+K[i].v;

     Dens = alpha * T[p3];
     F[p1] += K[i].result * Dens;
     Dens = alpha * T[p4];
     F[p2] += K[i].result * Dens;
     Dens = alpha * T[p1];
     F[p3] += K[i].result * Dens;
     Dens = alpha * T[p2];
     F[p4] += K[i].result * Dens;
   }

   timerF = omp_get_wtime();
}
