#ifndef OBTAIN_H
#define OBTAIN_H

#include "centros.h"

class Obtain
{
private:
       int M, M2, M3;
       int int_total;
       int method;
       int total_vec;
       double timerI, timerF;
       double cfactor; // Factor exchange 

       void method_A(double* T, FourCenter* K, double* F); // CALL FROM EXCITED STATES
       void method_B(double* T, FourCenter* K, double* F); // CALL FROM SCF

public:
       Obtain(int vec, int size, int dimension, int MET, bool PBE0); // constructor
       ~Obtain(); // Destructor

       void calculate(double* Tmat, FourCenter* Kmat, double* Fmat);
};

#endif

