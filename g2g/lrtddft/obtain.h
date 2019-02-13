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

       void method_A(double* T, FourCenter* K, double* F);

public:
       Obtain(int vec, int size, int dimension); // constructor
       ~Obtain(); // Destructor

       void calculate(double* Tmat, FourCenter* Kmat, double* Fmat);
};

#endif

