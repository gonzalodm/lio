#ifndef OBTAIN_H
#define OBTAIN_H

class Obtain
{
private:
       int M, M2, M3;
       int method;
       int total_vec;
       double timerI, timerF;

       void method_A(double* T, double* K, double* F);

public:
       Obtain(int vec, int size); // constructor
       ~Obtain(); // Destructor

       void calculate(double* Tmat, double* Kmat, double* Fmat);
};

#endif

