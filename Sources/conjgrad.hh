#pragma once
#include <vector>
// #include "ai.hh"

void multiplyDiag(std::vector<double> &y, std::vector<std::vector<double> > &A,std::vector<double > &x,
                        size_t Nx , size_t NxNy , size_t N_dof);

double MultiplyVV(std::vector<double>&a, std::vector<double>&b);


double NormV(std::vector<double>&x);



void conjGrad(std::vector<double> &x,
                    std::vector<std::vector<double> > &A,
                            std::vector<double> &b,
                            size_t Nx , size_t NxNy , size_t N_dof);
