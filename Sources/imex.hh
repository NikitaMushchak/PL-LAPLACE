#pragma once

#include "planar3D.hh"

void multiplyFiveDiagonalsFast(
    std::vector<double> &vector,
    std::vector<double> &result,
    const double coeff,
    const std::size_t Nx,
    const std::size_t Ny
);
double MultiplyVV(std::vector<double>&a, std::vector<double>&b);
double NormV(std::vector<double>&x);
void conjugateGradient(
    std::vector<double> &b,
    std::vector<double> &x,
    const double coeff,
    const std::size_t Nx,
    const std::size_t Ny
);