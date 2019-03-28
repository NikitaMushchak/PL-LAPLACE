#pragma once

#include "planar3D.hh"

/*!
 \brief Умножение матрицы на вектор
*/
inline void multiply(
    std::vector< std::vector<double> > &matrix,
    std::vector<double> &vector,
    std::vector<double> &result
);

/*!
 \brief Расчет давления в активных элементах
*/
// void calculatePressure(
//     std::vector<double> &pressure,
//     std::vector< std::vector<std::size_t> > &index,
//     std::vector< std::vector<std::size_t> > &activeElements,
//     std::vector< std::vector<double> > &influenceMatrix,
//     std::vector<double> &opening,
//     std::vector<double> &stress
// );
void calculatePressure(
    std::vector<double> &pressure,
    std::vector< std::vector<std::size_t> > &index,
    std::vector< std::vector<std::size_t> > &activeElements,
    std::vector< std::vector<double> > &influenceMatrix,
    std::vector<double> &opening,
    std::vector<double> &stress,
    std::vector<double> &T,
    std::vector<double> &b,
    std::vector<double> &r1,
    std::vector<double> &r2,
    std::vector<double> &p,
    std::vector<double> &A_p,
    size_t N_dof,
    size_t Nx,
    size_t Ny
);
/*!
 \brief Расчет скорости фронта
*/
void calculateFrontVelocity(
    std::vector< std::vector<double> > &velocities,
    std::vector<double> &leakOff,
    std::vector< std::vector<size_t> > &index,
    std::vector<Ribbon> &ribbons,
    std::vector< std::vector<double> > &distances,
    std::vector<double> &opening,
    double n
);

/*!
 \brief Расчет производной раскрытия по времени с учетом влияния пропанта и
 производной концентрации пропанта по времени
*/
void calculateOpeningAndConcentrationSpeeds(
    std::vector<double> &dWdt,
    std::vector<double> &dCdt,
    std::vector<double> &opening,
    std::vector<double> &concentration,
    std::vector<double> &pressure,
    std::vector<double> &leakOff,
    std::vector< std::vector<size_t> > &index,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector< std::vector<bool> > &elementIsActive,
    std::vector<double> &activationTime,
    double currentTime,
    double fluidDensity,
    double proppantDensity,
    double n
);
