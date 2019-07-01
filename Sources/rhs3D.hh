#pragma once

#include "planar3D.hh"

#if defined(USE_EIGEN)
    #include "eigen/Dense"
    #include "eigen/Sparse"
    #include "eigen/Core"
#endif

/*!
 \brief Умножение матрицы на вектор
*/
#if defined(USE_EIGEN)
inline void multiply(
    Eigen::MatrixXd &matrix,
    Eigen::VectorXd &vector,
    std::vector<double> &result
);
#else
inline void multiply(
    std::vector< std::vector<double> > &matrix,
    std::vector<double> &vector,
    std::vector<double> &result
);
#endif

/*!
 \brief Расчет давления в активных элементах
*/
void calculatePressure(
    std::vector<double> &pressure,
    std::vector< std::vector<std::size_t> > &index,
    std::vector< std::vector<std::size_t> > &activeElements,
    #if defined(USE_EIGEN)
    Eigen::MatrixXd &influenceMatrix,
    Eigen::VectorXd &opening,
    #else
    std::vector< std::vector<double> > &influenceMatrix,
    std::vector<double> &opening,
    #endif
    std::vector<double> &stress
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
    double fluidDensity,
    double proppantDensity,
    double proppantDiameter,
    double maximumConcentration,
    double currentTime,
    double n,
    double wn
);


/*!
 \brief Расчет производной раскрытия по времени с учетом влияния пропанта и
 производной концентрации пропанта по времени. Маркеры для проппантов
*/
void calculateOpeningAndConcentrationSpeedsDP(
    std::vector<double> &dWdt,
    std::vector<double> &opening,
    std::vector<double> &concentration,
    std::vector<double> &pressure,
    std::vector<double> &leakOff,
    std::vector< std::vector<size_t> > &index,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector< std::vector<bool> > &elementIsActive,
    std::vector< std::vector<double> > &markers,
    std::vector<double> &activationTime,
    std::vector<double> &proppantDensity,
    std::vector<double> &proppantDiameter,
    std::vector<double> &maximumConcentration,
    std::vector<double> &markerMass,
    std::vector<double> &markerVolume,
    double fluidDensity,
    double currentTime,
    double n,
    double wn
);

/*!
 \brief Расчет производной раскрытия по времени с учетом влияния пропанта и
 производной концентрации пропанта по времени. Маркеры для жидкостей
*/
void calculateOpeningAndConcentrationSpeedsDF(
    std::vector<double> &dWdt,
    std::vector<double> &dCdt,
    std::vector<double> &opening,
    std::vector<double> &concentration,
    std::vector<double> &pressure,
    std::vector<double> &leakOff,
    std::vector< std::vector<size_t> > &index,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector< std::vector<bool> > &elementIsActive,
    std::vector< std::vector<double> > &markers,
    std::vector<double> &activationTime,
    std::vector<double> &fluidViscosity,
    std::vector<double> &rheologyIndex,
    std::vector<double> &fluidDensity,
    double proppantDensity,
    double proppantDiameter,
    double maximumConcentration,
    double markerVolume,
    double currentTime,
    double n,
    double wn
);
