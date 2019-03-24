#pragma once

/*!
 \brief Определение растояния до фронта трещины
*/
double collectDistanceFromVelocity(
    const size_t i,
    const size_t j,
    double n,
    double opening,
    std::vector<double> &leakOff,
    std::vector<Ribbon> &ribbons,
    std::vector< std::vector<double> > &velocities
);
