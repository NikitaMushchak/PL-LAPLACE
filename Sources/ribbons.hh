#pragma once

#include "planar3D.hh"

/*!
 \brief Нахождение начальных граничных элементов
*/
std::vector<Ribbon> findRibbons(
    std::vector< std::vector<Cell> > &mesh,
    const std::vector< std::vector<size_t> > activeElements,
    std::vector< std::vector<double> > &distances,
    double testDistance,
    const double dMin,
    const double dMax
);

/*!
 \brief Нахождение новых граничных элементов
*/
void findNewRibbons(
    const size_t i,
    const size_t j,
    const double d,
    const double dMin,
    const double dCenter,
    double n,
    std::vector<double> &Wt,
    std::vector<double> &leakOff,
    std::vector< std::vector<Cell> > &mesh,
    std::vector< std::vector<size_t> > &index,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector< std::vector<bool> > &elementIsActive,
    std::vector<double> &activationTime,
    std::vector<Ribbon> &ribbons,
    std::vector<Ribbon> &ribbonsOld,
    std::vector< std::vector<double> > &distances,
    std::vector< std::vector<double> > &velocities,
    const double currentTime,
    const double opening
);