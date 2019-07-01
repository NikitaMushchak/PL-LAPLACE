#pragma once

#include <vector>

#include "planar3D.hh"

int completeMesh(
    std::size_t meshSize,
    std::vector<double> &x,
    std::vector<double> &y,
    std::size_t &i00,
    std::size_t &j00,
    std::vector< std::vector<Cell> > &mesh,
    std::vector<Ribbon> &ribbons,
    std::vector< std::vector<std::size_t> > &index,
    std::vector<double> &opening,
    std::vector<double> &pressure,
    std::vector<double> &concentration,
    std::vector< std::vector<double> > &distances,
    std::vector< std::vector<std::size_t> > &activeElements,
    std::vector< std::vector<bool> > &elementIsActive,
    std::vector<double> &activationTime,
    std::vector< std::vector<double> > &layers,
    std::vector<double> &flatYoungsModulus,
    std::vector<double> &stress,
    std::vector<double> &leakOff,
    std::vector<double> &toughness,
    const double wn,
    const double timeScale,
    const double mu,
    const double n,
    std::size_t &xSize,
    std::size_t &ySize,
    bool considerElasticModulusContrast,
    std::vector< std::vector<double> > &influenceMatrix
);

/// \todo aaaaa
int scaleMesh(
    std::vector<double> &x,
    std::vector<double> &y,
    std::vector< std::vector<Cell> > &mesh,
    const std::vector< std::vector<std::size_t> > &index,
    std::vector< std::vector<std::size_t> > &activeElements,
    std::vector< std::vector<bool> > &elementIsActive,
    std::vector<double> &activationTime,
    std::vector<Ribbon> &ribbons,
    std::vector< std::vector<double> > &distances,
    std::vector<double> &opening,
    std::vector<double> &concentration,
    const std::vector< std::vector<double> > &layers,
    std::vector<double> &flatYoungsModulus,
    std::vector<double> &stress,
    std::vector<double> &leakOff,
    std::vector<double> &toughness,
    const std::vector<double> &zeroVectorXY,
    const std::vector< std::vector<double> > &zeroMatrixXY,
    const std::size_t xSize,
    const std::size_t ySize,
    double &initialRadius,
    double &axMax,
    double &dMin1,
    double &dCenter1,
    double &dMax1,
    double &dMin2,
    double &dCenter2,
    double &dx,
    double &dy,
    double &dt,
    const double mu,
    const double n,
    const double wn,
    const double timeScale,
    const double T,
    const double T0,
    double &E,
    const bool considerElasticModulusContrast,
	#if defined(PROPPANT_MARKERS)
    std::vector< std::vector<double> > &markers,
    std::vector<double> &markerVolume,
	#endif
	int &meshScalingCounter
);