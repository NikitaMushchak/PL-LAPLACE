#pragma once

#include "io.hh"

/// \brief ������� �������� ������ � json �������
///
std::string ExportJson(
	std::vector<double> &Wk,
	double wn,
	std::vector<double> &pressure,
	std::vector<double>&concentration,
	std::vector<double> &x,
	std::vector<double> &y,
	double dx,
	double dy,
	size_t i00,
	size_t j00,
	std::vector< std::vector<size_t> > &index,
	double Time,
	double timeScale,
	double fluidEfficiency,
	double fluidDensity,
	double proppantDensity,
	double fluidInjection,
	double proppantInjection,
	double nominalStress,
	double Z_coordinate,
	std::string IdDesign,
	std::string IdStage
);

void ImportJSON(
	std::string J_String,
	std::string &IdDesign,
	std::string &IdStage,
	std::vector< std::vector<double> > &layers,
	std::vector< std::vector<double> > &injection,
	double &modelingTime,
	double &Emit_time,
	double &Z_coordinate
);

void ApproximateLayers(
	std::vector<std::vector<double> >&layers,
	double& cellSize
);
