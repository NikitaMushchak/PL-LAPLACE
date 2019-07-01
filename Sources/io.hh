#pragma once

#include <vector>

#include "planar3D.hh"

struct DLL_Param;

/// \todo aaaa
template <typename T>
void saveConcentrationData(
	const std::string filename,
	std::vector<T> &data,
	std::vector< std::vector<std::size_t> > &index,
	const double multiplier = 1.,
	const std::string comment = std::string()
) {
	std::ofstream output(filename + std::string("_m.txt"));

	if (!output.good()) {
		throw std::runtime_error(
			ai::string("Exception while saving the matrix into the file: ")
			+ filename
		);
	}

	if (std::string() != comment) {
		output << comment << std::endl;
	}

	for (int j = (int)index[0].size() - 1; j >= 0; --j) {
		for (std::size_t i = 0; i < index.size(); ++i) {
#if defined(PROPPANT_MARKERS)
			if (i00 == i) {
				output << std::setw(14) << 2. * multiplier * data[index[i][j]];
			}
			else {
				output << std::setw(14) << multiplier * data[index[i][j]];
			}
#else
			output << std::setw(14) << multiplier * data[index[i][j]];
#endif
		}

		output << std::endl;
	}

	output.close();
}

/// \todo aaaa
template<typename T>
void saveData(
	const std::string filename,
	std::vector<T> &data,
	std::vector< std::vector<std::size_t> > &index,
	const double multiplier = 1.,
	const std::string comment = std::string()
) {
	std::ofstream output(filename + std::string("_m.txt"));

	if (!output.good()) {
		throw std::runtime_error(
			ai::string("Exception while saving the matrix into the file: ")
			+ filename
		);
	}

	if (std::string() != comment) {
		output << comment << std::endl;
	}

	for (int j = (int)index[0].size() - 1; j >= 0; --j) {
		for (std::size_t i = 0; i < index.size(); ++i) {
			output << std::setw(14) << multiplier * data[index[i][j]];
		}

		output << std::endl;
	}

	output.close();
}

/// \todo aaaa
template<typename T>
void saveData(
	const std::string filename,
	std::vector< std::vector<T> > &data,
	const double multiplier = 1.,
	const std::string comment = std::string()
) {
	std::ofstream output(filename + std::string("_m.txt"));

	if (!output.good()) {
		throw std::runtime_error(
			ai::string("Exception while saving the matrix into the file: ")
			+ filename
		);
	}

	if (std::string() != comment) {
		output << comment << std::endl;
	}

	for (int j = (int)data[0].size() - 1; j >= 0; --j) {
		for (std::size_t i = 0; i < data.size(); ++i) {
			output << std::setw(14) << data[i][j];
		}

		output << std::endl;
	}

	output.close();
}


std::string ExportJson(
	std::vector<double> &Wk,
	double wn,
	std::vector<double> &pressure,
	std::vector<double>&concentration,
	std::vector< std::vector<double>>  &markers,	//Добавил информацию для вывода маркеров от светы
	std::vector<double> &markerVolume,				//
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
	std::vector<double> &fluidDensity,
	std::vector<double> &proppantDensity,
	double fluidInjection,
	double proppantInjection,
	double nominalStress,
	DLL_Param &DLL_Parametrs
);

void ImportJSON(
	DLL_Param &DLL_Parametrs,
	std::vector< std::vector<double> > &layers_new,
	std::vector< std::vector<double> > &injection
);

void ApproximateLayers(std::vector< std::vector<double> >& layers, double& cellSize);
