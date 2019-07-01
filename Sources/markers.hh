#pragma once

#include "planar3D.hh"

/*!
 \brief Расчет массы маркеров для проппантов
*/
void defineMarkerMass(
	std::vector< std::vector<double> > &injection,
	std::vector<double> &markerMass,
	std::vector<double> &markerVolume,
	std::vector<double> &proppantDensity,
	std::vector<double> &proppantDiameter,
	std::vector<double> &maximumConcentration,
	std::vector<size_t> &proppantType,
	double wn
);

/*!
 \brief Добавление новых маркеров для проппанта
*/
void addMarkers(
	std::vector< std::vector<double> > &markers,
	std::vector<double> &markerMass,
	std::vector<size_t> &proppantType,
	std::vector<double> &proppantDensity,
	double &concentration,
	double opening,
	double fluidInjection,
	double proppantInjection,
	double injectionIndex,
	double &injectedMass,
	double wn
);

/*!
 \brief Распределение маркеров по ячейкам
*/
void generateListOfMarkersInCells(
	std::vector< std::vector<size_t> > &markersInCells,
	std::vector< std::vector<double> > &markers,
	std::vector< std::vector<size_t> > &index
);

/*!
 \brief Расчет нового положения маркеров
*/
void markersTransport(
	std::vector< std::vector<size_t> > &markersInCells,
	std::vector< std::vector<double> > &markers,
	std::vector< std::vector<size_t> > &index,
	std::vector< std::vector<bool> > &elementIsActive,
	std::vector<double> &opening,
	std::vector<double> &pressure,
	std::vector<double> &concentration,
	std::vector<double> &proppantDensity,
	std::vector<double> &proppantDiameter,
	std::vector<double> &maximumConcentration,
	std::vector<double> &markerVolume,
	std::vector<double> &settlingVelocity,
	const double WPow,
	const double PPow,
	const double CPow,
	double fluidVelocityRight,
	double fluidVelocityBottom,
	double pressureIJ,
	double wn,
	size_t i,
	size_t j,
	double T
);

/*!
 \brief Перерасчет концентрации проппанта при удвоении сетки
*/
void recalcConcentation(
	std::vector< std::vector<double> > &markers,
	const std::vector< std::vector<size_t> > &index,
	std::vector<double> &concentration,
	std::vector<double> &opening,
	std::vector<double> &markerVolume,
	const double wn
);

/*!
\brief Вывод значений концентраций проппантов
*/
std::string outputConcentration(
	std::vector< std::vector<double> > &markers,
	std::vector< std::vector<size_t> > &index,
	std::vector<double> &concentration,
	std::vector<double> &opening,
	std::vector<double> &markerVolume,
	double wn,
	size_t i,
	size_t j,
	int &num_fluids,
	int &num_proppants
);


///Блок функций для маркеров для нескольких жидкостей (?)

/*!
 \brief Расчет объема маркера для жидкостей и заполнение начальной трещины маркерами
*/
void defineMarkerVolumeAndFillCrackWithMarkers(
	std::vector< std::vector<size_t> > &activeElements,
	std::vector< std::vector<size_t> > &index,
	std::vector< std::vector<double> > &markers,
	std::vector< size_t > &fluidType,
	std::vector<double> &opening,
	std::vector<double> &x,
	std::vector<double> &y,
	double &markerVolume
);

/*!
 \brief Добавление новых маркеров для жидкости
*/
void addMarkersFluid(
	std::vector< std::vector<double> > &markers,
	std::vector< size_t > &fluidType,
	double fluidInjection,
	double injectionIndex,
	double &injectedVolume,
	double markerVolume
);

/*!
 \brief Расчет эффективных вязкости, индекса течения и плотности жидкости в ячейке
*/
void calcEffectiveViscosity(
	std::vector< std::vector<size_t> > &markersInCells,
	std::vector< std::vector<size_t> > &index,
	std::vector< std::vector<double> > &markers,
	std::vector<double> &fluidViscosity,
	std::vector<double> &rheologyIndex,
	std::vector<double> &fluidDensity,
	double &effectiveViscosity,
	double &effectiveN,
	double &effectiveDensity,
	size_t i,
	size_t j
);

/*!
 \brief Расчет нового положения маркеров
*/
void markersFluidTransport(
	std::vector< std::vector<size_t> > &markersInCells,
	std::vector< std::vector<double> > &markers,
	std::vector< std::vector<size_t> > &index,
	std::vector< std::vector<bool> > &elementIsActive,
	std::vector<double> &opening,
	std::vector<double> &pressure,
	std::vector<double> &concentration,
	std::vector<double> &fluidViscosity,
	const double WPow,
	const double PPow,
	const double CPow,
	double fluidVelocityRight,
	double fluidVelocityBottom,
	double pressureIJ,
	double wn,
	double markerVolume,
	size_t i,
	size_t j
);

