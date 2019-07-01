#include "ailibrary/ai.hh"
#include "planar3D.hh"

/*!
 \details Функция рассчитывает массу маркеров для разных 
 проппантов на разных этапах закачки

 \param[in] injection - план закачки
 \param[in] markerMass - вектор масс маркеров
 \param[in] proppantDensity - вектор плотностей проппантов
 \param[in] proppantDiameter - вектор диаметров проппантов
 \param[in] maximumConcentration - вектор максимальных концентраций проппантов
 \param[in] proppantType - вектор типов проппантов
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
) {
	double injectionMass;
	for (std::size_t i = 0; i < proppantDensity.size(); ++i) {
		markerMass.push_back(maximumConcentration[i] * proppantDensity[i]
			* proppantDiameter[i] * dx * dy / 5.0 );
	}
	for (std::size_t i = 0; i < injection.size(); ++i) {
		if (injection[i][3] > epsilon) {
			injectionMass = injection[i][2] * wn * injection[i][3] * (injection[i][1]
				- injection[i][0]) / 2.;
			markerMass[proppantType[i]] = injectionMass / (std::ceil(injectionMass
				/ markerMass[proppantType[i]]));
			markerVolume.push_back(markerMass[proppantType[i]] / proppantDensity[proppantType[i]]);
		}
	}
}

/*!
 \details Функция добавляет в массив маркеров новые маркеры по ходу закачки

 \param[in] index - матрица индексов
 \param[in] markers - матрица положения маркеров
 \param[in] markerMass - вектор масс маркеров
 \param[in] proppantType - вектор типов проппантов
 \param[in] proppantDensity - вектор плотностей проппантов
 \param[in] opening - вектор раскрытий
 \param[in] concentration - вектор концентраций проппанта
 \param[in] fluidInjection - расход жидкости
 \param[in] proppantInjection - плотность закачки проппанта
 \param[in] injectionIndex - номер стадии закачки
 \param[in] injectedMass - масса закаченного проппанта
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
) {
	injectedMass += fluidInjection * proppantInjection * dt * wn / 2.;
	int nMarkToAdd = int(floor(injectedMass / markerMass[proppantType[injectionIndex - 1]]));
	if (nMarkToAdd > 0) {
		concentration += nMarkToAdd * markerMass[proppantType[injectionIndex - 1]]
			/ (proppantDensity[proppantType[injectionIndex - 1]] * opening * wn * dx * dy);

		for (size_t k = 0; k < nMarkToAdd; ++k) {
			std::vector< double > xCoordinate;
			std::vector< double > yCoordinate;
			ai::generateRandomVector(xCoordinate, 1, epsilon, 0.5 * dx - epsilon);
			ai::generateRandomVector(yCoordinate, 1, - 0.5 * dx + epsilon, 
				0.5 * dx - epsilon);
			markers.push_back(std::vector<double>{xCoordinate[0], yCoordinate[0], 
				double(proppantType[injectionIndex - 1])});
		}
		injectedMass -= nMarkToAdd * markerMass[proppantType[injectionIndex - 1]];
	}
}

/*!
 \details Функция заполняющая индексами маркеров, находящихся в соответствующих ячейках,
 вектор векторов markersInCells с длиной, равной количеству ячеек

 \param[in] markersInCells - матрица индексов маркеров в соответствующих ячейках
 \param[in] markers - матрица положения маркеров
 \param[in] index - матрица индексов
*/
void generateListOfMarkersInCells(
	std::vector< std::vector<size_t> > &markersInCells,
	std::vector< std::vector<double> > &markers,
	std::vector< std::vector<size_t> > &index
) {
	for (std::size_t i = 0; i < markers.size(); ++i) {
		int cell = index[int(round(markers[i][0] / dx + double(i00)))][int(round(markers[i][1]
			/ dx + double(j00)))];
		markersInCells[cell].resize(markersInCells[cell].size() + 1);
		markersInCells[cell][markersInCells[cell].size() - 1] = i;
	}
}

/*!
 \details Функция переноса маркеров

 \param[in] markersInCells - матрица индексов маркеров в соответствующих ячейках
 \param[in] markers - матрица положения маркеров
 \param[in] index - матрица индексов
 \param[in] elementIsActive - матрица активных элементов
 \param[in] opening - вектор раскрытий
 \param[in] pressure - вектор давлений
 \param[in] concentration - вектор концентраций проппанта
 \param[in] proppantDensity - вектор плотностей проппантов
 \param[in] proppantDiameter - вектор диаметров проппантов
 \param[in] maximumConcentration - вектор максимальных концентраций проппантов
 \param[in] markerVolume - вектор объемов маркеров
 \param[in] settlingVelocity - вектор скоростей оседания проппантов
 \param[in] WPow - степень для раскрытия в уравнении Пуассона
 \param[in] PPow - степень для давления в уравнении Пуассона
 \param[in] CPow - степень для концентрации в уравнении Пуассона
 \param[in] fluidFlowRight - поток жидкости справа
 \param[in] fluidFlowBottom - поток жидкости снизу
 \param[in] pressureIJ - давление в ячейке
 \param[in] wn - нормирующий множитель
 \param[in] i - индекс ячейки (x)
 \param[in] j - индекс ячейки (y)
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
) {
	double fluidVelocityLeft;
	double fluidVelocityTop;
	double pressureDrop;
	double openingOnTheBorder;
	double concentrationOnTheBorder;

	double alphaX;
	double alphaY;
	double velocityX;
	double velocityY;
	double xNew;
	double yNew;
	int iNew;
	int jNew;
	int markerIndex;
	
	fluidVelocityLeft = 0.;
	fluidVelocityTop = 0.;

	if (i == i00) {
		fluidVelocityLeft = -fluidVelocityRight;
	} else {
		if (elementIsActive[i - 1][j]) {
			pressureDrop = pressure[index[i - 1][j]] - pressureIJ;
			openingOnTheBorder = 0.5
				* (opening[index[i][j]] + opening[index[i - 1][j]]);
			concentrationOnTheBorder = 0.5
				* (concentration[index[i][j]] + concentration[index[i - 1][j]]);
			fluidVelocityLeft = -ai::sign(pressureDrop)
				* std::pow(openingOnTheBorder, WPow - 1.)
				* std::pow(std::abs(pressureDrop) / dx, PPow)
				/ std::pow(
					1. - concentrationOnTheBorder / 0.6, // / maximumconcentration,
					CPow
				);
		}
	}

	if (elementIsActive[i][j - 1]) {
		pressureDrop = pressure[index[i][j - 1]] - pressureIJ;
		openingOnTheBorder = 0.5
			* (opening[index[i][j]] + opening[index[i][j - 1]]);
		concentrationOnTheBorder = 0.5
			* (concentration[index[i][j]] + concentration[index[i][j - 1]]);
		fluidVelocityTop = -ai::sign(pressureDrop)
			* std::pow(openingOnTheBorder, WPow - 1.)
			* std::pow(std::abs(pressureDrop) / dy, PPow)
			/ std::pow(
				1. - concentrationOnTheBorder / 0.6,// / maximumConcentration,
				CPow
			);
	}

	for (std::size_t rInd = 0; rInd < markersInCells[index[i][j]].size(); ++rInd) {
		markerIndex = int(markers[markersInCells[index[i][j]][rInd]][2]);

		alphaX = std::abs(markers[markersInCells[index[i][j]][rInd]][0] - ((double)i - (double)i00 - 0.5) * dx) / dx;
		velocityX = (1. - alphaX) * fluidVelocityLeft + alphaX * fluidVelocityRight;
		xNew = markers[markersInCells[index[i][j]][rInd]][0] - dt * velocityX;

		alphaY = std::abs(markers[markersInCells[index[i][j]][rInd]][1] - ((double)j - (double)j00 - 0.5) * dx) / dx;
		velocityY = (1. - alphaY) * fluidVelocityTop + alphaY * fluidVelocityBottom;
		yNew = markers[markersInCells[index[i][j]][rInd]][1] - dt * (velocityY + std::pow((1.	//Учет оседания проппанта '-' между velocityY - std::pow означает оседание вниз.
			- concentration[index[i][j]] / maximumConcentration[markerIndex]), 5)
			* settlingVelocity[markerIndex]);
		//yNew = markers[rInd][1] - dt * velocityY;

		iNew = int(round(xNew / dx + i00));
		jNew = int(round(yNew / dy + j00));

		if ((iNew < 0) || (jNew < 0)) {
			saveData(
				ai::string("./Results/Concentration/")
				+ ai::prependNumber(T, 3),
				concentration,
				index
			);
			std::cout << "index = " << iNew << " " << jNew << std::endl;
			std::cout << "time = " << T << std::endl;

			std::cout << "i = " << i << " j= " << j << std::endl;
			std::cout << "/t /t " << concentration[index[i][j-1]]  << std::endl;
			std::cout << " " << /*concentration[index[i-1 ][j]]  << */concentration[index[i][j]] /*<< concentration[index[i+1][j]]*/<< std::endl;
			std::cout << "/t /t " << concentration[index[i][j+1]] << std::endl << std::endl;


			std::cout << "/t /t " << opening[index[i][j - 1]] << std::endl;
			std::cout << " " << opening[index[i - 1][j]] << opening[index[i][j]] << opening[index[i + 1][j]] << std::endl;
			std::cout << "/t /t " << opening[index[i][j + 1]] << std::endl;

		}

		if ((iNew != i) || (jNew != j)) {
			if ((concentration[index[iNew][jNew]] + markerVolume[markerIndex] /
				(opening[index[iNew][jNew]] * wn * dx * dy) < maximumConcentration[markerIndex])
				&& (opening[index[iNew][jNew]] * wn >= proppantDiameter[markerIndex])) {
				concentration[index[i][j]] -= markerVolume[markerIndex] / (opening[index[i][j]]
					* wn * dx * dy);
				concentration[index[iNew][jNew]] += markerVolume[markerIndex] /
					(opening[index[iNew][jNew]] * wn * dx * dy);
			}
			else {
				xNew = markers[markersInCells[index[i][j]][rInd]][0];
				yNew = markers[markersInCells[index[i][j]][rInd]][1];
			}
		}
		markers[markersInCells[index[i][j]][rInd]][0] = xNew;
		markers[markersInCells[index[i][j]][rInd]][1] = yNew;
	}
}

/*!
 \details Функция перерасчета концентрации проппанта при удвоении сетки

 \param[in] markers - матрица положения маркеров
 \param[in] index - матрица индексов
 \param[in] concentration - вектор концентраций проппанта
 \param[in] opening - вектор раскрытий
 \param[in] proppantDensity - вектор плотностей проппантов
 \param[in] markerMass - вектор масс маркеров
 \param[in] wn - нормирующий множитель
*/
void recalcConcentation(
	std::vector< std::vector<double> > &markers,
	const std::vector< std::vector<size_t> > &index,
	std::vector<double> &concentration,
	std::vector<double> &opening,
	std::vector<double> &markerVolume,
	const double wn
) {
	concentration.resize(opening.size());
	std::fill(concentration.begin(), concentration.end(), 0.);
	if (markers.size() != 0) {
		for (std::size_t i = 0; i < markers.size(); ++i) {
			int cell = index[int(round(markers[i][0] / dx + double(i00)))]
				[int(round(markers[i][1] / dx + double(j00)))];
			concentration[cell] += markerVolume[markers[i][2]] / (opening[cell] * wn * dx * dx);
		}
	}
}

/*!
\details Функция выводит значения концентраций проппантов

\param[in] markers - матрица положения маркеров
\param[in] index - матрица индексов
\param[in] concentration - вектор концентраций проппанта
\param[in] opening - вектор раскрытий
\param[in] markerVolume - вектор объемов маркеров
\param[in] wn - нормирующий множитель
\param[in] i - индекс ячейки (x)
\param[in] j - индекс ячейки (y)
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
) {
	std::ostringstream output;
	if (markers.size() != 0)
	{
		std::vector< std::vector<size_t> > markersInCells;
		markersInCells.resize(opening.size());
		generateListOfMarkersInCells(
			markersInCells,
			markers,
			index
		);
		std::vector< double > propConcentration;
		double check = 0.;
		propConcentration.resize(markerVolume.size());
		for (std::size_t ind = 0; ind < markersInCells[index[i][j]].size(); ++ind) {
			++propConcentration[int(markers[markersInCells[index[i][j]][ind]][2])];
		}

		num_fluids = 1; //Пока у нас одна жидкость выводим еденицу. Далее все будет зависить от размера вектора жидкости
		num_proppants = propConcentration.size();
		for (std::size_t ind = 0; ind < propConcentration.size(); ++ind) {
			propConcentration[ind] *= markerVolume[ind] / (opening[index[i][j]] * dx * dy * wn);
			check += propConcentration[ind];
			if (i == i00)														//Если столбец центральный, то удваиваем концентрацию проппанта
				output << "\t\t  " << 2 * propConcentration[ind] << ",\n";			//Вывод концентрации проппанта для центрального столбца
			else
				output << "\t\t  " << propConcentration[ind] << ",\n";			//Вывод концентрации проппанта
		}
		if (i == i00)														//Если столбец центральный, то удваиваем концентрацию проппанта
			output << "\t\t  " << 1 - 2 * concentration[index[i][j]] << "\n";           //Вывод концентрации жидкости
		else
			output << "\t\t  " << 1 - concentration[index[i][j]] << "\n";           //Вывод концентрации жидкости
	}
	else
	{
		//output << "\t\t  " << 0 << ",\n";				  //Вывод концентрации проппанта
		output << "\t\t  " << 1 << "\n";             //Вывод концентрации жидкости
	}

	return	output.str();
}


///Блок функций для маркеров для нескольких жидкостей (?)
/*!
\details Функция рассчитывает объем маркера для
жидкости и заполняет начальную трещину маркерами

\param[in] activeElements - массив активных элементов
\param[in] index - матрица индексов
\param[in] markers - матрица положения маркеров
\param[in] fluidType - вектор типов жидкостей
\param[in] opening - вектор раскрытий
\param[in] x - вектор x-координат узлов
\param[in] y - вектор y-координат узлов
\param[in] markerVolume  - объем маркера
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
) {
	markerVolume = ai::max(opening);
	for (std::size_t i = 0; i < activeElements.size(); ++i) {
		if (opening[index[activeElements[i][0]][activeElements[i][1]]] < markerVolume)
			markerVolume = opening[index[activeElements[i][0]][activeElements[i][1]]];
	}
	//markerVolume /= 10.;
	int nM;
	std::vector< double > xCoordinate;
	std::vector< double > yCoordinate;
	for (std::size_t i = 0; i < activeElements.size(); ++i) {
		nM = int(opening[index[activeElements[i][0]][activeElements[i][1]]] / markerVolume);
		for (std::size_t j = 0; j < nM; ++j) {
			if (activeElements[i][0] == 0) {
				ai::generateRandomVector(xCoordinate, 1, epsilon,
					x[activeElements[i][0]] + 0.5 * dx);
			}
			else {
				ai::generateRandomVector(xCoordinate, 1, x[activeElements[i][0]] -
					0.5 * dx, x[activeElements[i][0]] + 0.5 * dx);
			}
			ai::generateRandomVector(yCoordinate, 1, y[activeElements[i][1]] -
				0.5 * dx, y[activeElements[i][1]] + 0.5 * dx);
			markers.push_back(std::vector<double>{xCoordinate[0], yCoordinate[0], 
				double(fluidType[0])});
		}
	}
	markerVolume *= dx * dx;
}

/*!
 \details Функция добавляет в массив маркеров новые маркеры по ходу закачки

 \param[in] markers - матрица положения маркеров
 \param[in] fluidType - вектор типов жидкостей
 \param[in] fluidInjection - расход жидкости
 \param[in] injectionIndex - номер стадии закачки
 \param[in] injectedVolume - объем закаченной жидкости
 \param[in] markerVolume  - объем маркера
*/
void addMarkersFluid(
	std::vector< std::vector<double> > &markers,
	std::vector< size_t > &fluidType,
	double fluidInjection,
	double injectionIndex,
	double &injectedVolume,
	double markerVolume
) {
	injectedVolume += fluidInjection * dt;
	int nMarkToAdd = int(injectedVolume / (markerVolume));
	if (nMarkToAdd > 0) {
		for (size_t k = 0; k < nMarkToAdd; ++k) {
			std::vector< double > xCoordinate;
			std::vector< double > yCoordinate;
			ai::generateRandomVector(xCoordinate, 1, epsilon, 0.5 * dx - epsilon);
			ai::generateRandomVector(yCoordinate, 1, -0.5 * dx + epsilon,
				0.5 * dx - epsilon);
			markers.push_back(std::vector<double>{xCoordinate[0], yCoordinate[0],
				double(fluidType[injectionIndex - 1])});
		}
		injectedVolume -= nMarkToAdd * markerVolume;
	}
}

/*!
 \details Функция определяет эффективную вязкость, индекс течения и 
 плотность жидкости в ячейке

 \param[in] markersInCells - матрица индексов маркеров в соответствующих ячейках
 \param[in] index - матрица индексов
 \param[in] markers - матрица положения маркеров
 \param[in] fluidViscosity - вектор вязкостей жидкостей
 \param[in] rheologyIndex - вектор индексов течения жидкостей
 \param[in] fluidDensity - вектор плотностей жидкостей
 \param[in] effectiveViscosity - эффективная вязкость жидкости в ячейке
 \param[in] effectiveN - эффективный индекс течения жидкости в ячейке
 \param[in] effectiveDensity - эффективная плотность жидкости в ячейке
 \param[in] i - индекс ячейки (x)
 \param[in] j - индекс ячейки (y)
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
) {
	std::vector<double> cellM;
	cellM.resize(fluidViscosity.size());
	std::fill(cellM.begin(), cellM.end(), 0);

	effectiveViscosity = 0.;
	effectiveN = 0.;
	effectiveDensity = 0.;

	if (markersInCells[index[i][j]].size() == 0) {
		effectiveViscosity = fluidViscosity[0];
		effectiveN = rheologyIndex[0];
		effectiveDensity = fluidDensity[0];
	}
	else {
		for (std::size_t ind = 0; ind < markersInCells[index[i][j]].size(); ++ind) {
			++cellM[int(markers[markersInCells[index[i][j]][ind]][2])];
		}
		for (std::size_t ind = 0; ind < cellM.size(); ++ind) {
			effectiveViscosity += cellM[ind] * fluidViscosity[ind] 
				/ markersInCells[index[i][j]].size();
			effectiveN += cellM[ind] * rheologyIndex[ind]
				/ markersInCells[index[i][j]].size();
			effectiveDensity += cellM[ind] * fluidDensity[ind]
				/ markersInCells[index[i][j]].size();
		}
	}
}

/*!
 \details Функция переноса маркеров

 \param[in] markersInCells - матрица индексов маркеров в соответствующих ячейках
 \param[in] markers - матрица положения маркеров
 \param[in] index - матрица индексов
 \param[in] elementIsActive - матрица активных элементов
 \param[in] opening - вектор раскрытий
 \param[in] pressure - вектор давлений
 \param[in] concentration - вектор концентраций проппанта
 \param[in] WPow - степень для раскрытия в уравнении Пуассона
 \param[in] PPow - степень для давления в уравнении Пуассона
 \param[in] CPow - степень для концентрации в уравнении Пуассона
 \param[in] fluidFlowRight - поток жидкости справа
 \param[in] fluidFlowBottom - поток жидкости снизу
 \param[in] pressureIJ - давление в ячейке
 \param[in] wn - нормирующий множитель
 \param[in] markerVolume - объем маркера
 \param[in] i - индекс ячейки (x)
 \param[in] j - индекс ячейки (y)
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
) {
	double fluidVelocityLeft;
	double fluidVelocityTop;
	double pressureDrop;
	double openingOnTheBorder;
	double concentrationOnTheBorder;

	double alphaX;
	double alphaY;
	double velocityX;
	double velocityY;
	double xNew;
	double yNew;
	int iNew;
	int jNew;

	fluidVelocityLeft = 0.;
	fluidVelocityTop = 0.;

	if (i == i00) {
		fluidVelocityLeft = -fluidVelocityRight;
	}
	else {
		if (elementIsActive[i - 1][j]) {
			pressureDrop = pressure[index[i - 1][j]] - pressureIJ;
			openingOnTheBorder = 0.5
				* (opening[index[i][j]] + opening[index[i - 1][j]]);
			concentrationOnTheBorder = 0.5
				* (concentration[index[i][j]] + concentration[index[i - 1][j]]);
			fluidVelocityLeft = -ai::sign(pressureDrop)
				* std::pow(openingOnTheBorder, WPow - 1.)
				* std::pow(std::abs(pressureDrop) / dx, PPow)
				/ std::pow(
					1. - concentrationOnTheBorder / 0.6, // / maximumconcentration,
					CPow
				);
		}
	}

	if (elementIsActive[i][j - 1]) {
		pressureDrop = pressure[index[i][j - 1]] - pressureIJ;
		openingOnTheBorder = 0.5
			* (opening[index[i][j]] + opening[index[i][j - 1]]);
		concentrationOnTheBorder = 0.5
			* (concentration[index[i][j]] + concentration[index[i][j - 1]]);
		fluidVelocityTop = -ai::sign(pressureDrop)
			* std::pow(openingOnTheBorder, WPow - 1.)
			* std::pow(std::abs(pressureDrop) / dy, PPow)
			/ std::pow(
				1. - concentrationOnTheBorder / 0.6,// / maximumConcentration,
				CPow
			);
	}

	for (std::size_t rInd = 0; rInd < markersInCells[index[i][j]].size(); ++rInd) {

		alphaX = std::abs(markers[markersInCells[index[i][j]][rInd]][0]
			- ((double)i - (double)i00 - 0.5) * dx) / dx;
		velocityX = (1. - alphaX) * fluidVelocityLeft + alphaX * fluidVelocityRight;
		xNew = markers[markersInCells[index[i][j]][rInd]][0] - dt * velocityX
			/ fluidViscosity[int(markers[markersInCells[index[i][j]][rInd]][2])];

		alphaY = std::abs(markers[markersInCells[index[i][j]][rInd]][1]
			- ((double)j - (double)j00 - 0.5) * dx) / dx;
		velocityY = (1. - alphaY) * fluidVelocityTop + alphaY * fluidVelocityBottom;
		yNew = markers[markersInCells[index[i][j]][rInd]][1] - dt * velocityY
			/ fluidViscosity[int(markers[markersInCells[index[i][j]][rInd]][2])];

		iNew = int(round(xNew / dx + i00));
		jNew = int(round(yNew / dy + j00));

		if (((iNew != i) || (jNew != j)) && opening[index[iNew][jNew]] < markerVolume 
			/ (dx * dx)) {
			xNew = markers[markersInCells[index[i][j]][rInd]][0];
			yNew = markers[markersInCells[index[i][j]][rInd]][1];
		}
		markers[markersInCells[index[i][j]][rInd]][0] = xNew;
		markers[markersInCells[index[i][j]][rInd]][1] = yNew;
	}
}




