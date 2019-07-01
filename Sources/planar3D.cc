/*!
 \file planar3D.cc
 \brief файл с функцией planar3D - главной функцией расчета и вызова функций
 различных режимов
 \details Данный файл содержит в себе определения основных классов,
 используемых в программе
*/

#include <vector>
#include <iostream>
#include <algorithm>

#include "nlohmann/json.hpp"
#include "ailibrary/ai.hh"

 //#define BUILD_DLL

#if defined(USE_EIGEN)
    #include "eigen/Dense"
    #include "eigen/Sparse"
    #include "eigen/Core"
#endif

#include "io.hh"
#include "mesh.hh"
#include "rhs3D.hh"
#include "ribbons.hh"
#include "automodel.hh"
#include "initialData.hh"
#include "influenceMatrix.hh"
#include "findActiveElements.hh"

#if defined(IMEX)
    #include "imex.hh"
#endif

#if defined(PROPPANT_MARKERS) || defined(FLUID_MARKERS)
    #include "markers.hh"
#endif

#if defined(BUILD_DLL)
    #include "api.h"
    #include "api_callback.h"
#endif

#ifdef OPENMP
    #include <omp.h>
    int getNumberOfThreads(){
        int counter = 0;
        #pragma omp parallel
        {
            counter = omp_get_num_threads();
        }
        return counter;
    }
#else
    int getNumberOfThreads(){return 1;}
    void omp_set_num_threads(std::size_t numberOfThreads){};
#endif

#if defined(FLUID_MARKERS) && defined(PROPPANT_MARKERS)
    #error Cannot use markers for both fluid and proppant!
#endif

int regime;

double dx;
double dy;
double dt;

double fluidInjection;
double proppantInjection;

double E;
double alpha;
double Amu;
double epsilon;

size_t i00;
size_t j00;

/*!
 \details Функция возвращает название текущего режима распространения трещины
 \return Название режима
*/
std::string regimeName(){
    switch(regime){
        case LEAK_OFF:
            return "leak-off dominated";
        case TOUGHNESS:
            return "toughness dominated";
        case VISCOSITY:
            return "viscosity dominated";
        default:
            return "unknown";
    }
}

/*!
 \details Функция выводит текстовое лого программы
*/
void printLogo(){
    std::cout << " ____  _                       _____ ____     ____ _     "
        << "___ " << std::endl;
    std::cout << "|  _ \\| | __ _ _ __   __ _ _ _|___ /|  _ \\   / ___| |   "
        << "|_ _|" << std::endl;
    std::cout << "| |_) | |/ _` | '_ \\ / _` | '__||_ \\| | | | | |   | |    "
        << "| | " << std::endl;
    std::cout << "|  __/| | (_| | | | | (_| | |  ___) | |_| | | |___| |___ | "
        << "| " << std::endl;
    std::cout << "|_|   |_|\\__,_|_| |_|\\__,_|_| |____/|____/   "
        << "\\____|_____|___|" << std::endl;
    std::cout << std::endl << "Developed by REC \"Gazpromneft-Polytech\"."
        << std::endl << std::endl;
}


int calculateEverything(
    std::vector<double> &opening,
    std::vector<double> &pressure,
    std::vector<double> &concentration,
    #if defined(USE_EIGEN)
    Eigen::VectorXd &openingNew,
    #else
    std::vector<double> &openingNew,
    #endif
    std::vector< std::vector<double> > &distances,
    std::vector< std::vector<double> > &velocities,
    std::vector< std::vector<double> > &influenceMatrix,
    #if defined(USE_EIGEN)
    Eigen::MatrixXd &partialInfluenceMatrix,
    #else
    std::vector< std::vector<double> > &partialInfluenceMatrix,
    #endif
    std::vector<double> &zeroVectorXY,
    std::vector< std::vector<double> > &zeroMatrixXY,
    std::vector< std::vector<Cell> > &mesh,
    std::vector<double> x,
    std::vector<double> y,
    std::vector<Ribbon> &ribbons,
    std::vector< std::vector<std::size_t> > &index,
    std::vector< std::vector<std::size_t> > &activeElements,
    std::vector< std::vector<bool> > &elementIsActive,
    std::vector<double> &activationTime,
    std::vector< std::vector<double> > &injection,
    std::vector< std::vector<double> > &layers,
    std::vector<double> &stress,
    std::vector<double> &leakOff,
    std::vector<double> &toughness,
    std::vector<double> &flatYoungsModulus,
    std::vector<double> &proppantDensity,
    std::vector<double> &proppantDiameter,
    std::vector<double> &maximumConcentration,
    std::vector<size_t> &proppantType,
    std::vector<size_t> &fluidType,
    std::vector<double> &fluidDensity,
    std::vector<double> &fluidViscosity,
    std::vector<double> &rheologyIndex,
    double &n,
    double &mu,
    double &wn,
    const double timeScale,
    const double timeStep,
    std::vector< std::vector<double> > &fracture,
    double &length,
    double &height,
    double &initialRadius,
    double &axMax,
    int &meshScalingCounter,
    bool &meshIsNotExhausted,
    double &T,
    const double T0,
    std::size_t &injectionIndex,
    double &timeToChangeInjection,
    std::size_t stepToCheck,
    const double modelingTime,
    const bool considerElasticModulusContrast,
    const bool saveSteps,
    const bool runningFromGUI,
    #if defined(BUILD_DLL)
    const bool override,
    DLL_Param &DLL_Parametrs
    #else
    const bool override
    #endif
){
    double LastEmit = 0.;
    alpha = 2. / (n + 2.);
    double BAlpha = 0.25 * alpha * tan(0.5 * M_PI - M_PI * (1 - alpha));
    Amu = std::pow(BAlpha * (1 - alpha), -0.5 * alpha);

    std::size_t globalStep = (std::size_t) T / dt;
    std::size_t step = 0;
    double savedTime = T;
    std::size_t savingStep = 1;
    std::size_t stepToSave = (std::size_t) 1. / dt;
    std::size_t xSize = index.size();
    std::size_t ySize = index[0].size();

    double dMin1 = std::sqrt(std::pow(0.5 * dx, 2)
        + std::pow(1.5 * dy, 2));
    double dCenter1 = dx;
    double dMax1 = std::sqrt(std::pow(0.5 * dx, 2)
        + std::pow(0.5 * dy, 2));
    double dMin2 = std::sqrt(std::pow(1.5 * dx, 2)
        + std::pow(1.5 * dy, 2));
    double dCenter2 = std::sqrt(2) * dx;

    std::vector<double> dWdt;
    std::vector<double> dCdt;
    std::vector<double> dTdt;
    std::vector<double> savedDistances(ribbons.size(), 0.);
    // std::vector<double> temperature(activeElements.size(), 0.);

    for(size_t k = 0; k < ribbons.size(); ++k){
        const size_t i = ribbons[k].i;
        const size_t j = ribbons[k].j;

        savedDistances[k] = distances[i][j];
    }

    /////////////////////////////////////////////////////////////////////////////////////////
    /// Задание маркеров
    /////////////////////////////////////////////////////////////////////////////////////////
    /// Задание маркеров для проппантов
    #if defined(PROPPANT_MARKERS)
    std::vector< std::vector<double> > markers;
    std::vector<double> markerMass;
    std::vector<double> markerVolume;
    double injectedMass = 0.;
    defineMarkerMass(
        injection,
        markerMass,
        markerVolume,
        proppantDensity,
        proppantDiameter,
        maximumConcentration,
        proppantType,
        wn
    );
    #endif


#if !defined(PROPPANT_MARKERS)
	std::vector< std::vector<double> > markers;
	std::vector<double> markerMass;
	std::vector<double> markerVolume;
	double injectedMass = 0.;
	/*defineMarkerMass(
		injection,
		markerMass,
		markerVolume,
		proppantDensity,
		proppantDiameter,
		maximumConcentration,
		proppantType,
		wn
	);*/
#endif

    ///Задание маркеров для жидкостей
    #if defined(FLUID_MARKERS)
    std::vector< std::vector<double> > markersFluid;
    double markerFluidVolume;
    double injectedVolume = 0.;
    defineMarkerVolumeAndFillCrackWithMarkers(
        activeElements,
        index,
        markersFluid,
        fluidType,
        opening,
        x,
        y,
        markerFluidVolume
    );
    #endif
    /////////////////////////////////////////////////////////////////////////////////////////

    #if defined(IMEX)
        double imexD;
        double imexCoeff;

        std::vector<double> B(xSize * ySize);
        std::vector<double> imexOpening;
    #endif

    #if !defined(BUILD_DLL)
    std::cout << std::endl << "Starting calculations..." << std::endl;
    if(runningFromGUI){
        std::cout << "Progress: 0.0" << std::endl;
    }else{
        ai::showProgressBar(0.);
    }
    #endif

    #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
    std::vector<double> timeMeasurements(9, 0.);

    std::chrono::high_resolution_clock::time_point time0;
    std::chrono::high_resolution_clock::time_point time1;
    std::chrono::high_resolution_clock::time_point time2;
    std::chrono::high_resolution_clock::time_point time3;
    std::chrono::high_resolution_clock::time_point time4;
    std::chrono::high_resolution_clock::time_point time5;
    std::chrono::high_resolution_clock::time_point time6;
    std::chrono::high_resolution_clock::time_point time7;
    #endif

    auto startTime = ai::time();

    bool pauseCallback = false;	//Флаг того,  что на паузе мы послали кэлбак

    bool flag = true;

	std::cout << "start while" << std::endl;
	std::cout << "modelingTime: " << modelingTime << ">= T: " << T << "&& meshIsNotExhausted: " << meshIsNotExhausted << std::endl;

    double maxOpen  = 0.;
    double maxSpeed = 0.;

    while(modelingTime >= T && meshIsNotExhausted){
    maxOpen = wn*ai::max(opening);

    maxSpeed = ai::max(velocities)>0. ? ai::max(velocities) : 1.;


    #if defined(BUILD_DLL)
        //Обработчик нажатия на кнопку RUN в интерфейсе
        if (DLL_Parametrs.Dll_State == DLL_Parametrs.Running)
        {
            pauseCallback = false;
        }
        //Обработчик нажатия на кнопку PAUSE в интерфейсе
        if (DLL_Parametrs.Dll_State == DLL_Parametrs.Paused)
        {
            if (!pauseCallback)
            {
                double fluidEfficiency;
                double proppantEfficiency;

                /// Вычисляем эффективность
                calculateEfficiency(
                    fluidEfficiency,
                    proppantEfficiency,
                    injection,
                    opening,
                    wn,
                    concentration,
                    index,
                    T,
                    T0
                );

                if (std::isnan(fluidEfficiency)) {
                    fluidEfficiency = 0.;
                    proppantEfficiency = 0.;
                }

                /// Сохраняем сжимающие напряжения в центральном слое
                double nominalStress = 0.;
                for (std::size_t i = 0; i < layers.size(); ++i) {
                    if (0. > layers[i][0] * layers[i][1]) {
                        nominalStress = layers[i][2] * std::pow(10., -6);
                    }
                }

				DLL_Parametrs.J_String = ExportJson(
					opening,
					wn,
					pressure,
					concentration,
					markers,
					markerVolume,
					x,
					y,
					dx,
					dy,
					i00,
					j00,
					index,
					T,
					timeScale,
					fluidEfficiency,
					fluidDensity,
					proppantDensity,
					fluidInjection,
					proppantInjection,
					nominalStress,
					DLL_Parametrs
				);

                Solver::DataCallback dataCallback = Solver::Callback::GetDataCallback();
//              DLL_Parametrs.J_String = ExportJson(opening, wn, pressure, concentration, x, y, dx, dy, i00, j00, index, T, timeScale, fluidEfficiency, fluidDensity, proppantDensity, fluidInjection, proppantInjection, nominalStress, DLL_Parametrs).c_str();
                ai::saveLine("pause.json", DLL_Parametrs.J_String);

                dataCallback(DLL_Parametrs.J_String.c_str());


            }
            pauseCallback = true;
            continue;
        }
        //Обработчик нажатия на кнопку STOP в интерфейсе
        if (DLL_Parametrs.Dll_State == DLL_Parametrs.Idle)
        {
            ai::saveLine("break", DLL_Parametrs.J_String.c_str());
            break;
        }
        #endif
        // ai::printMarker(500);

        if(-1. == timeStep){
            #if defined(FLUID_MARKERS)
            dt = 0.05 * ai::min(fluidViscosity) * injection[0][5]
                / (E * std::pow(wn * ai::max(opening) / dx, 3));
            #else
            dt = 0.05 * mu / std::pow(maxOpen , 2. + n) /
                        std::pow(maxSpeed, 1. - n ) / (E  / std::pow(dx , 3));
            #endif
        }
        #if defined(IMEX)
        dt *= 1.7;
        imexD = 0.01 * 0.05 * std::pow(1000.*wn*ai::max(opening), 3);
        imexCoeff = (dx * dy / (imexD * dt) + 4.);
        #endif

        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        time0 = ai::time();
        #endif
        calculatePressure(
            pressure,
            index,
            activeElements,
            partialInfluenceMatrix,
            openingNew,
            stress
        );
        // ai::printMarker();
        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        timeMeasurements[0] += ai::duration(time0, "us");
        #endif

        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        auto time1 = ai::time();
        #endif
/////////////////////////////////////////////////////////////////
        #if defined(PROPPANT_MARKERS)
        /// Добавление маркеров для проппанта
        if (proppantInjection > epsilon) {
            addMarkers(
                markers,
                markerMass,
                proppantType,
                proppantDensity,
                concentration[index[i00][j00]],
                opening[index[i00][j00]],
                fluidInjection,
                proppantInjection,
                injectionIndex,
                injectedMass,
                wn
            );
        }
        #endif

        #if defined(FLUID_MARKERS)
        ///Добавление маркеров для жидкостей
        addMarkersFluid(
            markersFluid,
            fluidType,
            fluidInjection,
            injectionIndex,
            injectedVolume,
            markerFluidVolume
        );
        #endif
        /////////////////////////////////////////////////////////////////

        #if defined(PROPPANT_MARKERS)
        calculateOpeningAndConcentrationSpeedsDP(
            dWdt,
            opening,
            concentration,
            pressure,
            leakOff,
            index,
            activeElements,
            elementIsActive,
            markers,
            activationTime,
            proppantDensity,
            proppantDiameter,
            maximumConcentration,
            markerMass,
            markerVolume,
            fluidDensity[0],
            T,
            n,
            wn
        );
        #elif defined(FLUID_MARKERS)
        calculateOpeningAndConcentrationSpeedsDF(
            dWdt,
            dCdt,
            opening,
            concentration,
            pressure,
            leakOff,
            index,
            activeElements,
            elementIsActive,
            markersFluid,
            activationTime,
            fluidViscosity,
            rheologyIndex,
            fluidDensity,
            proppantDensity[0],
            proppantDiameter[0],
            maximumConcentration[0],
            markerFluidVolume,
            T,
            n,
            wn
        );
        #else
        // ai::saveVector("opening",opening );
        // ai::printMarker(999999);
        calculateOpeningAndConcentrationSpeeds(
            dWdt,
            dCdt,
            opening,
            concentration,
            pressure,
            leakOff,
            index,
            activeElements,
            elementIsActive,
            activationTime,
            fluidDensity[0],
            proppantDensity[0],
            proppantDiameter[0],
            maximumConcentration[0],
            T,
            n,
            wn
        );
        // ai::printMarker();
        #endif
        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        timeMeasurements[1] += ai::duration(time1, "us");
        #endif

        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        auto time2 = ai::time();
        #endif
        if(0 == regime){
            double maxDeltaDistance = epsilon;

            for(size_t k = 0; k < ribbons.size(); ++k){
                const size_t i = ribbons[k].i;
                const size_t j = ribbons[k].j;

                distances[i][j] = 0.25 *  sqrt(0.5 * M_PI) * E
                    * opening[index[i][j]] / toughness[j];

                if(epsilon > distances[i][j]){
                    distances[i][j] = epsilon;
                }else{
                    distances[i][j] *= distances[i][j];
                }

                maxDeltaDistance = ai::max(
                    distances[i][j] - savedDistances[k],
                    maxDeltaDistance
                );
            }

            if(epsilon < maxDeltaDistance){
                stepToCheck = std::round(
                    dx * (T - savedTime) / (20. * dt * maxDeltaDistance)
                );
            }else{
                stepToCheck = 2000;
            }
        }else{
            calculateFrontVelocity(
                velocities,
                leakOff,
                index,
                ribbons,
                distances,
                opening,
                n
            );

            for(size_t i = 0; i < xSize; ++i){
                for(size_t j = 0; j < ySize; ++j){
                    distances[i][j] += velocities[i][j] * dt;
                }
            }
        }
        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        timeMeasurements[2] += ai::duration(time2, "us");
        #endif

        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        auto time3 = ai::time();
        #endif
        #if defined(IMEX)
        std::fill(B.begin(), B.end(), 0.);
        for(std::size_t k = 0; k < activeElements.size(); ++k){
            std::size_t i = activeElements[k][0];
            std::size_t j = activeElements[k][1];

            double lapl = opening[index[i + 1][j]] + opening[index[i][j + 1]]
                + opening[index[i][j - 1]];

            if(i00 == i){
                lapl += opening[index[i + 1][j]];
            }else{
                lapl += opening[index[i - 1][j]];
            }

            B[index[i][j]] = opening[index[i][j]] * imexCoeff
                + dWdt[index[i][j]] * dx * dy / imexD - lapl;
        }
        conjugateGradient(B, imexOpening, imexCoeff, xSize, ySize);
        #endif
        for(std::size_t i = 0; i < opening.size(); ++i){
            #if defined(PROPPANT_MARKERS)
                concentration[i] = concentration[i] * opening[i];
            #else
                concentration[i] = concentration[i] * opening[i] + dCdt[i] * dt;
            #endif

            #if defined(IMEX)
                opening[i] = ai::max(imexOpening[i], 0.);
            #else
                opening[i] = opening[i] + dWdt[i] * dt;
                // opening[i] = ai::max(opening[i] + dWdt[i] * dt, 0.);
            #endif

            if(epsilon < opening[i]){
                concentration[i] /= opening[i];
            }else{
                if(epsilon < concentration[i]){
                    concentration[i] = 0.;
                    #if defined(DEBUG)
                    ai::printLine("Warning: check is required.");
                    #endif
                }
            }
        }

        #if defined(DEBUG) && !defined(BUILD_DLL)
        if(std::isnan(opening[index[i00][j00]])){
            ai::printLine(
                ai::string("Planar3D is dead. Time = ") + ai::string(T)
            );
            break;
        }
        #endif

        for(std::size_t k = 0; k < activeElements.size(); ++k){
            const std::size_t i = activeElements[k][0];
            const std::size_t j = activeElements[k][1];

            openingNew[k] = opening[index[i][j]];
        }
        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        timeMeasurements[3] += ai::duration(time3, "us");
        #endif



        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        auto time4 = ai::time();
        #endif
        if(step >= stepToCheck){
            #if !defined(BUILD_DLL) && defined(DEBUG)
            if(runningFromGUI){
                std::cout << "Progress: " << (T - T0) / (modelingTime - T0)
                    << std::endl;
            }else{
                ai::showProgressBar((T - T0) / (modelingTime - T0));
            }
            #endif
            ai::showProgressBar((T - T0) / (modelingTime - T0));
            step = 0;


            if(0 != regime){
                stepToCheck = std::round(
                    dx / (20. * dt * ai::max(velocities))
                );
            }

            const size_t savedSize = activeElements.size();

            std::vector<Ribbon> oldRibbons = ribbons;

            for(size_t k = 0; k < ribbons.size(); ++k){
                const size_t i = ribbons[k].i;
                const size_t j = ribbons[k].j;

                if(2 >= j || ySize - 2 <= j || xSize - 2 <= i){
                    /// При исчерпании сетки может сделать одно из трёх:
                    /// 1) закончить вычисления с предупреждением,
                    /// 2) изменить масштаб сетки,
                    /// 3) достроить сетку.

                    /// 1) Выход из цикла
                    // meshIsNotExhausted = false;
					// std::cout << "1. Mesh is NotExhausted. End calculate" << std::endl;
                    // break;

                    int returnCode;
                    // делаем удвоение при заканчивающейся сетки
                    meshScalingCounter = 1;

                    if(0 < meshScalingCounter){
                        /// 2) Масштабирование
						std::cout << "2. Scalling. Continue calculate" << std::endl;
						returnCode = scaleMesh(
                            x,
                            y,
                            mesh,
                            index,
                            activeElements,
                            elementIsActive,
                            activationTime,
                            ribbons,
                            distances,
                            opening,
                            concentration,
                            layers,
                            flatYoungsModulus,
                            stress,
                            leakOff,
                            toughness,
                            zeroVectorXY,
                            zeroMatrixXY,
                            xSize,
                            ySize,
                            initialRadius,
                            axMax,
                            dMin1,
                            dCenter1,
                            dMax1,
                            dMin2,
                            dCenter2,
                            dx,
                            dy,
                            dt,
                            mu,
                            n,
                            wn,
                            timeScale,
                            T,
                            T0,
                            E,
                            considerElasticModulusContrast,
                            #if defined(PROPPANT_MARKERS)
                            markers,
                            markerVolume,
                            #endif
                            meshScalingCounter
                        );
                    }else{
                        /// 3) Достроение сетки
                        /// Добавляем по 10 элементов вдоль каждой оси
						std::cout << "3. Adding 10 elements. Continue calculate" << std::endl;
						returnCode = completeMesh(
                            xSize + 10,
                            x,
                            y,
                            i00,
                            j00,
                            mesh,
                            ribbons,
                            index,
                            opening,
                            pressure,
                            concentration,
                            distances,
                            activeElements,
                            elementIsActive,
                            activationTime,
                            layers,
                            flatYoungsModulus,
                            stress,
                            leakOff,
                            toughness,
                            wn,
                            timeScale,
                            mu,
                            n,
                            xSize,
                            ySize,
                            considerElasticModulusContrast,
                            influenceMatrix
                        );
                    }
                    // oldRibbons = ribbons; todo ???
                    //
                    /// todo savedDistances должны быть в scaleMesh
                    /// и completeMesh!
                    //
                    if(0 != returnCode){
						std::cout << "return from point calculate1";
                        return returnCode;
                    }
                }

                const double d = distances[i][j];

                if(0 < i){
                    findNewRibbons(i - 1, j, d, dMin1, dCenter1, n, opening, leakOff, toughness, mesh,
                        index, activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                        velocities, T, opening[index[i - 1][j]]);

                    findNewRibbons(i - 1, j - 1, d, dMin2, dCenter2, n, opening, leakOff, toughness, mesh,
                        index, activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                        velocities, T, opening[index[i - 1][j - 1]]);

                    findNewRibbons(i - 1, j + 1, d, dMin2, dCenter2, n, opening, leakOff, toughness, mesh,
                        index, activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                        velocities, T, opening[index[i - 1][j + 1]]);
                }

                findNewRibbons(i + 1, j, d, dMin1, dCenter1, n, opening, leakOff, toughness, mesh, index,
                    activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                    velocities, T, opening[index[i + 1][j]]);

                findNewRibbons(i, j - 1, d, dMin1, dCenter1, n, opening, leakOff, toughness, mesh, index,
                    activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                    velocities, T, opening[index[i][j - 1]]);

                findNewRibbons(i, j + 1, d, dMin1, dCenter1, n, opening, leakOff, toughness, mesh, index,
                    activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                    velocities, T, opening[index[i][j + 1]]);

                findNewRibbons(i + 1, j - 1, d, dMin2, dCenter2, n, opening, leakOff, toughness, mesh,
                    index, activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                    velocities, T, opening[index[i + 1][j - 1]]);

                findNewRibbons(i + 1, j + 1, d, dMin2, dCenter2, n, opening, leakOff, toughness, mesh,
                    index, activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                    velocities, T, opening[index[i + 1][j + 1]]);

                if(dMin2 < d){
                    mesh[i][j].type = CHANNEL;
                }
            }

            for(size_t k = 0; k < ribbons.size(); ++k){
                const size_t i = ribbons[k].i;
                const size_t j = ribbons[k].j;

                if(RIBBON != mesh[i][j].type){
                    ribbons.erase(ribbons.begin() + k);
                    distances[i][j] = 0;
                }
            }

            if(savedSize != activeElements.size() && meshIsNotExhausted){
                buildPartialInfluenceMatrix(influenceMatrix, activeElements,
                    opening, openingNew, partialInfluenceMatrix, index
                );
            }

            // Сохраняем параметры трещины
            calculateCrackGeometry(mesh, distances, length, height);

            // ai::printMarker();

            fracture.push_back(
                std::vector<double>{
                    T,
                    1000 * wn * opening[index[i00][j00]],
                    pressure[index[i00][j00]],
                    length,
                    height,
                    length / height
                }
            );

            if(0 == regime){
                savedDistances.resize(ribbons.size());

                for(size_t k = 0; k < ribbons.size(); ++k){
                    const size_t i = ribbons[k].i;
                    const size_t j = ribbons[k].j;

                    savedDistances[k] = distances[i][j];
                }

                savedTime = T;
            }
        }
        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        timeMeasurements[4] += ai::duration(time4, "us");
        #endif

        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        auto time5 = ai::time();
        #endif


        #if defined(BUILD_DLL)
		////////////////////////////////////////////////////////////////////////////////////////////////
		//Обработка создания выгузки значений расчета по emit'у
		////////////////////////////////////////////////////////////////////////////////////////////////
		if (T * 60. - DLL_Parametrs.Emit_time >= LastEmit){
		std::cout << T * 60. - DLL_Parametrs.Emit_time << "Last emit " << LastEmit << std::endl;
			//Смещаем точку последнего имита
			LastEmit = LastEmit + DLL_Parametrs.Emit_time;
			double fluidEfficiency;
			double proppantEfficiency;

			/// Вычисляем эффективность
			calculateEfficiency(
				fluidEfficiency,
				proppantEfficiency,
				injection,
				opening,
				wn,
				concentration,
				index,
				T,
				T0
			);

			if (std::isnan(fluidEfficiency)) {
				fluidEfficiency = 0.;
				proppantEfficiency = 0.;
			}

			/// Сохраняем сжимающие напряжения в центральном слое
			double nominalStress = 0.;
			for (std::size_t i = 0; i < layers.size(); ++i) {
				if (0. > layers[i][0] * layers[i][1]) {
					nominalStress = layers[i][2] * std::pow(10., -6);
				}
			}

			DLL_Parametrs.J_String = ExportJson(
				opening,
				wn,
				pressure,
				concentration,
				markers,
				markerVolume,
				x,
				y,
				dx,
				dy,
				i00,
				j00,
				index,
				T,
				timeScale,
				fluidEfficiency,
				fluidDensity,
				proppantDensity,
				fluidInjection,
				proppantInjection,
				nominalStress,
				DLL_Parametrs
			).c_str();

			std::string filename = "emit_" + std::to_string(LastEmit) + ".json";
			ai::saveLine(filename, DLL_Parametrs.J_String);
			std::cout << "emit" << ai::string(LastEmit) << std::endl;

			Solver::DataCallback dataCallback = Solver::Callback::GetDataCallback();

			dataCallback(DLL_Parametrs.J_String.c_str());
		}
        #endif
		//       #else
		if (
         saveSteps //&& T > 340.74
					  && 0 == globalStep % stepToSave
			) {
			saveData(
				ai::string("./Results/Opening/")
				+ ai::prependNumber(savingStep, 3),
				opening,
				index,
				1000. * wn
			);
			saveData(
				ai::string("./Results/Pressure/")
				+ ai::prependNumber(savingStep, 3),
				pressure,
				index
			);
			saveData(
				ai::string("./Results/Concentration/")
				+ ai::prependNumber(savingStep, 3),
				concentration,
				index
			);

			++savingStep;
		}

#if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        timeMeasurements[5] += ai::duration(time5, "us");
        #endif

        T += dt;
        ++step;
        ++globalStep;

        // Меняем скорость закачки в соответсвии с планом

        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        auto time6 = ai::time();
        #endif
        if(T > timeToChangeInjection){
            // Начало экспериментального блока

            /// \todo убрать
            // Это вариант для мгновенной замены жидкости
            // Один из них должен быть закомментирован
            // n = injection[injectionIndex][4];
            // mu = injection[injectionIndex][5];
            // Это вариант для получения эффективной жидкости
            // Один из них должен быть закомментирован
            n = (
                    (double) injectionIndex * n + injection[injectionIndex][4]
                )
                / ((double) injectionIndex + 1.);
            mu = (
                    (double) injectionIndex * mu + injection[injectionIndex][5]
                )
                / ((double) injectionIndex + 1.);

            // Ниже идёт учёт измения параметров из-за новой жидкости.
            // Это только для режима вязкости!
            // Может требоваться пересчёт шага по времени!
            fluidInjection /= timeScale / (wn * 60.);
            for(std::size_t i = 0; i < stress.size(); ++i){
                stress[i] *= (wn * E);
            }
            for(std::size_t i = 0; i < leakOff.size(); ++i){
                leakOff[i] /= std::sqrt(timeScale) / wn;
            }
            for(std::size_t i = 0; i < toughness.size(); ++i){
                toughness[i] *= std::pow(
                    mu * E * std::pow(E / timeScale, n),
                    1. / (n + 2.)
                );
            }
            for(std::size_t i = 0; i < injection.size(); ++i){
                injection[i][2] /= timeScale / (wn * 60.);
            }
            wn = std::pow(
                mu * 2. * std::pow(2. * (2. * n + 1.) / n, 1. / n)
                / (E * std::pow(timeScale, n)), 1. / (n + 2.)
            );
            for(std::size_t i = 0; i < stress.size(); ++i){
                stress[i] /= (wn * E);
            }
            for(std::size_t i = 0; i < leakOff.size(); ++i){
                leakOff[i] *= std::sqrt(timeScale) / wn;
            }
            for(std::size_t i = 0; i < toughness.size(); ++i){
                toughness[i] /= std::pow(
                    mu * E * std::pow(E / timeScale, n),
                    1. / (n + 2.)
                );
            }
            for(std::size_t i = 0; i < injection.size(); ++i){
                injection[i][2] *= timeScale / (wn * 60.);
            }
            alpha = 2. / (n + 2.);
            BAlpha = 0.25 * alpha * tan(0.5 * M_PI - M_PI * (1 - alpha));
            Amu = std::pow(BAlpha * (1 - alpha), -0.5 * alpha);
            fluidInjection = injection[injectionIndex][2];
            proppantInjection = injection[injectionIndex][3];
            // Конец экспериментального блока

            ++injectionIndex;

            if(injection.size() > injectionIndex){
                timeToChangeInjection = injection[injectionIndex][0];
            }else{
                timeToChangeInjection = modelingTime + 10 * dt;
            }

           // #if !defined(BUILD_DLL)
            std::string line("Time = ");
            line += ai::string(injection[injectionIndex - 1][0])
                + ai::string(". Injection was changed.");
            ai::printLine(line);
            ai::printLine(
                ai::string("Current fluid: n = ") + ai::string(n)
                + ai::string(", mu = ") + ai::string(mu)
            );
            //#endif
        }
        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        timeMeasurements[6] += ai::duration(time6, "us");
        #endif

        // Масштабируем сетку при выполнении критерия

        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        auto time7 = ai::time();
        #endif

        if(
            0 < meshScalingCounter
            && length > 4 * initialRadius
            && height > 4 * initialRadius
            // && length >= 0.9 * dx *
        ){
            std::cout<<"Mesh scaling..";
            int returnCode = scaleMesh(
                x,
                y,
                mesh,
                index,
                activeElements,
                elementIsActive,
                activationTime,
                ribbons,
                distances,
                opening,
                concentration,
                layers,
                flatYoungsModulus,
                stress,
                leakOff,
                toughness,
                zeroVectorXY,
                zeroMatrixXY,
                xSize,
                ySize,
                initialRadius,
                axMax,
                dMin1,
                dCenter1,
                dMax1,
                dMin2,
                dCenter2,
                dx,
                dy,
                dt,
                mu,
                n,
                wn,
                timeScale,
                T,
                T0,
                E,
                considerElasticModulusContrast,
                #if defined(PROPPANT_MARKERS)
                markers,
                markerVolume,
                #endif
                meshScalingCounter
            );
            std::cout<<"at Time : "<<T<<std::endl;
            if(0 != returnCode){
				std::cout << "return from point calculate2";
				return returnCode;
            }
        }
        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        timeMeasurements[7] += ai::duration(time7, "us");
        #endif
    }

    auto finishTime = ai::time();

    #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
    timeMeasurements[8] = ai::duration(startTime, finishTime, "us");
    #endif

    #if !defined(BUILD_DLL)
    if(meshIsNotExhausted){
        if(runningFromGUI){
            std::cout << "Progress: 1.0" << std::endl;
        }else{
            ai::showProgressBar(1.);
        }

        std::cout << std::endl;
    }else{
        if(runningFromGUI){
            std::cout << "Progress: " << (T - T0) / (modelingTime - T0) << std::endl;
        }else{
            ai::showProgressBar((T - T0) / (modelingTime - T0));
        }

        std::cout << std::endl;
        std::cout << "Attention. Mesh is exhausted! "
            << "Cannot continue calculation." << std::endl;
        std::cout << "Finished at time = " << T << "." << std::endl;
    }

    std::cout << "Time used: " << ai::duration(startTime, finishTime, "s")
        << "s" << std::endl;
    #endif

    #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
    ai::saveVector("./Results/times", timeMeasurements);

    std::cout << "Measured times in us." << std::endl;
    #endif


	/////////////////////////////////////////////////////////////////////////////////////////
	//Выдаем окончательный вид трещины. Нужно при однородной среде и одной закачке (автомоделка)
	/////////////////////////////////////////////////////////////////////////////////////////
	double fluidEfficiency;
	double proppantEfficiency;

	/// Вычисляем эффективность
	calculateEfficiency(
		fluidEfficiency,
		proppantEfficiency,
		injection,
		opening,
		wn,
		concentration,
		index,
		T,
		T0
	);

	if (std::isnan(fluidEfficiency)) {
		fluidEfficiency = 0.;
		proppantEfficiency = 0.;
	}

	/// Сохраняем сжимающие напряжения в центральном слое
	double nominalStress = 0.;
	for (std::size_t i = 0; i < layers.size(); ++i) {
		if (0. > layers[i][0] * layers[i][1]) {
			nominalStress = layers[i][2] * std::pow(10., -6);
		}
	}

    // DLL_Param DLL_Parametrs;
    #if defined(BUILD_DLL)

	DLL_Parametrs.J_String = ExportJson(
		opening,
		wn,
		pressure,
		concentration,
		markers,
		markerVolume,
		x,
		y,
		dx,
		dy,
		i00,
		j00,
		index,
		T,
		timeScale,
		fluidEfficiency,
		fluidDensity,
		proppantDensity,
		fluidInjection,
		proppantInjection,
		nominalStress,
		DLL_Parametrs
	).c_str();

	std::string filename = "emit_" + std::to_string(LastEmit) + ".json";
	ai::saveLine(filename, DLL_Parametrs.J_String);
	std::cout << "emit" << ai::string(LastEmit) << std::endl;

	Solver::DataCallback dataCallback = Solver::Callback::GetDataCallback();

	dataCallback(DLL_Parametrs.J_String.c_str());
	/////////////////////////////////////////////////////////////////////////////////////////
	//Конец. Выдаем окончательный вид трещины. Нужно при однородной среде и одной закачке (автомоделка)
	/////////////////////////////////////////////////////////////////////////////////////////

    #endif

	std::cout << "End Time: " << T << std::endl;

	std::cout << "return from point calculate 3";

	return 0;
}

/*!
 \details Основная расчетная функция, считающая распространение трещины в
 заданном режиме по планарной модели и сохраняющая данные о результатах расчёта
 в виде текстовых файлов

 \param[in] modelingTime - время окончания расчета
 \param[in] timeStep - ограничение на максимальный шаг по времени в минутах
 \param[in] timeScale - масштаб времени (число реальных секунд в расчетном
 времени)
 \param[in] cellSize - длина стороны квадратной ячейки в метрах
 \param[in] meshSize - число ячеек в расчётной области вдоль горизонтальной оси
 \param[in] meshScalingCounter - сколько раз допускается масштабировать сетку
 \param[in] pathToLayersFile - путь (относительный или абсолютный) к
 расположению файла слоёв
 \param[in] pathToInjectionFile - путь (относительный или абсолютный) к
 расположению файла закачки
 \param[in] pathToImportFolder - путь (относительный или абсолютный) к
 расположению директории с данными для импорта
 \param[in] considerElasticModulusContrast - флаг подключения модуля учёта
 влияния слоистости по упругим модулям
 \param[in] saveSteps – флаг сохранения значений промежуточных результатов
 \param[in] runningFromGUI – флаг запуска программы через графический интерфейс
 \param[in] numberOfThreads – число параллельных потоков
 \return Код завершения
*/


int planar3D(
    #if defined(BUILD_DLL)
    DLL_Param &DLL_Parametrs,
    std::vector< std::vector<double> > &layers,
    std::vector< std::vector<double> > &injection,
    bool &initialized
    //Флаг, что инициализация (массштабирование) переменных прошла. от него нужно избавляться
    #else
    double modelingTime,
    double timeStep,
    const double timeScale,
    double cellSize,
    double meshSize,
    int meshScalingCounter,
    std::string pathToLayersFile,
    std::string pathToInjectionFile,
    const std::string pathToImportFolder,
    const bool considerElasticModulusContrast,
    const bool saveSteps,
    const bool runningFromGUI,
    const std::size_t numberOfThreads,
    const bool override = false
    #endif
){
    #if defined(BUILD_DLL)
    double modelingTime = DLL_Parametrs.modelingTime / 60.;
    double timeStep = -1.;
    double timeScale = 60.;
	double cellSize = -1.;
	double meshSize = 61.;
	int meshScalingCounter = 2;
	//Управление размером сетки из интерфейса
	if (DLL_Parametrs.IS_CLOSURE)
	{
		cellSize = DLL_Parametrs.dx;
		meshSize = DLL_Parametrs.nx;
		meshScalingCounter = (int)DLL_Parametrs.AZ;
	}
		std::cout << "cellSize=" << cellSize << "  meshSize=" << meshSize << "  Scaling:" << meshScalingCounter << std::endl;
    const bool considerElasticModulusContrast = false;
    const bool saveSteps = false;
    const bool runningFromGUI = false;
    const std::size_t numberOfThreads = 1;
    const bool override = true;
    #else
    bool initialized = false;

    if(!runningFromGUI){
        printLogo();
    }

    // Проверяем соответствие структуры рабочей директории
    if(!ai::folderExists("./Results")){
        std::cerr << "Cannot find './Results', create directory and restart."
            << std::endl;

        return 21;
    }

    if(saveSteps){
        if(!ai::folderExists("./Results/Concentration")){
            std::cerr << "Cannot find './Results/Concentration', create "
                << "directory and restart." << std::endl;

            return 21;
        }
        if(!ai::folderExists("./Results/Opening")){
            std::cerr << "Cannot find './Results/Opening', create directory "
                << "and restart." << std::endl;

            return 21;
        }
        if(!ai::folderExists("./Results/Pressure")){
            std::cerr << "Cannot find './Results/Pressure', create directory "
                << "and restart." << std::endl;

            return 21;
        }
    }

    if(
        std::string() != pathToImportFolder
        && !ai::folderExists(pathToImportFolder)
    ){
        std::cerr << "Cannot find '" << pathToImportFolder << "' to import "
            << "data, create directory and restart." << std::endl;

        return 21;
    }

    // Подгружаем входные данные

    // nlohmann::json importData;

    std::vector< std::vector<double> > layers;
    std::vector< std::vector<double> > injection;

    // pathToParametersFile
    // if(!importInitialData(pathToImportFolder, importData)){
    //     return 22;
    // }

    if(std::string() != pathToImportFolder){
        pathToLayersFile = pathToImportFolder
            + std::string("/") + pathToLayersFile;
        pathToInjectionFile = pathToImportFolder
            + std::string("/") + pathToInjectionFile;
    }

    // if(0 < importData.size()){
    //     modelingTime = importData["time"].get<double>();
    //     pathToLayersFile = std::string();
    //     pathToInjectionFile = std::string();
    // }

    if(
        !setInitialData(
            pathToLayersFile,
            layers,
            pathToInjectionFile,
            injection
        )
    ){
        return 22;
    }
    std::cout<<"layers:"<<std::endl;
    ai::printMatrix(layers);

    std::cout<<"injection:"<<std::endl;
    ai::printMatrix(injection);
    // Пересчитываем значения для слоёв в величины СИ

    for(std::size_t i = 0; i < layers.size(); ++i){
        layers[i][2] *= std::pow(10, 6);
        layers[i][3] *= std::pow(10, 9);
        layers[i][5] *= std::pow(10, -6);
    }
    #endif

    // Находим эффективный плоский модуль Юнга

    if(!recalculateElasticModulusContrast(
            layers,
            E,
            considerElasticModulusContrast
        )
    ){
        return 31;
    }

    // Проверяем параметры закачки, добавляем данные о свойствах жидкости

    if(!recalculateInjection(injection, modelingTime)){
        return 23;
    }

    std::size_t injectionIndex = 0;

    double timeToChangeInjection = modelingTime * 1.1;
    fluidInjection = injection[injectionIndex][2];
    proppantInjection = injection[injectionIndex][3];
    double n = injection[injectionIndex][4];
    double mu = injection[injectionIndex][5];
    ////////////////////////////////////////////////////////////////////////////

    std::vector<double> proppantDensity;
    std::vector<double> proppantDiameter;
    std::vector<double> maximumConcentration;
    std::vector<size_t> proppantType;
    std::vector<size_t> fluidType;
    std::vector<double> fluidDensity;
    std::vector<double> fluidViscosity;
    std::vector<double> rheologyIndex;
    size_t j = 0;
    #if !defined(PROPPANT_MARKERS)
    proppantDensity.push_back(2500.);
    proppantDiameter.push_back(0.002);
    maximumConcentration.push_back(0.585);
    #endif
    for (std::size_t i = 0; i < injection.size(); ++i) {
        fluidViscosity.push_back(injection[i][5] / injection[0][5]);
        rheologyIndex.push_back(injection[i][4]);
        fluidType.push_back(i);
        fluidDensity.push_back(1000.);
        if (injection[i][3] > epsilon) {
            proppantDensity.push_back(2500.);
            proppantDiameter.push_back(0.002);
            maximumConcentration.push_back(0.585);
            proppantType.push_back(j);
            ++j;
        }
        else {
            proppantType.push_back(0);
        }
    }
    ////////////////////////////////////////////////////////////////////////////

    ++injectionIndex;

    if(injection.size() > injectionIndex){
        timeToChangeInjection = injection[injectionIndex][0];
    }


    double T0 = 1.;
    double wn = std::pow(
        mu * 2. * std::pow(2. * (2. * n + 1.) / n, 1. / n)
        / (E * std::pow(timeScale, n)), 1. / (n + 2.)
    );
    const double gammaR = 1. / 3. * (1. + n / (n + 2.));

    #if !defined(BUILD_DLL)
    // Выводим входные параметры

    if(override){
        std::cout << std::endl;
        std::cout << "Ignoring warnings: override is set." << std::endl;
    }

    if(!runningFromGUI){
        std::cout << std::endl;
        std::cout << "Incoming parameters:" << std::endl
            << "  Q = " << fluidInjection << ", " << "n = " << n
            << ", mu = " << mu << ";" << std::endl
            << "  E\' = " << E << ";" << std::endl
            << "  time = " << modelingTime << ", time step = " << timeStep
            << ", time scale = " << timeScale << ";"
            << std::endl;
    }
    #endif

    // Масштабируем величины

    fluidInjection *= timeScale / (wn * 60.);

    // Строим автомодельное решение

    std::vector< std::vector<double> > modelSolution;

    double initialRadius;

    setModelSolution(n, initialRadius, modelSolution);

    // Определяем момент времени, вплоть до которого будет использоваться
    // автомодельное решение для радиальной плоской трещины в однородной среде
    // (диаметр трещины – 90% от мощности центрального слоя)

    double requiredRadius = -1.;

    for(std::size_t i = 0; i < layers.size(); ++i){
        if(0. > layers[i][0] * layers[i][1]){
            requiredRadius = 0.45 * std::abs(layers[i][0] - layers[i][1]);
        }
    }

    T0 = std::pow(
        requiredRadius / (initialRadius * std::pow(fluidInjection, 1. / 3.)),
        1. / gammaR
    );

    T0 = ai::min(T0, timeToChangeInjection);

    double T = T0;

    initialRadius *= std::pow(fluidInjection, 1. / 3.) * std::pow(T0, gammaR);

    // Пересчитываем автомодельное решение в соответсвии с параметрами закачки

    std::vector<double> zP;
    std::vector<double> openingAtTheStart;

    for(size_t i = 0; i < modelSolution.size(); ++i){
        modelSolution[i][1] *= std::pow(fluidInjection, 1. / 3.);
        openingAtTheStart.push_back(
            modelSolution[i][1] * std::pow(T0, (1. - 2. * n / (n + 2.)) / 3.)
        );
        zP.push_back(modelSolution[i][0]);
    }

    modelSolution.clear();

    const double initialVelocity = initialRadius * gammaR
        * std::pow(T0, gammaR - 1);

    // Автоматически задаём размер ячейки, если требуется

    if(-1. == cellSize){
        cellSize = initialRadius / 5. - epsilon;
    }

    #if !defined(BUILD_DLL)
    // Выводим входные параметры (продолжение)

    if(!runningFromGUI){
        std::cout << "  mesh = " << meshSize << ", cell = " << cellSize << ","
            << " scaling = " << meshScalingCounter << "." << std::endl;
        std::cout << "Initial regime: " << regimeName() << "." << std::endl;
    }
    #endif

    dt = 0.0001 * std::pow(5. / floor(initialRadius / cellSize), 2);

    if(-1. != timeStep && !override){
        timeStep *= 0.001;

        if(timeStep < dt){
            dt = timeStep;
        }else{
            #if !defined(BUILD_DLL)
            std::cerr << "Warning: time step is too big. It was changed to "
                << dt << "!"<< std::endl;
            #endif
        }
    }

    // Проверяем, что на радиус трещины приходится хотя бы пять элементов

    if(5 > floor(initialRadius / cellSize) && !override){
        cellSize = initialRadius / 5. - epsilon;

        #if !defined(BUILD_DLL)
        std::cerr << "Warning: cell size is too big. It was changed to "
            << cellSize << "!"<< std::endl;
    #endif
    }

    #if !defined(BUILD_DLL)
    std::cout << "Number of cells per initial crack: "
        << floor(initialRadius / cellSize) << "." << std::endl;
    #endif

    // Задаём расчётную область

    dx = cellSize;
    dy = dx;

    double axMax = meshSize * dx;

    std::vector<double> x;
    std::vector<double> y;

    for(double i = 0; i <= axMax + epsilon; i += dx){
        x.push_back(i);
    }
    for(double i = -axMax; i <= axMax + epsilon; i += dy){
        y.push_back(i);
    }

    std::size_t xSize = x.size();
    std::size_t ySize = y.size();

    std::vector< std::vector<Cell> > mesh(xSize);

    for(size_t i = 0; i < xSize; ++i){
        mesh[i].resize(ySize);

        for(size_t j = 0; j < ySize; ++j){
            mesh[i][j].setCoordinates(x[i], y[j]);
        }
    }

    // Задаём положение точечного источника закачки

    i00 = 0;
    j00 = floor(0.5 * ySize);

    #if !defined(BUILD_DLL)
    // Сохраняем параметры расчёта

    saveInitialData(
        "./Results/parameters",
        modelingTime,
        ySize,
        xSize,
        dy,
        dx,
        injection,
        layers
    );
    #endif

    // Пересчитываем контрасты по построенной расчётной сетке

    std::vector<double> stress;
    std::vector<double> flatYoungsModulus;
    std::vector<double> leakOff;
    std::vector<double> toughness;

    if(!recalculateStressContrast(layers, stress, y)){
        return 31;
    }

    if(!recalculateElasticModulusContrast(
            layers,
            E,
            flatYoungsModulus,
            y,
            considerElasticModulusContrast
        )
    ){
        return 31;
    }

    if(!recalculateLeakOffContrast(layers, leakOff, y)){
        return 31;
    }
    if(!recalculateToughnessContrast(layers, toughness, y)){
        return 31;
    }

    // Масштабируем величины (продолжение)

    for(std::size_t i = 0; i < stress.size(); ++i){
        stress[i] /= (wn * E);
    }
    for(std::size_t i = 0; i < leakOff.size(); ++i){
        leakOff[i] *= std::sqrt(timeScale) / wn;
    }
    for(std::size_t i = 0; i < toughness.size(); ++i){
        // toughness[i] *= std::pow(10., 6);
        toughness[i] /= std::pow(
            mu * E * std::pow(E / timeScale, n),
            1. / (n + 2.)
        );
    }

    // define PROPPANT_MARKERS
    for(std::size_t i = 0 ; i < injection.size();++i){
            if(injection[i][3] > 0.00001){
                // #define PROPPANT_MARKERS
                break;
            }
    }

    if(!initialized){
        for(std::size_t i = 0; i < injection.size(); ++i){
            injection[i][2] *= timeScale / (wn * 60.);
        }

        initialized = true;
    }

    std::vector<double> zeroVectorX(xSize, 0.);
    std::vector<double> zeroVectorY(ySize, 0.);
    std::vector<size_t> zeroSizeTVectorY(ySize, 0);
    std::vector<double> zeroVectorXY(xSize * ySize, 0.);
    std::vector< std::vector<double> > zeroMatrixXY(xSize, zeroVectorY);

    std::vector<double> opening = zeroVectorXY;
    std::vector< std::vector<size_t> > index(xSize, zeroSizeTVectorY);

    for(size_t i = 0; i < xSize; ++i){
        for(size_t j = 0; j < ySize; ++j){
            index[i][j] = i * ySize + j;
            opening[index[i][j]] = getInitialOpening(
                x[i],
                y[j],
                initialRadius,
                zP,
                openingAtTheStart
            );
        }
    }

    #if !defined(BUILD_DLL)
    std::cout << std::endl;
    std::cout << "Building influence matrix... ";
    #endif

    std::vector< std::vector<double> > influenceMatrix(xSize * ySize,
        zeroVectorXY);

    buildInfluenceMatrix(influenceMatrix, xSize, ySize);

    #if !defined(BUILD_DLL)
    std::cout << "OK." << std::endl;
    #endif

    if(considerElasticModulusContrast){
        #if !defined(BUILD_DLL)
        std::cout << "Applying elastic modulus contrast... ";
        #endif

        #if defined(TRUE_ELASTIC_CONTRAST)
        std::vector< std::vector<double> > additiveToInflunceMatrix;

        // Заполняю ноликами пока
        additiveToInflunceMatrix.resize(xSize * ySize);

        for(std::size_t i = 0; i < xSize * ySize; ++i){
            additiveToInflunceMatrix[i].resize(xSize * ySize);

            for(std::size_t j = 0; j < xSize * ySize; ++j){
                additiveToInflunceMatrix[i][j] = 0.;
            }
        }
        // Заполняю ноликами пока

        // Если не заполнять ноликами
        // void calculateInflunceAddition(
        //     std::vector< std::vector<double> > &additiveToInflunceMatrix,
        //     const std::vector<double> &flatYoungsModulus,
        //     const double xSize,
        //     const double ySize,
        //     const double dx,
        //     const double dy
        // );

        // calculateInflunceAddition(
        //     additiveToInflunceMatrix,
        //     flatYoungsModulus,
        //     xSize,
        //     ySize,
        //     dx,
        //     dy
        // );
        #endif


        /// Решение в виде кусочно-линейной функции
        for(std::size_t i = 0; i < xSize; ++i){
            for(std::size_t k = 0; k < ySize; ++k){
                for(std::size_t j = 0; j < influenceMatrix[0].size(); ++j){
                    influenceMatrix[j][index[i][k]] *= flatYoungsModulus[k];
                }
            }
        }

        #if defined(TRUE_ELASTIC_CONTRAST)
        /// Добавка  для точного решения
        for(std::size_t i = 0; i < influenceMatrix.size(); ++i){
            for(std::size_t j = 0; j < influenceMatrix[0].size(); ++j){
                influenceMatrix[i][j] += additiveToInflunceMatrix[i][j];
            }
        }

        additiveToInflunceMatrix.clear();
        #endif

        #if !defined(BUILD_DLL)
        std::cout << "OK." << std::endl;
        #endif

    }

    std::cout <<"cellSize:" << cellSize << " Mesh: " << mesh.size() << "x" << mesh[0].size() << "." << "Modeling time: " << modelingTime
        << std::endl;

    std::vector< std::vector<bool> > elementIsActive;

    std::vector< std::vector<std::size_t> > activeElements;

    findActiveElements(activeElements, elementIsActive, x, y, initialRadius);

    #if !defined(BUILD_DLL)
    std::cout << "Active elements: " << activeElements.size() << "."
        << std::endl;
    #endif

    // std::vector<std::vector<int> > ac( mesh.size() );
    // for(std::size_t i = 0; i< mesh.size(); ++i){
    //     ac[i].resize(mesh[0].size());
    // }
    //
    // for(size_t k = 0; k < activeElements.size(); ++k){
    //     const size_t i = activeElements[k][0];
    //     const size_t j = activeElements[k][1];
    //
    //     ac[i][j]=1;
    // }
    //
    // std::cout<<"ac:"<<std::endl;
    // ai::printMatrix(ac);
    std::cout<<" "<<std::endl;
    std::vector<double> activationTime;

    for(size_t k = 0; k < activeElements.size(); ++k){
        const size_t i = activeElements[k][0];
        const size_t j = activeElements[k][1];

        if(
            RIBBON != mesh[i][j].type
            && (i00 == i || elementIsActive[i - 1][j])
            && elementIsActive[i + 1][j]
            && elementIsActive[i][j - 1]
            && elementIsActive[i][j + 1]
        ){
            mesh[i][j].type = CHANNEL;
        }

        activationTime.push_back(T0 - dt);
    }

    std::vector< std::vector<double> > distances;
    std::vector<Ribbon> ribbons = findRibbons(
        mesh,
        activeElements,
        distances,
        initialRadius,
        std::sqrt(std::pow(0.5 * dx, 2) + std::pow(0.5 * dy, 2)),
        2. * dx
    );

    #if !defined(BUILD_DLL)
    std::cout << "Ribbons: " << ribbons.size() << "." << std::endl;
    #endif

    #if defined(USE_EIGEN)
    Eigen::VectorXd openingNew;
    Eigen::MatrixXd partialInfluenceMatrix;
    #else
    std::vector<double> openingNew;
    std::vector< std::vector<double> > partialInfluenceMatrix;
    #endif

    #if !defined(BUILD_DLL)
    std::cout << "Building partial influence matrix... ";
    #endif

    buildPartialInfluenceMatrix(
        influenceMatrix,
        activeElements,
        opening,
        openingNew,
        partialInfluenceMatrix,
        index
    );

    #if !defined(BUILD_DLL)
    std::cout << "OK." << std::endl;
    #endif

    double height = 0.;
    double length = 0.;

    std::vector< std::vector<double> > fracture;

    std::size_t stepToCheck = std::round(dx / (20. * dt * initialVelocity));

    bool meshIsNotExhausted = true;

    std::vector< std::vector<double> > velocities = zeroMatrixXY;

    std::vector<double> pressure = zeroVectorXY;
	std::vector<double> concentration = zeroVectorXY;

    omp_set_num_threads(numberOfThreads);

    #if !defined(BUILD_DLL)
    std::cout << std::endl << "Running threads... OK. Size: "
        << getNumberOfThreads() << "." << std::endl;
    #endif

	std::cout << "start calculate" << std::endl;
    int returnCode = calculateEverything(
        opening,
        pressure,
        concentration,
        openingNew,
        distances,
        velocities,
        influenceMatrix,
        partialInfluenceMatrix,
        zeroVectorXY,
        zeroMatrixXY,
        mesh,
        x,
        y,
        ribbons,
        index,
        activeElements,
        elementIsActive,
        activationTime,
        injection,
        layers,
        stress,
        leakOff,
        toughness,
        flatYoungsModulus,
        proppantDensity,
        proppantDiameter,
        maximumConcentration,
        proppantType,
        fluidType,
        fluidDensity,
        fluidViscosity,
        rheologyIndex,
        n,
        mu,
        wn,
        timeScale,
        timeStep,
        fracture,
        length,
        height,
        initialRadius,
        axMax,
        meshScalingCounter,
        meshIsNotExhausted,
        T,
        T0,
        injectionIndex,
        timeToChangeInjection,
        stepToCheck,
        modelingTime,
        considerElasticModulusContrast,
        saveSteps,
        runningFromGUI,
        #if defined(BUILD_DLL)
        override,
        DLL_Parametrs
        #else
        override
        #endif
    );

    if(0 != returnCode){
        return returnCode;
    }

    #if !defined(BUILD_DLL)
    std::cout << "Saving results... ";
    saveData("./Results/opening", opening, index, 1000. * wn);
    saveData("./Results/pressure", pressure, index);
    saveConcentrationData("./Results/concentration", concentration, index);
    saveData("./Results/distances", distances);

    calculateCrackGeometry(mesh, distances, length, height);

    fracture.push_back(
        std::vector<double>{
            T,
            1000 * wn * opening[index[i00][j00]],
            pressure[index[i00][j00]],
            length,
            height,
            length / height
        }
    );
    std::stringstream comment;
    comment << std::right << std::setw(14) << "Time [min]\t"
        << std::right << std::setw(14) << "Opening [mm]\t"
        << std::right << std::setw(14) << "Pressure [Pa]\t"
        << std::right << std::setw(14) << "Length [m]\t"
        << std::right << std::setw(14) << "Height [m]\t"
        << std::right << std::setw(14) << "Aspect ratio [n/d]\t";
    ai::saveMatrix("./Results/fracture", fracture, comment.str());
    comment.str(std::string());

    std::cout << "OK." << std::endl;
    #endif

    // Вычисляем эффективность

    double fluidEfficiency;
    double proppantEfficiency;

    calculateEfficiency(
        fluidEfficiency,
        proppantEfficiency,
        injection,
        opening,
        wn,
        concentration,
        index,
        T,
        T0
    );

    #if !defined(BUILD_DLL)
    std::cout << std::endl;
    std::cout << "Efficiency: " << fluidEfficiency << "% [fluid]";
    if(-1. < proppantEfficiency){
        std::cout << ", " << proppantEfficiency << "% [proppant]";
    }
    std::cout << "." << std::endl;
    #endif

 ////////////////////////////   #if defined(BUILD_DLL)
 ////////////////////////////   if(std::isnan(fluidEfficiency)){
 ////////////////////////////       fluidEfficiency = 0.;
 ////////////////////////////       proppantEfficiency = 0.;
 ////////////////////////////   }
 ////////////////////////////
 ////////////////////////////   /// Сохраняем сжимающие напряжения в центральном слое
 ////////////////////////////   double nominalStress = 0.;
 ////////////////////////////   for (std::size_t i = 0; i < layers.size(); ++i) {
 ////////////////////////////       if (0. > layers[i][0] * layers[i][1]) {
 ////////////////////////////           nominalStress = layers[i][2] * std::pow(10., -6);
 ////////////////////////////       }
 ////////////////////////////   }

	////////////////////////////DLL_Parametrs.J_String = ExportJson(
	////////////////////////////	opening,
	////////////////////////////	wn,
	////////////////////////////	pressure,
	////////////////////////////	concentration,
	////////////////////////////	markers,
	////////////////////////////	markerVolume,
	////////////////////////////	x,
	////////////////////////////	y,
	////////////////////////////	dx,
	////////////////////////////	dy,
	////////////////////////////	i00,
	////////////////////////////	j00,
	////////////////////////////	index,
	////////////////////////////	T,
	////////////////////////////	timeScale,
	////////////////////////////	fluidEfficiency,
	////////////////////////////	fluidDensity,
	////////////////////////////	proppantDensity,
	////////////////////////////	fluidInjection,
	////////////////////////////	proppantInjection,
	////////////////////////////	nominalStress,
	////////////////////////////	DLL_Parametrs
	////////////////////////////);

 ////////////////////////////   // ai::saveLine("rez_emit.json", DLL_Parametrs.J_String);
 ////////////////////////////
 ////////////////////////////   Solver::DataCallback dataCallback =
 ////////////////////////////        Solver::Callback::GetDataCallback();
 ////////////////////////////
 ////////////////////////////   dataCallback(DLL_Parametrs.J_String.c_str());
 ////////////////////////////   #endif

    return 0;
}

/*!
 \details Вычисление полной длины и высоты трещины по элементам на осях,
 проходящих через центральную ячейку

 \param[in] mesh - сеточная матрица
 \param[in] distances - матрица расстояний до фронта трещины
 \param[out] length - длина трещины
 \param[out] height - высота трещины
*/
void calculateCrackGeometry(
    std::vector< std::vector<Cell> > &mesh,
    std::vector< std::vector<double> > &distances,
    double &length,
    double &height
){
    length = 0.;
    height = 0.;
    double heightUp = 0.;
    double heightDown = 0.;

    const std::size_t xSize = distances.size();
    const std::size_t ySize = distances[0].size();

    for(size_t i = 0; i < xSize; ++i){
        if(RIBBON == mesh[i][j00].type){
            length = distances[i][j00]
                + std::sqrt(std::pow(mesh[i][j00].x, 2)
                + std::pow(mesh[i][j00].y, 2));
            length *= 2.;

            break;
        }
    }
    for(size_t j = 0; 2 * j < ySize; ++j){
        if(
            RIBBON == mesh[i00][j].type
            || RIBBON == mesh[i00][ySize - j - 1].type
        ){
            if(0. == heightUp && RIBBON == mesh[i00][j].type ){
                heightUp = distances[i00][j]
                    + std::sqrt(std::pow(mesh[i00][j].x, 2)
                    + std::pow(mesh[i00][j].y, 2));
            }

            if(0. == heightDown && RIBBON == mesh[i00][ySize - j - 1].type){
                heightDown = distances[i00][ySize - j - 1]
                    + std::sqrt(std::pow(mesh[i00][ySize - j - 1].x, 2)
                    + std::pow(mesh[i00][ySize - j - 1].y, 2));
            }

            if(0. != heightUp && 0 != heightDown){
                height = heightUp + heightDown;

                break;
            }
        }
    }
}

/*!
 \details Функция вычисляет отношения объемов жидкости и пропанта в трещине к
 их закаченным объёмам в процентах

 \param[out] fluidEfficiency - эффективность жидкости
 \param[out] proppantEfficiency - эффективность проппанта
 \param[in] injection - матрица, хранящая план закачки
 \param[in] opening - вектор раскрытий
 \param[in] wn - масштабирующий коэффициент
 \param[in] concentration - матрица концентраций проппанта
 \param[in] index - матрица индексов
 \param[in] T - текущее расчетное время
 \param[in] T0 - расчетное время построения автомодельного решения
*/
void calculateEfficiency(
    double &fluidEfficiency,
    double &proppantEfficiency,
    const std::vector< std::vector<double> > &injection,
    const std::vector<double> &opening,
    const double wn,
    const std::vector<double> &concentration,
    const std::vector< std::vector<std::size_t> > &index,
    const double T,
    const double T0
){
    const double xSize = index.size();
    const double ySize = index[0].size();

    double fluidVolumeOut = 0.;
    double fluidVolumeIn = 0.;
    double proppantVolumeOut = 0.;
    double proppantVolumeIn = 0.;

    for(std::size_t i = 0; i < injection.size(); ++i){
        double volumeOut = 0.;

        if(injection[i][1] < T){
            volumeOut = (injection[i][1] - injection[i][0]) * injection[i][2];
        }else{
            if(injection[i][0] < T){
                volumeOut = (T - injection[i][0]) * injection[i][2];
            }
        }

        fluidVolumeOut += volumeOut;
        proppantVolumeOut += volumeOut * injection[i][3];
    }

    /// Поправка на автомодельное решение

    proppantVolumeOut -= T0 * injection[0][2] * injection[0][3];

    fluidVolumeOut *= wn;
    proppantVolumeOut *= wn;

    for(std::size_t i = i00 + 1; i < xSize; ++i){
        for(std::size_t j = 0; j < ySize; ++j){
            fluidVolumeIn += opening[index[i][j]];
            proppantVolumeIn += concentration[index[i][j]]
                * opening[index[i][j]];
        }
    }

    fluidVolumeIn *= 2.;
    proppantVolumeIn *= 2.;

    for(std::size_t j = 0; j < ySize; ++j){
        fluidVolumeIn += opening[index[i00][j]];
        #if defined(PROPPANT_MARKERS)
        proppantVolumeIn += 2. * concentration[index[i00][j]]
            * opening[index[i00][j]];
        #else
        proppantVolumeIn += concentration[index[i00][j]]
            * opening[index[i00][j]];
        #endif
    }

    fluidVolumeIn *= wn * dx * dy;
    proppantVolumeIn *= wn * dx * dy * 2.65 * 1000;

    fluidEfficiency = 100 * fluidVolumeIn / fluidVolumeOut;
    if(0. < proppantVolumeOut){
        proppantEfficiency = 100 * proppantVolumeIn / proppantVolumeOut;
    }else{
        proppantEfficiency = -1.;
    }
}

/*!
 \details Функция сохраняет результаты расчетов в формате json

 \param[in] filename - имя сохраняемого файла(без расширения)
 \param[in] opening - вектор раскрытий
 \param[in] wn - масштабирующий коэффициент
 \param[in] pressure - вектор давлений
 \param[in] concentration - матрица концентраций пропанта
 \param[in] index - матрица индексов
 \param[in] Time - текущее расчетное время
 \param[in] timeScale - масштаб времени
 \param[in] fluidEfficiency - эффективность жидкости
 \param[in] fluidDensity - плотность жидкости
 \param[in] proppantEfficiency - эффективность пропанта
 \param[in] proppantDensity - плотность пропанта
 \param[in] nominalStress - сжимающее напряжение в центральном слое
*/
