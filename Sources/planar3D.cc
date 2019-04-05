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

#include "io.hh"
#include "mesh.hh"
#include "rhs3D.hh"
#include "ribbons.hh"
#include "automodel.hh"
#include "initialData.hh"
#include "influenceMatrix.hh"
#include "findActiveElements.hh"

 //#define BUILD_DLL

#if defined(BUILD_DLL)
    #include "dll/api.h"
    #include "dll/api_callback.h"
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
void createMatrixDiag(
                        std::vector<double>& A1,
                        std::vector<double>& A2,
                        std::vector<double>& A3,
                        std::vector<double>& A4,
                        std::vector<double>& A5,
                        size_t N_dof,
                        size_t Nx,
                        size_t Ny,
                        size_t Nz){

    size_t NxNy = Nx * Ny;
//     % create all diagonals
    for(size_t i = 0; i < N_dof ;++i){
        A2[i] = 1.;
        A3[i] = -6.;
        A4[i] = 1.;
    }
    for(size_t i = 0; i<N_dof - Nx; ++i){
        A1[i] = 1.;
        A5[i] = 1.;
    }
    // ai::printMarker();
    //% set some elements of A equal to zero using boundary conditions
   //     % zero flux at x=0 and x=Nx*dx

   size_t n;
   size_t p ;
   for(size_t j = 0 ; j < Ny; ++j){
       for(size_t k = 0 ; k < Nz; ++k){
           n = j * Nx + k * NxNy;
           A2[n] = 0.;
           A4[n] = 2.;

           p = n + Nx - 1;
           A4[p] = 0.;
           A3[p] = -5.;
       }
   }
// ai::printMarker();
   for(size_t i = 0; i < Nx; ++i ){
      for(size_t k = 0; k < Nz; ++k){
          n = i  + k * NxNy;
          if(n-Nx < N_dof-Nx){
            A1[n-Nx] = 0.;
          }
          A3[n] = -5.;

          p = n + (Ny-1)*Nx;
          if(p <= N_dof-Nx-1){
            A5[p] = 0.;
          }
          A3[p] = -5.;
      }
  }
// ai::printMarker();
  for(size_t i = 0; i<Nx ;++i){
      for(size_t j = 0; j < Ny; ++j){
          n = i  + j * Nx;

          if(n + (Nz-1)*NxNy < N_dof-Nx){
              A5[n + (Nz-1)*NxNy] = 0.;
          }
          A3[n + (Nz-1)*NxNy] = -5.;
      }
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
    std::vector<double> &openingNew,
    std::vector< std::vector<double> > &distances,
    std::vector< std::vector<double> > &velocities,
    std::vector <double>& A1, //пять диагоналей матрицы
    std::vector <double>& A2,
    std::vector <double>& A3,
    std::vector <double>& A4,
    std::vector <double>& A5,
    // std::vector< std::vector<double> > &influenceMatrix,
    // std::vector< std::vector<double> > &partialInfluenceMatrix,
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
    double &fluidDensity,
    double &proppantDensity,
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
    size_t &Nx,
    size_t &Ny,
    size_t &Nz,
    size_t &N_dof,
    #if defined(BUILD_DLL)
    const bool override,
    DLL_Param &DLL_Parametrs
    #else
    const bool override
    #endif
){
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
    std::vector<double> temperature(activeElements.size(), 0.);

    for(size_t k = 0; k < ribbons.size(); ++k){
        const size_t i = ribbons[k].i;
        const size_t j = ribbons[k].j;

        savedDistances[k] = distances[i][j];
    }

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

    std::vector<double> b(N_dof, 0.);  // right side of finite difference discretization of Laplace equation
    std::vector<double> r1(N_dof, 0.);
    std::vector<double> r2(N_dof, 0.);
    std::vector<double> p(N_dof, 0.);
    std::vector<double> A_p(N_dof, 0.);
    std::vector<double> T1(N_dof, 0.);  // unknowns in  finite difference discretization of Laplace equation (AT = b)


    auto startTime = ai::time();

	bool pauseCallback = false;	//Флаг того, что на паузе мы послали кэлбак
    size_t iter = 0;
    std::cout<<"meshflag = "<<meshIsNotExhausted<<std::endl;

    std::cout<<"modeling time  = "<<modelingTime<<"  T = "<<T<<std::endl;
    while(modelingTime >= T && meshIsNotExhausted){
        iter ++;

        if(1 == iter){
            meshIsNotExhausted = 0;
        }

std::cout<<"iter = "<<iter<<std::endl;
        #if defined(BUILD_DLL)
		if (DLL_Parametrs.Dll_State == DLL_Parametrs.Running)
		{
			pauseCallback = false;
		}
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

				ai::saveLine("rez_emit.json", DLL_Parametrs.J_String.c_str());

				Solver::DataCallback dataCallback =	Solver::Callback::GetDataCallback();

				dataCallback(ExportJson(opening, wn, pressure, concentration, x, y, dx, dy, i00, j00, index, T, timeScale, fluidEfficiency, fluidDensity, proppantDensity, fluidInjection, proppantInjection, nominalStress, DLL_Parametrs).c_str());

				dataCallback(DLL_Parametrs.J_String.c_str());


			}
			pauseCallback = true;
			continue;
		}
		if (DLL_Parametrs.Dll_State == DLL_Parametrs.Idle)
		{
			//ai::saveLine("break", DLL_Parametrs.J_String.c_str());
			//break;
		}
        #endif

        if(-1. == timeStep){
            dt = 0.05 * mu / (E * std::pow(wn * ai::max(opening) / dx, 3));
        }

        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        time0 = ai::time();
        #endif

        calculatePressure(
            pressure,
            index,
            activeElements,
            A1,
            A2,
            A3,
            A4,
            A5,
            openingNew,
            stress,
            T1,
            b,
            r1,
            r2,
            p,
            A_p,
            N_dof,
            Nx,
            Ny
        );
        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        timeMeasurements[0] += ai::duration(time0, "us");
        #endif

        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        auto time1 = ai::time();
        #endif
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
            T,
            fluidDensity,
            proppantDensity,
            n
        );
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
        for(std::size_t i = 0; i < opening.size(); ++i){
            concentration[i] = ai::max(
                concentration[i] * opening[i] + dCdt[i] * dt,
                0.
            );
            opening[i] = ai::max(opening[i] + dWdt[i] * dt, 0.);

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

            /// \todo проверить баланс масс
            if(0.585 - epsilon < concentration[i]){
                concentration[i] = 0.585 - epsilon;
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
            if(runningFromGUI){
                std::cout << "Progress: " << (T - T0) / (modelingTime - T0)
                    << std::endl;
            }else{
                ai::showProgressBar((T - T0) / (modelingTime - T0));
            }

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
                    meshIsNotExhausted = false;
                    break;

                    /// 2) Масштабирование
                    // int returnCode = scaleMesh(
                    //     x,
                    //     y,
                    //     mesh,
                    //     index,
                    //     activeElements,
                    //     elementIsActive,
                    //     activationTime,
                    //     ribbons,
                    //     distances,
                    //     opening,
                    //     concentration,
                    //     layers,
                    //     flatYoungsModulus,
                    //     stress,
                    //     leakOff,
                    //     toughness,
                    //     zeroVectorXY,
                    //     zeroMatrixXY,
                    //     xSize,
                    //     ySize,
                    //     initialRadius,
                    //     axMax,
                    //     dMin1,
                    //     dCenter1,
                    //     dMax1,
                    //     dMin2,
                    //     dCenter2,
                    //     dx,
                    //     dy,
                    //     dt,
                    //     mu,
                    //     n,
                    //     wn,
                    //     timeScale,
                    //     T,
                    //     T0,
                    //     E,
                    //     considerElasticModulusContrast,
                    //     meshScalingCounter
                    // );

                    /// 3) Достроение сетки
                    /// Добавляем по 10 элементов вдоль каждой оси

                    // int returnCode = completeMesh(
                    //     xSize + 10,
                    //     x,
                    //     y,
                    //     i00,
                    //     j00,
                    //     mesh,
                    //     index,
                    //     opening,
                    //     pressure,
                    //     concentration,
                    //     distances,
                    //     elementIsActive,
                    //     activationTime,
                    //     layers,
                    //     flatYoungsModulus,
                    //     stress,
                    //     leakOff,
                    //     toughness,
                    //     xSize,
                    //     ySize,
                    //     considerElasticModulusContrast,
                    //     influenceMatrix
                    // );
                    //
                    /// \todo savedDistances должны быть в scaleMesh
                    /// и completeMesh!
                    /// \todo изменить запись в ribbons
                    //
                    // if(0 != returnCode){
                    //     return returnCode;
                    // }
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
                buildPartialInfluenceMatrix( activeElements,
                    opening, openingNew, index
                );
            }

            // Сохраняем параметры трещины

            calculateCrackGeometry(mesh, distances, length, height);

            if((Nx - 4)* dx < length/2|| (Ny-4) * dy < height){

                Nx = std::ceil(length/2.)+4;
                Ny = std::ceil(height)+4;
                Nz = Nx;
                N_dof = Nx*Ny*Nz;
                // std::cout<<"Nx = "<<Nx<<"  Ny  = "<<Ny<<" Nz = "<<Nz<<" N_dof = "<<N_dof<<std::endl;

                size_t oldsize  = T1.size();

                int difsize = N_dof - oldsize;

                if(difsize > 0){
                    std::vector<double> A1(N_dof - Nx, 0.);
                    std::vector<double> A2(N_dof, 0.);
                    std::vector<double> A3(N_dof, 0.);
                    std::vector<double> A4(N_dof, 0.);
                    std::vector<double> A5(N_dof - Nx, 0.);


                    createMatrixDiag(A1,A2,A3,A4,A5, N_dof, Nx , Ny, Nz);

                    for(size_t i = 0 ; i< difsize; ++i){
                        T1.push_back(0.);
                        b.push_back(0.);
                        r1.push_back(0.);
                        r2.push_back(0.);
                        p.push_back(0.);
                        A_p.push_back(0.);
                    }
                }
                else{
                    if(difsize < 0){
                        std::vector<double> A1(N_dof - Nx, 0.);
                        std::vector<double> A2(N_dof, 0.);
                        std::vector<double> A3(N_dof, 0.);
                        std::vector<double> A4(N_dof, 0.);
                        std::vector<double> A5(N_dof - Nx, 0.);

                        createMatrixDiag(A1,A2,A3,A4,A5, N_dof, Nx , Ny, Nz);

                        T1.resize(N_dof);
                        b.resize(N_dof);
                        r1.resize(N_dof);
                        r2.resize(N_dof);
                        p.resize(N_dof);
                        A_p.resize(N_dof);
                    }
                }
            }


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
        if(
            (T * 60.) / DLL_Parametrs.Emit_time
            - (int) T * 60. / (int) DLL_Parametrs.Emit_time < dt
        ){
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

            if(std::isnan(fluidEfficiency)){
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

 //           ai::saveLine("rez_emit.json", DLL_Parametrs.J_String.c_str());

            Solver::DataCallback dataCallback = Solver::Callback::GetDataCallback();

            dataCallback(
                ExportJson(opening, wn, pressure, concentration, x, y, dx, dy, i00, j00, index, T, timeScale, fluidEfficiency, fluidDensity, proppantDensity, fluidInjection, proppantInjection, nominalStress, DLL_Parametrs).c_str());

            dataCallback(DLL_Parametrs.J_String.c_str());
        }
        #else
        if(saveSteps && 0 == globalStep % stepToSave){
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
        #endif
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

            #if !defined(BUILD_DLL)
            std::string line("Time = ");
            line += ai::string(injection[injectionIndex - 1][0])
                + ai::string(". Injection was changed.");
            ai::printLine(line);
            ai::printLine(
                ai::string("Current fluid: n = ") + ai::string(n)
                + ai::string(", mu = ") + ai::string(mu)
            );
            #endif
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
            && height > 2 * initialRadius
        ){
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
                meshScalingCounter
            );

            if(0 != returnCode){
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
    double modelingTime = DLL_Parametrs.modelingTime;
    double timeStep = -1.;
    double timeScale = 60.;
    double cellSize = -1.;
    double meshSize = 61.;
    int meshScalingCounter = 2;
    const bool considerElasticModulusContrast = false;
    const bool saveSteps = false;
    const bool runningFromGUI = false;
    const std::size_t numberOfThreads = 1;
    const bool override = false;
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
    double fluidDensity = 1000.;
    double proppantDensity = 2500.;

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
std::cout<<"TT = "<<T<<std::endl;
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

    // Масштабируем выличины (продолжение)

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

    std::vector<std::vector<double> >opep(xSize, zeroVectorY);

    std::vector<double> opening = zeroVectorXY;
    std::vector< std::vector<size_t> > index(xSize, zeroSizeTVectorY);
    // initialRadius = 10.;
    // double co;
    for(size_t i = 0; i < xSize; ++i){
        for(size_t j = 0; j < ySize; ++j){
            index[i][j] = i * ySize + j;

            // co = 0.5*cos(0.5 *M_PI* std::sqrt(x[i]*x[i]+y[j]*y[j])/(
            // xSize*0.8/dx
            // initialRadius
            // ) );
            opening[ index[i][j] ] = getInitialOpening(
                x[i],
                y[j],
                //xSize*0.8/dx,
                initialRadius,
                //10.,
                zP,
                openingAtTheStart);
            // )>epsilon?co:0.;


			opep[i][j] = abs(opening[index[i][j]]);//>epsilon?co:0.;
        }
    }
    ai::saveVector("iniopen",opening);
    ai::saveMatrix("opep",opep);

    #if !defined(BUILD_DLL)
    std::cout << std::endl;
    // std::cout << "Building influence matrix... ";
    #endif

    #if !defined(BUILD_DLL)
    // std::cout << "OK." << std::endl;
    #endif

    if(considerElasticModulusContrast){
        #if !defined(BUILD_DLL)
        std::cout << "Applying elastic modulus contrast... ";
        #endif

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

        // for(std::size_t i = 0; i < xSize; ++i){
        //     for(std::size_t k = 0; k < ySize; ++k){
        //         for(std::size_t j = 0; j < influenceMatrix[0].size(); ++j){
        //             influenceMatrix[index[i][k]][j] *= flatYoungsModulus[k];
        //         }
        //     }
        // }

        // for(std::size_t i = 0; i < influenceMatrix.size(); ++i){
        //     for(std::size_t j = 0; j < influenceMatrix[0].size(); ++j){
        //         influenceMatrix[i][j] += additiveToInflunceMatrix[i][j];
        //     }
        // }
        //
        // additiveToInflunceMatrix.clear();

    //     #if !defined(BUILD_DLL)
    //     std::cout << "OK." << std::endl;
    //     #endif
    //
     }

    std::cout << "Mesh: " << mesh.size() << "x" << mesh[0].size() << "."
        << std::endl;

    std::vector< std::vector<bool> > elementIsActive;

    std::vector< std::vector<std::size_t> > activeElements;

    findActiveElements(activeElements, elementIsActive, x, y, initialRadius);//////////////////////


    size_t Nx = ((int)std::ceil(initialRadius)+5)%2==0?(std::ceil(initialRadius)+6):(std::ceil(initialRadius)+5);
    size_t Ny = ((int)std::ceil(2*initialRadius)+12)%2==0?(std::ceil(2*initialRadius)+13):(std::ceil(2*initialRadius)+12);;
    size_t Nz = Nx;

    std::cout<<"Nx  = "<<Nx<<"  Ny= "<<Ny<<"  Nz = "<<Nz<<std::endl;

    size_t N_dof = Nx*Ny*Nz;

    std::vector<double> A1(N_dof - Nx, 0.);
    std::vector<double> A2(N_dof, 0.);
    std::vector<double> A3(N_dof, 0.);
    std::vector<double> A4(N_dof, 0.);
    std::vector<double> A5(N_dof - Nx, 0.);

    createMatrixDiag(   A1,
                        A2,
                        A3,
                        A4,
                        A5,
                        N_dof,
                        Nx,
                        Ny,
                        Nz);

    std::cout<<"Nx = "<<Nx<<" Ny = "<<Ny<<" Nz = "<<Nz<<std::endl;

    #if !defined(BUILD_DLL)
    std::cout << "Active elements: " << activeElements.size() << "."
        << std::endl;
    #endif

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

    std::vector<double> openingNew(activeElements.size() ,0.);
    // std::vector< std::vector<double> > partialInfluenceMatrix;

    #if !defined(BUILD_DLL)
    std::cout << "Building partial influence matrix... ";
    #endif

    buildPartialInfluenceMatrix(

        activeElements,
        opening,
        openingNew,

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

    int returnCode = calculateEverything(
        opening,
        pressure,
        concentration,
        openingNew,
        distances,
        velocities,
        A1,
        A2,
        A3,
        A4,
        A5,
        //influenceMatrix,
        //partialInfluenceMatrix,
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
        fluidDensity,
        proppantDensity,
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
        Nx,
        Ny,
        Nz,
        N_dof,
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
    saveData("./Results/concentration", concentration, index);
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

    #if defined(BUILD_DLL)
    if(std::isnan(fluidEfficiency)){
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

    // ai::saveLine("rez_emit.json", DLL_Parametrs.J_String);

    Solver::DataCallback dataCallback =
         Solver::Callback::GetDataCallback();

    dataCallback(DLL_Parametrs.J_String.c_str());
    #endif

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
        proppantVolumeIn += concentration[index[i00][j]]
            * opening[index[i00][j]];
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
