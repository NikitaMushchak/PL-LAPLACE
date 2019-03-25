/*!
 \file planar3D.cc
 \brief
 \details
файл с функцией planar3D - главной функцией расчета и вызова функций различных режимов
Данный файл содержит в себе определения основных
классов, используемых в демонстрационной программе
*/

#include <vector>
#include <iostream>
#include <algorithm>
#include    <math.h>

#include "nlohmann/json.hpp"
#include "ailibrary/ai.hh"


#include "io.hh"
#include "rhs3D.hh"
#include "ribbons.hh"
#include "automodel.hh"
#include "initialData.hh"
#include "createMatrix.hh"
#include "influenceMatrix.hh"
#include "findActiveElements.hh"

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
double Kic;
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

void createMatrixDiag(std::vector<std::vector<double> >& A,
                            size_t& N_dof, size_t& Nx, size_t& Ny, size_t& Nz){
    size_t NxNy = Nx * Ny;
    for(size_t i = 0; i < N_dof ;++i){
        A[i][0] = 1.;
        A[i][1] = 1.;
        A[i][2] = 1.;
        A[i][3] = -6.;
        A[i][4] = 1.;
        A[i][5] = 1.;
        A[i][6] = 1.;
    }
   double n;
   double p ;
   for(size_t j = 0 ; j < Ny; ++j){
       for(size_t k = 0 ; k < Nz; ++k){
           n = j * Nx + k * NxNy;
           A[n][2] = 0.;
           A[n][4] = 2.;

           p = n + Nx - 1;
           A[p][4] = 0.;
           A[p][3] = -5.;
       }
   }
   for(size_t i = 0; i < Nx; ++i ){
      for(size_t k = 0; k < Nz; ++k){
          n = i  + k * NxNy;
          A[n][1] = 0.;
          A[n][3] = -5.;

          p = n + (Ny-1)*Nx;
          A[p][5] = 0.;
          A[p][3] = -5.;
      }
  }

  for(size_t i = 0; i<Nx ;++i){
      for(size_t j = 0; j < Ny; ++j){
          n = i  + j * Nx;
          A[i + j*Nx][0] = 0.;

          A[n + (Nz-1)*NxNy][5] = 0.;
          A[n + (Nz-1)*NxNy][3] = -5.;
      }
  }
}

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
    const double wn,
    const double timeScale,
    const double T,
    const double T0,
    double &E,
    const bool considerElasticModulusContrast,
    int &meshScalingCounter
){
    --meshScalingCounter;

    ai::printLine(
        ai::string("Doubling the mesh at time = ") + ai::string(T)
    );

    axMax *= 2.;
    dx *= 2.;
    dy = dx;
    dt *= 1.2;

    x.clear();
    y.clear();
    for(double i = 0; i <= axMax + epsilon; i += dx){
        x.push_back(i);
    }
    for(double i = -axMax; i <= axMax + epsilon; i += dy){
        y.push_back(i);
    }

    if(xSize != x.size() || ySize != y.size()){
        std::cerr << "Cannot scale the mesh: wrong sizes"
            << std::endl;

        return 31;
    }

    for(size_t i = 0; i < xSize; ++i){
        for(size_t j = 0; j < ySize; ++j){
            mesh[i][j].setCoordinates(x[i], y[j]);
        }
    }
    std::vector<double> openingNew = zeroVectorXY;
    std::vector<double> concentrationNew = zeroVectorXY;
    std::vector< std::vector<double> > activationTimeTable = zeroMatrixXY;
    std::vector< std::vector<double> > distancesNew = zeroMatrixXY;

    for(size_t i = 0; i < xSize; ++i){
        for(size_t j = 0; j < ySize; ++j){
            activationTimeTable[i][j] = T0 - 10. * dt;
        }
    }
    for(std::size_t k = 0; k < activeElements.size(); ++k){
        const std::size_t i = activeElements[k][0];
        const std::size_t j = activeElements[k][1];

        activationTimeTable[i][j] = activationTime[k];
    }

    ribbons.clear();
    activeElements.clear();
    activationTime.clear();
    for(int j = 0; 2 * j < ySize - j00 - 1; ++j){
        for(int i = 0; 2 * i < xSize - 1; ++i){
            openingNew[index[i00 + i][j00 + j]] = opening[index[i00 + 2 * i][j00 + 2 * j]];
            openingNew[index[i00 + i][j00 - j]] = opening[index[i00 + 2 * i][j00 - 2 * j]];
            concentrationNew[index[i00 + i][j00 + j]] = concentration[index[i00 + 2 * i][j00 + 2 * j]];
            concentrationNew[index[i00 + i][j00 - j]] = concentration[index[i00 + 2 * i][j00 - 2 * j]];

            std::size_t k = i00 + 2 * i;
            std::size_t l = j00 + 2 * j;
            std::size_t m = j00 - 2 * j;

            double distance = 10 * dMin2;

            if(RIBBON == mesh[k][l].type){
                if(distances[k][l] >= dMax1){
                    distance = ai::min(distances[k][l], distance);
                }
            }

            l += 1;

            if(RIBBON == mesh[k][l].type){
                if(distances[k][l] >= dMax1){
                    distance = ai::min(distances[k][l] + dx, distance);
                }
            }

            l -= 2;

            if(RIBBON == mesh[k][l].type){
                if(distances[k][l] >= dMax1){
                    distance = ai::min(distances[k][l] + dx, distance);
                }
            }

            l += 1;
            k += 1;

            if(RIBBON == mesh[k][l].type){
                if(distances[k][l] >= dMax1){
                    distance = ai::min(distances[k][l] + dx, distance);
                }
            }

            l += 1;

            if(RIBBON == mesh[k][l].type){
                if(distances[k][l] >= dMax1){
                    distance = ai::min(distances[k][l] + sqrt(2.) * dx, distance);
                }
            }

            l -= 2;

            if(RIBBON == mesh[k][l].type){
                if(distances[k][l] >= dMax1){
                    distance = ai::min(distances[k][l] + sqrt(2.) * dx, distance);
                }
            }

            k -= 1;
            l += 1;

            if(i00 < k){
                k -= 1;

                if(RIBBON == mesh[k][l].type){
                    if(distances[k][l] >= dMax1){
                        distance = ai::min(distances[k][l] + dx, distance);
                    }
                }

                l += 1;

                if(RIBBON == mesh[k][l].type){
                    if(distances[k][l] >= dMax1){
                        distance = ai::min(distances[k][l] + sqrt(2.) * dx, distance);
                    }
                }

                l -= 2;

                if(RIBBON == mesh[k][l].type){
                    if(distances[k][l] >= dMax1){
                        distance = ai::min(distances[k][l] + sqrt(2.) * dx, distance);
                    }
                }

                l += 1;
                k += 1;
            }

            if(distance <= dMin2 && distance >= dMax1){
                distancesNew[i00 + i][j00 + j] = distance;
            }

            distance = 10 * dMin2;

            if(RIBBON == mesh[k][m].type){
                if(distances[k][m] >= dMax1){
                    distance = ai::min(distances[k][m], distance);
                }
            }

            m += 1;

            if(RIBBON == mesh[k][m].type){
                if(distances[k][m] >= dMax1){
                    distance = ai::min(distances[k][m] + dx, distance);
                }
            }

            m -= 2;

            if(RIBBON == mesh[k][m].type){
                if(distances[k][m] >= dMax1){
                    distance = ai::min(distances[k][m] + dx, distance);
                }
            }

            m += 1;
            k += 1;

            if(RIBBON == mesh[k][m].type){
                if(distances[k][m] >= dMax1){
                    distance = ai::min(distances[k][m] + dx, distance);
                }
            }

            m += 1;

            if(RIBBON == mesh[k][m].type){
                if(distances[k][m] >= dMax1){
                    distance = ai::min(distances[k][m] + sqrt(2.) * dx, distance);
                }
            }

            m -= 2;

            if(RIBBON == mesh[k][m].type){
                if(distances[k][m] >= dMax1){
                    distance = ai::min(distances[k][m] + sqrt(2.) * dx, distance);
                }
            }

            k -= 1;
            m += 1;

            if(i00 < k){
                k -= 1;

                if(RIBBON == mesh[k][m].type){
                    if(distances[k][m] >= dMax1){
                        distance = ai::min(distances[k][m] + dx, distance);
                    }
                }

                m += 1;

                if(RIBBON == mesh[k][m].type){
                    if(distances[k][m] >= dMax1){
                        distance = ai::min(distances[k][m] + sqrt(2.) * dx, distance);
                    }
                }

                m -= 2;

                if(RIBBON == mesh[k][m].type){
                    if(distances[k][m] >= dMax1){
                        distance = ai::min(distances[k][m] + sqrt(2.) * dx, distance);
                    }
                }

                m += 1;
                k += 1;
            }

            if(distance <= dMin2 && distance >= dMax1){
                distancesNew[i00 + i][j00 - j] = distance;
            }
        }
    }

    for(size_t i = 0; i < xSize; ++i){
        for(size_t j = 0; j < ySize; ++j){
            elementIsActive[i][j] = false;
            if(distancesNew[i][j] > epsilon){
                ribbons.push_back(Ribbon(i, j));
                mesh[i][j].type = RIBBON;
            }else{
                mesh[i][j].type = OUTSIDE;
            }
            if(epsilon < openingNew[index[i][j]] || RIBBON == mesh[i][j].type){
                activeElements.push_back(std::vector<size_t>{i, j});
                elementIsActive[i][j] = true;

                if(T0 - 9. * dt < activationTimeTable[i][j]){
                    activationTime.push_back(activationTimeTable[i][j]);
                }else{
                    activationTime.push_back(T);
                }
            }
        }
    }

    activationTimeTable.clear();

    std::cout << "Active elements: " << activeElements.size() << "."
        << std::endl;


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
    }
    initialRadius *= 2;
    std::cout << "Ribbons: " << ribbons.size() << "." << std::endl;

    opening = openingNew;
    distances = distancesNew;
    concentration = concentrationNew;

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

    // Масштабируем выличины (продолжение)

    for(std::size_t i = 0; i < stress.size(); ++i){
        stress[i] /= (wn * E);
    }
    for(std::size_t i = 0; i < leakOff.size(); ++i){
        leakOff[i] *= std::sqrt(timeScale) / wn;
    }

    dMin1 = std::sqrt(std::pow(0.5 * dx, 2)
        + std::pow(1.5 * dy, 2));
    dCenter1 = dx;
    dMax1 = std::sqrt(std::pow(0.5 * dx, 2)
        + std::pow(0.5 * dy, 2));
    dMin2 = std::sqrt(std::pow(1.5 * dx, 2)
        + std::pow(1.5 * dy, 2));
    dCenter2 = std::sqrt(2) * dx;

    return 0;
}

/// \todo aaaa
inline void saveData(
    const std::string filename,
    std::vector<double> &data,
    std::vector< std::vector<std::size_t> > &index,
    const double multiplier = 1.,
    const std::string comment = std::string()
){
    std::ofstream output(filename + std::string("_m.txt"));

    if(!output.good()){
        throw std::runtime_error(
            ai::string("Exception while saving the matrix into the file: ")
            + filename
        );
    }

    if(std::string() != comment){
        output << comment << std::endl;
    }

    for(int j = (int) index[0].size() - 1; j >= 0; --j){
        for(std::size_t i = 0; i < index.size(); ++i){
            output << std::setw(14) << multiplier * data[index[i][j]];
        }

        output << std::endl;
    }

    output.close();
}

/// \todo aaaa
template<typename T>
inline void saveData(
    const std::string filename,
    std::vector< std::vector<T> > &data,
    const double multiplier = 1.,
    const std::string comment = std::string()
){
    std::ofstream output(filename + std::string("_m.txt"));

    if(!output.good()){
        throw std::runtime_error(
            ai::string("Exception while saving the matrix into the file: ")
            + filename
        );
    }

    if(std::string() != comment){
        output << comment << std::endl;
    }

    for(int j = (int) data[0].size() - 1; j >= 0; --j){
        for(std::size_t i = 0; i < data.size(); ++i){
            output << std::setw(14) << data[i][j];
        }

        output << std::endl;
    }

    output.close();
}

int calculateEverything(
    std::vector<double> &opening,
    std::vector<double> &pressure,
    std::vector<double> &concentration,
    std::vector<double> &openingNew,
    std::vector< std::vector<double> > &distances,
    std::vector< std::vector<double> > &velocities,
    std::vector< std::vector<double> > &influenceMatrix,
    std::vector< std::vector<double> > &partialInfluenceMatrix,
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
    const bool override,
    size_t &N_dof,
    size_t &Nx,
    size_t &Ny,
    size_t &Nz
){


    double dx = x[1]-x[0];
    double dy = y[1]-y[0];
    std::cout<<"dx = "<<dx <<" dy = "<<dy<<std::endl;

    alpha = 2. / (n + 2.);
    double BAlpha = 0.25 * alpha * tan(0.5 * M_PI - M_PI * (1 - alpha));
    Amu = std::pow(BAlpha * (1 - alpha), -0.5 * alpha);

    std::size_t globalStep = (std::size_t) T / dt;
    std::size_t step = 0;
    double savedTime = T;
    std::size_t savingStep = 1;
    std::size_t stepToSave = (std::size_t) 1. / dt;
    const std::size_t xSize = index.size();
    const std::size_t ySize = index[0].size();

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
    std::vector<double> savedDistances(ribbons.size(), 0.);

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

    auto startTime = ai::time();

	// std::vector<double> val1;
	// std::vector<double> val2;

    while(modelingTime >= T && meshIsNotExhausted){
		// ai::printLine("1");
		//
        // meshIsNotExhausted=0;
		double maxOpen = wn*ai::max(opening);

        if(-1. == timeStep){
            dt = 0.06 * mu / (E * std::pow(maxOpen / dx, 3));
        }
        std::cout<<"active Elements = "<<activeElements.size()<<std::endl;
        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        time0 = ai::time();
        #endif
        calculatePressure(
            pressure,
            index,
            activeElements,
            partialInfluenceMatrix,
            openingNew,
            stress,
            N_dof,
            Nx,
            Ny
        );

		//ai::saveVector("openingnew",openingNew);
		//ai::saveVector("pressure",pressure);
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
        //ai::saveVector("dWdt", dWdt);
        //ai::saveVector("dCdt", dCdt);
        //ai::saveVector("calcope",opening);
        //std::cout<<"calcOpenAndConce"<<std::endl;
        ai::printMarker();//2
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
                    * opening[index[i][j]] / Kic;

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
        ai::printMarker();//3
        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        timeMeasurements[2] += ai::duration(time2, "us");
        #endif

        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        auto time3 = ai::time();
        #endif
        //dt= 0.0001;
        std::cout<<"dt = "<<std::fixed<<std::setprecision(7)<<dt<<std::endl;
        for(std::size_t i = 0; i < opening.size(); ++i){
            //std::cout<<"opening  cicle"<<std::endl;
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
        ai::saveVector("ateropen", opening);
        std::cout<<"epsilon()"<<std::endl;
        ai::printMarker();//4
        #if defined(DEBUG) && !defined(BUILD_DLL)
        if(std::isnan(opening[index[i00][j00]])){
            ai::printLine(
                ai::string("Planar3D is dead. Time = ") + ai::string(T)
            );
            break;
        }
        #endif
        std::cout<<"Opening new"<<std::endl;
        ai::saveVector("opeconce",opening);
        for(std::size_t k = 0; k < activeElements.size(); ++k){
            const std::size_t i = activeElements[k][0];
            const std::size_t j = activeElements[k][1];

            openingNew[k] = opening[index[i][j]];
        }
        ai::printMarker();//5
        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        timeMeasurements[3] += ai::duration(time3, "us");
        #endif

        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        auto time4 = ai::time();
        #endif

        ////// STEPTOCHECK \\\\\\
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
                    meshIsNotExhausted = false;
                    break;
                }

                const double d = distances[i][j];

                if(0 < i){
                    findNewRibbons(i - 1, j, d, dMin1, dCenter1, n, opening, leakOff, mesh,
                        index, activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                        velocities, T, opening[index[i - 1][j]]);

                    findNewRibbons(i - 1, j - 1, d, dMin2, dCenter2, n, opening, leakOff, mesh,
                        index, activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                        velocities, T, opening[index[i - 1][j - 1]]);

                    findNewRibbons(i - 1, j + 1, d, dMin2, dCenter2, n, opening, leakOff, mesh,
                        index, activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                        velocities, T, opening[index[i - 1][j + 1]]);
                }

                findNewRibbons(i + 1, j, d, dMin1, dCenter1, n, opening, leakOff, mesh, index,
                    activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                    velocities, T, opening[index[i + 1][j]]);

                findNewRibbons(i, j - 1, d, dMin1, dCenter1, n, opening, leakOff, mesh, index,
                    activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                    velocities, T, opening[index[i][j - 1]]);

                findNewRibbons(i, j + 1, d, dMin1, dCenter1, n, opening, leakOff, mesh, index,
                    activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                    velocities, T, opening[index[i][j + 1]]);

                findNewRibbons(i + 1, j - 1, d, dMin2, dCenter2, n, opening, leakOff, mesh,
                    index, activeElements, elementIsActive, activationTime, ribbons, oldRibbons, distances,
                    velocities, T, opening[index[i + 1][j - 1]]);

                findNewRibbons(i + 1, j + 1, d, dMin2, dCenter2, n, opening, leakOff, mesh,
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
            ai::saveVector("PREOPEN",opening );
            if(savedSize != activeElements.size() && 1/*meshIsNotExhausted*/){
                buildPartialInfluenceMatrix(influenceMatrix, activeElements,
                    opening, openingNew, partialInfluenceMatrix, index
                );
            }

            // Сохраняем параметры трещины

            calculateCrackGeometry(mesh, distances, length, height);

            std::cout<<"length = "<<length<<"   heigth = "<<height<<std::endl;

            if(Nx * dx < (3./2.)*length/2. || (Ny+1) * dy < (3./2.)*height){

                Nx = (3./2.)*std::ceil(length/2.);
                Ny = (3./2.)*std::ceil(height)-1;
                Nz = Nx;
                N_dof = Nx*Ny*Nz;
                std::cout<<"Nx = "<<Nx<<"  Ny  = "<<Ny<<" Nz = "<<Nz<<" N_dof = "<<N_dof<<std::endl;

                createMatrixDiag(partialInfluenceMatrix, N_dof, Nx , Ny, Nz);
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
            //break;
        }
        ai::printMarker();
        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        timeMeasurements[4] += ai::duration(time4, "us");
        #endif

        #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
        auto time5 = ai::time();
        #endif
        if(saveSteps && 0 == globalStep % stepToSave){
            #if defined(BUILD_DLL)
            /// \todo Сохраняться в json
            #else
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
            #endif

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
            for(std::size_t i = 0; i < injection.size(); ++i){
                injection[i][2] *= timeScale / (wn * 60.);
            }
            Kic /= std::pow(mu * E * std::pow(E / timeScale, n), 1. / (n + 2.));
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
            && length > 4.* initialRadius
            && height > 2.* initialRadius
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
    //}

    auto finishTime = ai::time();

    #if defined(MEASURE_TIME) && !defined(BUILD_DLL)
    timeMeasurements[8] = ai::duration(startTime, finishTime, "us");
    #endif

    //#if !defined(BUILD_DLL)
    #if defined(BUILD_DLL)
    if(meshIsNotExhausted=0){
        if(runningFromGUI){
            std::cout << "Progress: 1.0" << std::endl;
        }else{
            ai::showProgressBar(1.);
        }

        std::cout << " "<< std::endl;
    }else{
        if(runningFromGUI=0){
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

	// ai::saveVector("./Results/val1", val1);
	// ai::saveVector("./Results/val2", val2);

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
){

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

    nlohmann::json importData;

    std::vector< std::vector<double> > layers;
    std::vector< std::vector<double> > injection;

    if(!importInitialData(pathToImportFolder, importData)){
        return 22;
    }

    if(0 < importData.size()){
        modelingTime = importData["time"].get<double>();
        pathToLayersFile = std::string();
        pathToInjectionFile = std::string();
    }

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
    //ai::printMatrix(layers);
    //Интерполируем на планаровкую сетку
     //ApproximateLayers(layers);

    // Пересчитываем значения для слоёв в величины СИ

    for(std::size_t i = 0; i < layers.size(); ++i){
        layers[i][2] *= std::pow(10, 6);
        layers[i][3] *= std::pow(10, 9);
        layers[i][5] *= std::pow(10, -6);
    }

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

    size_t injectionIndex = 0;

    double timeToChangeInjection = modelingTime * 1.1;
    fluidInjection = injection[injectionIndex][2];
    proppantInjection = injection[injectionIndex][3];
    double n = injection[injectionIndex][4];
    double mu = injection[injectionIndex][5];
    double fluidDensity = 1000;
    double proppantDensity = 2500;

    ++injectionIndex;

    if(injection.size() > injectionIndex){
        timeToChangeInjection = injection[injectionIndex][0];
    }


    double T0 = 1.;
    double wn = std::pow(
       mu * 2. * std::pow(2. * (2. * n + 1.) / n, 1. / n)
       / (E * std::pow(timeScale, n)), 1. / (n + 2.)
    );
    const double gammaR =1. / 3. * (1. + n / (n + 2.));

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
            << "  E\' = " << E << ", Kic = " << Kic << ";" << std::endl
            << "  time = " << modelingTime << ", time step = " << timeStep
            << ", time scale = " << timeScale << ";"
            << std::endl;
    }

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

    if(0 < importData.size()){
        cellSize = round(
            initialRadius / importData["mesh"]["cell"]["length"].get<double>()
        );
        meshSize = floor(
            (importData["mesh"]["length"].get<double>() - 1.) / cellSize
        );
    }

    // Автоматически задаём размер ячейки, если требуется

    if(-1. == cellSize){
        //cellSize = initialRadius / 5. - epsilon;
    }

    // Выводим входные параметры (продолжение)

    if(!runningFromGUI){
        std::cout << "  mesh = " << meshSize << ", cell = " << cellSize << ","
            << " scaling = " << meshScalingCounter << "." << std::endl;
        std::cout << "Initial regime: " << regimeName() << "." << std::endl;
    }

    ///////// Костыль!!!!
    //initialRadius = 10.;
 /////////////////////
    dt = 0.0001 * std::pow(5. / floor(initialRadius / cellSize), 2);

    if(-1. != timeStep && !override){
        timeStep *= 0.001;

        if(timeStep < dt){
            dt = timeStep;
        }else{
            std::cerr << "Warning: time step is too big. It was changed to "
                << dt << "!"<< std::endl;
        }
    }

    // Проверяем, что на радиус трещины приходится хотя бы пять элементов

    if(5 > floor(initialRadius / cellSize) && !override){
        //cellSize = initialRadius / 5. - epsilon;

        std::cerr << "Warning: cell size is too big. It was changed to "
            << cellSize << "!"<< std::endl;
    }

    std::cout << "Number of cells per initial crack: "
        << floor(initialRadius / cellSize) << "." << std::endl;

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

    const size_t xSize = x.size();
    const size_t ySize = y.size();

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

    // Пересчитываем контрасты по построенной расчётной сетке

    std::vector<double> stress;
    std::vector<double> flatYoungsModulus;
    std::vector<double> leakOff;

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

    // Масштабируем выличины (продолжение)

    for(std::size_t i = 0; i < stress.size(); ++i){
        stress[i] /= (wn * E);
    }
    for(std::size_t i = 0; i < leakOff.size(); ++i){
        leakOff[i] *= std::sqrt(timeScale) / wn;
    }
    for(std::size_t i = 0; i < injection.size(); ++i){
        injection[i][2] *= timeScale / (wn * 60.);
    }
    Kic /= std::pow(mu * E * std::pow(E / timeScale, n), 1. / (n + 2.));

    std::vector<double> zeroVectorX(xSize, 0.);
    std::vector<double> zeroVectorY(ySize, 0.);
    std::vector<size_t> zeroSizeTVectorY(ySize, 0);
    std::vector<double> zeroVectorXY(xSize * ySize, 0.);
    std::vector< std::vector<double> > zeroMatrixXY(xSize, zeroVectorY);


	std::cout<<"xSize = "<<xSize<<std::endl;
	std::cout<<"dx = "<<dx<<std::endl;

    std::vector<double> opening = zeroVectorXY;
	std::vector<std::vector<double> > opep = zeroMatrixXY;
    std::vector< std::vector<size_t> > index(xSize, zeroSizeTVectorY);
    //double co;
    for(size_t i = 0; i < xSize; ++i){
        for(size_t j = 0; j < ySize; ++j){
            index[i][j] = i * ySize + j;

            // co = 0.5*cos(0.5 *M_PI* std::sqrt(x[i]*x[i]+y[j]*y[j])/(
            // //xSize*0.8/dx
            // 10.
            // ));
            opening[ index[i][j] ] = getInitialOpening(
                x[i],
                y[j],
                //xSize*0.8/dx,
                initialRadius,
                //10.,
                zP,
                openingAtTheStart
            );

            // if(x[i]*x[i]+y[j]*y[j] <= 0.8*0.8*xSize*xSize)
			//          opening[index[i][j]]=1;
			//opep[i][j] = abs(opening[index[i][j]])>epsilon?co:0.;
        }
    }
	//ai::printVector(opening);
	// ai::printMarker();
	ai::saveVector("iniopening", opening);
	ai::saveMatrix("iniopep",opep);
	//ai::printMarker();
    std::cout << std::endl;
    std::cout << "Building A matrix... ";

    // std::vector< std::vector<double> > influenceMatrix(xSize * ySize,
    //     zeroVectorXY);


   size_t Nx = xSize;

   size_t Ny = 2 * Nx - 1;

   size_t Nz = Nx;

   size_t N_dof = Nx * Ny * Nz;




   std::vector<std::vector<double> > influenceMatrix; // matrix corresponding to finite difference discretization of Laplace equation (AT = b)
   influenceMatrix.resize(N_dof);                     //. Nonzero elements are stored only.

   for(size_t i = 0; i < N_dof; ++i){
       influenceMatrix[i].resize(7);
   }

   // create matrix corresponding to finite difference discretization of Laplace equation (AT = b)

   createMatrixDiag(influenceMatrix, N_dof, Nx , Ny, Nz);



    std::cout << "OK." << std::endl;



    // if(considerElasticModulusContrast){
    //     // std::cout << "Applying elastic modulus contrast... ";
    //     //
    //     // std::vector< std::vector<double> > additiveToInflunceMatrix;
    //     //
    //     // // Заполняю ноликами пока
    //     // additiveToInflunceMatrix.resize(xSize * ySize);
    //     //
    //     // for(std::size_t i = 0; i < xSize * ySize; ++i){
    //     //     additiveToInflunceMatrix[i].resize(xSize * ySize);
    //     //
    //     //     for(std::size_t j = 0; j < xSize * ySize; ++j){
    //     //         additiveToInflunceMatrix[i][j] = 0.;
    //         }
    //     }
        // Заполняю ноликами пока

        // for(std::size_t i = 0; i < xSize; ++i){
        //     for(std::size_t k = 0; k < ySize; ++k){
        //         for(std::size_t j = 0; j < influenceMatrix[0].size(); ++j){
        //             influenceMatrix[index[i][k]][j] *= flatYoungsModulus[k];
        //         }
        //     }
        // }
        //
        // for(std::size_t i = 0; i < influenceMatrix.size(); ++i){
        //     for(std::size_t j = 0; j < influenceMatrix[0].size(); ++j){
        //         influenceMatrix[i][j] += additiveToInflunceMatrix[i][j];
        //     }
        // }
        //
        // additiveToInflunceMatrix.clear();

        //std::cout << "OK." << std::endl;

    //}
    std::cout<< "Initial Radius = "<<initialRadius <<std::endl;
    std::cout << "Mesh: " << mesh.size() << "x" << mesh[0].size() << "."
        << std::endl;

    std::vector< std::vector<bool> > elementIsActive;

    std::vector< std::vector<std::size_t> > activeElements;

    findActiveElements(activeElements, elementIsActive, x, y, initialRadius);

    std::cout << "Active elements: " << activeElements.size() << "."
        << std::endl;

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

    std::cout << "Ribbons: " << ribbons.size() << "." << std::endl;

    std::vector<double> openingNew = opening;
    std::vector< std::vector<double> > partialInfluenceMatrix ;



    ///////////////////////////
    partialInfluenceMatrix = influenceMatrix;
    ///////////////////////////
    //std::cout << "Building partial influence matrix... ";

    // buildPartialInfluenceMatrix(influenceMatrix, activeElements, opening, openingNew,
    //     partialInfluenceMatrix, index);

    //std::cout << "OK." << std::endl;

    double height = 0.;
    double length = 0.;

    std::vector< std::vector<double> > fracture;

    std::size_t stepToCheck = std::round(dx / (20. * dt * initialVelocity));

    bool meshIsNotExhausted = true;

    std::vector< std::vector<double> > velocities = zeroMatrixXY;

    std::vector<double> pressure = zeroVectorXY;
	std::vector<double> concentration = zeroVectorXY;

    omp_set_num_threads(numberOfThreads);

    std::cout << std::endl << "Running threads... OK. Size: "
        << getNumberOfThreads() << "." << std::endl;

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
        override,
        N_dof,
        Nx,
        Ny,
        Nz
    );

    if(0 != returnCode){
        return returnCode;
    }

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

    std::cout << std::endl;
    std::cout << "Efficiency: " << fluidEfficiency << "% [fluid]";
    if(-1. < proppantEfficiency){
        std::cout << ", " << proppantEfficiency << "% [proppant]";
    }
    std::cout << "." << std::endl;

    // Сохраняем сжимающие напряжения в центральном слое

    double nominalStress = 0.;

    for(std::size_t i = 0; i < layers.size(); ++i){
        if(0. > layers[i][0] * layers[i][1]){
            nominalStress = layers[i][2] * std::pow(10., -6);

            break;
        }
    }
    double Z_coordinate = layers[0][0];
    // Сохраняем JSON для интерфейса КиберГРП
    std::string js;

     js = ExportJson(
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
        Z_coordinate,
        "a",
        "b"
    );
    ai::saveLine("out.json", js);
    // SaveJson(
    //     "./Results/output",
    //     opening,
    //     wn,
    //     pressure,
    //     concentration,
    //     index,
    //     T,
    //     timeScale,
    //     fluidEfficiency,
    //     fluidDensity,
    //     proppantDensity,
    //     nominalStress
    // );

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
void SaveJson(
    const std::string filename,
    std::vector<double> &opening,
    double wn,
    std::vector<double> &pressure,
    std::vector<double> &concentration,
    std::vector< std::vector<size_t> > &index,
    double Time,
    double timeScale,
    double fluidEfficiency,
    double fluidDensity,
    double proppantDensity,
    double nominalStress
){
    // Число атмосфер в 1 МПа
    const double atmosphereCoefficient = 9.869;

    const double xSize = index.size();
    const double ySize = index[0].size();

    double azimuth = 0;
    int id = 0;
    int stage_id = 0;
    int num_fluids = 1;
    int num_proppants = 1;


    std::string extension(".json");
    std::ofstream output(filename+ ai::string(Time) + extension);

    if (!output.good()) {
        throw std::runtime_error(
            ai::string("Exception while saving the matrix into the file: ")
            + filename
        );
    }


    output << "{\n";
    output << std::setprecision(14) << "\"Results\": {\n";
    //Smin - минимальное значиние напряжений в интервале перфорации, [атм]
    output << "\t\"Smin\": " << nominalStress * atmosphereCoefficient << ",\n";
    output << "\t\"accumulated data\": {\n";
    //    fluid efficiency - эффективность жидкости в момент времени "time", [%]
    output << "\t  \"fluid efficiency\": " << fluidEfficiency << ",\n";
    //    net pressure - значение чистого давления на устье трещины в слое интеравала перфорации с напряжением "Smin", [атм]
    output << "\t  \"net pressure\": " << pressure[index[i00][j00]] * atmosphereCoefficient << ",\n";
    //    proppant concentration - конентрация пропанта в момент времени "time" на входе в трещину, [кг / куб.м.]
    output << "\t  \"proppant concentration\": " << proppantInjection << ",\n";
    //    rate - расход смеси в момент времени "time" на входе в трещину, [куб.м / мин]
    output << "\t  \"rate\": " << fluidInjection / (timeScale / (wn * 60.)) << ",\n";
    //    time - момент времени на который записан результат расчета
    output << "\t  \"time\": " << Time << "\n";
    output << "\t},\n"; //close accumulated data

    //////////////////////////////////////////////////////////////////////////////////////////
    //     geometry
    //////////////////////////////////////////////////////////////////////////////////////////
    output << "\t\"geometry\": {\n";

    //     branches - обозначение полукрыла трещины(0)
    output << "\t\t\"branches\": [\n\t\t {\n";
    //     cells - расчетные ячейки записанные в последовательности : сверху - вниз(j индексы), слева - направо(i индексы)
    output << "\t\t \"cells\": [\n";
    int activeElements = 0;                 //Число активных элементов (те элементы, где вектор раскрытия не нулевой)

    //////////////////////////////////////////////////////////////////////////////////////////
    //     Расчет числа активных элементов
    //////////////////////////////////////////////////////////////////////////////////////////
    for (int i = 0; i < xSize; i++)
    {
         for (int j = 0; j < ySize; j++)
         {
                if (opening[j*xSize + i] != 0) activeElements++; // если ячейка не раскрылась, то ее не выводим!!!
         }
    }

    int Index=0;                                                 //Индексная переменная - номер выводимого массива данных

    for (int i = 0; i < xSize; i++)
    {
         for (int j = 0; j < ySize; j++)
         {
                if (opening[j*xSize + i] == 0) continue; // если ячейка не раскрылась, то ее не выводим!!!
                Index++;
                output << "\t\t{ \n";
                //     azimuth - азимут расчетной ячейки в пространстве(задается 0)
                output << "\t\t  \"azimuth\": " << azimuth << ",\n";
                //     concentrations - объемная доля концентрации фаз, [д.ед.]
                 output << "\t\t  \"concentrations\": [\n ";
                       //for (int index=0; index< num_proppants; index++) // вывод нескольких проппантов (когда будет несколько проппантов)
                output << "\t\t  " << concentration[index[i][j]] << ",\n";                //Вывод концентрации проппанта
                output << "\t\t  " << 1-concentration[index[i][j]] << "\n";                                   //Вывод концентрации жидкости

                output << "\t\t ],\n"; //close concentrations
                                                     //  dx, dy, dz - размеры расчетной ячейки(dy - раскрытие), [м]
                output << "\t\t  \"dx\": " << dx << ",\n";
                output << "\t\t  \"dy\": " << opening[j*xSize+i] << ",\n";
                output << "\t\t  \"dz\": " << dy << ",\n";
                //     id - 0; stage id - 0
                output << "\t\t  \"i\": " << i << ",\n";
                output << "\t\t  \"id\": " << id << ",\n";
                output << "\t\t  \"j\": " << j << ",\n";
                output << "\t\t  \"stage id\": " << stage_id << ",\n";
                //     x, y, z - координаты центра расчетной ячеки в пространстве, [м]
                output << "\t\t  \"x\": " << i*dx << ",\n";
                output << "\t\t  \"y\": " << j* opening[j*xSize + i] << ",\n";
                output << "\t\t  \"z\": " << j* dy << "\n";
                if ( activeElements != Index) //i != xSize - 1 && j != ySize - 1  &&
                       output << "\t\t },\n";
                else{
                       output << "\t\t }\n";
                   }
         }
    }

    output << "\n\t\t]\n"; //close cells
    output << "\t\t }\n\t\t]\n"; //close branches
    output << "\t},\n"; //close geometry

    //////////////////////////////////////////////////////////////////////////////////////////
    //     grid - параметры расчетной области
    //////////////////////////////////////////////////////////////////////////////////////////
    output << "\t\"grid\": {\n";
    //     nx - общее количество расчетных ячеек по Х
    output << "\t\t  \"nx\": " << xSize << ",\n";
    //     nz - общее количество расчетных ячеек по Z
    output << "\t\t  \"nz\": " << ySize << "\n";
    output << "\t},\n"; //close grid

    //////////////////////////////////////////////////////////////////////////////////////////
    //     slurry - свойства смеси в закачке
    //////////////////////////////////////////////////////////////////////////////////////////
    output << "\t\"slurry\": {\n";

    output << "\t\t\"mass density of components\": [\n";
    //     mass density of components - плотность компонент(в той же последовательности, что и в "concentrations"), [кг / куб.м]
    output << "\t\t " << proppantDensity <<",\n";
    output << "\t\t " << fluidDensity <<"\n";

    output << "\t\t ],\n"; //close mass density of components

    output << "\t\t\"name of components\": [\n";
    //     name of components - названия компонент(в той же последовательности, что и в "concentrations")
    output << "\t\t \"Propant 1\",\n";
    output << "\t\t \"Fluid 1\"\n";
    output << "\t\t ],\n"; //close name of componentss

    //     number of fluids - общее количество флюидов
    output << "\t\t  \"number of fluids\": " << num_fluids << ",\n";
    //     number of proppants - общее количество пропантов
    output << "\t\t  \"number of proppants\": " << num_proppants << "\n";
    output << "\t}\n"; //close slurry

    output << "},\n"; //close Results

    //coordinates - координаты центра интервала перфорации, [м]
    output << "\"coordinates\": {\n";
    output << "\t\t  \"x\": " << i00 * dx << ",\n";
    output << "\t\t  \"y\": " << j00 * dy << ",\n";
    output << "\t\t  \"z\": " << 0 << "\n";

    output << "}\n";
    output << "}\n";


    output.close();
}
