/*!
 \file io.cc
 \brief
 \details
*/
#include <vector>

#include "planar3D.hh"
#include "initialData.hh"
#include "influenceMatrix.hh"

#include "ailibrary/ai.hh"

int completeMesh(
    std::size_t meshSize,
    std::vector<double> &x,
    std::vector<double> &y,
    std::size_t i00,
    std::size_t j00,
    std::vector< std::vector<Cell> > &mesh,
    std::vector< std::vector<std::size_t> > &index,
    std::vector<double> &opening,
    std::vector<double> &pressure,
    std::vector<double> &concentration,
    std::vector< std::vector<double> > &distances,
    std::vector< std::vector<bool> > &elementIsActive,
    std::vector<double> &activationTime,
    std::vector< std::vector<double> > &layers,
    std::vector<double> &flatYoungsModulus,
    std::vector<double> &stress,
    std::vector<double> &leakOff,
    std::vector<double> &toughness,
    std::size_t &xSize,
    std::size_t &ySize,
    bool considerElasticModulusContrast,
    std::vector< std::vector<double> > &influenceMatrix
){
    
    if(meshSize < x.size()){
        return 0;
    }
    double axMax = (meshSize - 1) * dx;
    for(double i = x[x.size() - 1] + dx; i <= axMax + epsilon; i += dx){
        x.push_back(i);
    }
    for(double i = y[y.size() - 1] + dy; i <= axMax + epsilon; i += dy){
        y.push_back(i);
    }
    double minimumY = y[0];
    std::vector<double>::iterator iterator = y.begin();
    for(double i = -axMax; i <= minimumY - dy + epsilon; i += dy){
        y.insert(iterator, i);
        iterator += 1;
    }

    std::size_t oldXSize = xSize;
    std::size_t oldYSize = ySize;
    
    xSize = x.size();
    ySize = y.size();
    i00 = 0;
    j00 = floor(0.5 * ySize);
    
    std::size_t maxIndex = 0;
    
    for(std::size_t i = 0; i < index.size(); ++i){
        for(std::size_t j = 0; j < index[0].size(); ++j){
            if(index[i][j] > maxIndex){
                maxIndex = index[i][j];
            }
        }
    }

    mesh.resize(xSize);
    index.resize(xSize);
    distances.resize(xSize);
    elementIsActive.resize(xSize);
    influenceMatrix.resize(xSize * ySize);

    std::size_t k = 1;
    std::size_t l = 0;

    for(size_t i = 0; i < oldXSize; ++i){
        mesh[i].resize(ySize);
        distances[i].resize(ySize);
        elementIsActive[i].resize(ySize);

        for(size_t j = oldYSize; j < ySize; ++j){
            distances[i][j] = 0.;
            elementIsActive[i][j] = false;
        }
        for(l = 0; l < (ySize - oldYSize) / 2; ++l){
            index[i].insert(index[i].begin() + l, maxIndex + k);
            
            ++k;
        }
        for(l += oldYSize; l < ySize; ++l){
            index[i].push_back(maxIndex + k);
            
            ++k;
        }
    }
    for(size_t i = oldXSize; i < xSize; ++i){
        mesh[i].resize(ySize);
        index[i].resize(ySize);
        distances[i].resize(ySize);
        elementIsActive[i].resize(ySize);

        for(size_t j = 0; j < ySize; ++j){
            distances[i][j] = 0.;
            index[i][j] = maxIndex + k;
            elementIsActive[i][j] = false;
            
            ++k;
        }
    }
    for(std::size_t i = 0; i < xSize; ++i){
        for(std::size_t j = 0; j < ySize; ++j){
            mesh[i][j].setCoordinates(x[i], y[j]);
        }
    }
    for(std::size_t i = 0; i < influenceMatrix.size(); ++i){
        influenceMatrix[i].resize(xSize * ySize);
    }
    buildInfluenceMatrix(influenceMatrix, xSize, ySize);
    
    ai::printMatrix(index);
    
    opening.resize(xSize * ySize);
    pressure.resize(xSize * ySize);
    concentration.resize(xSize * ySize);
    activationTime.resize(xSize * ySize);

    for(std::size_t i = oldXSize * oldYSize; i < xSize * ySize; ++i){
        opening[i] = 0.;
        pressure[i] = 0.;
        concentration[i] = 0.;
        activationTime[i] = 0.;
    }
    
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

    // for(std::size_t i = 0; i < stress.size(); ++i){
    //     stress[i] /= (wn * E);
    // }
    // for(std::size_t i = 0; i < leakOff.size(); ++i){
    //     leakOff[i] *= std::sqrt(timeScale) / wn;
    // }
    // for(std::size_t i = 0; i < toughness.size(); ++i){
    //     toughness[i] /= std::pow(
    //         mu * E * std::pow(E / timeScale, n),
    //         1. / (n + 2.)
    //     );
    // }

    ai::printLine(
        ai::string("Completing mesh... Added ")
        + ai::string(k - 1) + ai::string(" new cells.")
    );
    
    return 0;
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
        toughness[i] /= std::pow(
            mu * E * std::pow(E / timeScale, n),
            1. / (n + 2.)
        );
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