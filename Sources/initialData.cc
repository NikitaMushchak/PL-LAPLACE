#include <vector>
#include <iostream>

#include "nlohmann/json.hpp"
#include "ailibrary/ai.hh"

#include "planar3D.hh"
#include "linearOperator.hh"
#include "automodel.hh"

/*!
 \details Функция определяет начальное раскрытие в точке круглой трещины
 известного радиуса, проецируя автомодельное решение

 \param[in] x - координата центра ячейки по абсциссе
 \param[in] y - координата центра ячейки по ординате
 \param[in] radius - радиус трещины
 \param[in] zP - Вектор расчетной области для массштабирования автомодельного решения
 \param[in] W0 - Вектор начальных раскрытий трещины на момент окончания автомодельного решения
 \return Значение раскрытия
*/
double getInitialOpening(
    const double x,
    const double y,
    const double radius,
    const std::vector<double> &zP,
    const std::vector<double> &W0
){
    const double distance = std::sqrt(std::pow(x, 2) + std::pow(y, 2))
        / radius;

    double opening = 0;

    if(1 > distance){
        for(size_t i = 0; i < zP.size() - 1; ++i){
            if(distance >= zP[i] && distance < zP[i + 1]){
                opening = W0[i] * zP[i + 1] - W0[i + 1] * zP[i]
                    - distance * (W0[i] - W0[i + 1]);
                opening /= (zP[i + 1] - zP[i]);
            }
        }
    }

    return opening;
}

/*!
 \details Функция подгружает данные о пласте и план закачки из соответствующих
 текстовых файлов, проверяя их на совместимость

 \param[in] pathToLayersFile - путь (относительный или абсолютный) к
 расположению файла слоёв
 \param[out] layers - матрица, хранящая данные о слоях
 \param[in] pathToInjectionFile - путь (относительный или абсолютный) к
 расположению файла закачки
 \param[out] injection - матрица, хранящая план закачки
 \return Корректно ли отработала функция
*/
bool setInitialData(
    const std::string pathToLayersFile,
    std::vector< std::vector<double> > &layers,
    const std::string pathToInjectionFile,
    std::vector< std::vector<double> > &injection
){
    if(std::string() != pathToInjectionFile){
        try{
            std::cout << "Injection... ";

            ai::parseFileInMatrix(pathToInjectionFile, ' ', injection);

            if(2 > injection.size() || 6 > injection[0].size()){
                std::cerr << "Cannot open '" << pathToInjectionFile
                    << ", file is not appropriate." << std::endl;

                return false;
            }
            if(injection[0][0] < 1){
                std::cerr << "Cannot open '" << pathToInjectionFile << ", "
                << "standard version should be at least 1.0.0 "
                << "for this type of file." << std::endl;

                return false;
            }

            injection.erase(injection.begin());

            std::cout << "OK. Size: " << injection.size() << "." << std::endl;
        }catch(const std::exception &error){
            std::cout << "Fail!" << std::endl;

            std::cerr << "Cannot open '" << pathToInjectionFile << "', "
                << "provide file and restart." << std::endl;

            std::cerr << error.what() << std::endl;

            return false;
        }
    }else{
        std::cout << "Injection... No file." << std::endl;
    }

    if(std::string() != pathToLayersFile){
        try{
            std::cout << "Layers... ";

            ai::parseFileInMatrix(pathToLayersFile, ' ', layers);

            if(2 > layers.size() || 6 > layers[0].size()){
                std::cerr << "Cannot open '" << pathToLayersFile
                    << ", file is not appropriate." << std::endl;

                return false;
            }
            if(layers[0][0] < 1){
                std::cerr << "Cannot open '" << pathToInjectionFile << ", "
                    << "standard version should be at least 1.0.0 "
                    << "for this type of file." << std::endl;

                return false;
            }

            layers.erase(layers.begin());

            std::cout << "OK. Size: " << layers.size() << "." << std::endl;
        }catch(const std::exception &error){
            std::cout << "Fail!" << std::endl;

            std::cerr << "Cannot open '" << pathToLayersFile << "', "
                << "provide file and restart." << std::endl;

            std::cerr << error.what() << std::endl;

            return false;
        }
    }else{
        std::cout << "Layers... No file." << std::endl;
    }

    return true;
}

/*!
 \details Функция подгружает начальные данные из файла формата JSON

 \param[in] pathToImportFolder - путь (относительный или абсолютный) к
 расположению директории с данными для импорта
 \param[out] importData - json-объект с начальными данными
 \return Корректно ли отработала функция
*/
bool importInitialData(
    const std::string pathToImportFolder,
    nlohmann::json &importData
){
    std::cout << "Reading initial conditions from text files:" << std::endl;

    if(std::string() != pathToImportFolder){
        try{
            std::cout << "Parameters (JSON)... ";

            std::string importDataString = std::string();

            ai::parseFileIntoString(
                pathToImportFolder + std::string("/parameters.json"),
                importDataString
            );

            importData = nlohmann::json::parse(importDataString);

            std::cout << "OK. Size: " << importData.size() << "." << std::endl;
        }catch(const std::exception &error){
            std::cout << "Fail!" << std::endl;

            std::cerr << "Cannot open '" << pathToImportFolder
                << "/parameters.json', provide file and restart." << std::endl;

            std::cerr << error.what() << std::endl;

            return false;
        }
    }else{
        std::cout << "Parameters... No file." << std::endl;
    }

    return true;
}

/*!
 \details Функция модифицирует полученные данные из плана закачки с учетом
 общего времени моделирования

 \param[out] injection - матрица, хранящая план закачки
 \param[in] modelingTime - расчетное время
 \return Корректно ли отработала функция
*/
bool recalculateInjection(
    std::vector< std::vector<double> > &injection,
    const double modelingTime
){
    if(1 > injection.size()){
        std::cerr << "Format error in the injection file: matrix is empty."
            << std::endl;

        return false;
    }

    for(std::size_t i = 0; i < injection.size(); ++i){
        if(6 != injection[i].size()){
            std::cerr << "Format error in the injection file: matrix size "
                << "must be six." << std::endl;

            return false;
        }
    }

    for(size_t i = 0; i + 1 < injection.size(); ++i){
        if(injection[i][1] < injection[i + 1][0]){
            injection.insert(
                injection.begin() + i + 1,
                std::vector<double>{
                    injection[i][1],
                    injection[i + 1][0],
                    0.,
                    0.,
                    injection[i][4],
                    injection[i][5]
                }
            );
        }
    }

    if(modelingTime > injection[injection.size() - 1][1]){
        injection.push_back(
            std::vector<double>{
                injection[injection.size() - 1][1],
                modelingTime,
                0.,
                0.,
                injection[injection.size() - 1][4],
                injection[injection.size() - 1][5]
            }
        );
    }

    return true;
}

/*!
 \details Функция производит пересчет данных о сжимающих напряжениях в слоях
 с учетом заданной расчетной области и линейной интерполяции

 \param[in] layers - матрица, хранящая данные о слоях
 \param[out] stress - вектор учитываемых сжимающих напряжений
 \param[in] y - вектор координат центров ячеек по ординате
 \return Корректно ли отработала функция
*/
bool recalculateStressContrast(
    const std::vector< std::vector<double> > &layers,
    std::vector<double> &stress,
    const std::vector<double> &y
){
    std::vector<double> coordinates;
    std::vector<double> coefficients;

    double nominalStress = 0.;

    for(int i = layers.size() - 1; i >= 0; --i){
        coordinates.push_back(layers[i][0]);
        coordinates.push_back(layers[i][1] - 0.001);
        coefficients.push_back(layers[i][2]);
        coefficients.push_back(layers[i][2]);

        if(0. > layers[i][0] * layers[i][1]){
            nominalStress = layers[i][2];
        }
    }

    try{
        createLinearOperator(coordinates, coefficients);
    }catch(const std::exception &error){
        std::cerr << "Error during calculations of the linear operator for "
            << "stress contrast." << std::endl;

        std::cerr << error.what() << std::endl;

        return false;
    }

    stress.clear();

    for(std::size_t i = 0; i < y.size(); ++i){
        stress.push_back(
            calculateValueWithLinearOperator(
                y[i],
                coordinates,
                coefficients
            ) - nominalStress
        );
    }

    return true;
}

/*!
 \details Функция производит пересчет данных об упругих модулях в слоях
 в эффективное значение

 \param[in] layers - матрица, хранящая данные о слоях
 \param[out] flatYoungsModulus - эффективный плоский модуль Юнга
 \param[in] considerElasticModulusContrast - флаг подключения модуля учёта
 влияния слоистости по упругим модулям
 \return Корректно ли отработала функция
*/
bool recalculateElasticModulusContrast(
    const std::vector< std::vector<double> > &layers,
    double &flatYoungsModulus,
    const bool considerElasticModulusContrast
){
    flatYoungsModulus = 0.;

    if(considerElasticModulusContrast){
        for(std::size_t i = 0; i < layers.size(); ++i){
            if(0. > layers[i][0] * layers[i][1]){
                flatYoungsModulus = layers[i][3]
                    / (1. - std::pow(layers[i][4], 2));

                break;
            }
        }
    }else{
        for(std::size_t i = 0; i < layers.size(); ++i){
            flatYoungsModulus += std::abs(layers[i][1] - layers[i][0])
                * layers[i][3] / (1. - std::pow(layers[i][4], 2));
        }

        flatYoungsModulus /= std::abs(
            layers[0][1] - layers[layers.size() - 1][0]
        );
    }

   if(epsilon > flatYoungsModulus){
        std::cerr << "Error during calculations of the flat Young's modulus."
            << std::endl;

        return false;
    }

    return true;
}

/*!
 \details Функция производит пересчет данных об упругих модулях в слоях
 с учетом заданной расчетной области и линейной интерполяции

 \param[in] layers - матрица, хранящая данные о слоях
 \param[in] flatYoungsModulus - эффективный плоский модуль Юнга
 \param[out] flatYoungsModulusContrast - вектор учитываемых плоских модулей Юнга
 \param[in] y - вектор координат центров ячеек по ординате
 \param[in] considerElasticModulusContrast - флаг подключения модуля учёта
 влияния слоистости по упругим модулям
 \return Корректно ли отработала функция
*/
bool recalculateElasticModulusContrast(
    const std::vector< std::vector<double> > &layers,
    double &flatYoungsModulus,
    std::vector<double> &flatYoungsModulusContrast,
    const std::vector<double> &y,
    const bool considerElasticModulusContrast
){
    if(considerElasticModulusContrast){
        std::vector<double> coordinates;
        std::vector<double> coefficients;

        for(int i = layers.size() - 1; i >= 0; --i){
            double value = layers[i][3] / (1. - std::pow(layers[i][4], 2));

            coordinates.push_back(layers[i][0]);
            coordinates.push_back(layers[i][1] - 0.001);
            coefficients.push_back(value);
            coefficients.push_back(value);
        }

        try{
            createLinearOperator(coordinates, coefficients);
        }catch(const std::exception &error){
            std::cerr << "Error during calculations of the linear operator for "
                << "elastic modulus contrast." << std::endl;

            std::cerr << error.what() << std::endl;

            return false;
        }

        flatYoungsModulusContrast.clear();

        for(std::size_t i = 0; i < y.size(); ++i){
            flatYoungsModulusContrast.push_back(
                calculateValueWithLinearOperator(
                    y[i],
                    coordinates,
                    coefficients
                ) / flatYoungsModulus
            );
        }
    }

    return true;
}

/*!
 \details Функция производит пересчет данных об утечках в слоях
 с учетом заданной расчетной области и линейной интерполяции

 \param[in] layers - матрица, хранящая данные о слоях
 \param[out] leakOff - вектор учитываемых коэффициентов Картера
 \param[in] y - вектор координат центров ячеек по ординате
 \return Корректно ли отработала функция
*/
bool recalculateLeakOffContrast(
    const std::vector< std::vector<double> > &layers,
    std::vector<double> &leakOff,
    const std::vector<double> &y
){
    std::vector<double> coordinates;
    std::vector<double> coefficients;

    for(int i = layers.size() - 1; i >= 0; --i){
        coordinates.push_back(layers[i][0]);
        coordinates.push_back(layers[i][1] - 0.001);
        coefficients.push_back(layers[i][5]);
        coefficients.push_back(layers[i][5]);
    }

    try{
        createLinearOperator(coordinates, coefficients);
    }catch(const std::exception &error){
        std::cerr << "Error during calculations of the linear operator for "
            << "leak-off contrast." << std::endl;

        std::cerr << error.what() << std::endl;

        return false;
    }

    leakOff.clear();

    for(std::size_t i = 0; i < y.size(); ++i){
        leakOff.push_back(
            calculateValueWithLinearOperator(
                y[i],
                coordinates,
                coefficients
            )
        );
    }

    return true;
}

/*!
 \details Функция сохраняет начальные параметры в файл формата JSON по
 установленному шаблону

 \param[in] filename - имя выходного файла
 \param[in] modelingTime - расчетное время
 \param[in] meshHeight - высота сетки в ячейках
 \param[in] meshLength - ширина сетки в ячейках
 \param[in] cellHeight - высота ячейки в метрах
 \param[in] cellLength - ширина ячейка в метрах
 \param[in] injection - матрица, хранящая план закачки
 \param[in] layers - матрица, хранящая данные о слоях
 \param[in] tab - значение тектового разделителя (опционально)
*/
void saveInitialData(
    const std::string filename,
    const double modelingTime,
    const std::size_t meshHeight,
    const std::size_t meshLength,
    const double cellHeight,
    const double cellLength,
    std::vector< std::vector<double> > &injection,
    std::vector< std::vector<double> > &layers,
    const std::string tab = std::string("    ")
){
    std::ofstream output(filename + ai::string(".json"));

    if(!output.good()){
        throw std::runtime_error(
            ai::string("Exception while saving the data into the file: ")
            + filename
        );
    }

    std::size_t indentLevel = 0;

    auto indent = [&indentLevel, &tab]() -> std::string
    {
        std::string output;

        for(std::size_t i = 0; i < indentLevel; ++i){
            output += tab;
        }

        return output;
    };

    output << "{" << std::endl;
    ++indentLevel;

    output << indent() << "\"model\": \"Planar3D\"," << std::endl;
    output << indent() << "\"time\": " << modelingTime << "," << std::endl;

    output << indent() << "\"mesh\": {" << std::endl;
    ++indentLevel;
        output << indent() << "\"height\": " << meshHeight << "," << std::endl;
        output << indent() << "\"length\": " << meshLength << "," << std::endl;
        output << indent() << "\"cell\": {" << std::endl;
        ++indentLevel;
            output << indent() << "\"height\": " << cellHeight << ","
                << std::endl;
            output << indent() << "\"length\": " << cellLength
                << std::endl;
        --indentLevel;
        output << indent() << "}" << std::endl;
    --indentLevel;
    output << indent() << "}," << std::endl;

    output << indent() << "\"injection\": {" << std::endl;
    ++indentLevel;
        for(size_t i = 0; i < injection.size(); ++i){
            if(injection[i].size() != 6){
                throw std::runtime_error(
                    ai::string("Exception in size of the vector: injection")
                );
            }
            output << indent() << "\"" << ai::string(i) << "\": {"
                << std::endl;
            ++indentLevel;
                output << indent() << "\"start time\": " << injection[i][0]
                    << "," << std::endl;
                output << indent() << "\"stop time\": " << injection[i][1]
                    << "," << std::endl;
                output << indent() << "\"pumping rate\": "
                    << injection[i][2] << "," << std::endl;
                output << indent() << "\"bulk proppant\": "
                    << injection[i][3] << "," << std::endl;
                output << indent() << "\"rheology index\": " << injection[i][4]
                    << "," << std::endl;
                output << indent() << "\"dynamic viscosity\": "
                    << injection[i][5] << std::endl;
            --indentLevel;
            output << indent() << "}";
            if(injection.size() > i + 1){
                output << ",";
            }
            output << std::endl;
        }
    --indentLevel;
    output << indent() << "}," << std::endl;

    output << indent() << "\"layers\": {" << std::endl;
    ++indentLevel;
        for(size_t i = 0; i < layers.size(); ++i){
            if(layers[i].size() != 6){
                throw std::runtime_error(
                    ai::string("Exception in size of the vector: layers")
                );
            }

            output << indent() << "\"" << ai::string(i) << "\": {"
                << std::endl;
            ++indentLevel;
                output << indent() << "\"y1\": " << layers[i][1] << ","
                    << std::endl;
                output << indent() << "\"y2\": " << layers[i][0] << ","
                    << std::endl;
                output << indent() << "\"stress\": " << layers[i][2]
                    << "," << std::endl;
                output << indent() << "\"Young's modulus\": "
                    << layers[i][3] << "," << std::endl;
                output << indent() << "\"Poisson's ratio\": "
                    << layers[i][4] << "," << std::endl;
                output << indent() << "\"Carter's coefficient\": "
                    << layers[i][5] << std::endl;
            --indentLevel;
            output << indent() << "}";
            if(layers.size() > i + 1){
                output << ",";
            }
            output << std::endl;
        }
    --indentLevel;
    output << indent() << "}" << std::endl;

    output << "}";
}
