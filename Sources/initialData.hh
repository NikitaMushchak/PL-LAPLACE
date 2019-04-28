#pragma once

#include "planar3D.hh"
#include "nlohmann/json.hpp"
/*!
 \brief Нахождение начального раскрытия
*/
double getInitialOpening(
    const double x,
    const double y,
    const double radius,
    const std::vector<double> &zP,
    const std::vector<double> &W0
);

/*!
 \brief Загрузка начальных данных (пласт, план закачки)
*/
bool setInitialData(
    const std::string pathToLayersFile,
    std::vector< std::vector<double> > &layers,
    const std::string pathToInjectionFile,
    std::vector< std::vector<double> > &injection
);

/*!
 \brief Импорт начальных данных из файла формата JSON
*/
bool importInitialData(
    const std::string pathToImportFolder,
    nlohmann::json &importData
);

/*!
 \brief Пересчет плана закачки
*/
bool recalculateInjection(
    std::vector< std::vector<double> > &injection,
    const double modelingTime
);

/*!
 \brief Пересчет контраста напряжений
*/
bool recalculateStressContrast(
    const std::vector< std::vector<double> > &layers,
    std::vector<double> &stress,
    const std::vector<double> &y
);

/*!
 \brief Пересчет контраста упругих модулей в эффективное значение
*/
bool recalculateElasticModulusContrast(
    const std::vector< std::vector<double> > &layers,
    double &flatYoungsModulus,
    const bool considerElasticModulusContrast
);

/*!
 \brief Пересчет контраста упругих модулей
*/
bool recalculateElasticModulusContrast(
    const std::vector< std::vector<double> > &layers,
    double &flatYoungsModulus,
    std::vector<double> &flatYoungsModulusContrast,
    const std::vector<double> &y,
    const bool considerElasticModulusContrast
);

/*!
 \brief Пересчет контраста утечек
*/
bool recalculateLeakOffContrast(
    const std::vector< std::vector<double> > &layers,
    std::vector<double> &leakOff,
    const std::vector<double> &y
);

/*!
 \brief Сохранение начальных параметров в файл формата JSON
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
);
