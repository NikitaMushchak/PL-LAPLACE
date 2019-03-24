#pragma once

#include <cmath>
#include <string>
#include <vector>


#if !defined M_PI
    /*!
        \brief Число пи
        \details Число пи
    */
    #define M_PI 3.14159265358979323846264338327950288
#endif

/*!
 \brief Тип внешних для трещины сеточных элементов
 \details Тип внешних для трещины сеточных элементов
*/
#define OUTSIDE 0

/*!
 \brief Тип внутренних для трещины сеточных элементов, близких к фронту трещины
 \details Тип внутренних для трещины сеточных элементов, близких к фронту
 трещины
*/
#define RIBBON  1

/*!
 \brief Тип внутренних для трещины сеточных элементов
 \details Тип внутренних для трещины сеточных элементов
*/
#define CHANNEL 2

/*!
 \brief Режим доминирующих утечек
 \details Режим доминирующих утечек
*/
#define LEAK_OFF -1

/*!
 \brief Режим доминирующей трещиностойкости
 \details Режим доминирующей трещиностойкости
*/
#define TOUGHNESS 0

/*!
 \brief Режим доминирующей вязкости
 \details Режим доминирующей вязкости
*/
#define VISCOSITY 1

/*!
 \brief Текущий режим распространения
 \details Текущий режим распространения
*/
extern int regime;

/*!
 \brief Шаг по расстоянию (абсцисса)
 \details Шаг по расстоянию (абсцисса)
*/
extern double dx;

/*!
 \brief Шаг по расстоянию (ордината)
 \details Шаг по расстоянию (ордината)
*/
extern double dy;

/*!
 \brief Шаг по времени
 \details Шаг по времени
*/
extern double dt;

/*!
 \brief Текущее значение закачки жидкости
 \details Текущее значение закачки жидкости
*/
extern double fluidInjection;

/*!
 \brief Текущее значение закачки пропанта
 \details Текущее значение закачки пропанта
*/
extern double proppantInjection;

/*!
 \brief Эффективный плоский модуль Юнга
 \details Эффективный плоский модуль Юнга
*/
extern double E;

/*!
 \brief Эффективная трещиностойкость
 \details Эффективная трещиностойкость
*/
extern double Kic;

/*!
 \brief Значение амплитуды для асимптотического зонтика
 \details Значение амплитуды для асимптотического зонтика
*/
extern double Amu;

/*!
 \brief Степенной коэфициент для асимптотического зонтика
 \details Степенной коэфициент для асимптотического зонтика
*/
extern double alpha;

/*!
 \brief Пренебрежимо малая величина
 \details Пренебрежимо малая величина, используется для сравнений
*/
extern double epsilon;

/*!
 \brief Координата источника
 \details Координата источника по абсциссе
*/
extern size_t i00;

/*!
 \brief Координата источника
 \details Координата источника по ординате
*/
extern size_t j00;

/*!
 \brief Класс граничного элемента
 \details Класс, задающий граничный элемент (координаты – индексы ячейки) с
 оператором сравнения двух граничных элементов
*/
class Ribbon{
    public:
        /*!
            \brief Индекс элемента
            \details Индекс элемента по абсциссе
        */
        size_t i;
        /*!
            \brief Индекс элемента
            \details Индекс элемента по ординате
        */
        size_t j;

        /*!
            \brief Конструктов
            \details Конструктор принимает на вход значения индексов
            \param[in] i - индекс элемента по абсциссе
            \param[in] j - индекс элемента по ординате
        */
        Ribbon(const size_t i, const size_t j): i(i), j(j){};

        /*!
            \brief Сравнение двух элементов
            \details Метод проверяет два граничных элемента на совпадению их
            индексов
            \param[in] ribbon - граничный элемент для сравнения
            \return Совпадают ли два граничных элемента
        */
        bool operator == (const Ribbon ribbon) const{
            return this->i == ribbon.i && this->j == ribbon.j;
        }
};

/*!
 \brief Класс ячейки сетки
 \details Класс, задающий ячейку расчётной сетки (координаты центра, тип)
*/
class Cell{
    public:
        /*!
            \brief Координата ячейки
            \details Координата центра ячейки по абсциссе
        */
        double x = 0;
        /*!
            \brief Координата ячейки
            \details Координата центра ячейки по ординате
        */
        double y = 0;

        /*!
            \brief Тип элемента
            \details Тип элемента (внешний, граничный, внутренний)
        */
        int type = OUTSIDE;

        /*!
            \brief Задание координат ячейки
            \details Метод записывает значения координат центра ячейки в
            экземпляр класса
            \param[in] x - координата по абсциссе
            \param[in] y - координата по ординате
        */
        void setCoordinates(const double x, const double y){
            this->x = x;
            this->y = y;
        }
};
void ApproximateLayers(
    std::vector<std::vector<double> >&layers
    );
/*!
 \brief Название режима распространения
*/
std::string regimeName();

/*!
 \brief Вывод лого программы
*/
void printLogo();

/*!
 \brief Основная расчетная функция
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
    const bool override
);

/*!
 \brief Нахождение длины и высоты трещины
*/
void calculateCrackGeometry(
    std::vector< std::vector<Cell> > &mesh,
    std::vector< std::vector<double> > &distances,
    double &length,
    double &height
);

/*!
 \brief Вычисление эффективности жидкости и проппанта
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
);

/*!
 \brief Сохранение результатов в формате json
*/
void SaveJson(
    const std::string filename,
    std::vector<double> &Wk,
    double wn,
    std::vector<double> &pressure,
    std::vector<double>&concentration,
    std::vector< std::vector<size_t> > &index,
    double Time,
    double timeScale,
    double fluidEfficiency,
    double fluidDensity,
    double proppantDensity,
    double nominalStress
);
