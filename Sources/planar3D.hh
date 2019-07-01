#pragma once

#include <cmath>
#include <string>
#include <vector>

#if defined(USE_EIGEN)
    #include "eigen/Dense"
    #include "eigen/Sparse"
    #include "eigen/Core"
#endif

// #define BUILD_DLL

#include "io.hh"

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

class Fluid_Param
{
	//"0": "KEY",
	//"1" : "CODE",
	//"2" : "NAME",
	//"3" : "SPECIFIC_GRAVITY",		Относительная плотность[д.ед.]
	//"4" : "SHEAR_RATE",			Скорость сдвига[1 / с]
	//"5" : "SPECIFIC_HEAT",		Удельная теплоемкость[Дж / кг * К]
	//"6" : "RHEOLOGY_TABLE",
	//"7" : "FRICTION_TABLE"

public:
	std::string KEY = std::string();
	std::string CODE = std::string();
	//Вектор зависимости реология жидкости от температуры
	//[1] - температура
	//[2] - K [Па*cек^n] динамическая вязкость
	//[3] - n реология
	//[4] - Эффективная вязкость [Па*cек^n]
	std::vector <std::vector<double>> rheology;

	double SPECIFIC_GRAVITY = -1.;						//Относительная плотность [д.ед.]
	double SHARE_RATE = -1.;							//Скорость сдвига [1/с]
	double SPECIFIC_HEAT = -1.;							//Удельная теплоемкость [Дж/кг*К]



														//конструктор
	void setReology(
		std::string KEY,
		std::string CODE,
		double SPECIFIC_GRAVITY,
		double SHARE_RATE,
		double SPECIFIC_HEAT,
		std::vector <std::vector<double>> rheology

	) {
		this->KEY = KEY;
		this->CODE = CODE;
		this->SPECIFIC_GRAVITY = SPECIFIC_GRAVITY;
		this->SPECIFIC_HEAT = SPECIFIC_HEAT;
		this->rheology = rheology;
	};
};

class Proppant_Param
{
	//"1" : "CODE",
	//"4" : "DIAMETER",
	//"5" : "SPECIFIC_GRAVITY",
	//"6" : "BULK_GRAVITY",
	//"7" : "CRITICAL_DENSITY",
	//"8" : "SPECIFIC_HEAT",
	//"9" : "EMBEDMENT_RATIO"
	//"10" : "SETTLING_VELOSITY_FACTOR",
	//"11" : "DRIFT_VELOSITY_FACTOR",
public:
	std::string KEY = std::string();
	std::string CODE = std::string();
	double DIAMETER = -1.;								//Средний диаметр пропанта [мм]
	double SPECIFIC_GRAVITY = -1.;						//Относительная плотность [д.ед.]
	double BULK_GRAVITY = -1.;							//Относительная насыпная плотность [д.ед.]
	double CRITICAL_DENSITY = -1.;						//Критическая концентрация [д.ед.]
	double SPECIFIC_HEAT = -1.;							//Удельная теплоемкость [Дж/кг*К]
	double EMBEDMENT_RATIO = -1.;						//Ширина остаточной проводимости [д.ед.]
	double SETTLING_VELOSITY_FACTOR = -1.;				//Коэффициент скорости осаждения [д.ед.]
	double DRIFT_VELOSITY_FACTOR = -1.;					//Коэффициент скорости дрейфа [д.ед.]
	void setCoordinates(const std::string KEY,
		const std::string CODE,
		const double DIAMETER,
		const double SPECIFIC_GRAVITY,
		const double BULK_GRAVITY,
		const double CRITICAL_DENSITY,
		const double SPECIFIC_HEAT,
		const double EMBEDMENT_RATIO,
		const double SETTLING_VELOSITY_FACTOR,
		const double DRIFT_VELOSITY_FACTOR) {
		this->KEY = KEY;
		this->CODE = CODE;
		this->DIAMETER = DIAMETER;
		this->SPECIFIC_GRAVITY = SPECIFIC_GRAVITY;
		this->BULK_GRAVITY = BULK_GRAVITY;
		this->CRITICAL_DENSITY = CRITICAL_DENSITY;
		this->SPECIFIC_HEAT = SPECIFIC_HEAT;
		this->EMBEDMENT_RATIO = EMBEDMENT_RATIO;
		this->SETTLING_VELOSITY_FACTOR = SETTLING_VELOSITY_FACTOR;
		this->DRIFT_VELOSITY_FACTOR = DRIFT_VELOSITY_FACTOR;
	};
};

class Stages
{
public:

	int id_stage;			//номер стадии (от 0)
	int start_time;			//Время начала стадии
	int end_time;			//Время окончания стадии
	//std::string IdStage;	//id стадии (копия идентификатора переданного от МФТИ)
	int index_fluid;
	int index_proppant;
	//конструктор
	void setStages(
		int id_stage,
		int start_time,
		int end_time,
		int index_fluid,
		int index_proppant
	) {
		this->id_stage = id_stage;
		this->start_time = start_time;
		this->end_time = end_time;
		this->index_fluid = index_fluid;
		this->index_proppant = index_proppant;
	};
};


struct DLL_Param
{
	std::string J_String = std::string();
	std::string IdDesign = std::string();
	std::string IdStage = std::string();
	double modelingTime = -1.;
	double dx = -1.;				//	dx - размер расчетных ячеек по x
	double dz = -1.;				//	dz - размер расчетных ячеек по z
	double nx = -1.;				//	nx - количество ячеек по x
	double nz = -1.;				//	nz - количество ячеек по z
	double W_MIN = -1.;				//	W_MIN - минимальное раскрытие трещины при визуализации в GUI
	double AZ = -1.;				//	AZ - Азимут главных напряжений[град]
	double FLUID_TEMPERATURE = -1.;	//	FLUID_TEMPERATURE - температура флюида[град]
	double Emit_time = -1.;			//	emit - периодичность выгрузки результатов расчета в GUI[секунды]
	bool IS_CLOSURE = false;		//	IS_CLOSURE - вкл / выкл расчет закрытия трещины
	bool IS_EMBEDMENT = false;		//	IS_EMBEDMENT - вкл / выкл опции по учету вдавливания пропанта
	bool IS_GRAVITY = false;		//	IS_GRAVITY	- вкл / выкл опции по учету гравитации
	bool IS_SETTLING = false;		//	IS_SETTLING - вкл / выкл опции по учету осаждения пропанта
	bool IS_TEMPERATURE = false;	//	IS_TEMPERATURE - вкл / выкл опции по расчету температуры

	double Z_coordinate = -1.;
	enum State
	{
		Idle = 0,		//не считает
		Running = 1,	//В работе
		Paused = 2,		//Остановлен
	};
	int Dll_State = Idle;
	std::vector<Proppant_Param> ProppantParam;
	std::vector<Fluid_Param> FluidParam;
	std::vector<Stages> Stages;
};




/*!
 \brief Название режима распространения
*/
std::string regimeName();

/*!
 \brief Вывод лого программы
*/
void printLogo();

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
    #if defined(BUILD_DLL)
    const bool override,
    DLL_Param &DLL_Parametrs
    #else
    const bool override
    #endif
);

/*!
 \brief Основная расчетная функция
*/
int planar3D(
    #if defined(BUILD_DLL)
    DLL_Param &DLL_Parametrs,
    std::vector< std::vector<double> > &layers,
    std::vector< std::vector<double> > &injection,
    bool &initialized
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
    const bool override
    #endif
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
