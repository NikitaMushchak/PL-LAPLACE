/*!
 \file io.cc
 \brief
 \details
*/
#include <vector>
#include <fstream>

#include "planar3D.hh"

#include "ailibrary/ai.hh"
#include "nlohmann/json.hpp"


void SaveEverything(
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
	const bool override
)
{
	ai::saveVector("1/opening", opening);
	ai::saveVector("1/pressure", pressure);
	ai::saveVector("1/concentration", concentration);
	ai::saveVector("1/openingNew", openingNew);
	ai::saveVector("1/x", x);
	ai::saveVector("1/y", y);
	//ai::saveVector("ribbons",ribbons);
	ai::saveVector("1/activationTime", activationTime);
	ai::saveVector("1/stress", stress);
	ai::saveVector("1/leakOff", leakOff);
	ai::saveMatrix("1/distances", distances);
	ai::saveMatrix("1/velocities", velocities);
//	ai::saveVector("1/flatYoungsModulus", flatYoungsModulus);


//	ai::saveMatrix("1/influenceMatrix", influenceMatrix);
	ai::saveMatrix("1/partialInfluenceMatrix", partialInfluenceMatrix);
	//ai::saveMatrix("mesh",mesh);
	ai::saveMatrix("1/index", index);
	ai::saveMatrix("1/activeElements", activeElements);
	ai::saveMatrix("1/injection", injection);
	ai::saveMatrix("1/layers", layers);
	ai::saveMatrix("1/elementIsActive", elementIsActive);
	ai::saveMatrix("1/fracture", fracture);


	std::vector<double> save;

	save.push_back(fluidDensity);
	save.push_back(proppantDensity);
	save.push_back(n);
	save.push_back(mu);
	save.push_back(wn);
	save.push_back(timeScale);
	save.push_back(timeStep);
	save.push_back(length);
	save.push_back(height);
	save.push_back(initialRadius);
	save.push_back(axMax);
	save.push_back((double)meshScalingCounter);
	save.push_back(T);
	save.push_back(T0);
	save.push_back((double)injectionIndex);
	save.push_back(timeToChangeInjection);
	save.push_back((double)stepToCheck);
	ai::saveVector("1/Save", save);

}

/*!
\brief  Экспорт выходных данных в симулятор из расчетного модуля

\details Функция сохраняет выходные данные в строке  в формате json.

\param[in] Wk - вектор раскрытий
\param[in] wn - массштабный коэффициент
\param[in] pressure - матрица давлений
\param[in] concentration - матрица концентраций проппанта
\param[in] x
\param[in] y
\param[in] dx
\param[in] dy
\param[in] i00
\param[in] j00
\param[in] index - матрица активных элементов
\param[in] Time - Текущее время (модельное)
\param[in] timeScale -
\param[in] fluidEfficiency - эффективность жидкости
\param[in] fluidDensity - плотность жидкости
\param[in] proppantEfficiency - эффективность проппанта
\param[in] proppantDensity - плотность проппанта
\param[in] fluidInjection
\param[in] proppantInjection
\param[in] nominalStress - минимальное значиние напряжений в интервале перфорации, [атм]
\param[in] IdDesign
\param[in] IdStage
*/
std::string ExportJson(
	std::vector<double> &Wk,
	double wn,
	std::vector<double> &pressure,
	std::vector<double>&concentration,
	std::vector<double> &x,
	std::vector<double> &y,
	double dx,
	double dy,
	size_t i00,
	size_t j00,
	std::vector< std::vector<size_t> > &index,
	double Time,
	double timeScale,
	double fluidEfficiency,
	double fluidDensity,
	double proppantDensity,
	double fluidInjection,
	double proppantInjection,
	double nominalStress,
	DLL_Param &DLL_Parametrs
) {
	// Число атмосфер в 1 МПа
	const double atmosphereCoefficient = 9.869;
	//комментарий с лекции !!!!!!!!!!!!!!!
	const double xSize = x.size();
	const double ySize = y.size();

	//////////////////////////////////////////////////////////////////////////////////////////
	//Эти данные пока МФТИ просит заполнять нулями
	//////////////////////////////////////////////////////////////////////////////////////////
	double azimuth = 0;								//азимут расчетной ячейки в пространстве(задается 0)
	int id = 0;										//(задается 0)
	int stage_id = 0;								//(задается 0)
	int num_fluids = 1;								//number of fluids - общее количество жидкостей
	int num_proppants = 0;							//number of proppants - общее количество пропантов


	if (fluidEfficiency > 100)fluidEfficiency = 100.;


	//////////////////////////////////////////////////////////////////////////////////////////
	//Структура выходного файла json
	//////////////////////////////////////////////////////////////////////////////////////////
	//Results
	//	Smin - минимальное значиние напряжений в интервале перфорации, [атм]
	//	accumulated data
	//	fluid efficiency - эффективность жидкости в момент времени "time", [%]
	//	net pressure - значение чистого давления на устье трещины в слое интеравала перфорации с напряжением "Smin", [атм]
	//	proppant concentration - конентрация пропанта в момент времени "time" на входе в трещину, [кг / куб.м.]
	//	rate - расход смеси в момент времени "time" на входе в трещину, [куб.м / мин]
	//	time - момент времени на который записан результат расчета
	//	geometry
	//	branches - обозначение полукрыла трещины(0)
	//	cells - расчетные ячейки записанные в последовательности : сверху - вниз(j индексы), слева - направо(i индексы)
	//	azimuth - азимут расчетной ячейки в пространстве(задается 0)
	//	concentrations - объемная доля концентрации фаз, [д.ед.]
	//	dx, dy, dz - размеры расчетной ячейки(dy - раскрытие), [м]
	//	id - 0; stage id - 0
	//	x, y, z - координаты центра расчетной ячеки в пространстве, [м]
	//	grid - параметры расчетной области
	//	nx - общее количество расчетных ячеек по Х
	//	nz - общее количество расчетных ячеек по Z
	//	slurry - свойства смеси в закачке
	//	mass density of components - плотность компонент(в той же последовательности, что и в "concentrations"), [кг / куб.м]
	//	name of components - названия компонент(в той же последовательности, что и в "concentrations")
	//	number of fluids - общее количество флюидов
	//	number of proppants - общее количество пропантов
	//	coordinates - координаты центра интервала перфорации, [м]
	
	
	
	
	
	std::ostringstream output;

	output << "{\n";
	//Report
	output << "\"Report\":\n ";
	output << "{\n";
	output << "\t\"relativeMassBalanceError\": " << "0.0" << ",\n";	//нужно понять откуда эту цифру брать!!!
	output << "\t\"stageId\": \"" << DLL_Parametrs.IdStage.c_str() << "\"\n";
	output << "},\n"; //close Report
	//Ports
	output << "\t\"" << DLL_Parametrs.IdStage.c_str() << "\" :\n";	output << "{\n";
	output << "\"Ports\":\n [ \n";
	//////////////////////////////////////////////////////////////////////////////////////////
	//Results
	//////////////////////////////////////////////////////////////////////////////////////////
	{
	output << "{\n";
	output << std::setprecision(14) << "\"Results\": {\n";
	//Smin - минимальное значиние напряжений в интервале перфорации, [атм]
	output << "\t\"Smin\": " << nominalStress * atmosphereCoefficient << ",\n";
	output << "\t\"accumulated data\": {\n";
	//	fluid efficiency - эффективность жидкости в момент времени "time", [%]
	output << "\t  \"fluid efficiency\": " << fluidEfficiency << ",\n";
	//	net pressure - значение чистого давления на устье трещины в слое интеравала перфорации с напряжением "Smin", [атм]
	output << "\t  \"net pressure\": " << pressure[index[i00][j00]] * atmosphereCoefficient << ",\n";
	//	proppant concentration - конентрация пропанта в момент времени "time" на входе в трещину, [кг / куб.м.]
	output << "\t  \"proppant concentration\": " << proppantInjection << ",\n";
	//	rate - расход смеси в момент времени "time" на входе в трещину, [куб.м / мин]
	output << "\t  \"rate\": " << fluidInjection / (timeScale / (wn * 60.)) << ",\n";
	//	time - момент времени на который записан результат расчета
	output << "\t  \"time\": " << Time << "\n";
	output << "\t},\n"; //close accumulated data

	//////////////////////////////////////////////////////////////////////////////////////////
	//     geometry
	//     branches
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
			//if (Wk[j*xSize + i] != 0) 
				activeElements++; // если ячейка не раскрылась, то ее не выводим!!!
		}
	}

	int Index = 0;                                                 //Индексная перменная - номер выводимого массива данных

for (int i = 0; i < xSize; i++)
	{
	for (int j = 0; j < ySize; j++)
		{
			//*if (Wk[j*xSize + i] == 0) continue; // если ячейка не раскрылась, то ее не выводим!!!
			Index++;
			output << "\t\t{ \n";
			//     azimuth - азимут расчетной ячейки в пространстве(задается 0)
			output << "\t\t  \"azimuth\": " << azimuth << ",\n";
			//     concentrations - объемная доля концентрации фаз, [д.ед.]
			output << "\t\t  \"concentrations\": [\n ";
			//for (int index=0; index< num_proppants; index++) // вывод нескольких проппантов (когда будет несколько проппантов)
			output << "\t\t  " << concentration[index[i][j]] << ",\n";                //Вывод концентрации проппанта
			output << "\t\t  " << 1 - concentration[index[i][j]] << "\n";             //Вывод концентрации жидкости

			output << "\t\t ],\n"; //close concentrations
								   //  dx, dy, dz - размеры расчетной ячейки(dy - раскрытие), [м]
			output << "\t\t  \"dx\": " << dx << ",\n";
			//if (std::isnan(Wk[j*xSize + i]))
			//	Wk[j*xSize + i] = 100;
			output << "\t\t  \"dy\": " << Wk[i*ySize + j]*wn << ",\n"; // Wk[j*xSize + i]
			output << "\t\t  \"dz\": " << dy << ",\n";
			//     id - 0; stage id - 0
			output << "\t\t  \"i\": " << i << ",\n";
			output << "\t\t  \"id\": \"" << DLL_Parametrs.IdStage.c_str() << "\",\n";
			output << "\t\t  \"j\": " << j << ",\n";
			output << "\t\t  \"stage id\": " << stage_id << ",\n";
			//     x, y, z - координаты центра расчетной ячеки в пространстве, [м]
			output << "\t\t  \"x\": " << i * dx << ",\n";
			output << "\t\t  \"y\": " << j * Wk[j*xSize + i] << ",\n";  //0
			output << "\t\t  \"z\": " << DLL_Parametrs.Z_coordinate + j * dy << "\n";				// Z_coordinate + смещение

			if (activeElements != Index) //i != xSize - 1 && j != ySize - 1  &&
				output << "\t\t },\n";
			else {	//Если выводим последний элемент массива активных элементов - закрываем JSON массив
				output << "\t\t }\n";
			}
		}
	}

	output << "\n\t\t]\n"; //close cells
	
		//////////////////////////////////////////////////////////////////////////////////////////
		//close branches
		//close geometry
		//////////////////////////////////////////////////////////////////////////////////////////
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
		//close grid
		//////////////////////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////////////////////
		//     slurry - свойства смеси в закачке
		//////////////////////////////////////////////////////////////////////////////////////////
		output << "\t\"slurry\": {\n";

		output << "\t\t\"mass density of components\": [\n";
		//     mass density of components - плотность компонент(в той же последовательности, что и в "concentrations"), [кг / куб.м]
		output << "\t\t " << proppantDensity << ",\n";
		output << "\t\t " << fluidDensity << "\n";

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
		output << "\t}\n"; 
		//////////////////////////////////////////////////////////////////////////////////////////
		//close slurry
		//////////////////////////////////////////////////////////////////////////////////////////

	output << "},\n"; 
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	//close Results
	//////////////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////////////
	//coordinates - координаты центра интервала перфорации, [м]
	//////////////////////////////////////////////////////////////////////////////////////////
	output << "\"coordinates\": {\n";
	output << "\t\t  \"x\": " << 0 << ",\n";
	output << "\t\t  \"y\": " << 0 << ",\n";
	output << "\t\t  \"z\": " << DLL_Parametrs.Z_coordinate << "\n";	//Нужно передать сюда координату центра J_D_S_Ports[index]["md"].get<double>()
	output << "},\n";
	output << "\t\t\"id\": \"" << DLL_Parametrs.IdDesign << "\"\n";

	output << "}\n ";
	output << "]\n ";
	//Закрыли Ports
	output << "}\n }\n";
	//Закрыли все

	return	output.str();
	//////////////////////////////////////////////////////////////////////////////////////////
	//End function ExportJson
	//////////////////////////////////////////////////////////////////////////////////////////
}

void ImportJSON(
	DLL_Param &DLL_Parametrs,
	std::vector< std::vector<double> > &layers_new,
	std::vector< std::vector<double> > &injection
)
{
	using json = nlohmann::json;
	std::string J_String = DLL_Parametrs.J_String;
	//	ai::parseFileIntoString("../Planar Challenge - Test 12_init.json", J_String);
	json J_IN_DATA = json::parse(J_String.c_str());

	std::vector<double> Time;
	std::vector<double> TimeStart;
	std::vector<double> TimeEnd;


	//	std::vector< std::vector <double> > layers;
	//Вектор layers:
	//{L_start - middle, L_end - middle, stress, young, poisson, carter * 1000000.}
	//[0] - начало слоя (смещение относительно центрального стоя)
	//[1] - конец слоя (смещение относительно центрального стоя)
	//[2] - горизонтальное напряжение в слое
	//[3] - Модуль Юнга
	//[4] - Коэффициент пуассона
	//[5] - Коэффициент утечек по Картеру

	//struct LAYER {
	//	double	start,
	//			end,
	//			stress,
	//			young,
	//			poisson,
	//			carter;
	//};
	//std::vector< LAYER > layers;
	std::vector< std::vector<double>> layers;
	//int index_middle;					//Индекс центрального слоя в векторе layers
	//LAYER layer;
	//std::vector< std::vector<double> > injection;
	std::vector<double> reologies;
	std::vector<double> viscosities;
	std::vector<double> prop;
	std::vector<double> fluid;
	std::vector<std::string> flui;


	std::vector<double> timefl;


	double maxstress = 0;
	double minstress = 99999999;

	DLL_Parametrs.Emit_time = J_IN_DATA["Design"]["Settings"]["emit"].get<double>();
	DLL_Parametrs.IdDesign = J_IN_DATA["Design"]["id"].get<std::string>();
	json  J_D_Stages = J_IN_DATA["Design"]["Stages"];

	for (int index = 0; index < J_D_Stages.size(); index++)
	{
		//информация о слоях ReservoirFormation - Таблица с планшетом литологии
		//"0": "ZONE_NAME",
		//"1" : "SYMBOL",
		//"2" : "TVD_TOP",				TVD_TOP - кровля слоя [м]
		//"3" : "TVD_BOT",				TVD_BOT - подошва слоя [м]
		//"4" : "MD_TOP",
		//"5" : "MD_BOT",
		//"6" : "H_LAYER",				H_LAYER - толщина слоя [м]
		//"7" : "H_LOSS",				H_LOSS - проницаемая толщина слоя [м]
		//"8" : "GRAD_STRESS_MIN",
		//"9" : "STRESS_MIN"			STRESS_MIN - минимальное горизонтальное напряжение в слое [Па]
		//"10" : "YOUNG_MODULUS",		YOUNG_MODULUS - Модуль Юнга [Па]
		//"11" : "POISSON_RATIO",		POISSON_RATIO - Коэффициент Пуассона [д.ед.]
		//"12" : "TOUGHNESS",			TOUGHNESS - Коэффициент трещиностойкости [атм * м^0.5]
		//"13" : "FLUID_LOSS",			FLUID_LOSS - Коэффициент утечек по Картеру [м/сек^0.5]
		//"14" : "SPURT_LOSS",			SPRUT_LOSS - Мгновенные утечки [м^3/м^2]
		//"15" : "GRAD_PORE_PRESSURE",
		//"16" : "PORE_PRESSURE",
		//"17" : "POROSITY",
		//"18" : "PERMEABILITY",
		//"19" : "FORMATION_COMPRESSIBILITY",
		//"20" : "TOTAL_COMPRESSIBILITY",
		//"21" : "SPECIFIC_GRAVITY",
		//"22" : "OIL_SATURATION",
		//"23" : "WATER_SATURATION",
		//"24" : "GAS_SATURATION",
		//"25" : "TEMPERATURE",
		//"26" : "SPECIFIC_HEAT",
		//"27" : "HEAT_CONDUCTIVITY",
		//"28" : "LIMESTONE_SATURATION",
		//"29" : "DOLOMITE_SATURATION",


		//		std::cout << "Stage " << index << ":" << std::endl;
		DLL_Parametrs.IdStage = J_D_Stages[index]["id"].get<std::string>();
		json  J_D_S_Ports = J_D_Stages[index]["Ports"];
		std::vector<std::vector<json>> FracturePumpingSchedule_data = J_D_Stages[index]["FracturePumpingSchedule"]["data"];
		std::string modelingTime_s = FracturePumpingSchedule_data[0][1];
		//		sscanf_s(modelingTime_s.c_str(), "%lf", &DLL_Parametrs.modelingTime);
		DLL_Parametrs.modelingTime = std::stof(modelingTime_s);		//Время окончания закачки из JSON

		for (int index = 0; index < J_D_S_Ports.size(); index++)
		{
			//			std::cout << "Port " << index << ":" << std::endl;
			std::vector<std::vector<json>> Ports_DATA = J_D_S_Ports[index]["ReservoirFormation"]["data"];

			//Информация о послойной геологии модели.
			double middle = J_D_S_Ports[index]["md"].get<double>();		//координата центра центррального слоя
			DLL_Parametrs.Z_coordinate = middle;										//координата центра перфорации (для отправки в выходной json)

			double L_start;
			double L_end;
			double stress;	//Горизонтальное напряжение в слое
			double young;	//Модуль Юнга
			double poisson;	//Коэффициент пуассона
			double carter;	//Коэффициент утечек по Картеру
							//Временные переменный для создания трехслойной схемы
			double	L_first;
			bool flag1 = false;	//Флаг исполнения одного действия
			bool flag2 = false;	//Флаг исполнения одного действия

			for (int index = 0; index < Ports_DATA.size(); index++)
			{
				//Глубины начала и окончания index слоя
				std::string L_start_s = Ports_DATA[index][2];
				L_start = std::stof(L_start_s);
				if (index == 0)L_first = L_start;
				std::string L_end_s = Ports_DATA[index][3];
				L_end = std::stof(L_end_s);

				//Горизонтальное напряжение в слое
				std::string stress_s = Ports_DATA[index][9];
				stress = std::stof(stress_s);
				//Модуль Юнга
				std::string young_s = Ports_DATA[index][10];
				young = std::stof(young_s);
				//Коэффициент пуассона
				std::string poisson_s = Ports_DATA[index][11];
				poisson = std::stof(poisson_s);
				//Коэффициент утечек по Картеру
				std::string carter_s = Ports_DATA[index][13];
				carter = std::stof(carter_s);

				if ((L_end - middle) == 0)	//Проверяем случай, когда центр перфорацмм попадает на границу между слоями
				{
					index++;

					std::string L_end_s = Ports_DATA[index][3];
					L_end = std::stof(L_end_s);
					double	stress_temp = 0.;
					double	young_temp = 0.;
					double	poisson_temp = 0.;
					double	carter_temp = 0.;
					//Горизонтальное напряжение в слое
					std::string stress_s = Ports_DATA[index][9];
					stress_temp = std::stof(stress_s);
					//Модуль Юнга
					std::string young_s = Ports_DATA[index][10];
					young_temp = std::stof(young_s);
					//Коэффициент пуассона
					std::string poisson_s = Ports_DATA[index][11];
					poisson_temp = std::stof(poisson_s);
					//Коэффициент утечек по Картеру
					std::string carter_s = Ports_DATA[index][13];
					carter_temp = std::stof(carter_s);
					//Усредняем соседние слои в соответствии с значениями толщин слоев.
					stress = (stress*fabs(middle - L_start) + stress_temp * fabs(middle - L_end)) / fabs(L_end - L_start);
					young = (young*fabs(middle - L_start) + young_temp * fabs(middle - L_end)) / fabs(L_end - L_start);
					poisson = (poisson*fabs(middle - L_start) + poisson_temp * fabs(middle - L_end)) / fabs(L_end - L_start);
					carter = (carter*fabs(middle - L_start) + carter_temp * fabs(middle - L_end)) / fabs(L_end - L_start);

				}
				layers.push_back(std::vector<double>{L_start - middle, L_end - middle, stress*std::pow(10, -6), young*std::pow(10, -9), poisson, carter * 1000000.});
			}

		}
		//информация о режимах закачки
		//"0": "STAGE_TYPE",
		//"1" : "TIME",
		//"2" : "SLURRY_RATE",
		//"3" : "SLURRY_VOLUME",
		//"4" : "FLUID_VOLUME",
		//"5" : "FLUID_TYPE",
		//"6" : "PROPANT_TYPE",
		//"7" : "PROPANT_CONCENTRATION_START",
		//"8" : "PROPANT_CONCENTRATION_END",
		//"9" : "ACID_TYPE"
		//"10" : "ACID_CONCENTRATION",
		//"11" : "ACID_EQUILIBRIUM_CONCENTRATION",
		//"12" : "DIFFUSION_COEFFICIENT_MULTIPLIER",
		//"13" : "PROPANT_MASS",
		//"14" : "TOTAL_FLUID_VOLUME",
		//"15" : "TOTAL_PROPANT_MASS",
		//"16" : "TOTAL_SLURRY_VOLUME",
		//"17" : "TOTAL_PUMP_TIME",
		//"18" : "PROPANT_VOLUME",
		//"19" : "FLUID_MASS",
		//"20" : "FLUID_RATE",
		//"21" : "PROPANT_RATE",
		//	std::cout << "time " << index << ":" << std::endl;
		std::vector<std::vector<json>> Pumping_DATA = J_D_Stages[index]["UserPumpingSchedule"]["data"];		//["UserPumpingSchedule"]["data"];

		Time.push_back(0.);
		//	TimeStart.push_back(0.);
		int start_stage_time = 0;
		for (size_t index_PD = 0; index_PD < Pumping_DATA.size(); index_PD++)
		{
			//вычисление параметров стадий закачки
			//Время стадии
			std::string fluid_type = Pumping_DATA[index_PD][5];
			//Тип жидкости. (дальше по этому типу жидкости мы будем добавлять временные промежутки для конкретной жидкости)
			std::string end_time_s = Pumping_DATA[index_PD][1];
			int end_stage_time = std::stoi(end_time_s);				//продолжительность текущей стадии
																	//			std::cout << end_stage_time << "     ";
																	//Расход жидкости за стадию, куб. М
			std::string FluodVol_s = Pumping_DATA[index_PD][2];
			double FluodVol = std::stof(FluodVol_s);
			FluodVol = FluodVol * 60; //;/ end_time;					//перерасчет в м3/мин
									  //Масса проппанта за стадию, кг
			std::string PropMass_s = Pumping_DATA[index_PD][7]; ///Pumping_DATA[index][13];	//концентрация проппанта на старте промежутка стадии
			double PropMass_start = std::stof(PropMass_s);
			std::string PropMass_e = Pumping_DATA[index_PD][8]; ///Pumping_DATA[index][13]; //концентрация пропанта в конце промежутка стадии
			double PropMass_end = std::stoi(PropMass_e);


			//Цикл поиска параметров жидкости в базе жидкостей
			std::vector<std::vector<json>> F_Base = J_IN_DATA["FluidBase"]["data"];
			// std::cout<<"F_Base.size() "<<F_Base.size()<<std::endl;
			for (int index = 0; index < F_Base.size(); index++)	//Цикл перебора жидкостей
			{
				//Модуль для нестационарной реалогии жидкости
				std::string F_Type = F_Base[index][0];
				if (fluid_type == F_Type)
				{
					//Для конкретной жидкости получили двумерный вектор ее параметров.
					std::vector<std::vector<json>> F_Base_data = F_Base[index][6]["data"];
					//std::cout<<"F_Base_data size = "<< F_Base_data.size()<<std::endl;
					int start_time = 0;
					for (int index1 = 0; index1 < F_Base_data.size(); index1++)
					{
						//double steptime;// = 0.;

						std::vector<std::vector<json> > F_operator_time = F_Base_data[0][0]["data"];
						int stages = F_Base_data[0][0]["rowCount"].get<int>();
						size_t j = 0;
						size_t k = 0;
						for (int index2 = 0; index2 < stages; index2++)
						{
							///////////////////////////////////////////////////////////////////////////
							//Времена изменения характеристик для конкретной жидкости
							//Добавляем в вектор общих времен, который потом будем сортировать
							std::string end_time_s = F_operator_time[index2][0];
							int end_time = std::stoi(end_time_s);
							//Зависимость вязкости от времени для конкретной жидкости
							std::string viscosity_s = F_operator_time[index2][1];
							double viscosity = std::stof(viscosity_s);
							//Зависимость реологии от времени для конкретной жидкости
							std::string reology_s = F_operator_time[index2][2];
							double reology = std::stof(reology_s);
							//Заполнение вектора закачек.
							//start_stage_time						
							//Это интересующий нас промежуток
							//start_stage_time+ end_stage_time
							////////////////////////////////////////////////////////////////////////////////////////////////////////////////////						
							//							|																		|
							//					start_stage_time												start_stage_time + end_stage_time													
							//							|																		|
							//							|																		|
							//			1				|								2										|				3
							//							|																		|
							//							|																		|
							////////////////////////////////////////////////////////////////////////////////////////////////////////////////////						

							if (start_time <= start_stage_time + end_stage_time && start_time + end_time >= start_stage_time && start_time != end_time)
								if (start_time < start_stage_time && start_time + end_time > start_stage_time && start_time + end_time < start_stage_time + end_stage_time) // мы находимся началом в 1 интервале, концом во 2 интервале
								{
									injection.push_back(std::vector<double> {start_stage_time / 60., (start_time + end_time) / 60.,
										FluodVol,
										PropMass_end,			//нужно уточнить что сюда класть
										reology,
										viscosity});
								}
								else if (start_time >= start_stage_time && start_time + end_time <= start_stage_time + end_stage_time)// мы находимся в 2 интервале и начало и конец
								{
									injection.push_back(std::vector<double> {start_time / 60., (start_time + end_time) / 60.,
										FluodVol,
										PropMass_end,			//нужно уточнить что сюда класть
										reology,
										viscosity});
								}
								else if (start_time >= start_stage_time && end_time > start_stage_time + end_stage_time) // мы находимся началом в 2 интервале, концом во 3 интервале
								{
									injection.push_back(std::vector<double> {start_time / 60., (start_stage_time + end_stage_time) / 60.,
										FluodVol,
										PropMass_end,			//нужно уточнить что сюда класть
										reology,
										viscosity});
								}
								else // мы находимся началом в 1 интервале, концом во 3 интервале
								{
									injection.push_back(std::vector<double> {start_stage_time / 60., (start_stage_time + end_stage_time) / 60.,
										FluodVol,
										PropMass_end,			//нужно уточнить что сюда класть
										reology,
										viscosity});
								}
							//перемещаемся на следующую стадию реологии жидкости
							start_time = start_time + end_time;
						}
					}
					break; //выходим из цикла перебора жидкостей так как мы уже нужную жидкость нашли
				}
			}


			start_stage_time = start_stage_time + end_stage_time;		//перемещаемся на следующую стадию
		}

	}

	//Цикл заполнения структуры параметров еврейкости в базе еврейкостей
	std::vector<std::vector<json>> F_Base = J_IN_DATA["FluidBase"]["data"];

	Fluid_Param fluidparam;

	//Считываем параментры всех еврейкостей и заполняем вектор класса для хранения данных
	for (int index = 0; index < F_Base.size(); index++)	//Цикл перебора еврейкостей
	{
		std::string KEY = F_Base[index][0];
		std::string CODE = F_Base[index][1];
		std::string SPECIFIC_GRAVITY_s = F_Base[index][3];
		double SPECIFIC_GRAVITY = std::stof(SPECIFIC_GRAVITY_s);
		std::string SHARE_RATE_s = F_Base[index][4];
		double SHARE_RATE = std::stof(SHARE_RATE_s);
		std::string SPECIFIC_HEAT_s = F_Base[index][5];
		double SPECIFIC_HEAT = std::stof(SPECIFIC_HEAT_s);
		//Получаем значение температуры для каждого слоя
		std::vector<std::vector<json>> F_Operator = F_Base[index][6]["data"];

		std::vector <std::vector<double>> rheology_temp; //Вектор хранения реологии одной еврейкости

		for (int index1 = 0; index1 < F_Operator.size(); index1++)	//Цикл перебора температурных параметров еврейкостей
		{
			std::string  temperature_s = F_Operator[index1][0]["value"];
			double temperature = std::stof(temperature_s);
			std::vector<std::vector<json>> T_operator = F_Operator[index1][0]["data"];
			std::string D_Viscosity_s = T_operator[0][1];
			std::string Rheology_s = T_operator[0][2];
			std::string Eff_Viscosit_s = T_operator[0][3];
			double D_Viscosity = std::stof(D_Viscosity_s);
			double Rheology = std::stof(Rheology_s);
			double Eff_Viscosit = std::stof(Eff_Viscosit_s);
			rheology_temp.push_back({ temperature, D_Viscosity, Rheology, Eff_Viscosit });
		}
		fluidparam.setReology(KEY, CODE, SPECIFIC_GRAVITY, SHARE_RATE, SPECIFIC_HEAT, rheology_temp);
		DLL_Parametrs.FluidParam.push_back(fluidparam);
	}


	//Цикл заполнения структуры параметров проппанта в базе проппантов
	std::vector<std::vector<json>> P_Base = J_IN_DATA["ProppantBase"]["data"];
	Proppant_Param proppantparam;

	//Считываем параментры всех проппантов и заполняем вектор класса для хранения данных
	for (int index = 0; index < P_Base.size(); index++)	//Цикл перебора проппантов
	{
		std::string KEY = P_Base[index][0];
		std::string CODE = P_Base[index][1];
		std::string  DIAMETER_s = P_Base[index][4];
		double DIAMETER = std::stof(DIAMETER_s);
		std::string SPECIFIC_GRAVITY_s = P_Base[index][5];
		double SPECIFIC_GRAVITY = std::stof(SPECIFIC_GRAVITY_s);
		std::string BULK_GRAVITY_s = P_Base[index][6];
		double BULK_GRAVITY = std::stof(BULK_GRAVITY_s);
		std::string CRITICAL_DENSITY_s = P_Base[index][7];
		double CRITICAL_DENSITY = std::stof(CRITICAL_DENSITY_s);
		std::string SPECIFIC_HEAT_s = P_Base[index][8];
		double SPECIFIC_HEAT = std::stof(SPECIFIC_HEAT_s);
		std::string EMBEDMENT_RATIO_s = P_Base[index][9];
		double EMBEDMENT_RATIO = std::stof(EMBEDMENT_RATIO_s);
		std::string SETTLING_VELOSITY_FACTOR_s = P_Base[index][10];
		double SETTLING_VELOSITY_FACTOR = std::stof(SETTLING_VELOSITY_FACTOR_s);
		std::string DRIFT_VELOSITY_FACTOR_s = P_Base[index][11];
		double DRIFT_VELOSITY_FACTOR = std::stof(DRIFT_VELOSITY_FACTOR_s);

		//Создаем экземпляр класса для хранения параметров проппанта и пеердачи их в вектор
		proppantparam.setCoordinates(KEY, CODE, DIAMETER, SPECIFIC_GRAVITY, BULK_GRAVITY, CRITICAL_DENSITY, SPECIFIC_HEAT, EMBEDMENT_RATIO, SETTLING_VELOSITY_FACTOR, DRIFT_VELOSITY_FACTOR);
		DLL_Parametrs.ProppantParam.push_back(proppantparam);
	}



	for (int index = layers.size() - 1; index >= 0; index--)
		layers_new.push_back(layers[index]);
	//////////////////////////////////////////////////////////////////////////////////////////
	//End function ImportJSON
	//////////////////////////////////////////////////////////////////////////////////////////
	return;
}






//////////////////////////////////////////////////////////////////////////////////////////
//Function ApproximateLayers
//////////////////////////////////////////////////////////////////////////////////////////
/*!
\brief  Преобразование выходных данных о реологии

\details Функция производит усреднение и размытие слоев.

\param[in,out] layers - матрица послойной реологий
*/
void ApproximateLayers(std::vector< std::vector<double> > &layers) {
	const double epsilon = pow(10., -8);

	std::size_t centralIndex;

	double H;
	double dx;

	/// Меняем последовательность записи слоев (сортировка от большего к меньшему)
	std::vector< std::vector<double> > la;
	for (int i = layers.size() - 1; i >= 0; --i) {
		la.push_back(layers[i]);
	}
	layers = la;
	la.clear();

	/// Определяем центральный слой
	for (std::size_t i = 0; i < layers.size(); ++i) {
		if (layers[i][0] * layers[i][1] <= epsilon) {
			H = layers[i][1] - layers[i][0];
			centralIndex = i;

			break;
		}
	}

	double sigma = layers[centralIndex][2];

	std::vector<std::size_t> sloi;

	double h;
	double hnew;

	std::size_t it = centralIndex;

	sloi.push_back(it);

	h = H;
	sigma *= H;

	/// Убеждаемся, что центральный слой будет хотя бы 11.2 метра
	if (11.2 >= H) {
		--it;

		while (5.6 > std::abs(layers[it][1])) {
			hnew = layers[it][1] - layers[it][0];

			sigma += hnew * layers[it][2];
			h += hnew;

			sloi.push_back(it);

			--it;
		}

		it = centralIndex + 1;

		while (5.6 > layers[it][0]) {
			hnew = layers[it][1] - layers[it][0];

			sigma += hnew * layers[it][2];
			h += hnew;

			sloi.push_back(it);

			++it;
		}

		H = h;
	}

	sigma /= H;

	dx = 0.09 * (H - epsilon);

	/// Переписываем слои в новую матрицу    
	std::size_t min = ai::min(sloi);
	std::size_t max = ai::max(sloi);

	centralIndex = min;

	std::vector<std::vector<double> > layers1;

	for (std::size_t i = 0; i < layers.size(); ++i) {
		if (i == min) {
			layers1.push_back(
				std::vector<double>{
				layers[min][0],
					layers[max][1],
					sigma,
					layers[i][3],
					layers[i][4],
					layers[i][5],
					layers[i][6]
			}
			);

			i = max;
		}
		else {
			layers1.push_back(layers[i]);
		}
	}

	/// Симметризуем центральный слой
	const double deltaH = layers[max][1] - 0.5 * H;

	for (std::size_t i = 0; i < layers1.size(); ++i) {
		layers1[i][0] -= deltaH;
		layers1[i][1] -= deltaH;
	}

	std::vector< std::vector<double> > mesh;

	mesh.push_back(layers1[centralIndex]);

	double heightStep = dx;
	double currentHeight = layers1[centralIndex][0];

	for (int i = (int)centralIndex - 1; i >= 0; --i) {
		mesh.push_back(std::vector<double>{
			currentHeight - heightStep, currentHeight, 0., 0., 0., 0., 0.
		});
		currentHeight -= heightStep;
	}

	currentHeight = layers1[centralIndex][1];

	for (std::size_t i = centralIndex + 1; i < layers1.size(); ++i) {
		mesh.push_back(std::vector<double>{
			currentHeight, currentHeight + heightStep, 0., 0., 0., 0., 0.
		});
		currentHeight += heightStep;
	}

	std::sort(
		mesh.begin(),
		mesh.end(),
		[](std::vector<double> vec1, std::vector<double> vec2)->bool {
		return vec1[0] < vec2[0];
	}
	);
	std::cout << "layers1: " << std::endl;
	ai::print(layers1);

	for (std::size_t j = 0; j < mesh.size(); ++j) {
		std::size_t indexStart = 0;
		std::size_t indexFinish = 0;

		double heightStart = mesh[j][0];
		double heightFinish = mesh[j][1];
		std::cout << " " << std::endl;
		for (std::size_t i = 0; i < layers1.size(); ++i) {
			if (layers1[i][0] <= heightStart && heightStart <= layers1[i][1]) {
				indexStart = i;
			}
			if (layers1[i][0] <= heightFinish && heightFinish <= layers1[i][1]) {
				indexFinish = i;
			}
		}
		if (indexStart == indexFinish) {
			double sigmaValue = layers1[indexStart][2]
				* (heightFinish - heightStart) +
				(dx - heightFinish + heightStart)*layers1[indexStart + 1][2];

			double young = layers1[indexStart][3] * (heightFinish - heightStart)
				+ (dx - heightFinish + heightStart)*layers1[indexStart + 1][3];

			double poisson = layers1[indexStart][4]
				* (heightFinish - heightStart) +
				(dx - heightFinish + heightStart)*layers1[indexStart + 1][4];


			double carter = layers1[indexStart][5]
				* (heightFinish - heightStart) +
				(dx - heightFinish + heightStart)*layers1[indexStart + 1][5];

			double kin = layers1[indexStart][6]
				* (heightFinish - heightStart) +
				(dx - heightFinish + heightStart)*layers1[indexStart + 1][6];


			mesh[j][2] = sigmaValue / dx;

			mesh[j][3] = young / dx;

			mesh[j][4] = poisson / dx;

			mesh[j][5] = carter / dx;

			mesh[j][6] = kin / dx;


		}
		else {
			std::cout << "indexStart = " << indexStart << "  indexFinish = " << indexFinish << std::endl;
			double sigmaValue = (layers1[indexStart][1] - heightStart)
				* layers1[indexStart][2];

			double young = (layers1[indexStart][1] - heightStart)
				* layers1[indexStart][3];

			double poisson = (layers1[indexStart][1] - heightStart)
				* layers1[indexStart][4];

			double carter = (layers1[indexStart][1] - heightStart)
				* layers1[indexStart][5];

			double kin = (layers1[indexStart][1] - heightStart)
				* layers1[indexStart][6];


			std::cout << heightStart << "/" << heightFinish << " " << sigmaValue << " ixj " << indexStart << " " << indexFinish << std::endl;

			for (std::size_t i = indexStart + 1; i < indexFinish; ++i) {
				sigmaValue += std::abs(layers1[i][1] - layers1[i][0]) * layers1[i][2];
				young += std::abs(layers1[i][1] - layers1[i][0]) * layers1[i][3];
				poisson += std::abs(layers1[i][1] - layers1[i][0]) * layers1[i][4];
				carter += std::abs(layers1[i][1] - layers1[i][0]) * layers1[i][5];
				kin += std::abs(layers1[i][1] - layers1[i][0]) * layers1[i][6];
			}

			std::cout << heightStart << "/" << heightFinish << " " << sigmaValue << " ixj " << indexStart << " " << indexFinish << std::endl;

			sigmaValue += (heightFinish - layers1[indexFinish][0])
				* layers1[indexFinish][2];

			young += (heightFinish - layers1[indexFinish][0])
				* layers1[indexFinish][3];

			poisson += (heightFinish - layers1[indexFinish][0])
				* layers1[indexFinish][4];

			carter += (heightFinish - layers1[indexFinish][0])
				* layers1[indexFinish][5];

			kin += (heightFinish - layers1[indexFinish][0])
				* layers1[indexFinish][6];
			std::cout << heightStart << "/" << heightFinish << " " << sigmaValue << " ixj " << indexStart << " " << indexFinish << std::endl;

			sigmaValue /= heightFinish - heightStart;

			young /= heightFinish - heightStart;

			poisson /= heightFinish - heightStart;

			carter /= heightFinish - heightStart;

			kin /= heightFinish - heightStart;

			mesh[j][2] = sigmaValue;

			mesh[j][3] = young;

			mesh[j][4] = poisson;

			mesh[j][5] = carter;

			mesh[j][6] = kin;



			std::cout << heightStart << "/" << heightFinish << " " << sigmaValue << " ixj " << indexStart << " " << indexFinish << std::endl;
		}
	}
	std::cout << "mesh: " << std::endl;

	std::sort(
		mesh.begin(),
		mesh.end(),
		[](std::vector<double> vec1, std::vector<double> vec2)->bool {
		return vec1[0] > vec2[0];
	}
	);
	ai::print(mesh);
	ai::save("./m", mesh);

	layers = mesh;

	return;
	//Обсчтываем все сверху продуктивного

	//////////////////////////////////////////////////////////////////////////////////////////
	//End function ApproximateLayers
	//////////////////////////////////////////////////////////////////////////////////////////
}
