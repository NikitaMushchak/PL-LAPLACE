#include <vector>
#include <iostream>
#include <algorithm>

#include "nlohmann/json.hpp"
#include "ailibrary/ai.hh"


#include "io.hh"

/*!
\brief  ������� �������� ������ � ��������� �� ���������� ������

\details ������� ��������� �������� ������ � ������  � ������� json.

\param[in] Wk - ������ ���������
\param[in] wn - ����������� �����������
\param[in] pressure - ������� ��������
\param[in] concentration - ������� ������������ ���������
\param[in] x
\param[in] y
\param[in] dx
\param[in] dy
\param[in] i00
\param[in] j00
\param[in] index - ������� �������� ���������
\param[in] Time - ������� ����� (���������)
\param[in] timeScale -
\param[in] fluidEfficiency - ������������� ��������
\param[in] fluidDensity - ��������� ��������
\param[in] proppantEfficiency - ������������� ���������
\param[in] proppantDensity - ��������� ���������
\param[in] fluidInjection
\param[in] proppantInjection
\param[in] nominalStress - ����������� �������� ���������� � ��������� ����������, [���]
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
	double Z_coordinate,
	std::string IdDesign,
	std::string IdStage
) {
	// ����� �������� � 1 ���
	const double atmosphereCoefficient = 9.869;

	const double xSize = x.size();
	const double ySize = y.size();

	//////////////////////////////////////////////////////////////////////////////////////////
	//��� ������ ���� ���� ������ ��������� ������
	//////////////////////////////////////////////////////////////////////////////////////////
	double azimuth = 0;								//������ ��������� ������ � ������������(�������� 0)
	int id = 0;										//(�������� 0)
	int stage_id = 0;								//(�������� 0)
	int num_fluids = 1;								//number of fluids - ����� ���������� ���������
	int num_proppants = 0;							//number of proppants - ����� ���������� ���������


	if (fluidEfficiency > 100)fluidEfficiency = 100.;


	//////////////////////////////////////////////////////////////////////////////////////////
	//��������� ��������� ����� json
	//////////////////////////////////////////////////////////////////////////////////////////
	//Results
	//	Smin - ����������� �������� ���������� � ��������� ����������, [���]
	//	accumulated data
	//	fluid efficiency - ������������� �������� � ������ ������� "time", [%]
	//	net pressure - �������� ������� �������� �� ����� ������� � ���� ���������� ���������� � ����������� "Smin", [���]
	//	proppant concentration - ����������� �������� � ������ ������� "time" �� ����� � �������, [�� / ���.�.]
	//	rate - ������ ����� � ������ ������� "time" �� ����� � �������, [���.� / ���]
	//	time - ������ ������� �� ������� ������� ��������� �������
	//	geometry
	//	branches - ����������� ��������� �������(0)
	//	cells - ��������� ������ ���������� � ������������������ : ������ - ����(j �������), ����� - �������(i �������)
	//	azimuth - ������ ��������� ������ � ������������(�������� 0)
	//	concentrations - �������� ���� ������������ ���, [�.��.]
	//	dx, dy, dz - ������� ��������� ������(dy - ���������), [�]
	//	id - 0; stage id - 0
	//	x, y, z - ���������� ������ ��������� ����� � ������������, [�]
	//	grid - ��������� ��������� �������
	//	nx - ����� ���������� ��������� ����� �� �
	//	nz - ����� ���������� ��������� ����� �� Z
	//	slurry - �������� ����� � �������
	//	mass density of components - ��������� ���������(� ��� �� ������������������, ��� � � "concentrations"), [�� / ���.�]
	//	name of components - �������� ���������(� ��� �� ������������������, ��� � � "concentrations")
	//	number of fluids - ����� ���������� �������
	//	number of proppants - ����� ���������� ���������
	//	coordinates - ���������� ������ ��������� ����������, [�]





	std::ostringstream output;

	output << "{\n";
	//Report
	output << "\"Report\":\n ";
	output << "{\n";
	output << "\t\"relativeMassBalanceError\": " << "0.0" << ",\n";	//����� ������ ������ ��� ����� �����!!!
	output << "\t\"stageId\": \"" << IdStage.c_str() << "\"\n";
	output << "},\n"; //close Report
	//Ports
	output << "\t\"" << IdStage.c_str() << "\" :\n";	output << "{\n";
	output << "\"Ports\":\n [ \n";
	//////////////////////////////////////////////////////////////////////////////////////////
	//Results
	//////////////////////////////////////////////////////////////////////////////////////////
	{
	output << "{\n";
	output << std::setprecision(14) << "\"Results\": {\n";
	//Smin - ����������� �������� ���������� � ��������� ����������, [���]
	output << "\t\"Smin\": " << nominalStress * atmosphereCoefficient << ",\n";
	output << "\t\"accumulated data\": {\n";
	//	fluid efficiency - ������������� �������� � ������ ������� "time", [%]
	output << "\t  \"fluid efficiency\": " << fluidEfficiency << ",\n";
	//	net pressure - �������� ������� �������� �� ����� ������� � ���� ���������� ���������� � ����������� "Smin", [���]
	output << "\t  \"net pressure\": " << pressure[index[i00][j00]] * atmosphereCoefficient << ",\n";
	//	proppant concentration - ����������� �������� � ������ ������� "time" �� ����� � �������, [�� / ���.�.]
	output << "\t  \"proppant concentration\": " << proppantInjection << ",\n";
	//	rate - ������ ����� � ������ ������� "time" �� ����� � �������, [���.� / ���]
	output << "\t  \"rate\": " << fluidInjection / (timeScale / (wn * 60.)) << ",\n";
	//	time - ������ ������� �� ������� ������� ��������� �������
	output << "\t  \"time\": " << Time << "\n";
	output << "\t},\n"; //close accumulated data

	//////////////////////////////////////////////////////////////////////////////////////////
	//     geometry
	//     branches
	//////////////////////////////////////////////////////////////////////////////////////////
	output << "\t\"geometry\": {\n";

	//     branches - ����������� ��������� �������(0)
	output << "\t\t\"branches\": [\n\t\t {\n";
	//     cells - ��������� ������ ���������� � ������������������ : ������ - ����(j �������), ����� - �������(i �������)
	output << "\t\t \"cells\": [\n";
	int activeElements = 0;                 //����� �������� ��������� (�� ��������, ��� ������ ��������� �� �������)

											//////////////////////////////////////////////////////////////////////////////////////////
											//     ������ ����� �������� ���������
											//////////////////////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < xSize; i++)
	{
		for (int j = 0; j < ySize; j++)
		{
			//if (Wk[j*xSize + i] != 0)
				activeElements++; // ���� ������ �� ����������, �� �� �� �������!!!
		}
	}

	int Index = 0;                                                 //��������� ��������� - ����� ���������� ������� ������

for (int i = 0; i < xSize; i++)
	{
	for (int j = 0; j < ySize; j++)
		{
			//*if (Wk[j*xSize + i] == 0) continue; // ���� ������ �� ����������, �� �� �� �������!!!
			Index++;
			output << "\t\t{ \n";
			//     azimuth - ������ ��������� ������ � ������������(�������� 0)
			output << "\t\t  \"azimuth\": " << azimuth << ",\n";
			//     concentrations - �������� ���� ������������ ���, [�.��.]
			output << "\t\t  \"concentrations\": [\n ";
			//for (int index=0; index< num_proppants; index++) // ����� ���������� ���������� (����� ����� ��������� ����������)
			output << "\t\t  " << concentration[index[i][j]] << ",\n";                //����� ������������ ���������
			output << "\t\t  " << 1 - concentration[index[i][j]] << "\n";             //����� ������������ ��������

			output << "\t\t ],\n"; //close concentrations
								   //  dx, dy, dz - ������� ��������� ������(dy - ���������), [�]
			output << "\t\t  \"dx\": " << dx << ",\n";
			//if (std::isnan(Wk[j*xSize + i]))
			//	Wk[j*xSize + i] = 100;
			output << "\t\t  \"dy\": " << Wk[i*ySize + j]*wn << ",\n"; // Wk[j*xSize + i]
			output << "\t\t  \"dz\": " << dy << ",\n";
			//     id - 0; stage id - 0
			output << "\t\t  \"i\": " << i << ",\n";
			output << "\t\t  \"id\": \"" << IdStage.c_str() << "\",\n";
			output << "\t\t  \"j\": " << j << ",\n";
			output << "\t\t  \"stage id\": " << stage_id << ",\n";
			//     x, y, z - ���������� ������ ��������� ����� � ������������, [�]
			output << "\t\t  \"x\": " << i * dx << ",\n";
			output << "\t\t  \"y\": " << j * Wk[j*xSize + i] << ",\n";
			output << "\t\t  \"z\": " << j * dy << "\n";

			if (activeElements != Index) //i != xSize - 1 && j != ySize - 1  &&
				output << "\t\t },\n";
			else {	//���� ������� ��������� ������� ������� �������� ��������� - ��������� JSON ������
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
		//     grid - ��������� ��������� �������
		//////////////////////////////////////////////////////////////////////////////////////////
		output << "\t\"grid\": {\n";
		//     nx - ����� ���������� ��������� ����� �� �
		output << "\t\t  \"nx\": " << xSize << ",\n";
		//     nz - ����� ���������� ��������� ����� �� Z
		output << "\t\t  \"nz\": " << ySize << "\n";
		output << "\t},\n"; //close grid
		//////////////////////////////////////////////////////////////////////////////////////////
		//close grid
		//////////////////////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////////////////////
		//     slurry - �������� ����� � �������
		//////////////////////////////////////////////////////////////////////////////////////////
		output << "\t\"slurry\": {\n";

		output << "\t\t\"mass density of components\": [\n";
		//     mass density of components - ��������� ���������(� ��� �� ������������������, ��� � � "concentrations"), [�� / ���.�]
		output << "\t\t " << proppantDensity << ",\n";
		output << "\t\t " << fluidDensity << "\n";

		output << "\t\t ],\n"; //close mass density of components

		output << "\t\t\"name of components\": [\n";
		//     name of components - �������� ���������(� ��� �� ������������������, ��� � � "concentrations")
		output << "\t\t \"Propant 1\",\n";
		output << "\t\t \"Fluid 1\"\n";
		output << "\t\t ],\n"; //close name of componentss

							   //     number of fluids - ����� ���������� �������
		output << "\t\t  \"number of fluids\": " << num_fluids << ",\n";
		//     number of proppants - ����� ���������� ���������
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
	//coordinates - ���������� ������ ��������� ����������, [�]
	//////////////////////////////////////////////////////////////////////////////////////////
	output << "\"coordinates\": {\n";
	output << "\t\t  \"x\": " << 0 << ",\n";
	output << "\t\t  \"y\": " << 0 << ",\n";
	output << "\t\t  \"z\": " << Z_coordinate << "\n";	//����� �������� ���� ���������� ������ J_D_S_Ports[index]["md"].get<double>()
	output << "},\n";
	output << "\t\t\"id\": \"" << IdDesign << "\"\n";

	output << "}\n ";
	output << "]\n ";
	//������� Ports
	output << "}\n }\n";
	//������� ���

	return	output.str();
	//////////////////////////////////////////////////////////////////////////////////////////
	//End function ExportJson
	//////////////////////////////////////////////////////////////////////////////////////////
}




void ImportJSON(
	std::string J_String,
	std::string &IdDesign,
	std::string &IdStage,
	std::vector< std::vector<double> > &layers_new,
	std::vector< std::vector<double> > &injection,
	double &modelingTime,
	double &Emit_time,
	double &Z_coordinate
)
{
	using json = nlohmann::json;

	//std::string J_String;
	//	ai::parseFileIntoString("../Planar Challenge - Test 12_init.json", J_String);
	json J_IN_DATA = json::parse(J_String);

	std::vector<double> Time;
	std::vector<double> TimeStart;
	std::vector<double> TimeEnd;


	//	std::vector< std::vector <double> > layers;
	//������ layers:
	//{L_start - middle, L_end - middle, stress, young, poisson, carter * 1000000.}
	//[0] - ������ ���� (�������� ������������ ������������ ����)
	//[1] - ����� ���� (�������� ������������ ������������ ����)
	//[2] - �������������� ���������� � ����
	//[3] - ������ ����
	//[4] - ����������� ��������
	//[5] - ����������� ������ �� �������


	//std::vector< LAYER > layers;
	std::vector< std::vector<double>> layers;
	//int index_middle;					//������ ������������ ���� � ������� layers
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

	Emit_time = J_IN_DATA["Design"]["Settings"]["emit"].get<double>();
	IdDesign = J_IN_DATA["Design"]["id"].get<std::string>();
	json  J_D_Stages = J_IN_DATA["Design"]["Stages"];

	for (int index = 0; index < J_D_Stages.size(); index++)
	{
		//���������� � ����� ReservoirFormation - ������� � ��������� ���������
		//"0": "ZONE_NAME",
		//"1" : "SYMBOL",
		//"2" : "TVD_TOP",				TVD_TOP - ������ ���� [�]
		//"3" : "TVD_BOT",				TVD_BOT - ������� ���� [�]
		//"4" : "MD_TOP",
		//"5" : "MD_BOT",
		//"6" : "H_LAYER",				H_LAYER - ������� ���� [�]
		//"7" : "H_LOSS",				H_LOSS - ����������� ������� ���� [�]
		//"8" : "GRAD_STRESS_MIN",
		//"9" : "STRESS_MIN"			STRESS_MIN - ����������� �������������� ���������� � ���� [��]
		//"10" : "YOUNG_MODULUS",		YOUNG_MODULUS - ������ ���� [��]
		//"11" : "POISSON_RATIO",		POISSON_RATIO - ����������� �������� [�.��.]
		//"12" : "TOUGHNESS",			TOUGHNESS - ����������� ���������������� [��� * �^0.5]
		//"13" : "FLUID_LOSS",			FLUID_LOSS - ����������� ������ �� ������� [�/���^0.5]
		//"14" : "SPURT_LOSS",			SPRUT_LOSS - ���������� ������ [�^3/�^2]
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
		IdStage = J_D_Stages[index]["id"].get<std::string>();
		json  J_D_S_Ports = J_D_Stages[index]["Ports"];
		for (int index = 0; index < J_D_S_Ports.size(); index++)
		{
			//			std::cout << "Port " << index << ":" << std::endl;
			std::vector<std::vector<json>> Ports_DATA = J_D_S_Ports[index]["ReservoirFormation"]["data"];

			//���������� � ��������� �������� ������.
			double middle = J_D_S_Ports[index]["md"].get<double>();		//���������� ������ ������������� ����
			Z_coordinate = middle;										//���������� ������ ���������� (��� �������� � �������� json)

			double L_start;
			double L_end;
			double stress;	//�������������� ���������� � ����
			double young;	//������ ����
			double poisson;	//����������� ��������
			double carter;	//����������� ������ �� �������
							//��������� ���������� ��� �������� ����������� �����
			double	L_first;
			bool flag1 = false;	//���� ���������� ������ ��������
			bool flag2 = false;	//���� ���������� ������ ��������

			for (int index = 0; index < Ports_DATA.size(); index++)
			{
				//������� ������ � ��������� index ����
				std::string L_start_s = Ports_DATA[index][2];
				sscanf(L_start_s.c_str(), "%lf", &L_start);
				if (index == 0)L_first = L_start;
				std::string L_end_s = Ports_DATA[index][3];
				sscanf(L_end_s.c_str(), "%lf", &L_end);

				//�������������� ���������� � ����
				std::string stress_s = Ports_DATA[index][9];
				sscanf(stress_s.c_str(), "%lf", &stress);
				//������ ����
				std::string young_s = Ports_DATA[index][10];
				sscanf(young_s.c_str(), "%lf", &young);
				//����������� ��������
				std::string poisson_s = Ports_DATA[index][11];
				sscanf(poisson_s.c_str(), "%lf", &poisson);
				//����������� ������ �� �������
				std::string carter_s = Ports_DATA[index][13];
				sscanf(carter_s.c_str(), "%lf", &carter);

				if ((L_end - middle) == 0)	//��������� ������, ����� ����� ���������� �������� �� ������� ����� ������
				{
					index++;

					std::string L_end_s = Ports_DATA[index][3];
					sscanf(L_end_s.c_str(), "%lf", &L_end);
					double	stress_temp = 0.;
					double	young_temp = 0.;
					double	poisson_temp = 0.;
					double	carter_temp = 0.;
					//�������������� ���������� � ����
					std::string stress_s = Ports_DATA[index][9];
					sscanf(stress_s.c_str(), "%lf", &stress_temp);
					//������ ����
					std::string young_s = Ports_DATA[index][10];
					sscanf(young_s.c_str(), "%lf", &young_temp);
					//����������� ��������
					std::string poisson_s = Ports_DATA[index][11];
					sscanf(poisson_s.c_str(), "%lf", &poisson_temp);
					//����������� ������ �� �������
					std::string carter_s = Ports_DATA[index][13];
					sscanf(carter_s.c_str(), "%lf", &carter_temp);
					//��������� �������� ���� � ������������ � ���������� ������ �����.
					stress	= (stress*fabs(middle - L_start)+ stress_temp*fabs(middle - L_end)) / fabs(L_end - L_start);
					young	= (young*fabs(middle - L_start) + young_temp * fabs(middle - L_end)) / fabs(L_end - L_start);
					poisson = (poisson*fabs(middle - L_start) + poisson_temp * fabs(middle - L_end)) / fabs(L_end - L_start);
					carter	= (carter*fabs(middle - L_start) + carter_temp * fabs(middle - L_end)) / fabs(L_end - L_start);

				}
				layers.push_back(std::vector<double>{L_start - middle, L_end - middle, stress*std::pow(10, -6), young*std::pow(10, -9), poisson, carter * 1000000.});
			}

		}
		//���������� � ������� �������
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
		TimeStart.push_back(0.);

		for (size_t y = 0; y < Pumping_DATA.size(); ++y)
		{
			flui.push_back(Pumping_DATA[y][5]);
		}
		// std::cout << "flui size = " << flui.size() << std::endl;
		// ai::printVector(flui);
		/*double tmp_T1;
		for (int i = flui.size() - 1; i > 0; i--)
		for (int j = 0; j<(i - 1); j++)
		{

		if (flui[j] == flui[j + 1])
		flui.erase(flui.begin() + j + 1);
		}*/
		//� ������� ���������
		/*unique(flui.begin(), flui.begin() + flui.size() - 1);
		flui.resize(flui.begin() + flui.size() - 1 - flui.begin());*/

		for (int i = 0; i < flui.size();) {
			if (std::find(flui.begin(), flui.begin() + i, flui[i]) != flui.begin() + i)
				flui.erase(flui.begin() + i);
			else
				i++;
		}
		// std::cout << "flui size = " << flui.size() << std::endl;



		double start_time = 0;
		//std::cout << "Pumping_DATA.size()"<<Pumping_DATA.size()<<std::endl;
		for (int index = 0; index < Pumping_DATA.size(); index++)
		{
			///////////////////////////////////////////////////////////////////
			//	��� ������� ���� �������� ���������� ��������� ������ ������ ��������� ����� ������� (["UserPumpingSchedule"]["data"])
			//  ��� ������ ��� ����������� ��� ����, ����� ������� ������������ ���� ������� ��������� �������� ������� � �������� ��������� � ����������� �� �������.
			///////////////////////////////////////////////////////////////////
			TimeStart.clear();
			TimeEnd.clear();
			timefl.clear();
			Time.clear();
			viscosities.clear();
			reologies.clear();
			//flui.clear();
			double steptime = 0.;
			//������ �� ��������
			for (int index = 0; index < Pumping_DATA.size(); index++)
			{

				std::string end_time_s = Pumping_DATA[index][1];
				double end_time;
				sscanf(end_time_s.c_str(), "%lf", &end_time);
				TimeStart.push_back(steptime);
				steptime = (steptime + end_time);
				Time.push_back(steptime);
				TimeEnd.push_back(steptime);

				//������ �������� ���-�� �������� � ��������
				//for (size_t y = 0; y < Pumping_DATA.size(); ++y)
				//{
				//	 flui.push_back(Pumping_DATA[y][5]);
				//}
				//std::cout << "flui size = " << flui.size() << std::endl;
				//ai::printVector(flui);
				//double tmp_T1;
				//for (int i = flui.size() - 1; i>0; i--)
				//	for (int j = 0; j<(i - 1); j++)
				//	{
				//
				//		if (flui[j] == flui[j + 1])
				//			flui.erase(flui.begin() + j + 1);
				//	}
				////� ������� ���������
				////unique(flui.begin(), flui.begin() + flui.size() - 1);
				////flui.resize(flui.begin() + flui.size() - 1 - flui.begin());

				//std::cout << "flui size = " << flui.size()<< std::endl;


				//��������� ����. ��������
				std::string fl = Pumping_DATA[index][5];
				//std::cout << "fl = " << fl << std::endl;
				std::vector<std::vector<json>> F_Base1 = J_IN_DATA["FluidBase"]["data"];
				//std::cout << "F_Base1.data() = " << F_Base1.size() << std::endl;
				//int stages1 = F_Base_data1[0][0]["rowCount"].get<int>();
				//��������� ��������
				if (flui.size() > 1) {
					for (size_t j = 0; j < F_Base1.size(); ++j)
					{
						std::string fbt = F_Base1[j][0];
						// std::cout << "fbt = " << fbt << std::endl;
						if (fl == fbt)
						{
							//std::cout << "IN" << std::endl;;
							std::vector<std::vector<json>> F_Base_data1 = F_Base1[j][6]["data"];

							std::vector<std::vector<json> > F_operator_time1 = F_Base_data1[0][0]["data"];
							std::string viscosity_s = F_operator_time1[0][1];

							//std::cout << "viscosity_s = " << viscosity_s << std::endl;
							double viscosity;
							sscanf(viscosity_s.c_str(), "%lf", &viscosity);

							viscosities.push_back(viscosity);

							std::string reology_s = F_operator_time1[0][2];
							double reology;
							sscanf(reology_s.c_str(), "%lf", &reology);
							reologies.push_back(reology);

						}
					}
				}
			}
			///////////////////////////////////////////////////////////////////

			// std::cout << " pre viscosities" << std::endl;
			// if (viscosities.size()>0)
			// ai::printVector(viscosities);

			//���������� ���������� ������ �������
			std::string end_time_s = Pumping_DATA[index][1];
			double end_time;
			sscanf(end_time_s.c_str(), "%lf", &end_time);
			//������ �������� �� ������, ���. �
			std::string FluodVol_s = Pumping_DATA[index][2];
			double FluodVol;
			sscanf(FluodVol_s.c_str(), "%lf", &FluodVol);
			FluodVol = FluodVol * 60; //;/ end_time;					//���������� � �3/���
			fluid.push_back(FluodVol);
			//����� ��������� �� ������, �� ����!!!!!!!!!!!!
			std::string PropMass_s = Pumping_DATA[index][7]; ///Pumping_DATA[index][13];
			double PropMass;
			sscanf(PropMass_s.c_str(), "%lf", &PropMass);

			std::string PropMass_s1 = Pumping_DATA[index][8]; ///Pumping_DATA[index][13];
			double PropMass1;
			sscanf(PropMass_s1.c_str(), "%lf", &PropMass1);

			prop.push_back((PropMass + PropMass1) / 2.);
		}

		//���� ������ ���������� ��������� � ���� ����������
		std::string proppant_type = Pumping_DATA[index][6];
		std::vector<std::vector<json>> P_Base = J_IN_DATA["ProppantBase"]["data"];
		//for (int index = 0; index < P_Base.size(); index++)
		//{
		//	if (proppant_type == P_Base[index][0])
		//		// std::cout << "P_type: " << P_Base[index][1] << " ";
		//}

		//���� ������ ���������� �������� � ���� ���������
		std::string fluid_type = Pumping_DATA[index][5];
		std::vector<std::vector<json>> F_Base = J_IN_DATA["FluidBase"]["data"];
		// std::cout<<"F_Base.size() "<<F_Base.size()<<std::endl;
		for (int index = 0; index < F_Base.size(); index++)
		{



			//������ ��� �������������� �������� ��������
			if (fluid_type == F_Base[index][0])
			{
				//std::cout << "F_type: " << F_Base[index][1] << " \n";
				//��� ���������� �������� �������� ��������� ������ �� ����������.
				std::vector<std::vector<json>> F_Base_data = F_Base[index][6]["data"];
				//std::cout<<"F_Base_data size = "<< F_Base_data.size()<<std::endl;
				for (int index1 = 0; index1 < F_Base_data.size(); index1++)
				{
					double steptime;// = 0.;

					std::vector<std::vector<json> > F_operator_time = F_Base_data[0][0]["data"];
					int stages = F_Base_data[0][0]["rowCount"].get<int>();
					size_t j = 0;
					size_t k = 0;
					for (int index2 = 0; index2 < stages; index2++)
					{
						///////////////////////////////////////////////////////////////////////////
						//������� ��������� ������������� ��� ���������� ��������
						//��������� � ������ ����� ������, ������� ����� ����� �����������
						std::string end_time_s = F_operator_time[index2][0];
						double end_time;
						sscanf(end_time_s.c_str(), "%lf", &end_time);
						steptime = end_time * 60.;
						//std::cout << "Time: " << end_time << " \n";
						Time.push_back(steptime);

						//std::cout << " cicle flui size()" << flui.size() << std::endl;

						if (flui.size() == 1) {
							//����������� �������� �� ������� ��� ���������� ��������
							std::string viscosity_s = F_operator_time[index2][1];
							double viscosity;
							sscanf(viscosity_s.c_str(), "%lf", &viscosity);

							viscosities.push_back(viscosity);
							//����������� �������� �� ������� ��� ���������� ��������
							std::string reology_s = F_operator_time[index2][2];
							double reology;
							sscanf(reology_s.c_str(), "%lf", &reology);

							reologies.push_back(reology);
							if (steptime > 0.)
								for (size_t i = index2 - 1 + j; i < Pumping_DATA.size(); ++i)
								{
									// std::cout << "steptime = " << steptime << std::endl;
									// std::cout << "TimeEnd" << i << " = " << TimeEnd[i] << std::endl;
									if (steptime > TimeEnd[i])
									{
										j++;
										/*viscosities.push_back(viscosity);
										reologies.push_back(reology);*/
										viscosities.insert(viscosities.begin() + i, viscosities[i]);
										reologies.insert(reologies.begin() + i, reologies[i]);

										//std::cout << "I = " << i << std::endl;
									}

								}
							if (steptime > 0.)
								for (size_t i = index2; i < Pumping_DATA.size(); i++)
								{
									if (steptime < TimeEnd[i] && steptime > TimeStart[i])
									{
										fluid.insert(fluid.begin() + i + k, fluid[i + k]);
										k++;
									}
								}
						}
					}
				}
			}
		}

		// std::cout<<"All times"<<std::endl;
		// ai::printVector(Time);
		double tmp_T;
		for (int index_1 = Time.size() - 1; index_1>0; index_1--)
			for (int index = 0; index<(index_1 - 1); index++)
			{
				if (Time[index] > Time[index + 1])
				{
					tmp_T = Time[index + 1];
					Time[index + 1] = Time[index];
					Time[index] = tmp_T;
				}
				if (Time[index] == Time[index + 1])
					Time.erase(Time.begin() + index + 1);
			}
		//� ������� ���������
		unique(Time.begin(), Time.begin() + Time.size() - 1);
		Time.resize(Time.begin() + Time.size() - 1 - Time.begin());
		//���, ������ Time �������� ������ ���������� �����
	}
	// std::cout << "time size " << Time.size() << std::endl;
	//ai::printVector(Time);

	//start_time = start_time+end_time;


	// std::cout << "Startime" << std::endl;
	// ai::printVector(TimeStart);

	// std::cout << "Time End" << std::endl;
	// ai::printVector(TimeEnd);

	std::vector<double > ti;
	//ti.push_back(0.);

	for (size_t i = 0; i < Time.size(); i++) ti.push_back(Time[i] / 60.);

	// std::cout << "inj" << std::endl;

	// ai::printVector(fluid);


	// std::cout << "ti" << std::endl;

	// ai::printVector(ti);


	// std::cout << " new ti" << std::endl;

	// ai::printVector(ti);


	// std::cout << "viscosities" << std::endl;
	// if (viscosities.size()>0)
	// ai::printVector(viscosities);

	//	std::cout << "time fluid" << std::endl;
	//	ai::printVector(timefl);


	// std::cout << "realogies" << std::endl;
	// if (reologies.size()>0)
	// ai::printVector(reologies);


	for (int index = 0; index < Time.size() - 1; index++)
		injection.push_back(std::vector<double> {ti[index], ti[index + 1],
			fluid[index]<0.00001 ? 0. : fluid[index], prop[index]<0.00001 ? 0. : prop[index],
			reologies[index]<0.000001 ? reologies[index - 1] : reologies[index],
			viscosities[index]<0.00001 ? viscosities[index - 1] : viscosities[index] });


	modelingTime = injection[injection.size() - 1][1];
	//	std::cout << "injection" << std::endl;
	//	ai::printMatrix(injection);

	//std::cout << "layers" << std::endl;
	//ai::printMatrix(layers_new);


	//for (int index = layers.size() - 1; index > 0; index--)
	//{
	//	maxstress = ai::max(maxstress, layers[index].stress);
	//	minstress = ai::min(minstress, layers[index].stress);
	//	std::cout << layers[index].start << "  " << layers[index].end << "  " << layers[index].stress << "  " << layers[index].poisson << "  " << layers[index].young << "  " << layers[index].carter << std::endl;
	//}
	//	std::vector<std::vector<double >> layers_new;
	for (int index = layers.size() - 1; index >= 0; index--)
		layers_new.push_back(layers[index]);
	//layers = layers_new;

	//////////////////////////////////////////////////////////////////////////////////////////
	//End function ImportJSON
	//////////////////////////////////////////////////////////////////////////////////////////
	return;
}



/*!
\brief  �������������� �������� ������ � ��������

\details ������� ���������� ���������� � �������� �����.

\param[in,out] layers - ������� ��������� ��������
*/
void ApproximateLayers(std::vector<std::vector<double> >&layers) {

	std::size_t num;
	double H;
	double dx;


	//=======================================================================
	std::vector<std::vector<double> >la;
		for (size_t i = layers.size(); i > 0; ) {
		--i;
		la.push_back(std::vector<double>{layers[i][0], layers[i][1], layers[i][2], layers[i][3], layers[i][4], layers[i][5]});
	}
	layers = la;
	la.clear();
//=======================================================================

	for (std::size_t i = 1; i < layers.size(); ++i) {
		if (layers[i][0] * layers[i][1] <= 0.) {
			H = layers[i][1] - layers[i][0];
			num = i;
			break;
		}
	}
	dx = 0.09 * H;
	double sigma = layers[num][2];

	std::vector<std::size_t> sloi;

	double h;
	double hnew;
	if (H <= 11.2) {

		std::size_t it;
		it = num;

		sloi.push_back(it);

		h = layers[it][1] - layers[it][0];
		it--;


		while (11.2 / 2. > std::abs(layers[it][1])) {

			hnew = layers[it][1] - layers[it][0];

			sigma = (sigma * h + hnew * layers[it][2]) / (h + hnew);
			h = h + hnew;

			sloi.push_back(it);
			it--;
		}

		it = num + 1;

		while (11.2 / 2. > layers[it][0]) {

			hnew = layers[it][1] - layers[it][0];

			sigma = (sigma * h + hnew * layers[it][2]) / (h + hnew);
			h += hnew;

			sloi.push_back(it);
			it++;

		}
		H = 11.2;
		dx = 0.09 * H;
	}

	if (H > 11.2) {

		std::size_t it;
		it = num;

		sloi.push_back(it);

		h = layers[it][1] - layers[it][0];
		it--;

		hnew = layers[it][1] - layers[it][0];

		sigma = (sigma * h + hnew * layers[it][2]) / (h + hnew);
		h = h + hnew;

		sloi.push_back(it);
		it--;

		it = num + 1;

		hnew = layers[it][1] - layers[it][0];

		sigma = layers[num][2];
		h += hnew;

		sloi.push_back(it);
		it++;
	}
	std::size_t min = ai::min(sloi);
	std::size_t max = ai::max(sloi);

	std::vector<std::vector<double> > layers1;

	for (size_t i = 0; i < layers.size(); i++) {
		layers1.push_back(std::vector<double>{layers[i][0], layers[i][1], layers[i][2], layers[i][3], layers[i][4], layers[i][5]});
		if (i == min) {
			layers1.push_back(std::vector<double>{-H / 2., H / 2., sigma, layers[i][3], layers[i][4], layers[i][5]});
			i = max - 1;
		}

	}
	// ������������ ����������� �����


	for (size_t i = 0; i < layers1.size(); i++) {
		if (layers1[i][0] * layers1[i][1] < 0.) num = i;
	}

	for (size_t i = 0; i < layers1.size() - 1; ++i) {
		if (layers1[i][1] > layers1[i + 1][0]) {

			if (i != num)
				layers1[i][1] = layers1[i + 1][0];

			if (i == num)
				layers1[i + 1][0] = layers1[i][1];

		}

	}
	std::vector<std::vector<double> > mesh;
	//���������� ��� ������ �������������

	//��� ������������ ����
	mesh.push_back(std::vector<double>{layers1[num][0], layers1[num][1],

		layers1[num][2], layers1[num][3], layers1[num][4], layers1[num][5]});
	double up = layers1[num][0];

	while (std::abs(layers1[0][0]) > std::abs(up)) {
		mesh.push_back(std::vector<double>{up - dx, up, 0, 0, 0, 0});
		up -= dx;
	}



	std::vector<std::vector<double> > mesh1;
	for (size_t i = mesh.size() - 1; i>0; --i) {
		mesh1.push_back(std::vector<double>{mesh[i][0], mesh[i][1], mesh[i][2], mesh[i][3], mesh[i][4], mesh[i][5] });
	}
	mesh1.push_back(std::vector<double>{layers1[num][0], layers1[num][1],

		layers1[num][2], layers1[num][3], layers1[num][4], layers1[num][5]});




	up = layers1[num][1];

	size_t f = layers1.size() - 1;

	while (std::abs(layers1[f][1]) > std::abs(up)) {

		mesh1.push_back(std::vector<double>{up, up + dx, 0, 0, 0, 0});
		up += dx;
	}
	mesh.clear();



	double eps = pow(10, -8);

	size_t middle;
	for (size_t i = 0; i < mesh1.size(); ++i) {
		if (mesh1[i][0] * mesh1[i][1] < 0.) {
			middle = i;
			break;
		}
	}

	double rastlayer;
	double rastmesh;
	double dl;
	double hadd;


	size_t j = middle - 1;


	for (size_t i = num - 1; i >= 0; --i) {  //���� �� ����� layers1 i - ��������

		h = layers1[i][1] - layers1[i][0];     //j - �������� �� mesh1

		dl = std::ceil(h / dx);
		if ((int)dl == 1) {


			mesh1[j][2] += h * layers1[i][2];
			mesh1[j][3] += h * layers1[i][3];
			mesh1[j][4] += h * layers1[i][4];
			mesh1[j][5] += h * layers1[i][5];
			if (i >= 0) {

				i--;
				mesh1[j][2] += (dx - h)*layers1[i][2];
				mesh1[j][3] += (dx - h)*layers1[i][3];
				mesh1[j][4] += (dx - h)*layers1[i][4];
				mesh1[j][5] += (dx - h)*layers1[i][5];

				mesh1[j][2] /= dx;
				mesh1[j][3] /= dx;
				mesh1[j][4] /= dx;
				mesh1[j][5] /= dx;
			}
			i++;
			j--;
		}

		if ((int)dl > 1) {
			int iter = 0;

			double dh;

			while (iter < (int)dl) {
				if (j >= 0) {

				    mesh1[j][2] = layers1[i][2];
					mesh1[j][3] = layers1[i][3];
					mesh1[j][4] = layers1[i][4];
					mesh1[j][5] = layers1[i][5];
					if (j == 0) { break; }

					j--;

				}
				iter++;


			}

		}
		if (i == 0) { break; }
	}

	j = middle + 1;


	for (size_t i = num + 1; i < layers1.size(); ++i) {  //���� �� ����� layers1 i - ��������
														 // std::cout<<" "<<std::endl;
														 // std::cout<<"iteration = "<<i<<std::endl;
		h = layers1[i][1] - layers1[i][0];     //j - �������� �� mesh1
											   // std::cout<<"h curent layer = "<<h<<std::endl;
		dl = std::floor(h / dx + 0.51);
		if ((int)dl == 1) {

			mesh1[j][2] += h * layers1[i][2];
			mesh1[j][3] += h * layers1[i][3];
			mesh1[j][4] += h * layers1[i][4];
			mesh1[j][5] += h * layers1[i][5];

			if (i < layers1.size()) {
				i++;
				mesh1[j][2] += (dx - h)*layers1[i][2];
				mesh1[j][3] += (dx - h)*layers1[i][3];
				mesh1[j][4] += (dx - h)*layers1[i][4];
				mesh1[j][5] += (dx - h)*layers1[i][5];
				}
			i--;
			j++;
		}

		if ((int)dl > 1) {
			int iter = 0;
			h = 0;
			while (iter <= (int)dl) {
				if (j < mesh1.size()) {
					mesh1[j][2] = layers1[i][2];
					mesh1[j][3] = layers1[i][3];
					mesh1[j][4] = layers1[i][4];
					mesh1[j][5] = layers1[i][5];

					if (j == mesh1.size() - 1) { break; }

					j++;
				}
				iter++;
			}
			if (i == layers1.size() - 1) { break; }
			if (j == mesh1.size() - 1) { break; }

		}
	}

	std::vector<std::vector<double> > mesh2;
	double start, end;

	for (size_t i = 0; i < mesh1.size() - 1; ++i) {

		sigma = mesh1[i][2];
		end = mesh1[i][1];
		start = mesh1[i][0];
		while (mesh1[i][2] == mesh1[i + 1][2]) {
			sigma = mesh1[i + 1][2];
			end = mesh1[i + 1][1];
			i++;

			if (i == mesh1.size() - 1) { break; }

		}
		mesh2.push_back(std::vector<double>{start, end, sigma, mesh1[i][3], mesh1[i][4], mesh1[i][5]});


	}

	for (size_t i = 0; i < mesh2.size(); ++i) {
		if (mesh2[i][1] * mesh2[i][0]<0.) {
			middle = i;
			break;
		}
	}

	for (size_t i = middle; i>1; --i) {
		if (mesh2[i - 1][2] - mesh2[i][2]>6.) {
			mesh2[i - 1][2] = 6. + mesh[i][2];
			}
	}

	for (size_t i = middle; i<mesh2.size() - 1; ++i) {
		if (mesh2[i + 1][2] - mesh2[i][2]>6.) {
			mesh2[i + 1][2] = 6. + mesh[i][2];
		}
	}
	std::vector<std::vector<double> > layerss;


	for (size_t i = mesh2.size(); i > 0; ) {
		--i;
		layerss.push_back(std::vector<double>{mesh2[i][0], mesh2[i][1], mesh2[i][2], mesh2[i][3], mesh2[i][4], mesh2[i][5]});
	}
	layers = layerss;
//////////////////////////////////////////////////////////////////////////////////////////
//End function ApproximateLayers
//////////////////////////////////////////////////////////////////////////////////////////
}
