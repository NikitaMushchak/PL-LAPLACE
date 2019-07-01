//#include "../jsoncpp/include/json/worker.h"
//#include "../jsoncpp/include/json/api_callback.h"
#include "worker.h"
#include "../io.hh"
#include "api_callback.h"
#include "../planar3D.hh"
#include "ailibrary/ai.hh"

#include <algorithm>
#include <chrono>
#include <string>
#include <ctime>
#include <cmath>

#include <cstdio>
#include <windows.h>
#include <tlhelp32.h>



#include <windows.h>

using namespace Solver;


worker::worker()
  : play(false)
  , pause(false)
  , stop(false)
{
  plotCurves = 1;
}

worker::~worker()
{
  Stop();
}

void worker::NotifyFinished()
{
  StateCallback stateCallback = Callback::GetStateCallback();
  if (stateCallback)
  {
    stateCallback(Idle);
  }
}

void worker::NotifyPaused()
{
  StateCallback stateCallback = Callback::GetStateCallback();
  if (stateCallback)
  {
    stateCallback(Paused);
  }
}

void worker::NotifyRunning()
{
  StateCallback stateCallback = Callback::GetStateCallback();
  if (stateCallback)
  {
    stateCallback(Running);
  }
}

void worker::NotifyResult()
{
  DataCallback dataCallback = Callback::GetDataCallback();

  if (dataCallback)
  {
	FILE *stream;
	std::string output;
	errno_t err;
	err = fopen_s(&stream, "data_rez.json", "w+");
	if (err == 0)
	{
//		printf("The file 'data2' was opened\n");
	}
	else
	{
//		printf("The file 'data2' was not opened\n");
	}

	fprintf(stream,"%s", DLL_Parametrs.J_String.c_str());
	err = fclose(stream);


    dataCallback(DLL_Parametrs.J_String.c_str());//возврат в гуи
  }
}

Result worker::Init(const char* rawData)
{
	std::cout << "start Init" << std::endl;
	std::string strJson(rawData);

	ai::saveLine("inputjson", strJson.c_str());

	//  std::string mystr;
	// ai::parseFileIntoString("init.json", strJson);

	DLL_Parametrs.J_String = strJson.c_str();

	std::cout << "start parce" << std::endl;
	ImportJSON(
		DLL_Parametrs,
		layers,
		injection
	);

	ai::saveMatrix("layers_in", layers);
	std::cout << "end parce" << std::endl;

	ApproximateLayers(layers, DLL_Parametrs.dx);	//перерасчет слоев и интерполирование на планаровскую сетку

	ai::saveMatrix("layers_appr", layers);
	// Пересчитываем значения для слоёв в величины СИ
	for (std::size_t i = 0; i < layers.size(); ++i) {
		layers[i][2] *= std::pow(10, 6);
		layers[i][3] *= std::pow(10, 9);
		layers[i][5] *= std::pow(10, -6);
	}


	 //ai::saveMatrix("layers_planar", layers);
	//
	ai::saveMatrix("injection", injection);

	counter = 0;

	dt = t / plotCurves;
	std::cout << "end init" << std::endl;

	return Success;
}

void worker::DoWork()
{
	std::cout << "start work" << std::endl;

  std::srand(std::time(nullptr));

  // Reset here
  //resultObject.reset(new Json::Value);
  counter = 0;
  stop = false;
  double RET_MASSIV[1000] = { 0 };
  
  NotifyRunning();
//  Calculate(); //Тут вставляем расчет планара

#define BUILD_DLL
std::cout << "start planar" << std::endl;
  int a = planar3D(
	  DLL_Parametrs,
	  layers,
	  injection,
	  initialized
  );
  std::cout << "end planar. Return code =" << a << std::endl;

  initialized = true;
  
  NotifyResult();      
  NotifyFinished();

}

//#define BUILD_DLL

//void worker::Calculate()
//{
//	//Здесь осуществляем запуск планара с параметрами
//	 #define BUILD_DLL
//	int a = planar3D(
//		DLL_Parametrs,
//		layers,
//		injection,
//		initialized
//	);
//	initialized = true;
//}

void worker::Play()
{
	std::cout << "start play" << std::endl;
  mutex.lock();
  play = true;
  pause = false;
  stop = false;
  mutex.unlock();
  pauseManager.notify_all();
  DLL_Parametrs.Dll_State = DLL_Parametrs.Running;
  NotifyRunning();
  std::cout << "end play" << std::endl;

}

void worker::Pause()
{
  mutex.lock();
  play = false;
  pause = true;
  stop = false;
  mutex.unlock();
  DLL_Parametrs.Dll_State = DLL_Parametrs.Paused;
  NotifyPaused();
}

void worker::Stop()
{
  mutex.lock();
  play = false;
  pause = false;
  stop = true;
  mutex.unlock();
  DLL_Parametrs.Dll_State = DLL_Parametrs.Idle;
  pauseManager.notify_all();
}
