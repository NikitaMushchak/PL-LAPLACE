#pragma once

#include <vector>
#include <thread>
#include <condition_variable>
//#include "../jsoncpp/include/json/value.h"
//#include "../jsoncpp/include/json/reader.h"
#include <memory>

#include "../planar3D.hh"

#include "api.h"

namespace Solver
{
}
using namespace Solver;


class worker
{
public:
  worker();
  ~worker();

  Result Init(const char* rawData);
//  void read(const Json::Value &json);
  bool IsRunning() { return play; }
  bool IsPaused() { return pause; }
 // void Calculate();
 // void write(Json::Value &json) const;

  void NotifyFinished();
  void NotifyPaused();
  void NotifyRunning();
  void NotifyResult();

public:
  void DoWork();
  void Play();
  void Pause();
  void Stop();

private:
  bool play, pause, stop;

  std::condition_variable pauseManager;
  ///// Planar variables
  std::mutex mutex;
  DLL_Param DLL_Parametrs;

  std::vector< std::vector<double> > layers;
  std::vector< std::vector<double> > injection;
  std::string IdDesign;
  std::string IdStage;
  double modelingTime = -1;
  double Emit_time = -1;
  double Z_coordinate;
  std::string JSONstring;		
  bool initialized=false;
  /////End Planar variables


  //-----settings variables
  double K;
  double n;
  double E;
  double mu;
  double t;
  double q;
  double C;
  double dx;
  double dh;

  int n_time_nodes;
  int nx;
  double toughness;
  double leakOffCoeff;
  std::string modelType;

  // ----calculation variables
  double Vin;
  double a;
  double b;
  double visc;
  double td;
  double dt;
  double l;
  double w0;
  double Vf;
  double gamma;

  int plotCurves;
  int counter;

  bool start;

  std::vector<double> w;
  std::vector<double> L;

//  std::unique_ptr<Json::Value> resultObject;
};
