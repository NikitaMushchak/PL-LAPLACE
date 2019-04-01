#pragma once

#ifdef SOLVER_EXPORTING_API
   #define SOLVER_API extern "C" __declspec(dllexport)
#else  
   #define SOLVER_API extern "C" __declspec(dllimport)
#endif

namespace Solver
{
  /*
    Solver API call result 
  */
  enum Result
  {
    Success = 0,
    Error   = 1
  };

  /*
    Solver state
  */
  enum State
  {
    Idle      = 0,
    Running   = 1,
    Paused    = 2,
  };

  /*
    State callback function
    @param state - new solver state
  */
  typedef void (__cdecl *StateCallback)(State state);

  /*
    @param callbackFn - state callback function
    @return Result    - one of the Solver::Result values
  */
  SOLVER_API Result __cdecl SetStateCallback(StateCallback stateCallbackFn);

  /*
    Data callback function
    @param data       - zero terminated string representing a JSON object containing solver result data
  */
  typedef void (__cdecl *DataCallback)(const char* data);

  /*
    @param dataCallbackFn - data callback function
  */
  SOLVER_API Result __cdecl SetDataCallback(DataCallback dataCallbackFn);

  /*
    @param parameters - zero terminated string representing a JSON object containing data for seeding
    @return Result    - one of the Solver::Result values
  */
  SOLVER_API Result __cdecl SetParameters(const char* parameters);

  /*
    Begin or resume paused calculation
        
    @return Result    - one of the Solver::Result values
  */
  SOLVER_API Result __cdecl Run();
  
  /*
    Pause calculation request

    @return Result    - one of the Solver::Result values
  */
  SOLVER_API Result __cdecl Pause();
  
  /*
    Stop calculation request

    @return Result    - one of the Solver::Result values
  */
  SOLVER_API Result __cdecl Stop();

  /*
    Wait for Solver::Idle State
  */
  SOLVER_API Result __cdecl Wait();

  /*
    Shutdown
  */
  SOLVER_API Result __cdecl Shutdown();
}
