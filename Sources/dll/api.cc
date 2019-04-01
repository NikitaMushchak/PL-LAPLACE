#include "api.h"
#include "worker.h"
//#include "api.h"
//#include "worker.h"
#include <memory>
#include <chrono>
#include <thread>

namespace Solver
{
  namespace 
  {
    static StateCallback StateCallbackFn__ = 0;
    static DataCallback  DataCallbackFn__  = 0;

    std::unique_ptr<worker> Worker__;
    std::unique_ptr<std::thread> WorkerThread__;
    std::mutex WorkerMutex__;
  }

  namespace Callback
  {
    StateCallback GetStateCallback()
    {
      return StateCallbackFn__;
    }

    DataCallback GetDataCallback()
    {
      return DataCallbackFn__;
    }
  }

  Result SetStateCallback(StateCallback stateCallbackFn)
  {
    StateCallbackFn__ = stateCallbackFn;
    return Success;
  }

  Result SetDataCallback(DataCallback dataCallbackFn)
  {
    DataCallbackFn__ = dataCallbackFn;
    return Success;
  }

  Result SetParameters(const char* parameters)
  {
    std::lock_guard<std::mutex> lock(WorkerMutex__);
    if (!Worker__.get())
    {
      Worker__.reset(new worker);
    }

    Worker__->Init(parameters);
    return Success;
  }

  Result Run()
  {
    std::lock_guard<std::mutex> lock(WorkerMutex__);
    if (!Worker__.get())
      return Error;

    if (Worker__->IsPaused())
    {
      Worker__->Play();
      return Success;
    }

    if (WorkerThread__.get())
      WorkerThread__->join();
    
    WorkerThread__.reset(new std::thread(&worker::DoWork, Worker__.get()));
    return Success;
  }
  
  Result Pause()
  {
    std::lock_guard<std::mutex> lock(WorkerMutex__);
    if (!Worker__.get())
      return Error;

    Worker__->Pause();
    return Success;
  }
  
  Result Stop()
  {
    std::lock_guard<std::mutex> lock(WorkerMutex__);
    if (!Worker__.get())
      return Error;

    Worker__->Stop();
    return Success;
  }

  Result Wait()
  {
    std::lock_guard<std::mutex> lock(WorkerMutex__);
    if (WorkerThread__.get())
    {
      WorkerThread__->join();
      WorkerThread__.reset(0);
    }
    return Success;
  }

  Result Shutdown()
  {
    std::lock_guard<std::mutex> lock(WorkerMutex__);
    if (Worker__.get())
      Worker__->Stop();

    if (WorkerThread__.get())
      WorkerThread__->join();

    Worker__.reset(0);
    WorkerThread__.reset(0);
    return Success;
  }
}
