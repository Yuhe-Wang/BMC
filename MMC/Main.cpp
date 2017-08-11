#include "Penelope.h"

// Begin: essential for command line call
bool b_abort = false;
bool b_thread_active = false;
LogProvider Log;

void executeJob(const char* configFileName, MPS& configMacro);

#ifdef WIN32
#define DEXPORT __declspec(dllexport)
#else
#define DEXPORT __attribute__ ((visibility ("default")))
#endif

extern "C" DEXPORT void startSimulation(const char* configFileName, MPS& configMacro, bool bWait) //launch a thread to do the simulation
{
    b_thread_active = true;
    std::thread thd(executeJob, configFileName, std::ref(configMacro));
    if (bWait) thd.join();
    else thd.detach();
}

extern "C" DEXPORT void stopSimulation()
{
    if (b_thread_active) b_abort = true; //only if the thread is on
    b_thread_active = false;
}
// End: essential for command line call

void exitApp(const char *inf)
{
    Log("fatal error: %s", inf);
    Log.flush(); //flush the log file's buffer

    printf("\nPress enter key to exit...");
    getchar();
    exit(-1);
}