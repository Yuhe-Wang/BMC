#include "../Tools/Tools.h"

// the simulation dll/so at least provides the following two functions
typedef void(*StartSimulation)(const char*, MPS&, bool);
typedef void(*StopSimulation)();

StopSimulation stopSimulation = NULL;

//Let's put some global variables here
EmailNotify emailNotify;

#ifdef WIN32
BOOL WINAPI HandlerRoutine(_In_ DWORD dwCtrlType)
{
	stopSimulation();
	getchar();
	return TRUE;
}
#else //LINUX
void HandlerRoutine(int sig)
{
	stopSimulation();
	printf("press any key to exit...");
	getchar();
	return;
}
#endif

int main(int argc, char *argv[])
{
#ifdef WIN32
	WORD wVersionRequested;
	WSADATA wsaData;
	int err;
	wVersionRequested = MAKEWORD(2, 2);
	err = WSAStartup(wVersionRequested, &wsaData);
	if (err != 0)
	{
		printf("WSAStartup failed with error: %d\n", err);
		return 1;
	}

	//add handler when closing this program
	SetConsoleCtrlHandler(HandlerRoutine, TRUE);
#else // LINUX
	signal(SIGINT, HandlerRoutine);
#endif

	ConfigFile jobList;

	//detect if we have job list config file
	bool useJobList = jobList.parse("jobList.py");
	if (!useJobList) exitApp("Cannot open jobList.py");

	ConfigFile *pcf = NULL;
	pcf = jobList.getBlock("email");
	if (pcf) emailNotify.init(pcf);

	//do the micro replacement
	MPS configMacro;
	pcf = jobList.getBlock("ConfigVariables");
	if (pcf)
	{
		pcf->copyMps(configMacro);
		for (unsigned int i = 0; i < configMacro.size(); ++i)
		{
			string d("$");
			configMacro[i].first = d + configMacro[i].first + "$";
		}
		jobList.macroReplace(configMacro);
	}

	//load the dll/so first
	string engine = "gZeus"; // default engine name
	if (!jobList.getValue("MC engine", engine)) // If the config file doesn't specify the dll/so module
	{	//we use the name of this executable file name as the module name
		string mc(argv[0]);
		size_t islash = mc.find_last_of("\\/");
		size_t idot = mc.find_last_of(".");
		if (idot == string::npos) mc = mc.substr(islash + 1, 256);
		else mc = mc.substr(islash + 1, idot - islash - 1);
		engine = mc;
	}
	string BatchDir = engine;
	BatchDir += FileSeparator;
	BatchDir += "BatchMode";
#ifdef WIN32
	engine += ".dll";
	HINSTANCE hlib = LoadLibrary(engine.c_str());
#else
	engine += ".so";
	void* hlib = dlopen(engine.c_str(), RTLD_LAZY);
    if (hlib == NULL) Log("Failed to open the so file: %s", dlerror());
#endif
	if (hlib == NULL)
	{
		string msg = "Cannot load the Monte Carlo engine from ";
		msg += engine;
		exitApp(msg.c_str());
	}

	StartSimulation startSimulation = (StartSimulation)getLibAddress(hlib, "startSimulation");
	if (startSimulation == NULL) exitApp("cannot load function startSimulation");
	stopSimulation = (StopSimulation)getLibAddress(hlib, "stopSimulation");
	if (stopSimulation == NULL) exitApp("cannot load function stopSimulation");

	string batchMode;
	jobList.getValue("BatchMode", batchMode);
	if (0 == batchMode.compare("yes")) //search job folders in directory BatchMode
	{
		fs::directory_iterator end_itr;
		for (fs::directory_iterator ifile(BatchDir); ifile != end_itr; ++ifile)
		{
			if (fs::is_directory(ifile->status())) //check the directory
			{
				string dir = ifile->path().string();
				printf("Running the job defined in directory: %s", dir.c_str());

				string result = findFirstFileInDir(dir, ".*\\.skip");
				if (!result.empty())
				{
					printf("Will skip the job defined in directory: %s", dir.c_str());
					continue;
				}

				result = findFirstFileInDir(dir, ".*\\.py");
				if (result.empty())
				{
					printf("Will skip the job defined in directory: %s", dir.c_str());
					continue;
				}
				else startSimulation(result.c_str(), configMacro, true);
			}
		}

	}
	else //use job config defined in this file
	{
		//execute the job configurations defined in this file in sequence
		ConfigFile *jobConfig = NULL;
		jobList.resetSearchIndex(); //make sure we can find jobs in sequence
		int job = 0;
		while (true)
		{
			jobConfig = jobList.getBlock("job", true);
			if (NULL == jobConfig) break;
			++job;
			printf("\n/*>>>>>>>>>>>>>>>>>>>>>>>>>>> Starting job %d >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/\n\n", job);

			int NRepeat = 1;
			jobConfig->getValue("repeat times", NRepeat);
			if (NRepeat < 0) NRepeat = 0; //skip this job
			for (int nr = 0; nr < NRepeat; ++nr)
			{
				if (NRepeat > 1) printf("\n/*>>>>>>>>>>>>>>>>>>>>>>>>>>> Repeat %dth time >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/", nr + 1);
				string cmd;
				jobConfig->resetSearchIndex();
				while (jobConfig->getValue("pre command", cmd, true)) { system(cmd.c_str()); }
				if (!jobConfig->getValue("config file", cmd)) exitApp("cannot find config file for this job");
				//execute this job by config file
				startSimulation(cmd.c_str(), configMacro, true);
				jobConfig->resetSearchIndex();
				while (jobConfig->getValue("post command", cmd, true)) { system(cmd.c_str()); }
			}
			printf("\n/*<<<<<<<<<<<<<<<<<<<<<<<<<<< Finished job %d <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/\n\n", job);
			double restTime = 0;
			if (jobConfig->getValue("rest time", restTime) && restTime != 0)
			{
				printf("The job execution will sleep %f minutes to rest the CPU,\n Please wait...\n\n", restTime);
				//pauseSeconds(unsigned int(restTime * 60));
				std::this_thread::sleep_for(std::chrono::seconds(long(restTime * 60)));
			}
		}
	}

#ifdef WIN32
	WSACleanup();
#endif

	
	if (emailNotify.doNotify()) emailNotify.send(true, NULL);

	string exitConfirm;
	jobList.getValue("confirm to exit", exitConfirm);
	if (exitConfirm.compare("yes") == 0)
	{
		printf("\n\n Press enter to exit...");
		getchar();
	}

	if (hlib) freeLib(hlib);

	return 0;
}


void exitApp(const char *inf)
{
	printf("fatal error: %s", inf);
	if (emailNotify.doNotify()) emailNotify.send(false); //notify the user that the job was aborted abnormally
	printf("\nPress enter key to exit...");
	getchar();
	exit(-1);
}
