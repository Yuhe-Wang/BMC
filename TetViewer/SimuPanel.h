#ifndef _SIMUPANEL_H_
#define _SIMUPANEL_H_

#include "PreCompile.h"
#include "GLView.h"

wxDECLARE_EVENT(SIMULATION_FINISHED, wxCommandEvent); //give the reference of the customized event const

struct OneJob
{
	OneJob()
	{
		skip = false;
		percent = 0;
		timeUsed = 0;
		timeLeft = 0;
		speed = 0;
		uncertainty = 0;
		checked = false;
	}
	string configFile; // BatchMode/JobName/*.py
	bool skip;
	double percent;
	double timeUsed;
	double timeLeft;
	double speed;
	double uncertainty;

	bool checked; //auxiliary variable when updating the job list
	string getJobName()
	{
		size_t pos1 = configFile.find_last_of("\\/");
		size_t pos2 = configFile.find_first_of("\\/");
		return configFile.substr(pos2 + 1, pos1 - pos2 - 1);
	}
	string getJobDir()
	{
		size_t pos1 = configFile.find_last_of("\\/");
		return configFile.substr(0, pos1 + 1);
	}
	void removeSkipFile()
	{
		skip = false;
		size_t pos1 = configFile.find_last_of("\\/");
		string dir = configFile.substr(0, pos1 + 1);
		dir += "txt.skip";
		wxRemoveFile(dir.c_str());
	}
	bool addSkipFile()
	{
		size_t pos1 = configFile.find_last_of("\\/");
		string dir = configFile.substr(0, pos1 + 1);
		dir += "txt.skip";
		FILE* fp = fopen(dir.c_str(), "w");
		if (NULL == fp) return false;
		fprintf(fp, "If this file exist in a job folder, this job will be skipped");
		fclose(fp);
		skip = true;
		return true;
	}
};

//definitions of interface functions of the MC program
typedef void(*StartSimulation)(const char*, MPS&, bool);
typedef void(*SimpleCall)();
typedef void(*SetProgressCallBack)(ProgressCallBack);
typedef void(*SetPeekDoseCallBack)(PeekDoseCallBack);
typedef void(*SetJobFinishedCallBack)(JobFinishedCallBack);
typedef void(*SetLogCallBack)(LogCallBack);

class SimuPanel : public wxPanel
{
public:
	SimuPanel(wxWindow* parent);
	void OnSimuStart(wxCommandEvent& event);
	void OnSimuStop(wxCommandEvent& event);
	void OnUpdateJobList(wxCommandEvent& event);
	void OnMoveJobUp(wxCommandEvent& event);
	void OnMoveJobDown(wxCommandEvent& event);
	void OnEditJobListPY(wxMouseEvent& event);
	void OnShowOption(wxGridEvent& event);
	void OnJobSelected(wxGridEvent& event);
	void OnJobDClick(wxGridEvent& event);
	void OnPopupClick(wxCommandEvent &event);
	void OnOneJobFinished(bool abortJob);
	void OnTimer(wxTimerEvent& event);
	void OnPeekDose(wxCommandEvent& event);
	void OnCheckShowLog(wxCommandEvent& event);
	void OnCheckAutoPeek(wxCommandEvent& event);
	void OnSetPeekByDivision(wxCommandEvent &event);
	void OnSetPeekByTime(wxCommandEvent &event);

	void progressUpdate(double percent, double speed, double time, double time_left, double uncertainty);
	void setTimer();
	void getBinaryFILE(const char* config, BinaryFile& BF);
	bool isRunning(){ return _running; }
public:
	wxGrid* _grid;
	list<OneJob> jobList;
	MPS configMacro;
	wxTextCtrl* _logText;
private:
	bool getJobList(); // Vector<OneJob> joblist; 
	void showJobStatus();
	
	int _selectedRow;
	int _cJob; // job that will be processed(if not running) or being processed (if running)
	bool _running;
	bool _autoPeek;
	int _loadedEngine;// the index engine that has been loaded
	int _peekWay; // 0 means by division, 1 means by time intervals
	int _NPeekDivision;
	int _iLastPeek;
	int _PeekTimeSpan; //unit seconds
	wxDynamicLibrary wxDL;
	wxTimer _peekTimer;
	//function pointer of the calculation engine
	StartSimulation startSimulation;
	SimpleCall peekDose;
	SimpleCall stopSimulation;
	SetProgressCallBack setProgressCallBack;
	SetPeekDoseCallBack setPeekDoseCallBack;
	SetJobFinishedCallBack setJobFinishedCallBack;
	SetLogCallBack setLogCallBack;
	SetLogCallBack setExitCallBack;

	DECLARE_EVENT_TABLE()
};

#endif
