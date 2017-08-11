#include "PreCompile.h"
#include "SimuPanel.h"

#define ID_CHANGE_SKIP 2011
#define ID_RECALCULATE 2012
#define ID_REMOVE_JOB 2013
#define ID_MOVEUP_JOB 2014
#define ID_MOVEDOWN_JOB 2015
#define ID_OPEN_JOB_IN_FOLDER 2016
#define ID_OPEN_JOB_CONFIG_FILE 2017

wxDEFINE_EVENT(SIMULATION_FINISHED, wxCommandEvent); //define customized event type

wxString time2print(double sec)
{
	wxString str;
	int time = int(sec);
	int hour = time / 3600;
	time = time % 3600;
	int min = time / 60;
	time = time % 60;
	if (hour == 0) str.Printf("%d:%d", min, time);
	else str.Printf("%d:%d:%d", hour, min, time);
	return str;
}

OneJob* listAt(list<OneJob>& jobList, int i) // no safe check here
{
	int j = 0;
	for (list<OneJob>::iterator it = jobList.begin(); it != jobList.end(); ++it)
	{
		if (j == i) return &(*it);
		++j;
	}
	return NULL;
}

void listRemoveAt(list<OneJob>& jobList, int i) // no safe check here
{
	int j = 0;
	for (list<OneJob>::iterator it = jobList.begin(); it != jobList.end(); ++it)
	{
		if (j == i)
		{
			jobList.erase(it);
			return;
		}
		++j;
	}
}

void listSwap(list<OneJob>& jobList, int i, int j)
{
	int k = 0;
	list<OneJob>::iterator it = jobList.begin();
	for (; it != jobList.end(); ++it)
	{
		if (k == i) break;
		++k;
	}
	k = 0;
	list<OneJob>::iterator jt = jobList.begin();
	for (; jt != jobList.end(); ++jt)
	{
		if (k == j) break;
		++k;
	}
	OneJob temp(*it);
	*it = *jt;
	*jt = temp;
}

// Here are the global variables that used by these call back functions
SimuPanel* gSimu = NULL;

void LogCB(const char* logLine)
{
	wxTextCtrl& logText = *(gSimu->_logText);
	logText << logLine << "\n";
}

void progressCB(double percent, double speed, double time, double time_left, double uncertainty)
{
	gSimu->progressUpdate(percent, speed, time, time_left, uncertainty);
}
void peekDoseCB(BinaryFile& BF)
{
	ggl->showBF(BF);
}
void jobFinishedCB(bool abortJob, BinaryFile& BF)
{
	ggl->showBF(BF);
	gSimu->OnOneJobFinished(abortJob);
}

void exitCB(const char* exitMessage)
{
	wxMessageBox(exitMessage, "Fatal Error");
}

BEGIN_EVENT_TABLE(SimuPanel, wxPanel)
EVT_BUTTON(XRCID("m_button_simuStart"), SimuPanel::OnSimuStart)
EVT_BUTTON(XRCID("m_button_simuStop"), SimuPanel::OnSimuStop)
EVT_BUTTON(XRCID("m_button_updateJobList"), SimuPanel::OnUpdateJobList)
EVT_BUTTON(XRCID("m_button_move_job_up"), SimuPanel::OnMoveJobUp)
EVT_BUTTON(XRCID("m_button_move_job_down"), SimuPanel::OnMoveJobDown)
EVT_BUTTON(XRCID("m_button_peek_dose"), SimuPanel::OnPeekDose)
EVT_TIMER(-1, SimuPanel::OnTimer)
EVT_CHECKBOX(XRCID("m_checkBox_showLog"), SimuPanel::OnCheckShowLog)
EVT_CHECKBOX(XRCID("m_checkBox_autoPeek"), SimuPanel::OnCheckAutoPeek)
EVT_RADIOBUTTON(XRCID("m_radioBtn_peekByDivision"), SimuPanel::OnSetPeekByDivision)
EVT_RADIOBUTTON(XRCID("m_radioBtn_peekByTime"), SimuPanel::OnSetPeekByTime)
END_EVENT_TABLE()

SimuPanel::SimuPanel(wxWindow* parent) :_cJob(0), _selectedRow(-1), _loadedEngine(-1), _running(false), _peekWay(0),
_NPeekDivision(20), _iLastPeek(0), _PeekTimeSpan(5), _autoPeek(true)
{
	wxXmlResource::Get()->LoadPanel(this, parent, "Panel_simulation");
	_grid = (wxGrid*)FindWindowById(XRCID("m_grid_tasks"));
	//_grid->SetScrollbars(10, 10, 50, 50);
	_grid->CreateGrid(10, 7);
	_grid->EnableEditing(false);//set it read-only
	//_grid->SetSelectionMode(wxGrid::wxGridSelectRows); // row selection only
	double DPIScale = wxBitmapScale::getDPIScale();
	_grid->SetColLabelSize(20*DPIScale);
	_grid->SetRowLabelSize(40*DPIScale);
	
	_grid->SetColLabelValue(0, wxT(" s "));
	_grid->SetColLabelValue(1, wxT(" config "));
	_grid->SetColLabelValue(2, wxT("percent"));
	_grid->SetColLabelValue(3, wxT(" time used"));
	_grid->SetColLabelValue(4, wxT(" time left"));
	_grid->SetColLabelValue(5, wxT(" speed"));
	_grid->SetColLabelValue(6, wxT(" uncertainty"));
	for (int i = 0; i < 7; ++i) _grid->AutoSizeColLabelSize(i);
	_grid->Bind(wxEVT_GRID_CELL_RIGHT_CLICK, &SimuPanel::OnShowOption, this);
	_grid->Bind(wxEVT_GRID_SELECT_CELL, &SimuPanel::OnJobSelected, this);
	_grid->Bind(wxEVT_GRID_CELL_LEFT_DCLICK, &SimuPanel::OnJobDClick, this);
	
	if(!getJobList()) return;
	showJobStatus();
	
	wxChoice* pEngineChoice = (wxChoice*)FindWindowById(XRCID("m_choice_MC_engine"));
	ConfigFile ext_config;
	vector<string> engineNames;
	if (ext_config.parse("Settings.py") && ext_config.getValue("MC engines", engineNames) && engineNames.size() > 0)
	{
		int nEngine = engineNames.size();
		for (int i = 0; i < nEngine; ++i) pEngineChoice->Append(engineNames[i].c_str());
	}
	else
	{
		pEngineChoice->Append("gZeus");
		pEngineChoice->Append("gMeshMM");
	}
	pEngineChoice->SetSelection(0);

	wxStaticText* jobListLabel = (wxStaticText*)FindWindowById(XRCID("m_staticText_jobList"));
	jobListLabel->Bind(wxEVT_LEFT_DCLICK, &SimuPanel::OnEditJobListPY, this);

	_logText = (wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_log"));
	_peekTimer.SetOwner(GetEventHandler());
	gSimu = this;
}

void SimuPanel::OnSimuStart(wxCommandEvent& event)
{
	if (jobList.size() == 0) return;
	//load the engine dll/so first
	wxChoice* pEngineChoice = (wxChoice*)FindWindowById(XRCID("m_choice_MC_engine"));
	int iMC = pEngineChoice->GetSelection();
	if (iMC == wxNOT_FOUND) iMC = 0; // default engine
	if (iMC != _loadedEngine)
	{
		if (_loadedEngine >= 0) { wxDL.Unload(); } //must unload the engine first
		//now load the engine
		wxString engineName = pEngineChoice->GetString(iMC);
		engineName += wxDynamicLibrary::GetDllExt();
		if (!wxDL.Load(engineName))
		{
			//wxMessageBox("Cannot load the designated engine dll/so. Please check if it exists");
			return;
		}
		//get the pointers of those functions
		startSimulation = (StartSimulation)wxDL.GetSymbol("startSimulation");
		peekDose = (SimpleCall) wxDL.GetSymbol("peekDose");
		stopSimulation = (SimpleCall)wxDL.GetSymbol("stopSimulation");
		setProgressCallBack = (SetProgressCallBack)wxDL.GetSymbol("setProgressCallBack");
		setPeekDoseCallBack = (SetPeekDoseCallBack)wxDL.GetSymbol("setPeekDoseCallBack");
		setJobFinishedCallBack = (SetJobFinishedCallBack)wxDL.GetSymbol("setJobFinishedCallBack");
		setLogCallBack = (SetLogCallBack)wxDL.GetSymbol("setLogCallBack");
		setExitCallBack = (SetLogCallBack)wxDL.GetSymbol("setExitCallBack");
		if (NULL == startSimulation || NULL == peekDose || NULL == stopSimulation ||
			NULL == setProgressCallBack || NULL == setPeekDoseCallBack || NULL == setJobFinishedCallBack || NULL == setLogCallBack)
		{
			wxMessageBox("Cannot load the correct symbol designated engine dll/so.");
			return;
		}

		//set the call-back functions
		setProgressCallBack(progressCB);
		setPeekDoseCallBack(peekDoseCB);
		setJobFinishedCallBack(jobFinishedCB);
		setLogCallBack(LogCB);
		setExitCallBack(exitCB);
	}

	if (_cJob >= jobList.size())
	{
		int ans = wxMessageBox("The job queue has been finished. Yes to recalculate this queue, No to do nothing",
			"Want to recalculate?", wxYES_NO);
		if(ans == wxNO) return; //nothing to do
		_cJob = 0;
	}

	OneJob* job = NULL;
	while (true)
	{
		job = listAt(jobList, _cJob);
		if (!job->skip) break;
		++_cJob;
	}

	//disable the start button
	FindWindowById(XRCID("m_button_simuStart"))->Disable();
	//this function launches another thread to do the real job, and return immediately
	startSimulation(job->configFile.c_str(),configMacro, false);
	_running = true;
	
	if (_autoPeek && _peekWay == 1) // peek dose by time interval
	{
		setTimer();
	}
	showJobStatus();
}

void SimuPanel::OnOneJobFinished(bool abortJob)
{
	_iLastPeek = 0;
	if (abortJob)
	{
		FindWindowById(XRCID("m_button_simuStart"))->Enable();
		_running = false;
		if (_peekWay == 1) _peekTimer.Stop();
		showJobStatus();

		wxCommandEvent event_notify(SIMULATION_FINISHED);
		event_notify.SetInt(1); //abort
		wxPostEvent(ggl->frame, event_notify);
		return;
	}

	if (_cJob == jobList.size() - 1) //finished the last job
	{
		FindWindowById(XRCID("m_button_simuStart"))->Enable();
		++_cJob;
		_running = false;
		if (_peekWay == 1) _peekTimer.Stop();
		showJobStatus();

		wxCommandEvent event_notify(SIMULATION_FINISHED);
		event_notify.SetInt(0); //finished normally
		wxPostEvent(ggl->frame, event_notify);

		return; //nothing to do
	}

	//still has jobs left
	OneJob* job = NULL;
	do
	{
		_cJob += 1;
		job = listAt(jobList, _cJob);
	} while (job->skip);
	showJobStatus();
	//this function launches another thread to do the real job, and return immediately
	startSimulation(job->configFile.c_str(), configMacro, false);
}

void SimuPanel::OnSimuStop(wxCommandEvent& event)
{
	// this function label the simulation thread to exit, but it may takes a while
	if(_running) stopSimulation();
}

void SimuPanel::OnTimer(wxTimerEvent& event)
{
	if (_running) peekDose();
}

void SimuPanel::OnPeekDose(wxCommandEvent& event) //send peed dose command to the dll/so engine
{
	if(_running) peekDose();
}

void SimuPanel::OnUpdateJobList(wxCommandEvent& event)
{
	getJobList();
}

void SimuPanel::OnMoveJobUp(wxCommandEvent& event)
{
	if (_selectedRow < 0) return;
	if (_selectedRow == 0 || _selectedRow == _cJob) return; // cannot move up
	if (_running && _selectedRow - 1 == _cJob) return;
	listSwap(jobList, _selectedRow, _selectedRow - 1);
	_selectedRow--;
	if (_selectedRow < 0) _selectedRow = 0;
	showJobStatus();
}

void SimuPanel::OnMoveJobDown(wxCommandEvent& event)
{
	if (_selectedRow < 0) return;
	if (_selectedRow == jobList.size() - 1 || _selectedRow + 1 == _cJob) return; // cannot move down
	if (_running && _selectedRow == _cJob) return;
	listSwap(jobList, _selectedRow, _selectedRow + 1);
	_selectedRow++;
	if (_selectedRow >= jobList.size()) _selectedRow = jobList.size() - 1;
	showJobStatus();
}

void SimuPanel::OnCheckShowLog(wxCommandEvent& event)
{
	wxCheckBox* checkLog = (wxCheckBox*)FindWindowById(XRCID("m_checkBox_showLog"));
	wxWindow* logText = FindWindowById(XRCID("m_textCtrl_log"));
	if (checkLog->IsChecked()) logText->Show();
	else logText->Hide();
}

void SimuPanel::OnCheckAutoPeek(wxCommandEvent& event)
{
	_autoPeek = event.IsChecked();
	if (_autoPeek) //from manual peek to auto peek
	{
		wxCommandEvent evtx;
		if (_peekWay == 0) OnSetPeekByDivision(evtx);
		else OnSetPeekByTime(evtx);
	}
	else //from auto peek to manual peek
	{
		if (_peekTimer.IsRunning()) _peekTimer.Stop();
	}
	
}

void SimuPanel::OnShowOption(wxGridEvent& event)
{
	int row = event.GetRow();
	int col = event.GetCol();
	_grid->GoToCell(row, col);
	if (row < jobList.size())
	{
		_selectedRow = row;
		wxMenu mnu;
		if(listAt(jobList,row)->skip) mnu.Append(ID_CHANGE_SKIP, "Unskip");
		else mnu.Append(ID_CHANGE_SKIP, "Skip");
		mnu.Append(ID_RECALCULATE, "Recalculate");
		mnu.Append(ID_REMOVE_JOB, "Remove");
		mnu.Append(ID_MOVEUP_JOB, "Move up");
		mnu.Append(ID_MOVEDOWN_JOB, "Move down");
		mnu.Append(ID_OPEN_JOB_IN_FOLDER, "Open in folder");
		mnu.Append(ID_OPEN_JOB_CONFIG_FILE, "Edit config file");
		mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&SimuPanel::OnPopupClick, NULL, this);
		PopupMenu(&mnu);
	}
}

void SimuPanel::OnJobSelected(wxGridEvent& event)
{
	_selectedRow = event.GetRow();
	if (_selectedRow >= jobList.size()) _selectedRow = jobList.size() - 1;
	showJobStatus();
}

void SimuPanel::OnJobDClick(wxGridEvent& event)
{
	int row = event.GetRow();
	int col = event.GetCol();
	if (row < jobList.size())
	{
		OneJob* job = listAt(jobList, row);
		if (col == 0) //open the folder 
		{
			if (_selectedRow == _cJob&&_running) return; // should not apply to the running job
			job = listAt(jobList, _selectedRow);
			if (job->skip) //from "skip" to "unskip"
			{
				job->removeSkipFile(); //first need to delete the txt.skip file
				if (_selectedRow < _cJob && _running) //has been processed, so need to move it to the queue end
				{
					OneJob onejob(*job);
					listRemoveAt(jobList, _selectedRow);
					jobList.push_back(onejob);
					--_cJob; // the current jobIndex is reduced by 1
					_selectedRow = jobList.size() - 1;
				}
				//else don't need to change the sequence

			}
			else // from "unskip" to "skip"
			{
				job->addSkipFile(); //first need to add the txt.skip file
			}
			showJobStatus();
		}
		else if (col == 1) //open the folder 
		{
			wxString p("file:");
			wxLaunchDefaultBrowser(p + job->getJobDir().c_str());
		}
		else if (col == 2) // edit the config file
		{
			wxMimeTypesManager manager;
			wxFileType *filetype = manager.GetFileTypeFromExtension(wxT("py"));
			wxString command = filetype->GetOpenCommand(job->configFile);
			wxExecute(command);
			delete filetype;
		}
		else if (col == 3) //show the phantom
		{
// 			BinaryFile BF;
// 			getBinaryFILE(job->configFile.c_str(), BF);
// 			ggl->showBF(BF);
			ggl->showCustomizedPhantom(job->configFile.c_str());
		}
	}
}

void SimuPanel::OnPopupClick(wxCommandEvent &event)
{
	OneJob* job = NULL;
	switch (event.GetId())
	{
	case ID_CHANGE_SKIP:
		if (_selectedRow == _cJob&&_running) break; // should not apply to the running job
		job = listAt(jobList, _selectedRow);
		if (job->skip) //from "skip" to "unskip"
		{
			job->removeSkipFile(); //first need to delete the txt.skip file
			if (_selectedRow < _cJob && _running) //has been processed, so need to move it to the queue end
			{
				OneJob onejob(*job);
				listRemoveAt(jobList, _selectedRow);
				jobList.push_back(onejob);
				--_cJob; // the current jobIndex is reduced by 1
				_selectedRow = jobList.size() - 1;
			}
			//else don't need to change the sequence
			
		}
		else // from "unskip" to "skip"
		{
			job->addSkipFile(); //first need to add the txt.skip file
		}
		break;

	case ID_RECALCULATE:
	{
		if (_selectedRow == _cJob) break; // should not apply to the running job
		job = listAt(jobList, _selectedRow);
		if (job->skip) job->removeSkipFile();
		if (_selectedRow < _cJob) //has been processed, so need to move it to the queue end
		{
			OneJob onejob(*job);
			listRemoveAt(jobList, _selectedRow);
			jobList.push_back(onejob);
			--_cJob; // the current jobIndex is reduced by 1
		}
	}
		break;

	case ID_REMOVE_JOB:
		if (_selectedRow == _cJob)
		{
			wxMessageBox("Cannot remove the running job. Please stop the simulation first");
			return;
		}
		{
			job = listAt(jobList, _selectedRow);
			string cmd("trash ");
			cmd += job->getJobDir();
			system(cmd.c_str());// delete the folder to recycle bin
			listRemoveAt(jobList, _selectedRow);
			if (_selectedRow < _cJob) --_cJob; //has been processed, so the current jobIndex is reduced by 1
			_selectedRow--;
			if (_selectedRow < 0) _selectedRow = 0;
		}
		break;

	case ID_MOVEUP_JOB:
		if (_selectedRow == 0 || _selectedRow == _cJob) break; // cannot move up
		if (_running && _selectedRow - 1 == _cJob) break;
		listSwap(jobList, _selectedRow, _selectedRow - 1);
		_selectedRow--;
		if (_selectedRow < 0) _selectedRow = 0;
		break;

	case ID_MOVEDOWN_JOB:
		if (_selectedRow == jobList.size() - 1 || _selectedRow + 1 == _cJob) break; // cannot move down
		if (_running && _selectedRow == _cJob) break;
		listSwap(jobList, _selectedRow, _selectedRow + 1);
		_selectedRow++;
		if (_selectedRow >= jobList.size()) _selectedRow = jobList.size() - 1;
		break;

	case ID_OPEN_JOB_IN_FOLDER:
	{
		wxString p("file:");
		job = listAt(jobList, _selectedRow);
		wxLaunchDefaultBrowser(p + job->getJobDir().c_str());
	}
		break;
	case ID_OPEN_JOB_CONFIG_FILE:
	{
		job = listAt(jobList, _selectedRow);
		wxMimeTypesManager manager;
		wxFileType *filetype = manager.GetFileTypeFromExtension(wxT("py"));
		wxString command = filetype->GetOpenCommand(job->configFile);
		wxExecute(command);
		delete filetype;
	}
		break;
	}

	showJobStatus();
}

void SimuPanel::OnEditJobListPY(wxMouseEvent& event)
{
	wxMimeTypesManager manager;
	wxFileType *filetype = manager.GetFileTypeFromExtension(wxT("py"));
	wxString command = filetype->GetOpenCommand(wxT("jobList.py"));
	wxExecute(command);
	delete filetype;
}

void SimuPanel::OnSetPeekByDivision(wxCommandEvent &event)
{
	_peekWay = 0;
	//need to update _NPeekDivision
	wxTextCtrl* pt = (wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_peekDivision"));
	wxString str = pt->GetLineText(0);
	long n;
	if (!str.ToLong(&n) || n <= 0)
	{
		wxNotificationMessage nmsg;
		nmsg.SetMessage("Invalid division number for peeking dose. Will use default 20 instead.");
		nmsg.Show(1);
		pt->SetLabelText("20");
		n = 20;
	}
	_NPeekDivision = n;
	_iLastPeek = 0; //force to peek once immediately
	if (_peekTimer.IsRunning()) //need changing from way 1 to 0
	{
		_peekTimer.Stop();
	}
}

void SimuPanel::OnSetPeekByTime(wxCommandEvent &event)
{
	_peekWay = 1;
	if (!_peekTimer.IsRunning() && _running && _autoPeek) //need to start the timer
	{
		setTimer();
	}
}

bool isJobUnchecked(const OneJob& job)
{
	return !job.checked;
}

bool SimuPanel::getJobList()
{
	ConfigFile jobs;
	if (!jobs.parse("jobList.py"))
	{
		//wxMessageBox("Corrupted jobList.py! Please repair it by reinstalling this software.");
		return false; //wrong config
	}

	//do the micro replacement
	ConfigFile* pcf = jobs.getBlock("ConfigVariables");
	if (pcf)
	{
		pcf->copyMps(configMacro);
		for (unsigned int i = 0; i < configMacro.size(); ++i)
		{
			string d("$");
			configMacro[i].first = d + configMacro[i].first + "$";
		}
		jobs.macroReplace(configMacro);
	}

	//mark all the existing job as unchecked
	for (list<OneJob>::iterator it = jobList.begin(); it != jobList.end(); ++it) it->checked = false;

	fs::directory_iterator end_itr;
	extern wxString configDir;
	for (fs::directory_iterator ifile((configDir + FileSeparator + "BatchMode").ToStdString()); ifile != end_itr; ++ifile)
	{
		if (fs::is_directory(ifile->status())) //check the directory
		{
			string dir = ifile->path().string();
			string result = findFirstFileInDir(dir, ".*\\.py");
			if (result.empty()) continue; //skip this directory
			string strskip = findFirstFileInDir(dir, ".*\\.skip");
			bool skip = false;
			if (!strskip.empty()) skip = true;
			bool found = false;
			//search the existing jobList to see if there's any match;
			for (list<OneJob>::iterator it = jobList.begin(); it != jobList.end(); ++it)
			{
				if (it->configFile.compare(result.c_str()) == 0 && it->skip == skip)
				{
					it->checked = true;
					found = true;
				}
			}
			if (!found) //not found, add to the list end
			{
				OneJob job;
				job.checked = true;
				job.configFile = result;
				job.skip = skip;
				jobList.push_back(job);
			}
		}
	}

	//delete all those unchecked job
	jobList.remove_if(isJobUnchecked);
	return true;
}

void SimuPanel::showJobStatus()
{
	int ncRow = _grid->GetNumberRows();
	int njob = int(jobList.size());
	if (njob > ncRow) _grid->AppendRows(njob - ncRow);
	int ijob = 0;
	for (list<OneJob>::iterator it = jobList.begin(); it != jobList.end(); ++it)
	{
		wxColour sColour = it->skip ? wxColour(128, 128, 128) : wxColour(207, 239, 232);
// 		if (ijob == _cJob)
// 		{
// 			if(_running) sColour = wxColour(0, 255, 0); //indicate the running job
// 			else sColour = wxColour(255, 255, 0); //indicate the job which is ready to run
// 		}
		_grid->SetCellBackgroundColour(ijob, 0, sColour);
		if (ijob == _selectedRow)
		{
			_grid->SetCellBackgroundColour(ijob, 1, wxColour(207, 239, 232));
			_grid->SetCellBackgroundColour(ijob, 2, wxColour(207, 239, 232));
			_grid->SetCellBackgroundColour(ijob, 3, wxColour(207, 239, 232));
			_grid->SetCellBackgroundColour(ijob, 4, wxColour(207, 239, 232));
			_grid->SetCellBackgroundColour(ijob, 5, wxColour(207, 239, 232));
			_grid->SetCellBackgroundColour(ijob, 6, wxColour(207, 239, 232));
		}
		else
		{
			_grid->SetCellBackgroundColour(ijob, 1, wxColour(255, 255, 255));
			_grid->SetCellBackgroundColour(ijob, 2, wxColour(255, 255, 255));
			_grid->SetCellBackgroundColour(ijob, 3, wxColour(255, 255, 255));
			_grid->SetCellBackgroundColour(ijob, 4, wxColour(255, 255, 255));
			_grid->SetCellBackgroundColour(ijob, 5, wxColour(255, 255, 255));
			_grid->SetCellBackgroundColour(ijob, 6, wxColour(255, 255, 255));
		}
		if (ijob == _cJob)
		{
			if (_running) sColour = wxColour(0, 255, 0); //indicate the running job
			else sColour = wxColour(255, 255, 0); //indicate the job which is ready to run
			_grid->SetCellBackgroundColour(ijob, 1, sColour);
		}
		_grid->SetCellValue(ijob, 1, it->getJobName());

		wxString str;
		if (it->percent == 100) _grid->SetCellBackgroundColour(ijob, 2, wxColour(0, 255, 0));
		str.Printf("%.1f", it->percent);
		_grid->SetCellAlignment(ijob, 2, wxALIGN_RIGHT, wxALIGN_BOTTOM);
		_grid->SetCellValue(ijob, 2, str);

		_grid->SetCellAlignment(ijob, 3, wxALIGN_RIGHT, wxALIGN_BOTTOM);
		if(it->timeUsed) _grid->SetCellValue(ijob, 3, time2print(it->timeUsed)); // time used
		_grid->SetCellAlignment(ijob, 4, wxALIGN_RIGHT, wxALIGN_BOTTOM);
		if (it->timeLeft) _grid->SetCellValue(ijob, 4, time2print(it->timeLeft)); // time left
		_grid->SetCellAlignment(ijob, 5, wxALIGN_RIGHT, wxALIGN_BOTTOM);
		if (it->speed)
		{
			str.Printf("%d", int(it->speed));
			_grid->SetCellValue(ijob, 5, str);
		}
		_grid->SetCellAlignment(ijob, 6, wxALIGN_RIGHT, wxALIGN_BOTTOM);
		if (it->uncertainty)
		{
			str.Printf("%.1f", it->uncertainty);
			_grid->SetCellValue(ijob, 6, str);
		}

		++ijob;
	}
	_grid->AutoSizeColumns();
}

//called by call-back function
void SimuPanel::progressUpdate(double percent, double speed, double time, double time_left, double uncertainty)
{
	OneJob* job = listAt(jobList, _cJob);
	job->percent = percent;
	job->speed = speed;
	job->timeUsed = time;
	job->timeLeft = time_left;
	job->uncertainty = uncertainty;

	wxString str;
	//print out those information
	str.Printf("%.1f", percent);
	_grid->SetCellValue(_cJob, 2, str);

	
	_grid->SetCellValue(_cJob, 3, time2print(time)); // time used

	_grid->SetCellValue(_cJob, 4, time2print(time_left)); // time left

	str.Printf("%d", int(speed));
	_grid->SetCellValue(_cJob, 5, str);

	if (uncertainty != 0)
	{
		str.Printf("%.1f", uncertainty);
		_grid->SetCellValue(_cJob, 2, str);
	}
	
	if (_autoPeek && _peekWay == 0)  // peek by division
	{
		int inow = int(percent*(_NPeekDivision + 1) / 100);
		if (inow > _iLastPeek)
		{
			_iLastPeek = inow;
			peekDose();
		}
	}
}

void SimuPanel::setTimer()
{
	wxTextCtrl* pt = (wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_peekInterval"));
	wxString str = pt->GetLineText(0);
	double dt = 5;
	if (!str.ToDouble(&dt) || dt <= 0)
	{
		//wxMessageBox("Invalid time interval for peeking dose. Will use default 5 seconds instead.");
		wxNotificationMessage nmsg;
		nmsg.SetMessage("Invalid time interval for peeking dose. Will use default 5 seconds instead.");
		nmsg.Show(1);
		pt->SetLabelText("5");
		dt = 5;
	}
	//set the timer
	_peekTimer.Start(dt * 1000);
}

void SimuPanel::getBinaryFILE(const char* config, BinaryFile& BF) // get the phantom
{
	//must have I/O data
	int NX, NY, NZ; //voxel number
	double DX, DY, DZ; // voxel size, unit cm
	ArrayMgr<SFloat> ph; //relative density matrix; input

	//optional I/O data
	ArrayMgr<short> matid; //material id matrix; default 0
	int nBatch; // default 0
	double Hist; // default 0
	double Bs; //magnetic field strength; default 0
	double Bx, By, Bz; //unit magnetic field direction; default (0, 0, -1)
	double xo, yo, zo; // position of the cuboid corner in isocenter coordinate system; determine where's the isocenter
	double COPX, COPY, COPZ;//DICOM coordinates info; determine where's the patient coordinate system's origin
	double prescriptionDose; // used in DoseViewer; default maxDose
	int treatmentFraction; // used in DoseViewer; default 1
	double trimDensityThreshold; // default 0
	int rngSeed;

	//auxiliary data, which can be regenerated by the dose file
	double rf; //radius factor= 100/C/B
	int uniform;
	double LX, LY, LZ; // side length Lx=DX*NX, not need to output

	//read NX,NY,NZ; DX,DY,DZ(LX, LY, LZ); density matrix
	ConfigFile gcf(config);
	ConfigFile* cf = gcf.getBlock("PHANTOM");

	bool b_fromFile = true;
	string fname;
	if (!cf->getValue("phantomName", fname)) b_fromFile = false; //no phantom file specified

	cf->getValue("trim dose threshold", trimDensityThreshold);
	if (trimDensityThreshold < 0) trimDensityThreshold = 0;

	double isoPX = 0, isoPY = 0, isoPZ = 0;//iso center position in patient coordinate
	//first try to load from the file
	if (b_fromFile)
	{
		BinaryFile BF;
		if (BF.read(fname.c_str())) return;//try to load infomation from *.phtm file
		else //try to load RED file then
		{
			VIEWRAY_FORMAT VF;
			if (VF.read(fname.c_str()))
			{
				NX = VF.nx;
				NY = VF.ny;
				NZ = VF.nz;
				DX = VF.dx;
				DY = VF.dy;
				DZ = VF.dz;
				COPX = VF.offset_x - DX / 2;
				COPY = VF.offset_y - DY / 2;
				COPZ = VF.offset_z - DZ / 2;
				ph = VF.m;
				matid.resize(NX, NY, NZ);
				matid = 278; //default water
			}
			else
			{
				Log("Warning: Cannot load phantom from file. Will try to use customized phantom in config file");
				b_fromFile = false;
			}
		}
	}

	if (!b_fromFile)//need to get regular phantom from user's definition. Patient's coordinate center is the cuboid center
	{
		vector<double> dc;
		cf->getValue("DX DY DZ", dc);
		DX = dc[0]; DY = dc[1]; DZ = dc[2];

		cf->resetSearchIndex();
		cf->getValue("matrix", dc, true);
		NX = int(dc[0]); NY = int(dc[1]); NZ = int(dc[2]);
		if (NX <= 0 || NY <= 0 || NZ <= 0) 
		{
			wxMessageBox("Incorrect voxel dimensions!"); 
			return;
		}
		//the origin of the patient coordinates lies at the center of the cuboid,
		COPX = -NX*DX / 2;
		COPY = -NY*DY / 2;
		COPZ = -NZ*DZ / 2;
		isoPX = -dc[3];
		isoPY = -dc[4];
		isoPZ = -dc[5];
		xo = COPX - isoPX;
		yo = COPY - isoPY;
		zo = COPZ - isoPZ;
		ph.resize(NX, NY, NZ);
		ph = (dc.size() >= 7 ? SFloat(dc[6]) : SFloat(1)); //1 means use the default density of the corresponding material
		matid.resize(NX, NY, NZ);
		matid = (dc.size() >= 8 ? short(dc[7]) : short(278));

		while (cf->getValue("matrix", dc, true)) //fill inside with different density
		{
			int snx = int(dc[0]), sny = int(dc[1]), snz = int(dc[2]);
			if (snx > NX || sny > NY || snz > NZ)
			{
				wxMessageBox("Incorrect voxel dimensions!");
				return;
			}
			//need to find the origin index of the new cuboid
			int cix = lround((dc[3] - snx*DX*0.5 - xo) / DX);
			int ciy = lround((dc[4] - sny*DY*0.5 - yo) / DY);
			int ciz = lround((dc[5] - snz*DZ*0.5 - zo) / DZ);

			int pix, piy, piz;
			for (int i = 0; i < snx; ++i)
				for (int j = 0; j < sny; ++j)
					for (int k = 0; k < snz; ++k)
					{
						pix = i + cix;
						piy = j + ciy;
						piz = k + ciz;
						if (0 <= pix && pix < NX && 0 <= piy && piy < NY && 0 <= piz && piz < NZ)
						{
							ph.a(pix, piy, piz) = (dc.size() >= 7 ? SFloat(dc[6]) : SFloat(1)); //1 means use the default density of the corresponding material
							matid.a(pix, piy, piz) = (dc.size() >= 8 ? short(dc[7]) : short(278)); //default water
						}
					}
		}

		//may adjust the isoPX, isoPY, isoPZ
		double SSD = 0;
		if (cf->getValue("SSD", SSD)) isoPY = SAD - SSD - NY*DY / 2;
		vector<double> centerXZ;
		if (cf->getValue("centerXZ", centerXZ))
		{
			isoPX = -centerXZ[0];
			isoPZ = -centerXZ[1];
		}

		//different from the definition in phantom.h, we need to convert from relative density to absolute density
		PENELOPE_MAT penMat;
		penMat.load();
		int LEN = NX*NY*NZ;
		for (int i = 0; i < LEN; ++i)
		{
			ph.a(i) = ph.a(i)*penMat.density(matid.a(i));
		}
	}

	//it's correct regardless b_fromFile
	xo = COPX - isoPX;
	yo = COPY - isoPY;
	zo = COPZ - isoPZ;

	LX = NX*DX; LY = NY*DY; LZ = NZ*DZ;

	//overwrite couch density into above phantom
	string couchFile;
	if (cf->getValue("couchConfig", couchFile)) //has the density information
	{
		FILE* fp = fopen(couchFile.c_str(), "rb");
		if (NULL == fp){
			wxMessageBox("Cannot open the couch density file you assigned"); return;
		}
		int CNX, CNY, CNZ;
		float CDX, CDY, CDZ;
		float offset_x, offset_y, offset_z;
		fread(&CNX, sizeof(int), 1, fp);
		fread(&CNY, sizeof(int), 1, fp);
		fread(&CNZ, sizeof(int), 1, fp);
		fread(&CDX, sizeof(float), 1, fp);
		fread(&CDY, sizeof(float), 1, fp);
		fread(&CDZ, sizeof(float), 1, fp);
		if (fabs(CDX - DX) > 1e-5 || fabs(CDY - DY) > 1e-5 || fabs(CDZ - DZ) > 1e-5) {
			wxMessageBox("Inconsistent couch voxel size to the phantom"); return;
		}
		fread(&offset_x, sizeof(float), 1, fp);
		fread(&offset_y, sizeof(float), 1, fp);
		fread(&offset_z, sizeof(float), 1, fp);
		ArrayMgr<SFloat> cMat(CNX, CNY, CNZ);
		float elem = 0;
		for (int k = 0; k < CNZ; ++k)
			for (int j = 0; j < CNY; ++j)
				for (int i = 0; i < CNX; ++i)
				{
					fread(&elem, sizeof(float), 1, fp);
					cMat.a(i, j, k) = SFloat(elem);
				}
		fclose(fp);

		//find the couch's location
		vector<double> pos;
		cf->getValue("couchXYZ", pos);
		if (pos.size() != 3) {
			wxMessageBox("Incorrect couch position!"); return;
		}
		//fit the couch density matrix to above phantom

		//first locate the index of couch's origin
		int cix = lround((pos[0] - CNX*CDX / 2 - xo) / DX);
		int ciy = lround((pos[1] - yo) / DY);
		int ciz = lround((pos[2] - CNZ*CDZ / 2 - zo) / DZ);
		//now the indx(couch)+cix == indx(phantom), and so forth
		//we may need to extend z direction of couch to fit above phantom
		int pix, piy, piz;
		if (ciz > 0) //need a face copy
		{
			for (int piz = 0; piz < ciz; ++piz)
			{
				for (int ix = 0; ix < CNX; ++ix)
					for (int iy = 0; iy < CNY; ++iy)
					{
						pix = ix + cix;
						piy = iy + ciy;
						if (pix >= 0 && pix < NX &&piy >= 0 && piy < NY)
							ph.a(pix, piy, piz) = cMat.a(ix, iy, 0);
					}
			}
		}
		if (ciz + CNZ < NZ) //need another face copy
		{
			for (int piz = ciz + CNZ; piz < NZ; ++piz)
			{
				for (int ix = 0; ix < CNX; ++ix)
					for (int iy = 0; iy < CNY; ++iy)
					{
						pix = ix + cix;
						piy = iy + ciy;
						if (pix >= 0 && pix < NX &&piy >= 0 && piy < NY)
							ph.a(pix, piy, piz) = cMat.a(ix, iy, CNZ - 1);
					}
			}
		}
		//do the major replacement
		for (int ix = 0; ix < CNX; ++ix)
			for (int iy = 0; iy < CNY; ++iy)
				for (int iz = 0; iz < CNZ; ++iz)
				{
					pix = ix + cix;
					piy = iy + ciy;
					piz = iz + ciz;
					if (pix >= 0 && pix < NX &&piy >= 0 && piy < NY && piz >= 0 && piz < NZ)
						ph.a(pix, piy, piz) = cMat.a(ix, iy, iz);
				}
	}

	//may change the couch density
	double couchDensityFactor = 1.0;
	if (cf->getValue("couch density factor", couchDensityFactor))
	{
		int ty = NY - 45;//scan for the position of ViewRa's nominal couch
		for (; ty >= 0; --ty) //the couch's height takes about 48 voxels
		{
			if (fabs(ph.a(NX / 2, ty, NZ / 2) - 1.819) < 1e-4 && fabs(ph.a(NX / 2, ty + 1, NZ / 2) - 0.0778) < 1e-4) break;
		}
		if (ty != -1) //we got the couch top position
		{
			for (int iy = ty; iy <= ty + 48 && iy < NY; ++iy) //the couch's height takes about 48 voxels
			{
				for (int ix = 0; ix < NX; ++ix)
					for (int iz = 0; iz < NZ; ++iz)
					{
						ph.a(ix, iy, iz) *= SFloat(couchDensityFactor);
					}
			}
		}
		else Log("Warning: failed to modify the couch's density. Make sure this couch exist.");
	}
	else Log("Warning: cannot get couch density factor in config file; will not modify the couch density!");

	//read configuration of the magnetic field
	cf->getValue("magnetic field strength", Bs);
	if (Bs <= 0) rf = 0;
	else rf = 100 / (lightSpeed*Bs);

	vector<double> dir;
	cf->getValue("magnetic field direction", dir);
	if (dir.size() == 2) //give two polar angles: theta and phi
	{
		dir[0] *= PI / 180;
		dir[1] *= PI / 180;
		Bx = sin(dir[0])*cos(dir[1]);
		if (fabs(Bx) < 1e-8) Bx = 0;
		By = sin(dir[0])*sin(dir[1]);
		if (fabs(By) < 1e-8) By = 0;
		Bz = cos(dir[0]);
		if (fabs(Bz) < 1e-8) Bz = 0;
	}
	else if (dir.size() == 3) //give unit vector
	{
		if (fabs(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2] - 1.0) > 1e-8) { wxMessageBox("Wrong unit vector!"); return; }
		Bx = dir[0];
		By = dir[1];
		Bz = dir[2];
	}
	else { wxMessageBox("Wrong magnetic field direction parameters!"); return; }

	int LEN = NX*NY*NZ;
	BF.push("NX", &NX, sizeof(int));
	BF.push("NY", &NY, sizeof(int));
	BF.push("NZ", &NZ, sizeof(int));
	BF.push("DX", &DX, sizeof(double));
	BF.push("DY", &DY, sizeof(double));
	BF.push("DZ", &DZ, sizeof(double));
	BF.push("phantom", ph.getP(), sizeof(SFloat)*LEN);
	BF.push("matid", matid.getP(), sizeof(short)*LEN);

	BF.push("Hist", &Hist, sizeof(double));
	BF.push("nBatch", &nBatch, sizeof(int));
	BF.push("Bs", &Bs, sizeof(double));
	BF.push("Bx", &Bx, sizeof(double));
	BF.push("By", &By, sizeof(double));
	BF.push("Bz", &Bz, sizeof(double));
	BF.push("xo", &xo, sizeof(double)); //tell DoseViewer where is the ISO center
	BF.push("yo", &yo, sizeof(double));
	BF.push("zo", &zo, sizeof(double));
	BF.push("COPX", &COPX, sizeof(double));//keep original DICOM coordinate(may load extra DICOM data later)
	BF.push("COPY", &COPY, sizeof(double));
	BF.push("COPZ", &COPZ, sizeof(double));

	BF.push("prescriptionDose", &prescriptionDose, sizeof(double));
	BF.push("treatmentFraction", &treatmentFraction, sizeof(int));
	BF.push("trimDensityThreshold", &trimDensityThreshold, sizeof(double));
	BF.push("rngSeed", &rngSeed, sizeof(int)); // used for continuous simulation

	//optional output; just for viewer, don't read it in function readDose()
	BF.push("uniform", &uniform, sizeof(int));
}