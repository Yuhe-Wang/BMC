#include "PreCompile.h"
#include "DVFrame.h"

float pixelSize()
{
	static float pixsz = 0;
	if (pixsz != 0) return pixsz;
	int wmm = 0;
	wxDisplaySizeMM(&wmm, NULL);
	int wpix = 0;
	wxDisplaySize(&wpix, NULL);
	pixsz = (float(wmm)) / wpix;
	return pixsz;
}
double wxBitmapScale::DPI_Scale = 1;

BEGIN_EVENT_TABLE(DVFrame, wxFrame)
EVT_MENU(XRCID("m_menuItem_exit"), DVFrame::OnExit)
EVT_MENU(XRCID("m_menuItem_about"), DVFrame::OnAbout)
EVT_TOOL(XRCID("m_tool_exit"), DVFrame::OnExit)
EVT_CLOSE(DVFrame::OnClose)
EVT_COMMAND(wxID_ANY, SIMULATION_FINISHED, DVFrame::OnSimuFinished)
END_EVENT_TABLE()

DVFrame::DVFrame(wxWindow* parent) :gl(NULL), pSimu(NULL), bMarkClose(false)
{
	wxXmlResource::Get()->LoadFrame(this, NULL, _T("Frame"));
	//add panels to the notebook manager
	wxNotebook* tabs = (wxNotebook*)FindWindowById(XRCID("m_notebook_panels"));
// 	wxPanel* controls = wxXmlResource::Get()->LoadPanel(tabs, _T("Panel_controls"));
// 	wxPanel* dose = wxXmlResource::Get()->LoadPanel(tabs, _T("Panel_dose"));
// 	wxPanel* contours = wxXmlResource::Get()->LoadPanel(tabs, _T("Panel_contours"));
	wxBitmapScale::setDPIScale(this);

	const char* panelName[] = { "Panel_controls", "Panel_dose", "Panel_contours", "Panel_external" };
	const char* showName[] = { "controls", "dose", "contours", "external" };
	int NPanel = sizeof(panelName) / sizeof(char*);
	for (int i = 0; i < NPanel; ++i)
	{
		wxScrolledWindow* sw = new wxScrolledWindow(tabs);
		wxBoxSizer* bs = new wxBoxSizer(wxHORIZONTAL);
		wxPanel* pl = wxXmlResource::Get()->LoadPanel(sw, panelName[i]);
		bs->Add(pl, 0, wxEXPAND, 5);
		sw->SetSizer(bs);
		// sw->Hide();
		sw->SetScrollbars(0, 20, 0, 50);
		if (i == 0) tabs->AddPage(sw, showName[i], true);
		else tabs->AddPage(sw, showName[i], false);
	}

	wxScrolledWindow* sw = new wxScrolledWindow(tabs);
	wxBoxSizer* bs = new wxBoxSizer(wxHORIZONTAL);
	pSimu = new SimuPanel(sw);
	bs->Add(pSimu, 0, wxEXPAND, 5);
	sw->SetSizer(bs);
	sw->SetScrollbars(0, 20, 0, 50);
	tabs->AddPage(sw, "simu", false);

// 	wxSize sz = tabs->GetSize();
// 	sz.y += 50 * wxBitmapScale::getDPIScale();
// 	tabs->SetSize(sz);
	
	wxSplitterWindow* pw = (wxSplitterWindow*)FindWindowById(XRCID("m_splitter"));
	wxWindow* pold = FindWindowById(XRCID("m_panel_gl"));
	gl = new GLView(pw,this);
	pw->ReplaceWindow(pold,gl);
	pw->SetSashGravity(1.0);
	pold->Hide();
	delete pold;

	int spos = pw->GetSashPosition();
	spos -= 24 * wxBitmapScale::getDPIScale();// the width of scrollbar is 20 and I set the distance 24 to make it look better
	pw->SetSashPosition(spos);
	/**********************create the tools manually***********************************/
	wxString cwd = wxGetCwd();
	wxSetWorkingDirectory(wxPathOnly(wxStandardPaths::Get().GetExecutablePath()) + FileSeparator + MODULENAME);

	wxToolBar* m_toolBar = this->CreateToolBar(wxTB_HORIZONTAL, wxID_ANY);
	//wxToolBar* m_toolBar = GetToolBar();
	m_toolBar->SetBackgroundColour(wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOW));
	const int NP4Item = 3; // how many parameters each item will have
	const char* tools[] =
	{
		//name, resource, hint
		"m_tool_open", "res/open.bmp", "open CT /dose/phantom",
		"m_tool_openDir", "res/openDir.bmp", "",
		"", "", "", //separator
		"m_tool_openContours","res/contours.bmp", "",
		"m_tool_openRED", "res/RED.bmp", "",
		"m_tool_beams", "res/beams.bmp", "",
		"m_tool_importDose", "res/importDose.bmp", "",
		"", "", "", //separator
		"m_tool_save", "res/save.bmp", "",
		"m_tool_saveSlice", "res/saveSlice.bmp", "",
		"m_tool_adjustDose", "res/adjust.bmp", "",
		"m_tool_pdd", "res/pdd.bmp", "",
		"m_tool_doseProfile", "res/doseProfile.bmp", "",
		"m_tool_zoom_in", "res/zoomin.bmp", "",
		"m_tool_zoom_out", "res/zoomout.bmp", "",
		"m_tool_resetSize", "res/resetSize.bmp", "",
		"m_tool_left", "res/left.bmp", "",
		"m_tool_right", "res/right.bmp", "",
		"m_tool_up", "res/up.bmp", "",
		"m_tool_down", "res/down.bmp", "",
		"m_tool_rotate_left", "res/rotateL.bmp", "",
		"m_tool_rotate_right", "res/rotateR.bmp", "",
		"m_tool_edit_data", "res/editData.bmp", "",
		"m_tool_drawBox", "res/drawBox.bmp", "",
		"m_tool_profile_horizontal", "res/h_ruler.bmp", "",
		"m_tool_profile_vertical", "res/v_ruler.bmp", "",
		"m_tool_profile_tilted", "res/t_ruler.bmp", "",
		"m_tool_exit", "res/exit.bmp", "",
	};
	int NEntry = sizeof(tools) / sizeof(char*) / NP4Item;

	
	for (int i = 0; i < NEntry; ++i)
	{
		if (strcmp(tools[NP4Item*i], "") == 0) m_toolBar->AddSeparator();
		else m_toolBar->AddTool(XRCID(tools[NP4Item*i]), wxEmptyString,
			wxBitmapScale(tools[NP4Item*i + 1]), wxNullBitmap, wxITEM_NORMAL, tools[NP4Item*i + 1], wxEmptyString, NULL);
	}
	m_toolBar->Realize();
	wxSetWorkingDirectory(cwd);
}

bool DVFrame::ProcessEvent(wxEvent& event)
{
	bool success = wxFrame::ProcessEvent(event);//try to process by this class first	

	if (!success && event.IsCommandEvent() && !event.IsKindOf(CLASSINFO(wxChildFocusEvent)) &&
		!event.IsKindOf(CLASSINFO(wxContextMenuEvent)) && gl != NULL)
	{
		event.StopPropagation();//not going to pass to parent again
		success = gl->GetEventHandler()->ProcessEvent(event);//pass the event to wxCanvasGL
	}
	return success;
}

void DVFrame::OnExit(wxCommandEvent& event)
{
	Close();
}

void DVFrame::OnAbout(wxCommandEvent& event)
{
	wxAboutDialogInfo aboutInfo;
	aboutInfo.SetName("Dose Viewer");
	aboutInfo.SetVersion("1.0");
	aboutInfo.SetDescription(_("This app generates phantom file for Monte Carlo Radiation dose calculation\nfrom CT and view the phantom plus the dose distribution in 3 dimensions!"));
	aboutInfo.SetCopyright("(C) 2014 WUSTL Medical School");
	aboutInfo.SetWebSite("https://www.linkedin.com/profile/view?id=308598340&trk=nav_responsive_tab_profile");
	aboutInfo.AddDeveloper("Yuhe Wang,   yuhewang.ustc@gmail.com");
	wxAboutBox(aboutInfo);
}

void DVFrame::OnClose(wxCloseEvent& event)
{
	if(!pSimu->isRunning()) Destroy();
	else
	{
		event.Veto(); //don't destroy it now
		bMarkClose = true;
		wxCommandEvent evtx;
		pSimu->OnSimuStop(evtx);
	}
}

void DVFrame::OnSimuFinished(wxCommandEvent& event)
{
	if(bMarkClose) Destroy();
}
