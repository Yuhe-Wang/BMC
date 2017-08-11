#include "PreCompile.h"
#include "DVFrame.h"

//the default 
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
END_EVENT_TABLE()

DVFrame::DVFrame(wxWindow* parent) :gl(NULL)
{
	wxXmlResource::Get()->LoadFrame(this, NULL, _T("Frame"));
	//add panels to the notebook manager
	wxNotebook* tabs = (wxNotebook*)FindWindowById(XRCID("m_notebook_panels"));
	wxPanel* controls = wxXmlResource::Get()->LoadPanel(tabs, _T("Panel_controls"));
	wxPanel* dose = wxXmlResource::Get()->LoadPanel(tabs, _T("Panel_dose"));

	tabs->AddPage(controls, _T("controls"), true);
	tabs->AddPage(dose, _T("dose"), false);

	wxSplitterWindow* pw = (wxSplitterWindow*)FindWindowById(XRCID("m_splitter"));
	wxWindow* pold = FindWindowById(XRCID("m_panel_gl"));
	gl = new GLView(pw,this);
	pw->ReplaceWindow(pold,gl);
	pw->SetSashGravity(1.0);
	pold->Hide();
	delete pold;

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
		"m_tool_open",	"res/open.bmp",	"open CT /dose/phantom",
		"m_tool_saveSlice", "res/saveSlice.bmp", "",
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
		"m_tool_profile_horizontal", "res/h_ruler.bmp", "",
		"m_tool_profile_vertical", "res/v_ruler.bmp", "",
		"m_tool_profile_tilted", "res/t_ruler.bmp", "",
		"m_tool_exit", "res/exit.bmp", "",
	};
	int NEntry = sizeof(tools) /sizeof(char*)/ NP4Item;

	wxBitmapScale::setDPIScale(this);
	wxToolBarToolBase* toolItem = NULL;
	for (int i = 0; i < NEntry; ++i)
	{
		toolItem = m_toolBar->AddTool(XRCID(tools[NP4Item*i]), wxEmptyString,
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

// 	static wxEvent* lastEvent= NULL;
// 	if (&event == lastEvent) return false; //already processed
// 	if (event.IsCommandEvent() && !event.IsKindOf(CLASSINFO(wxChildFocusEvent)) &&
// 		!event.IsKindOf(CLASSINFO(wxContextMenuEvent)) && gl != NULL) //need to transfer to gl view first
// 	{
// 		lastEvent = &event;
// 		bool success = gl->GetEventHandler()->ProcessEvent(event);
// 		if (!success) success=wxFrame::ProcessEvent(event);
// 		lastEvent = NULL;
// 		return success;
// 	}
// 	else return wxFrame::ProcessEvent(event);
}

void DVFrame::OnExit(wxCommandEvent& event)
{
	Destroy();
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
