#include "PreCompile.h"
#include "DVFrame.h"
//#define LOAD_COMPILED_RESOURCE 1

#ifdef LOAD_COMPILED_RESOURCE
#include "res.h"
#endif

// Define a new application type
class MyApp : public wxApp
{
public:
	bool OnInit();
};

IMPLEMENT_APP(MyApp)

wxString configDir;

// `Main program' equivalent, creating windows and returning main app frame
bool MyApp::OnInit()
{	
#ifdef WIN32
	SetProcessDPIAware();
#endif
	wxXmlResource::Get()->InitAllHandlers();

	wxString cwd = wxGetCwd();
	configDir = wxPathOnly(wxStandardPaths::Get().GetExecutablePath())+ FileSeparator + MODULENAME;
	wxSetWorkingDirectory(wxPathOnly(wxStandardPaths::Get().GetExecutablePath())+ FileSeparator + MODULENAME);
	
#ifdef LOAD_COMPILED_RESOURCE
	InitXmlResource(); //load resource from compiled file
#else
	if (!wxXmlResource::Get()->Load(_T("form.xrc"))) return false; //load GUI resource from xrc file
#endif

	wxSetWorkingDirectory(cwd);

	DVFrame* frame = new DVFrame((wxWindow*)NULL);
	if (wxApp::argc > 1)
	{
		frame->gl->processCmdLine(wxApp::argc, wxApp::argv);
	}
	else //
	{
		((wxNotebook*)frame->FindWindowById(XRCID("m_notebook_panels")))->SetSelection(4);
	}
	int screenX = wxSystemSettings::GetMetric(wxSYS_SCREEN_X);
	int screenY = wxSystemSettings::GetMetric(wxSYS_SCREEN_Y);
	frame->SetSize(screenX*0.8, screenY*0.9);
	frame->Center();

	frame->Maximize();
	frame->Show(true);
	return true;
}

void exitApp(const char inf[])
{
	wxString s;
	s.Printf(wxT("error: %s\n\n\nApplication will exit...\n"), inf);
	wxMessageBox(s);
	exit(-1);
}

