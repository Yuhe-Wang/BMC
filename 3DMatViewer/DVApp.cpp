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

// `Main program' equivalent, creating windows and returning main app frame
bool MyApp::OnInit()
{	
#ifdef WIN32
	SetProcessDPIAware();
#endif
	wxXmlResource::Get()->InitAllHandlers();

	wxString cwd = wxGetCwd();
	wxSetWorkingDirectory(wxPathOnly(wxStandardPaths::Get().GetExecutablePath()) + FileSeparator + MODULENAME);
	
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
	
	frame->Maximize();
	wxSplitterWindow* pw = (wxSplitterWindow*)frame->FindWindowById(XRCID("m_splitter"));
	wxSize sz = frame->GetClientSize();
	//pw->SetSashPosition(sz.GetWidth()*0.85);
	frame->Show(true);

	//wxSplitterWindow* pw = (wxSplitterWindow*)(frame->FindWindowById(XRCID("m_splitter")));
	//wxSize s = pw->GetSize();
	//pw->SetSashPosition(s.GetWidth()*0.8);
	if (!wxSetWorkingDirectory(cwd)) exitApp("cannot restore the working directory!");
	return true;
}

void exitApp(const char inf[])
{
	wxString s;
	s.Printf(wxT("error: %s\n\n\nApplication will exit...\n"), inf);
	wxMessageBox(s);
	exit(-1);
}

