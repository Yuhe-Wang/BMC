#ifndef _GRAPHGL_H_
#define _GRAPHGL_H_

#include "PreCompile.h"
#include "View3D.h"
//shut up the warning caused by the third party library MATHGL
#pragma warning(push)
#pragma warning(disable: 4800 190 5)
#include "mgl2/mgl.h"
#include "mgl2/opengl.h"
#pragma warning(pop)

mglData MG(ArrayMgr<double>& am);

//this is the base class of drawing various graphs. It provides basic zooming, translating, etc functions
class GraphGL : public wxGLCanvas
{
public:
	GraphGL(wxWindow* parent);
	virtual ~GraphGL();

	//event functions
	void mouseMoved(wxMouseEvent& event);
	void mouseDown(wxMouseEvent& event);
	void mouseLeftDouble(wxMouseEvent& event);
	void mouseWheelMoved(wxMouseEvent& event);
	void mouseReleased(wxMouseEvent& event);
	virtual void rightClick(wxMouseEvent& event);
	void mouseEnterWindow(wxMouseEvent& event);
	void mouseLeftWindow(wxMouseEvent& event);
	void keyPressed(wxKeyEvent& event);
	void keyReleased(wxKeyEvent& event);
	void resized(wxSizeEvent& event);
	void OnPaint(wxPaintEvent& event);

	void OnSave(wxCommandEvent& event);
	virtual void OnExport(wxCommandEvent& event);
	void OnNormalize(wxCommandEvent& event);
	void OnZoomIn(wxCommandEvent& event);
	void OnZoomOut(wxCommandEvent& event);
	void OnResetSize(wxCommandEvent& event);
	void OnMoveLeft(wxCommandEvent& event);
	void OnMoveRight(wxCommandEvent& event);
	void OnMoveUp(wxCommandEvent& event);
	void OnMoveDown(wxCommandEvent& event);

	void render();
	void bindContext()
	{
		SetCurrent(*m_context); //make sure this window is controlling the openGL context. It's very important
	}
public: // public data members
	mglCanvasGL* gl;
	mglGraph* gr; // pointer we use to send plot command
protected:
	wxGLContext* m_context;
	wxWindow* parent; // it should be the frame
	
	wxPoint pMouse; //current clicked position in y-down coordinate system
	FloatRect picRect;//the actual drawing region in pixel unit
	wxRect glRect;//previous window rectangle in pixel
	FloatRect picRectBackup;//backup current picRect when you click the window

	double phi, theta; //view angle
protected: //functions
	FloatRect getBestFitRect(wxRect& rect);
	
	void zoom(double r);
	virtual void Plot(){};

	DECLARE_EVENT_TABLE()
};

class FigureFrame : public wxFrame //this class manage the class GraphGL or its sub class
{
public:
	FigureFrame(wxWindow* parent):gg(NULL)
	{
		wxXmlResource::Get()->LoadFrame(this, parent, _T("Graph")); // This is equivalent to call create(...)
		double DPIScale = parent->GetContentScaleFactor();
		SetClientSize(640*DPIScale, 480*DPIScale); //default size
		SetBackgroundStyle(wxBG_STYLE_CUSTOM);
		//need to create the tool bar manually to fit to right DPI scaling
		wxString cwd = wxGetCwd();
		wxSetWorkingDirectory(wxPathOnly(wxStandardPaths::Get().GetExecutablePath()) + FileSeparator + MODULENAME);
		wxToolBar* m_toolBar = this->CreateToolBar(wxTB_HORIZONTAL, wxID_ANY);
		m_toolBar->SetBackgroundColour(wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOW));
		const int NP4Item = 3; // how many parameters each item will have
		const char* tools[] =
		{
			//name, resource, hint
			"m_tool_save", "res/save.bmp", "save",
			"m_tool_export", "res/origin.bmp", "",
			"m_tool_zoom_in", "res/zoomin.bmp", "",
			"m_tool_zoom_out", "res/zoomout.bmp", "",
			"m_tool_resetSize", "res/resetSize.bmp", "",
			"m_tool_left", "res/left.bmp", "",
			"m_tool_right", "res/right.bmp", "",
			"m_tool_up", "res/up.bmp", "",
			"m_tool_down", "res/down.bmp", "",
			"m_tool_rotate_left", "res/rotateL.bmp", "",
			"m_tool_rotate_right", "res/rotateR.bmp", "",
			"m_tool_exit", "res/exit.bmp", "",
		};
		int NEntry = sizeof(tools) / sizeof(char*) / NP4Item;

		wxBitmapScale::setDPIScale(this);
		wxToolBarToolBase* toolItem = NULL;
		for (int i = 0; i < NEntry; ++i)
		{
			toolItem = m_toolBar->AddTool(XRCID(tools[NP4Item*i]), wxEmptyString,
				wxBitmapScale(tools[NP4Item*i + 1]), wxNullBitmap, wxITEM_NORMAL, tools[NP4Item*i + 1], wxEmptyString, NULL);
		}
		m_toolBar->Realize();
		wxSetWorkingDirectory(cwd);

		SetFocus();
		//Raise();
		//create the wxGLCanvas window
		gg = new GraphGL(this);
		Show(true);
	}
	void ShowTop() //make sure this window pops up to front
	{
		SetFocus();
	}
	bool ProcessEvent(wxEvent& event)
	{
		if (event.IsCommandEvent() && !event.IsKindOf(CLASSINFO(wxChildFocusEvent)) &&
			!event.IsKindOf(CLASSINFO(wxContextMenuEvent)) && gg != NULL)
		{
			event.StopPropagation();//not going to pass to parent again
			bool success = wxFrame::ProcessEvent(event);//try to process by this class first
			if (!success) success = gg->GetEventHandler()->ProcessEvent(event);//pass the event to wxCanvasGL
			return success;
		}
		else return wxFrame::ProcessEvent(event);
	}

// 	bool Destroy() // may not need it anymore.
// 	{
// 		//notify the parent that this sub window is no longer valid
// 		wxCommandEvent event_notify(VIEW3D_NOTIFY);
// 		event_notify.SetInt(VIEW3D_GRAPH_DESTROYED);
// 		event_notify.SetEventObject((wxObject*)this);
// 		wxPostEvent(parent, event_notify);
// 		return wxFrame::Destroy();
// 	}
	void OnExit(wxCommandEvent& event) // may not need it anymore
	{
		Destroy();
	}
	void OnAbout(wxCommandEvent& event)
	{
		wxAboutDialogInfo aboutInfo;
		aboutInfo.SetName("Graph Plotting Module of Dose Viewer");
		aboutInfo.SetVersion("1.0");
		aboutInfo.SetDescription(_("This module generates 2D or 3D figures for dose Viewer\nIt's based on the C++ plotting library MathGL\nThanks to the author Alexey Balakin!"));
		aboutInfo.SetCopyright("(C) 2014 WUSTL Medical School");
		aboutInfo.SetWebSite("https://www.linkedin.com/profile/view?id=308598340&trk=nav_responsive_tab_profile");
		aboutInfo.AddDeveloper("Yuhe Wang,   yuhewang.ustc@gmail.com");
		wxAboutBox(aboutInfo);
	}

public:
	GraphGL *gg; //the openGL window it's managing
	DECLARE_EVENT_TABLE()
};

class FIGURE // it is designed like MATLAB's plot system 
{
public:
	FIGURE(FigureFrame* ff) // just wrap the frame
	{
		this->ff = ff;
		b_holdOn = false;
		g = ff->gg->gr;
		//ff->gg->bindContext(); // make sure all commands is issued to the correct context
	}
	static FIGURE Figure(wxWindow* parent, bool b_create_new = true) // create a new figure
	{
		FigureFrame* ff = NULL;
		if (!b_create_new) // check if last figure is valid
		{
			while (figureList.size() > 0)
			{
				long wid = figureList.back();
				ff = (FigureFrame*)parent->FindWindowById(wid);
				if (ff != NULL) break;
				else figureList.remove(wid);
			}
		}
		if (ff == NULL) // need to create a new window
		{
			ff = new FigureFrame(parent);
			figureList.push_back(ff->GetId());
		}
		return FIGURE(ff);
	}

	void Plot(const mglDataA &x, const mglDataA &y, const char *pen = "", const char *opt = "")
	{
		if (!b_holdOn)
		{
			ff->gg->gl->Clf();
			ff->gg->gl->InPlot(0, 1, 0, 1, false);
			_xmax = x.Maximal();
			_xmin = x.Minimal();
			_ymax = y.Maximal();
			_ymin = y.Minimal();
		}
		else //need to calculate the over all range
		{
			double xmax = x.Maximal();
			double xmin = x.Minimal();
			double ymax = y.Maximal();
			double ymin = y.Minimal();
			_xmax = max(_xmax, xmax);
			_xmin = min(_xmin, xmin);
			_ymax = max(_ymax, ymax);
			_ymin = min(_ymin, ymin);
		}
		
		double dx = _xmax - _xmin;
		double dy = _ymax - _ymin;
		double kx = 0.05, ky = 0.05;
		
		g->SetRanges(_xmin - kx*dx, _xmax + kx*dx, _ymin - ky*dy, _ymax + ky*dy);
		g->Box();
		g->Plot(x, y, pen, opt);
		g->Axis();
		g->Grid("xy", "Y;");
		ff->ShowTop();
		render();
	}
	void hold(bool on = true)
	{
		b_holdOn = on;
	}
	void render()
	{
		ff->gg->render();
	}
public:
	mglGraph* g; // pointer we use to send other plot command
private:
	static list<long> figureList;
	FigureFrame* ff;
	bool b_holdOn; // default true

	double _xmin, _xmax;
	double _ymin, _ymax;
	double _zmin, _zmax;
};
#endif


