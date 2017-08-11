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

class FigureFrame : public wxFrame //this class manage the class GraphGL or its sub class
{
public:
	static FigureFrame* Frame(wxWindow* parent, bool b_create_new = true, int figureType = FIGURE_2D) // create a new frame
	{
		FigureFrame* ff = NULL;
		if (!b_create_new) // check if last figure is valid
		{
			for (list<long>::iterator it = frameList.begin(); it != frameList.end();)
			{
				long wid = *it;
				ff = (FigureFrame*)parent->FindWindowById(wid);
				if (ff == NULL)
				{
					it = frameList.erase(it);
				}
				else
				{
					if (ff->getFigureType() == figureType) return ff;
					else
					{
						ff = NULL;
						++it;
					}
				}
			}
		}
		if (ff == NULL) // need to create a new window
		{
			ff = new FigureFrame(parent);
			ff->setFigureType(figureType);
			long wid = ff->GetId();
			frameList.push_front(wid);
		}
		return ff;
	}

	bool hasSubWindow()
	{
		if (gg != NULL) return true;
		else return false;
	}
	void setSubWindow(wxWindow* subw)
	{
		if (gg != NULL) gg->Destroy();
		gg = subw;
	}
	wxWindow* getSubWindow(){ return gg; }
	int getFigureType(){ return figureType; }
	void setFigureType(int ft){ figureType = ft; }
	FigureFrame(wxWindow* parent) :gg(NULL), figureType(0)
	{
		wxXmlResource::Get()->LoadFrame(this, parent, _T("Graph")); // This is equivalent to call create(...)
		// LoadFrame() will always assign the same window ID, we need to create another unique one manually
		double DPIScale = parent->GetContentScaleFactor();
		SetClientSize(640 * DPIScale, 480 * DPIScale); //default size
		SetBackgroundStyle(wxBG_STYLE_CUSTOM);

		/**********************create the toolbar manually***********************************/
		wxString cwd = wxGetCwd();
		wxSetWorkingDirectory(wxPathOnly(wxStandardPaths::Get().GetExecutablePath()) + FileSeparator + MODULENAME);
		wxToolBar* m_toolBar = this->CreateToolBar(wxTB_HORIZONTAL, wxID_ANY);
		//wxToolBar* m_toolBar = GetToolBar();
		m_toolBar->SetBackgroundColour(wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOW));
		const int NP4Item = 3; // how many parameters each item will have
		const char* tools[] =
		{
			//name, resource, hint		
			"m_tool_save", "res/save.bmp", "",
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
			"m_tool_drawBox", "res/drawBox.bmp", "",
			"m_tool_profile_horizontal", "res/h_ruler.bmp", "",
			"m_tool_profile_vertical", "res/v_ruler.bmp", "",
			"m_tool_profile_tilted", "res/t_ruler.bmp", "",
			"m_tool_exit", "res/exit.bmp", "",
		};
		int NEntry = sizeof(tools) / sizeof(char*) / NP4Item;

		wxBitmapScale::setDPIScale(this);
		for (int i = 0; i < NEntry; ++i)
		{
			if (strcmp(tools[NP4Item*i], "") == 0) m_toolBar->AddSeparator();
			else m_toolBar->AddTool(XRCID(tools[NP4Item*i]), wxEmptyString,
				wxBitmapScale(tools[NP4Item*i + 1]), wxNullBitmap, wxITEM_NORMAL, tools[NP4Item*i + 1], wxEmptyString, NULL);
		}
		m_toolBar->Realize();
		wxSetWorkingDirectory(cwd);
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
	wxWindow *gg; //the sub window it's managing
	int figureType;
	static list<long> frameList;
	DECLARE_EVENT_TABLE()
};

//this is the base class of drawing various graphs. It provides basic zooming, translating, etc functions
class FIGURE : public wxGLCanvas
{
public:
	static FIGURE& Figure(wxWindow* parent, bool b_create_new = true) // create a new figure
	{
		FigureFrame* ff = FigureFrame::Frame(parent, b_create_new,FIGURE_2D);
		FIGURE* gw = (FIGURE*)(ff->getSubWindow());
		if (gw == NULL)
		{
			gw = new FIGURE(ff);
			ff->setSubWindow(gw);
		}
		// gw->bindGLContext(); // it may not create a new window, so we need to rebind the context for existed window
		gw->hold(false);
		ff->Show(true);//show the figure at this moment will reduce the blank window flashing
		return *gw;
	}

	FIGURE(wxWindow* parent);
	virtual ~FIGURE();

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

	void renders();
	void bindGLContext()
	{
		SetCurrent(*m_context); //make sure this window is controlling the openGL context. It's very important
	}

	void Plot(const mglDataA &x, const mglDataA &y, const char *pen = "", const char *opt = "")
	{
		gl->Clf();
		gl->InPlot(0, 1, 0, 1, false);
		_xmax = x.Maximal();
		_xmin = x.Minimal();
		_ymax = y.Maximal();
		_ymin = y.Minimal();
		
		double dx = _xmax - _xmin;
		double dy = _ymax - _ymin;
		double kx = 0.05, ky = 0.05;

		g->SetRanges(_xmin - kx*dx, _xmax + kx*dx, _ymin - ky*dy, _ymax + ky*dy);
		
		g->Plot(x, y, pen, opt);
		
		SetFocus();// show top
		renders();
	}
	void Plot(const mglDataA &x, const mglDataA &y, const mglDataA &y2, const char *pen1 = "", const char *pen2 = "")
	{

		gl->Clf();
		gl->InPlot(0, 1, 0, 1, false);
		_xmax = x.Maximal();
		_xmin = x.Minimal();
		_ymax = y.Maximal();
		_ymin = y.Minimal();

		double ymax = y2.Maximal();
		double ymin = y2.Minimal();
		_ymax = max(_ymax, ymax);
		_ymin = min(_ymin, ymin);

		double dx = _xmax - _xmin;
		double dy = _ymax - _ymin;
		double kx = 0.05, ky = 0.05;

		g->SetRanges(_xmin - kx*dx, _xmax + kx*dx, _ymin - ky*dy, _ymax + ky*dy);

		g->Plot(x, y, pen1);
		g->Plot(x, y2, pen2);

		SetFocus();// show top
		renders();
	}

	void hold(bool on = true)
	{
		b_holdOn = on;
	}
public: // public data members
	mglCanvasGL* gl;
	mglGraph* g; // pointer we use to send plot command

protected:
	wxGLContext* m_context;
	wxWindow* parent; // it should be the frame
	
	wxPoint pMouse; //current clicked position in y-down coordinate system
	FloatRect picRect;//the actual drawing region in pixel unit
	wxRect glRect;//previous window rectangle in pixel
	FloatRect picRectBackup;//backup current picRect when you click the window

	double phi, theta; //view angle

	bool b_holdOn; // default true
	double _xmin, _xmax;
	double _ymin, _ymax;
	double _zmin, _zmax;

protected: //functions
	FloatRect getBestFitRect(wxRect& rect);
	void zoom(double r);

	DECLARE_EVENT_TABLE()
};
#endif


