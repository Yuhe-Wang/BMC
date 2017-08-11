#include "PreCompile.h"
#include "Graph.h"

list<long> FigureFrame::frameList;

BEGIN_EVENT_TABLE(FigureFrame, wxFrame)
EVT_MENU(XRCID("m_menuItem_exit"), FigureFrame::OnExit)
EVT_MENU(XRCID("m_menuItem_about"), FigureFrame::OnAbout)
EVT_TOOL(XRCID("m_tool_exit"), FigureFrame::OnExit)
END_EVENT_TABLE()

//----------------------------------class FIGURE---------------------------//
BEGIN_EVENT_TABLE(FIGURE, wxGLCanvas)
//begin mouse event
EVT_MOTION(FIGURE::mouseMoved)
EVT_LEFT_DOWN(FIGURE::mouseDown)
EVT_LEFT_DCLICK(FIGURE::mouseLeftDouble)
EVT_LEFT_UP(FIGURE::mouseReleased)
EVT_RIGHT_DOWN(FIGURE::rightClick)
EVT_LEAVE_WINDOW(FIGURE::mouseLeftWindow)
EVT_ENTER_WINDOW(FIGURE::mouseEnterWindow)
EVT_MOUSEWHEEL(FIGURE::mouseWheelMoved)
//end mouse event
EVT_KEY_DOWN(FIGURE::keyPressed)
EVT_KEY_UP(FIGURE::keyReleased)
EVT_SIZE(FIGURE::resized)
EVT_PAINT(FIGURE::OnPaint)

//begin command event
EVT_MENU(XRCID("m_menuItem_zoom_in"), FIGURE::OnZoomIn)
EVT_MENU(XRCID("m_menuItem_zoom_out"), FIGURE::OnZoomOut)
EVT_TOOL(XRCID("m_tool_save"), FIGURE::OnSave)
EVT_TOOL(XRCID("m_tool_export"), FIGURE::OnExport)
EVT_TOOL(XRCID("m_tool_zoom_in"), FIGURE::OnZoomIn)
EVT_TOOL(XRCID("m_tool_zoom_out"), FIGURE::OnZoomOut)
EVT_TOOL(XRCID("m_tool_resetSize"), FIGURE::OnResetSize)
EVT_TOOL(XRCID("m_tool_left"), FIGURE::OnMoveLeft)
EVT_TOOL(XRCID("m_tool_right"), FIGURE::OnMoveRight)
EVT_TOOL(XRCID("m_tool_up"), FIGURE::OnMoveUp)
EVT_TOOL(XRCID("m_tool_down"), FIGURE::OnMoveDown)
//end command event
END_EVENT_TABLE()

int sample(mglGraph *g)
{

	mglData x(100);
	mglData y(100);
	for (int i = 0; i < 100; ++i)
	{
		x.a[i] = 0.5*cos(i*PI / 18);
		y.a[i] = sin(i*PI / 18);
	}
	
	//g->SetOrigin(0, 0, 0);
	//g->Box();
	//g->Plot(y);

	//mglData y;  mgls_prepare1d(&y); 
// 	g->SetOrigin(0, 0, 0);
// 	g->SubPlot(2, 2, 0, ""); 
	g->Title("Plot plot (default)");
	g->SetRanges(0, 10, 0, 4);
	g->Aspect(NAN, NAN);
	//g->Box();  
	g->Plot(x,y);
	g->Axis();
	g->Grid();
	

// 	g->SubPlot(2, 2, 2, "");  g->Title("'!' style; 'rgb' palette");
// 	g->Box();  g->Plot(y, "o!rgb");
// 
// 	g->SubPlot(2, 2, 3, "");  g->Title("just markers");
// 	g->Box();  g->Plot(y, " +");
// 
// 	g->SubPlot(2, 2, 1); g->Title("3d variant");
// 	g->Rotate(50, 60);  g->Box();
// 	mglData yc(30), xc(30), z(30);  z.Modify("2*x-1");
// 	yc.Modify("sin(pi*(2*x-1))"); xc.Modify("cos(pi*2*x-pi)");
// 	g->Plot(xc, yc, z, "rs");

// 	g->SubPlot(2, 2, 0, "");
// 	g->Putsw(mglPoint(0, 1), L"Text can be in ASCII and in Unicode");
// 	g->Puts(mglPoint(0, 0.6), "It can be \\wire{wire}, \\big{big} or #r{colored}");
// 	g->Puts(mglPoint(0, 0.2), "One can change style in string: "
// 	 	"\\b{bold}, \\i{italic, \\b{both}}");
// 	g->Puts(mglPoint(0, -0.2), "Easy to \\a{overline} or "
// 	 	"\\u{underline}");
// 	g->Puts(mglPoint(0, -0.6), "Easy to change indexes ^{up} _{down} @{center}");
// 	g->Puts(mglPoint(0, -1), "It parse TeX: \\int \\alpha \\cdot "
// 	 	"\\sqrt3{sin(\\pi x)^2 + \\gamma_{i_k}} dx");
// 	 
// 	g->SubPlot(2, 2, 1, "");
// 	g->Puts(mglPoint(0, 0.5), "\\sqrt{\\frac{\\alpha^{\\gamma^2}+\\overset 1{\\big\\infty}}{\\sqrt3{2+b}}}", "@", -4);
// 	g->Puts(mglPoint(0, -0.5), "Text can be printed\non several lines");
// 	 
// 	g->SubPlot(2, 2, 2, "");
// 	mglData y;  mgls_prepare1d(&y);
// 	g->Box();  g->Plot(y.SubData(-1, 0));
// 	g->Text(y, "This is very very long string drawn along a curve", ":k");
// 	g->Text(y, "Another string drawn under a curve", "T:r");
// 	 
// 	g->SubPlot(2, 2, 3, "");
// 	g->Line(mglPoint(-1, -1), mglPoint(1, -1), "rA");
// 	g->Puts(mglPoint(0, -1), mglPoint(1, -1), "Horizontal");
// 	g->Line(mglPoint(-1, -1), mglPoint(1, 1), "rA");
// 	g->Puts(mglPoint(0, 0), mglPoint(1, 1), "At angle", "@");
// 	g->Line(mglPoint(-1, -1), mglPoint(-1, 1), "rA");
// 	g->Puts(mglPoint(-1, 0), mglPoint(-1, 1), "Vertical");
	return 0;
}
inline bool empty(mglData& gdata)
{
	return gdata.GetNN() == 1 && gdata.a[0] == 0;
}

FIGURE::FIGURE(wxWindow* parent) :gl(NULL), g(NULL),
wxGLCanvas(parent, wxID_ANY, wxGLCanvas_args, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE)
{
	//SetBackgroundStyle(wxBG_STYLE_CUSTOM);
	m_context = new wxGLContext(this);
	// SetCurrent(*m_context);
	this->parent = parent; //need a frame to manage this window

	glRect = parent->GetClientRect();
	glRect.x = glRect.y = 0;
	picRect = glRect;
	
	gl = new mglCanvasGL(glRect.width, glRect.height);
	//gl->SetQuality(3);
	//gl->AddLight(0, mglPoint(0, 0, 3), false);
	//gl->InPlot(0, 1, 0, 1, false);
	g = new mglGraph((mglBase *)gl);

	theta = 0;
	phi = 270;
	b_holdOn = false;
}

FIGURE::~FIGURE() 
{ 
	delete m_context;
	if (g) delete g;
	if (gl) delete gl;
}

//Begin event functions
void FIGURE::mouseMoved(wxMouseEvent& event)
{
	if (event.LeftIsDown() && event.ControlDown())
	{
		wxPoint shift = event.GetPosition() - pMouse;
		picRect.x = picRectBackup.x + shift.x;
		picRect.y = picRectBackup.y + shift.y;
		renders();
	}
}
void FIGURE::mouseDown(wxMouseEvent& event)
{
	SetCurrent(*m_context);
	SetFocus(); //to receive keyboard event
	pMouse = event.GetPosition(); //store mouse position for future use
	picRectBackup = picRect; //backup the picture region for future use
}
void FIGURE::mouseLeftDouble(wxMouseEvent& event)
{

}
void FIGURE::mouseWheelMoved(wxMouseEvent& event)
{
	if (event.ControlDown()) //enlarge or shrink the image
	{
		if (event.GetWheelRotation() > 0) zoom(1.1);
		else zoom(0.9);
		renders();
		return;
	}
}
void FIGURE::mouseReleased(wxMouseEvent& event)
{

}
void FIGURE::rightClick(wxMouseEvent& event)
{
	SetCurrent(*m_context);
	SetFocus(); //to recieve keyboard event
}
void FIGURE::mouseEnterWindow(wxMouseEvent& event)
{
	SetCursor(wxCursor(wxCURSOR_CROSS));
}
void FIGURE::mouseLeftWindow(wxMouseEvent& event)
{

}
void FIGURE::keyPressed(wxKeyEvent& event)
{
	switch (event.GetKeyCode())
	{
	case 'W':
		theta += 5;
		if (theta > 180) theta = 180;
		break;
	case 'S':
		theta -= 5;
		if (theta < 0) theta = 0;
		break;
	case 'A':
		phi += 5;
		if (phi > 360) phi -= 360;
		break;
	case 'D':
		phi -= 5;
		if (phi < 0) phi += 360;
		break;
	default:
		break;
	}
	renders();
	wxStatusBar* sb = (wxStatusBar*)(parent->FindWindowByName("m_statusBar_graph"));
	wxString info;
	info.Printf("(theta, phi)=(%d,%d)", (int)theta, (int)phi);
	sb->SetStatusText(info);
}
void FIGURE::keyReleased(wxKeyEvent& event)
{
}
void FIGURE::resized(wxSizeEvent& event)
{
	wxRect newGLRect = event.GetSize(); //new size of the house
	//resize picRect if it exist
	double cx = glRect.width*0.5;
	double cy = glRect.height*0.5;
	double new_cx = newGLRect.width*0.5;
	double new_cy = newGLRect.height*0.5;

	FloatRect bf = getBestFitRect(glRect);
	FloatRect newbf = getBestFitRect(newGLRect);
	double sf = newbf.width / bf.width;
	picRect.x = new_cx + sf*(picRect.x - cx);
	picRect.y = new_cy + sf*(picRect.y - cy);
	picRect.width *= sf;
	picRect.height *= sf;

	glRect = newGLRect;//update the window rectangle
}
void FIGURE::OnPaint(wxPaintEvent& event)
{
	SetCurrent(*m_context);
	renders();
}

void FIGURE::OnSave(wxCommandEvent& event)
{
	wxString caption = wxT("Input a file name");
	wxString defaultFilename = wxT("image.svg");
	wxString defaultDir = wxEmptyString;
	enum {SVG, EPS, PRC};
	wxString wildcard = wxT("svg file(*.svg)|*.svg|eps file(*.eps)|*.eps|prc file(*.prc)|*.prc");
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_SAVE);
	if (dialog.ShowModal() == wxID_OK)
	{
		int format = dialog.GetFilterIndex();
		wxString path = dialog.GetPath();
		if (SVG == format) g->WriteSVG(path.c_str());
		else if (EPS == format) g->WriteEPS(path.c_str());
		else if (PRC == format) g->WritePRC(path.c_str());
		else return;
		wxMessageBox("File is saved successfully!");
	}
}
void FIGURE::OnExport(wxCommandEvent& event)//export data to plot better in Origin
{
	wxMessageBox("you should implement it in your derivative class of FIGURE");
}

void FIGURE::OnZoomIn(wxCommandEvent& event)
{
	zoom(1.1);
	renders();
}
void FIGURE::OnZoomOut(wxCommandEvent& event)
{
	zoom(0.9);
	renders();
}
void FIGURE::OnResetSize(wxCommandEvent& event)
{
	picRect = getBestFitRect(glRect);
	renders();
}
void FIGURE::OnMoveLeft(wxCommandEvent& event)
{
	picRect.x -= 10;
	renders();
}
void FIGURE::OnMoveRight(wxCommandEvent& event)
{
	picRect.x += 10;
	renders();
}
void FIGURE::OnMoveUp(wxCommandEvent& event)
{
	picRect.y -= 10;
	renders();
}
void FIGURE::OnMoveDown(wxCommandEvent& event)
{
	picRect.y += 10;
	renders();
}
//End event functions

FloatRect FIGURE::getBestFitRect(wxRect& rect)
{
	FloatRect bf;
	double imageRatio = double(gl->GetHeight()) / double(gl->GetWidth());
	double windowRatio = double(rect.height) / double(rect.width);
	if (windowRatio < imageRatio)
	{
		bf.x = (rect.width - rect.height / imageRatio) / 2;
		bf.y = 0;

		bf.height = rect.height;
		bf.width = rect.height / imageRatio;
	}
	else
	{
		bf.x = 0;
		bf.y = (rect.height - rect.width*imageRatio) / 2;

		bf.height = rect.width*imageRatio;
		bf.width = rect.width;
	}
	return bf;
}

void FIGURE::renders()
{
	if (!IsShown()) return;

	SetCurrent(*m_context); //make sure this window is controlling the openGL context. It's very important
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_POLYGON_SMOOTH);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	glViewport(0, 0, glRect.width, glRect.height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(0, glRect.width, 0, glRect.height, -1e4, 1e4); 
	gl->gl_clf(); // clear the back ground with default white color
	//move to the right origin
	glTranslatef(picRect.x, glRect.height - picRect.y - picRect.height, 0);
	//fit the size of picRect
	float s = float(picRect.width) / float(gl->GetWidth());
	glScalef(s, s, s);

// 	glTranslatef(gl->GetWidth() / 2.0f, gl->GetHeight() / 2.0f, 0);
// 	//set the view angle
// 	double R = 1000;
// 	double ks1 = sin(theta*PI / 180.0), ks2 = sin(phi*PI / 180.0), kc1 = cos(theta*PI / 180.0), kc2 = cos(phi*PI / 180.0);
// 	gluLookAt(R*ks1*kc2, R*ks1*ks2, R*kc1, 0.0f, 0.0f, 0.0f, -kc1*kc2, -kc1*ks2, ks1);
// 	glTranslatef(-gl->GetWidth() / 2.0f, -gl->GetHeight() / 2.0f, 0);
	//gl->Rotate(theta, phi);
	gl->Finish();

	//draw the boundary of this picture, switch to 2D mode
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, glRect.width, glRect.height, 0, -1e4, 1e4);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glPushMatrix();
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glTranslatef(picRect.x, picRect.y, 0.0f);
	glColor3f(0.0f, 1.0f, 1.0f);
	glEnable(GL_LINE_STIPPLE);
	glLineStipple(1, 0x5555);
	glLineWidth(1);
	glBegin(GL_LINE_LOOP);
	glVertex2f(0, 0);
	glVertex2f(picRect.width, 0);
	glVertex2f(picRect.width, picRect.height);
	glVertex2f(0, picRect.height);
	glEnd();
	glDisable(GL_LINE_STIPPLE);
	glPopAttrib();
	glPopMatrix();

	glPopAttrib();
	glFlush();
	SwapBuffers();
}
void FIGURE::zoom(double r)
{
	picRect.width = r*picRect.width;
	picRect.height = r*picRect.height;
	double cx = 0.5*glRect.width;
	double cy = 0.5*glRect.height;
	picRect.x = r*(picRect.x - cx) + cx;
	picRect.y = r*(picRect.y - cy) + cy;
}

mglData MG(ArrayMgr<double>& am)
{
	mglData mdata;
	int size[4];
	int ndim = am.getDim(size);
	if (ndim <= 0 || ndim > 3) return mdata;
	if (ndim < 3) size[2] = 1;
	if (ndim < 2) size[1] = 1;
	mdata.Link(am.getP(), size[0], size[1], size[2]);
	return mdata;
}


