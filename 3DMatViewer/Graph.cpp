#include "PreCompile.h"
#include "Graph.h"

BEGIN_EVENT_TABLE(FigureFrame, wxFrame)
EVT_MENU(XRCID("m_menuItem_exit"), FigureFrame::OnExit)
EVT_MENU(XRCID("m_menuItem_about"), FigureFrame::OnAbout)
EVT_TOOL(XRCID("m_tool_exit"), FigureFrame::OnExit)
END_EVENT_TABLE()

//----------------------------------class GraphGL-------------------------//
BEGIN_EVENT_TABLE(GraphGL, wxGLCanvas)
//begin mouse event
EVT_MOTION(GraphGL::mouseMoved)
EVT_LEFT_DOWN(GraphGL::mouseDown)
EVT_LEFT_DCLICK(GraphGL::mouseLeftDouble)
EVT_LEFT_UP(GraphGL::mouseReleased)
EVT_RIGHT_DOWN(GraphGL::rightClick)
EVT_LEAVE_WINDOW(GraphGL::mouseLeftWindow)
EVT_ENTER_WINDOW(GraphGL::mouseEnterWindow)
EVT_MOUSEWHEEL(GraphGL::mouseWheelMoved)
//end mouse event
EVT_KEY_DOWN(GraphGL::keyPressed)
EVT_KEY_UP(GraphGL::keyReleased)
EVT_SIZE(GraphGL::resized)
EVT_PAINT(GraphGL::OnPaint)

//begin command event
EVT_MENU(XRCID("m_menuItem_zoom_in"), GraphGL::OnZoomIn)
EVT_MENU(XRCID("m_menuItem_zoom_out"), GraphGL::OnZoomOut)
EVT_TOOL(XRCID("m_tool_save"), GraphGL::OnSave)
EVT_TOOL(XRCID("m_tool_export"), GraphGL::OnExport)
EVT_TOOL(XRCID("m_tool_zoom_in"), GraphGL::OnZoomIn)
EVT_TOOL(XRCID("m_tool_zoom_out"), GraphGL::OnZoomOut)
EVT_TOOL(XRCID("m_tool_resetSize"), GraphGL::OnResetSize)
EVT_TOOL(XRCID("m_tool_left"), GraphGL::OnMoveLeft)
EVT_TOOL(XRCID("m_tool_right"), GraphGL::OnMoveRight)
EVT_TOOL(XRCID("m_tool_up"), GraphGL::OnMoveUp)
EVT_TOOL(XRCID("m_tool_down"), GraphGL::OnMoveDown)
//end command event
END_EVENT_TABLE()

int sample(mglGraph *gr)
{

	mglData x(100);
	mglData y(100);
	for (int i = 0; i < 100; ++i)
	{
		x.a[i] = 0.5*cos(i*PI / 18);
		y.a[i] = sin(i*PI / 18);
	}
	
	//gr->SetOrigin(0, 0, 0);
	//gr->Box();
	//gr->Plot(y);

	//mglData y;  mgls_prepare1d(&y); 
// 	gr->SetOrigin(0, 0, 0);
// 	gr->SubPlot(2, 2, 0, ""); 
	gr->Title("Plot plot (default)");
	gr->SetRanges(0, 10, 0, 4);
	gr->Aspect(NAN, NAN);
	//gr->Box();  
	gr->Plot(x,y);
	gr->Axis();
	gr->Grid();
	

// 	gr->SubPlot(2, 2, 2, "");  gr->Title("'!' style; 'rgb' palette");
// 	gr->Box();  gr->Plot(y, "o!rgb");
// 
// 	gr->SubPlot(2, 2, 3, "");  gr->Title("just markers");
// 	gr->Box();  gr->Plot(y, " +");
// 
// 	gr->SubPlot(2, 2, 1); gr->Title("3d variant");
// 	gr->Rotate(50, 60);  gr->Box();
// 	mglData yc(30), xc(30), z(30);  z.Modify("2*x-1");
// 	yc.Modify("sin(pi*(2*x-1))"); xc.Modify("cos(pi*2*x-pi)");
// 	gr->Plot(xc, yc, z, "rs");

// 	gr->SubPlot(2, 2, 0, "");
// 	gr->Putsw(mglPoint(0, 1), L"Text can be in ASCII and in Unicode");
// 	gr->Puts(mglPoint(0, 0.6), "It can be \\wire{wire}, \\big{big} or #r{colored}");
// 	gr->Puts(mglPoint(0, 0.2), "One can change style in string: "
// 	 	"\\b{bold}, \\i{italic, \\b{both}}");
// 	gr->Puts(mglPoint(0, -0.2), "Easy to \\a{overline} or "
// 	 	"\\u{underline}");
// 	gr->Puts(mglPoint(0, -0.6), "Easy to change indexes ^{up} _{down} @{center}");
// 	gr->Puts(mglPoint(0, -1), "It parse TeX: \\int \\alpha \\cdot "
// 	 	"\\sqrt3{sin(\\pi x)^2 + \\gamma_{i_k}} dx");
// 	 
// 	gr->SubPlot(2, 2, 1, "");
// 	gr->Puts(mglPoint(0, 0.5), "\\sqrt{\\frac{\\alpha^{\\gamma^2}+\\overset 1{\\big\\infty}}{\\sqrt3{2+b}}}", "@", -4);
// 	gr->Puts(mglPoint(0, -0.5), "Text can be printed\non several lines");
// 	 
// 	gr->SubPlot(2, 2, 2, "");
// 	mglData y;  mgls_prepare1d(&y);
// 	gr->Box();  gr->Plot(y.SubData(-1, 0));
// 	gr->Text(y, "This is very very long string drawn along a curve", ":k");
// 	gr->Text(y, "Another string drawn under a curve", "T:r");
// 	 
// 	gr->SubPlot(2, 2, 3, "");
// 	gr->Line(mglPoint(-1, -1), mglPoint(1, -1), "rA");
// 	gr->Puts(mglPoint(0, -1), mglPoint(1, -1), "Horizontal");
// 	gr->Line(mglPoint(-1, -1), mglPoint(1, 1), "rA");
// 	gr->Puts(mglPoint(0, 0), mglPoint(1, 1), "At angle", "@");
// 	gr->Line(mglPoint(-1, -1), mglPoint(-1, 1), "rA");
// 	gr->Puts(mglPoint(-1, 0), mglPoint(-1, 1), "Vertical");
	return 0;
}
inline bool empty(mglData& gdata)
{
	return gdata.GetNN() == 1 && gdata.a[0] == 0;
}

GraphGL::GraphGL(wxWindow* parent) :gl(NULL), gr(NULL),
wxGLCanvas(parent, wxID_ANY, wxGLCanvas_args, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE)
{
	//SetBackgroundStyle(wxBG_STYLE_CUSTOM);
	m_context = new wxGLContext(this);
	// SetCurrent(*m_context);
	this->parent = parent; //need a frame to manage this window

	glRect = parent->GetClientRect();
	glRect.x = glRect.y = 0;
	picRect = glRect;
	theta = 0;
	phi = 270;
	gl = new mglCanvasGL(glRect.width, glRect.height);
	//gl->SetQuality(3);
	//gl->AddLight(0, mglPoint(0, 0, 3), false);
	//gl->InPlot(0, 1, 0, 1, false);
	gr = new mglGraph((mglBase *)gl);

}

GraphGL::~GraphGL() 
{ 
	delete m_context;
	if (gr) delete gr;
	if (gl) delete gl;
}

//Begin event functions
void GraphGL::mouseMoved(wxMouseEvent& event)
{
	if (event.LeftIsDown() && event.ControlDown())
	{
		wxPoint shift = event.GetPosition() - pMouse;
		picRect.x = picRectBackup.x + shift.x;
		picRect.y = picRectBackup.y + shift.y;
		render();
	}
}
void GraphGL::mouseDown(wxMouseEvent& event)
{
	SetFocus(); //to recieve keyboard event
	pMouse = event.GetPosition(); //store mouse position for future use
	picRectBackup = picRect; //backup the picture region for future use
}
void GraphGL::mouseLeftDouble(wxMouseEvent& event)
{

}
void GraphGL::mouseWheelMoved(wxMouseEvent& event)
{
	if (event.ControlDown()) //enlarge or shrink the image
	{
		if (event.GetWheelRotation() > 0) zoom(1.1);
		else zoom(0.9);
		render();
		return;
	}
}
void GraphGL::mouseReleased(wxMouseEvent& event)
{

}
void GraphGL::rightClick(wxMouseEvent& event)
{
	SetFocus(); //to recieve keyboard event
}
void GraphGL::mouseEnterWindow(wxMouseEvent& event)
{
	SetCursor(wxCursor(wxCURSOR_CROSS));
}
void GraphGL::mouseLeftWindow(wxMouseEvent& event)
{

}
void GraphGL::keyPressed(wxKeyEvent& event)
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
	Plot();
	wxStatusBar* sb = (wxStatusBar*)(parent->FindWindowByName("m_statusBar_graph"));
	wxString info;
	info.Printf("(theta, phi)=(%d,%d)", (int)theta, (int)phi);
	sb->SetStatusText(info);
}
void GraphGL::keyReleased(wxKeyEvent& event)
{
}
void GraphGL::resized(wxSizeEvent& event)
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
void GraphGL::OnPaint(wxPaintEvent& event)
{
	render();
}

void GraphGL::OnSave(wxCommandEvent& event)
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
		if (SVG == format) gr->WriteSVG(path.c_str());
		else if (EPS == format) gr->WriteEPS(path.c_str());
		else if (PRC == format) gr->WritePRC(path.c_str());
		else return;
		wxMessageBox("File is saved successfully!");
	}
}
void GraphGL::OnExport(wxCommandEvent& event)//export data to plot better in Origin
{
	wxMessageBox("you should implement it in your derivative class of GraphGL");
}

void GraphGL::OnZoomIn(wxCommandEvent& event)
{
	zoom(1.1);
	render();
}
void GraphGL::OnZoomOut(wxCommandEvent& event)
{
	zoom(0.9);
	render();
}
void GraphGL::OnResetSize(wxCommandEvent& event)
{
	picRect = getBestFitRect(glRect);
	render();
}
void GraphGL::OnMoveLeft(wxCommandEvent& event)
{
	picRect.x -= 10;
	render();
}
void GraphGL::OnMoveRight(wxCommandEvent& event)
{
	picRect.x += 10;
	render();
}
void GraphGL::OnMoveUp(wxCommandEvent& event)
{
	picRect.y -= 10;
	render();
}
void GraphGL::OnMoveDown(wxCommandEvent& event)
{
	picRect.y += 10;
	render();
}
//End event functions

FloatRect GraphGL::getBestFitRect(wxRect& rect)
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

void GraphGL::render()
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
void GraphGL::zoom(double r)
{
	picRect.width = r*picRect.width;
	picRect.height = r*picRect.height;
	double cx = 0.5*glRect.width;
	double cy = 0.5*glRect.height;
	picRect.x = r*(picRect.x - cx) + cx;
	picRect.y = r*(picRect.y - cy) + cy;
}

list<long> FIGURE::figureList;

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


