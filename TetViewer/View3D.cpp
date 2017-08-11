#include "PreCompile.h"
#include "View3D.h"
#include "Graph.h"
#include <unordered_map>

wxDEFINE_EVENT(VIEW3D_NOTIFY, wxCommandEvent); //define customized event type

struct FaceKey // the key value for hash table
{
	FaceKey(){};
	FaceKey(int i0, int i1, int i2){ iv0 = i0; iv1 = i1; iv2 = i2; sort(); }
	void set(int i0, int i1, int i2){ iv0 = i0; iv1 = i1; iv2 = i2; sort(); }
	int iv0, iv1, iv2;
	void sort()
	{
		if (iv0 > iv1) swapData<int>(iv0, iv1);
		if (iv1 > iv2) swapData<int>(iv1, iv2);
		if (iv0 > iv1) swapData<int>(iv0, iv1);
	}
	string operator () () const
	{
		char ch[48];
		sprintf(ch, "%d,%d,%d", iv0, iv1, iv2);
		return string(ch);
	}
	bool operator < (const FaceKey& rk) const
	{
		if (iv0 < rk.iv0) return true;
		else if (iv0 == rk.iv0)
		{
			if (iv1 < rk.iv1) return true;
			else if (iv1 == rk.iv1)
			{
				if (iv2 < rk.iv2) return true;
				else return false;
			}
			else return false;
		}
		else return false;
	}
};

struct FaceValue
{
	FaceValue() :idt0(-1), idt1(-1){};
	void set(int idt, int idf)
	{
		if (idt0 == -1)
		{
			idt0 = idt;
			idf0 = idf;
			return;
		}
		else if (idt1 == -1)
		{
			idt1 = idt;
			idf1 = idf;
			return;
		}
		else
		{
			printf("error: this face belongs to more than two tetrahedrals.\n");
			getchar();
			exit(-1);
		}
	}
	int idt0;
	int idf0; // which face of this tetrahedral

	int idt1;
	int idf1; // which face of this tetrahedral
};

const int TempColor[90][3] =
{
	5, 0, 255,
	4, 0, 255,
	3, 0, 255,
	2, 0, 255,
	1, 0, 255,
	0, 0, 255,
	0, 2, 255,
	0, 18, 255,
	0, 34, 255,
	0, 50, 255,
	0, 68, 255,
	0, 84, 255,
	0, 100, 255,
	0, 116, 255,
	0, 132, 255,
	0, 148, 255,
	0, 164, 255,
	0, 180, 255,
	0, 196, 255,
	0, 212, 255,
	0, 228, 255,
	0, 255, 244,
	0, 255, 208,
	0, 255, 168,
	0, 255, 131,
	0, 255, 92,
	0, 255, 54,
	0, 255, 16,
	23, 255, 0,
	62, 255, 0,
	101, 255, 0,
	138, 255, 0,
	176, 255, 0,
	215, 255, 0,
	253, 255, 0,
	255, 250, 0,
	255, 240, 0,
	255, 230, 0,
	255, 220, 0,
	255, 210, 0,
	255, 200, 0,
	255, 190, 0,
	255, 180, 0,
	255, 170, 0,
	255, 160, 0,
	255, 150, 0,
	255, 140, 0,
	255, 130, 0,
	255, 120, 0,
	255, 110, 0,
	255, 100, 0,
	255, 90, 0,
	255, 80, 0,
	255, 70, 0,
	255, 60, 0,
	255, 50, 0,
	255, 40, 0,
	255, 30, 0,
	255, 20, 0,
	255, 10, 0,
	255, 0, 0,
	255, 0, 16,
	255, 0, 32,
	255, 0, 48,
	255, 0, 64,
	255, 0, 80,
	255, 0, 96,
	255, 0, 112,
	255, 0, 128,
	255, 0, 144,
	255, 0, 160,
	255, 0, 176,
	255, 0, 192,
	255, 0, 208,
	255, 0, 224,
	255, 0, 240,
	255, 1, 240,
	255, 2, 240,
	255, 3, 240,
	255, 4, 240,
	255, 5, 240,
	255, 6, 240,
	255, 7, 240,
	255, 8, 240,
	255, 9, 240,
	255, 10, 240,
	255, 11, 240,
	255, 12, 240,
	255, 13, 240,
	255, 14, 240,
};

const float MatColor[7][3] =
{
	0, 0.4470, 0.7410,
	0.8500, 0.3250, 0.0980,
	0.9290, 0.6940, 0.1250,
	0.4940, 0.1840, 0.5560,
	0.4660, 0.6740, 0.1880,
	0.3010, 0.7450, 0.9330,
	0.6350, 0.0780, 0.1840
};

void colorGradient(float& r, float& g, float&b, float kc, int scheme)
{
	if (kc < 0) return;
	if (0 == scheme)
	{
		if (kc > 1) kc = 1;
		r = min(max(1.5 - 4 * fabs(kc - 0.75), 0.0), 1.0);
		g = min(max(1.5 - 4 * fabs(kc - 0.5), 0.0), 1.0);
		b = min(max(1.5 - 4 * fabs(kc - 0.25), 0.0), 1.0);
	}
	else if (1 == scheme)
	{
		float x = kc * 89;//get a number between 0~89
		int ix = (int)x;
		if (ix >= 89)
		{
			r = TempColor[89][0] / 255.0f;
			g = TempColor[89][1] / 255.0f;
			b = TempColor[89][2] / 255.0f;
			return;
		}
		x -= ix; //the decimal part
		r = (TempColor[ix][0] * (1 - x) + TempColor[ix + 1][0] * x) / 255.0f;
		g = (TempColor[ix][1] * (1 - x) + TempColor[ix + 1][1] * x) / 255.0f;
		b = (TempColor[ix][2] * (1 - x) + TempColor[ix + 1][2] * x) / 255.0f;
	}
}

//----------------------------------------------------------------------------------------------------------//
BEGIN_EVENT_TABLE(View3D, wxGLCanvas)
EVT_MOTION(View3D::mouseMoved)
EVT_LEFT_DOWN(View3D::mouseDown)
EVT_LEFT_DCLICK(View3D::mouseLeftDouble)
EVT_LEFT_UP(View3D::mouseReleased)
EVT_RIGHT_DOWN(View3D::rightClick)
EVT_LEAVE_WINDOW(View3D::mouseLeftWindow)
EVT_ENTER_WINDOW(View3D::mouseEnterWindow)
EVT_SIZE(View3D::resized)
EVT_KEY_DOWN(View3D::keyPressed)
EVT_KEY_UP(View3D::keyReleased)
EVT_MOUSEWHEEL(View3D::mouseWheelMoved)
EVT_PAINT(View3D::OnPaint)
END_EVENT_TABLE()

void View3D::mouseMoved(wxMouseEvent& event) // it may generate events of mouse moving in/out sub-windows
{
	static int lastW = -1;
	wxPoint pos = event.GetPosition();
	int wid = getSubCoordinates(pos);
	event.SetPosition(pos);
	if (lastW != activeID) //last position is outside the window
	{
		if (wid == activeID) //get in this window now
		{
			mouseEnterWindow(wid, event);
			mouseMoved(wid, event);
		}
		//else still moving out side, so do nothing
	}
	else // lastW == activeID, last position is in the active window
	{
		if (wid != activeID) mouseEnterWindow(activeID, event); // get out of this window now, and stop receiving the message
		else mouseMoved(wid, event); // keep moving in this window
	}
	lastW = wid; // always update the window id
}
void View3D::mouseDown(wxMouseEvent& event)
{
	SetFocus(); //to receive keyboard event
	SetCurrent(*m_context); // bind the right openGL context
	wxPoint pos = event.GetPosition();
	int wid = getSubCoordinates(pos);
	if (wid < 0) return;
	if (wid != activeID) //need to notify the active window has changed
	{
		// notify the window that line sampling is finished. This message may be used to draw a sample graph
		wxCommandEvent event_notify(VIEW3D_NOTIFY);
		event_notify.SetInt(VIEW3D_ACTIVE_VIEW_CHANGED);
		wxPostEvent(this, event_notify); // to View3D
	}
	activeID = wid; // update the active window
	event.SetPosition(pos);
	mouseDown(wid, event);
}
void View3D::mouseLeftDouble(wxMouseEvent& event)
{
	SetFocus(); //to receive keyboard event
	SetCurrent(*m_context); // bind the right openGL context
	wxPoint pos = event.GetPosition();
	int wid = getSubCoordinates(pos);
	if (wid < 0) return;
	activeID = wid; // update the active window
	event.SetPosition(pos);
	mouseLeftDouble(wid, event);
}
void View3D::mouseWheelMoved(wxMouseEvent& event)
{
	wxPoint pos = event.GetPosition();
	int wid = getSubCoordinates(pos);
	if (wid < 0) return;
	event.SetPosition(pos);
	mouseWheelMoved(activeID, event); //only send message to the activeID window even if the mouse isn't in that window now
}
void View3D::mouseReleased(wxMouseEvent& event)
{
	wxPoint pos = event.GetPosition();
	int wid = getSubCoordinates(pos);
	if (wid < 0) return;
	event.SetPosition(pos);
	mouseReleased(activeID, event); //only send message to the activeID window even if the mouse isn't in that window now
}
void View3D::rightClick(wxMouseEvent& event)
{
	SetFocus(); //to receive keyboard event
	SetCurrent(*m_context); // bind the right openGL context
	wxPoint pos = event.GetPosition();
	int wid = getSubCoordinates(pos);
	if (wid < 0) return;
	if (wid != activeID) //need to notify the active window has changed
	{
		// notify the window that line sampling is finished. This message may be used to draw a sample graph
		wxCommandEvent event_notify(VIEW3D_NOTIFY);
		event_notify.SetInt(VIEW3D_ACTIVE_VIEW_CHANGED);
		wxPostEvent(this, event_notify); // to View3D
	}
	activeID = wid; // update the active window
	event.SetPosition(pos);
	rightClick(wid, event);
}
void View3D::mouseEnterWindow(wxMouseEvent& event)
{
	//change the cursor shape to + in order to pick certain pixel
	SetCursor(wxCursor(wxCURSOR_CROSS));
}
void View3D::mouseLeftWindow(wxMouseEvent& event)
{
	wxPoint pos = event.GetPosition();
	wxMouseEvent event_leftUp(wxEVT_LEFT_UP); // simulate mouse left button up event
	event_leftUp.SetPosition(pos);
	wxPostEvent(this, event_leftUp); // to View3D itself
}
void View3D::keyPressed(wxKeyEvent& event)
{
	keyPressed(activeID, event);
}
void View3D::keyReleased(wxKeyEvent& event) 
{
	keyPressed(activeID, event);
}
void View3D::resized(wxSizeEvent& event)
{
	if(b_windowCreated) calcSubWinRect(event.GetSize()); // recalculate all the sub-window positions
}
void View3D::OnPaint(wxPaintEvent& evt)
{
	InitGL();
	renders();
}

void View3D::InitGL()
{
	if (m_context == NULL)
	{
		m_context = new wxGLContext(this);
		SetCurrent(*m_context);
		//do initialization of openGL
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black Background
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		glEnable(GL_POINT_SMOOTH);
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		//glEnable(GL_POLYGON_SMOOTH);
		//glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
		glLineWidth(2.0f);
		glEnable(GL_SCISSOR_TEST); // enable the scissor test to make sure the sub-windows' contents won't overlap

		glText.init(); // to print openGL text. call init() only once before printing openGL text
		gluQuad = gluNewQuadric(); // only call once
		for (int i = 0; i < 3; ++i)
		{
			pv[i].initViewAngle(i);
			glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
			glGenTextures(1, &pv[i].gltexture); //create texture for this view
		}
	}
}

View3D::View3D(wxWindow* parent):
wxGLCanvas(parent, wxID_ANY, wxGLCanvas_args, wxDefaultPosition, wxSize(300, 300), wxFULL_REPAINT_ON_RESIZE)
{
	m_context = NULL;
	wxClientDC dc(this);
	dc.SetAxisOrientation(true, false);
	// To avoid flashing on MSW
	SetBackgroundStyle(wxBG_STYLE_CUSTOM);

	NX = NY = NZ = 0; // must be specified from external
	DX = DY = DZ = 0.1f; // default 0.1 cm

	// mat and bgImage must have one presented, and bgPainter is optional
	mat = NULL; // a 3D matrix that will be shown
	bgImage = NULL; // a mono-color background. Usually from CT or MRI
	bgPainter = NULL; // a colorful background that may override the color generated from bgImage

	// related to the matrix
	b_showMat = true;   // default true
	b_globalISOLines = true; // default true, show the isolines in different views with the same level standard
	NISOLevel = 9; // default 9, how many level should we divide the matrix range
	stdColor = 0.9; //default = 0.9, portion of the color gradient that defines the standard color
	stdValue = 0; // default = max of matrix.This isoline will be rendered with standard color
	lineWidth = 2; // default 2, line width of the isolines, range 1-32 in openGL
	matAlpha = 0.2; // default 0.2, transparent level of the matrix wash, range 0.0-1.0
	colorScheme = 0; // default 0, choose one color scheme to render the matrix

	// related to the background image
	b_showBGImg = true; // default true
	b_globalContrast = true; // default true. All sub-views have the same brightness standard
	brightness = 1.0; // default 1.0, brightness level of the background image, range 0.0-1.0
	contrast_L = 0.0; // default 0.0, low cut-off portion when showing background image, range 0.0-1.0
	contrast_R = 1.0; // default 1.0, high cut-off portion when showing background image, range 0.0-1.0

	b_showBGPainter = true; // default true
	b_showBox = true;  // default true
	b_showCrossLine = true; // default true, whether to show crossing line. It should be false when maximizing one of the sub-windows
	b_showMarkPoint = false; // default false

	b_maxGLWindow = false; //default false, if the view is maximized
	b_windowCreated = false; // default false. Used to deal with size event correctly. set it to true when calling initShow()
	//data related to the sub-windows
	activeID = 0; //current active window's id, default 0 (the first sub window), range 0-3
	theta = THETA_INITIAL;
	phi = PHI_INITIAL; // default theta = 40, phi = 50; used in the 3D rendering

	// data related to drawing something. Note at a time, there's only one drawing process
	drawingType = DRAW_NOTHING; // default 0 (DRAW_NOTHING = 0), drawing nothing
	drawingStatus = 0; // default 0, before drawing anything. Different drawing types will have different meanings.
	b_mouseLeftDown = false;  // default false. whether the left mouse button is pressed down.
	b_showSampledLine = true; // default true. Draw the sampled line after mouse is released in corresponding mode

	sliceAlpha = 0.1; // default 0.1
	sliceColor[0] = sliceColor[1] = sliceColor[2] = 0.1; // default (0.1,0.1,0.1)
	b_showTriangleEdges = true;

	//private member initialization
	g_BGInit = g_MATInit = false;
}

View3D::~View3D()
{
	delete m_context; // will it delete the texture resource? and gluQuad?
}
void View3D::zoomIn()
{
	subw[activeID].picRect.zoom(1.1f);
	drawingType = DRAW_NOTHING;
	renders();
}
void View3D::zoomOut()
{
	subw[activeID].picRect.zoom(0.9f);
	drawingType = DRAW_NOTHING;
	renders();
}
void View3D::resetView()
{
	//b_maxGLWindow = false;
	for (int i = 0; i < 4; ++i) subw[i].picRect.setEmpty();
	for (int i = 0; i < 3; ++i) // init the plan view parameters
	{
		pv[i].initViewAngle(i);
	}
	theta = THETA_INITIAL;
	phi = PHI_INITIAL;
	calcSubWinRect(GetSize());
	renders();
}
void View3D::moveLeft()
{
	subw[activeID].picRect.shift(-10, 0);
	renders();
}
void View3D::moveRight()
{
	subw[activeID].picRect.shift(10, 0);
	renders();
}
void View3D::moveUp()
{
	subw[activeID].picRect.shift(0, 10);
	renders();
}
void View3D::moveDown()
{
	subw[activeID].picRect.shift(0, -10);
	renders();
}
void View3D::activeLeft()
{
	if (activeID < 1) return;
	activeID--;
	renders();
}
void View3D::activeRight()
{
	if (activeID > 2) return;
	activeID++;
	renders();
}
void View3D::rotateLeft()
{
	rotate(1);
	renders();
}
void View3D::rotateRight()
{
	rotate(-1);
	renders();
}

void View3D::renders()
{
	if (!IsShown()) return;
	InitGL();
	SetCurrent(*m_context); //This line is very important because multiple openGL contexts may be created
	wxSize sz = GetSize();
	glViewport(0, 0, sz.x, sz.y);
	glScissor(0, 0, sz.x, sz.y);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f); // Black Background
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if (b_maxGLWindow) render(activeID);
	else
	{
		for (int i = 0; i < 4; ++i)	render(i);
	}

	glFlush();
	SwapBuffers();
}

void View3D::show() //prepare slice to show, and call render to do the real work
{
	bgColor.resize(NTET, 3);
	bgColor = 0;
	if (bgImage != NULL) // calculate the grey level of the phantom
	{
		SFloat minv, maxv;
		if (g_BGInit)
		{
			maxv = g_BGMax;
			minv = g_BGMin;
		}
		else
		{
			bgImage->getMaxMin(maxv, minv);
			g_BGInit = true;
			g_BGMax = maxv;
			g_BGMin = minv;
		}
		double d = maxv - minv;
		double VL = d*contrast_L + minv;
		double VR = d*contrast_R + minv;
		for (int i = 0; i < NTET; ++i)
		{
			if (bgImage->a(i) < VL)
			{
				for (char j = 0; j<3; ++j) bgColor.a(i, j) = 0;
			}
			else if (bgImage->a(i) > VR)
			{
				for (char j = 0; j < 3; ++j) bgColor.a(i, j) = 1;
			}
			else
			{
				for (char j = 0; j < 3; ++j)
				{
					if (VR != VL) bgColor.a(i, j) = (bgImage->a(i) - VL) / (VR - VL);
					else bgColor.a(i, j) = 0.5;
				}
			}
		}
	}
	if (mat != NULL)
	{
		SFloat maxv = 0, minv = 0;
		if (g_MATInit)
		{
			maxv = g_MatMax;
			minv = g_MatMin;
		}
		else
		{
			mat->getMaxMin(maxv, minv); //need performance improvement ???
			g_MATInit = true;
			g_MatMax = maxv;
			g_MatMin = minv;
		}
		SFloat dz = (stdValue - minv) / (NISOLevel + 1);
		NISOLine = int((maxv - minv) / dz);
		if (fabs(maxv - minv - NISOLine*dz) < 1e-9) --NISOLine;
		if (NISOLine == 0)
		{
			wxMessageBox("Cannot generate iso-lines because the two dose levels are tremendously different!\n Please provide similar dose levels^_^");
			exit(0);
		}
		SFloat z0 = minv + dz;
		float r = 0, g = 0, b = 0;
		for (int i = 0; i < NTET; ++i)
		{
			if (mat->a(i) - z0 > 0)
			{
				double kc = ((mat->a(i) - z0) / dz + 1) / double(NISOLevel + 1)*stdColor;
				colorGradient(r, g, b, kc, colorScheme);

				bgColor.a(i, 0) = bgColor.a(i, 0) *(1 - matAlpha) + r*matAlpha;
				bgColor.a(i, 1) = bgColor.a(i, 1) *(1 - matAlpha) + g*matAlpha;
				bgColor.a(i, 2) = bgColor.a(i, 2) *(1 - matAlpha) + b*matAlpha;
			}
		}
	}
	//Let's check bgPainter to override some pixels
	if (bgPainter != NULL)
	{
		for (int i = 0; i < NTET; ++i)
		{
			if (!bgPainter->a(i).isBlack())
			{
				bgColor.a(i, 0) = bgPainter->a(i).r / 255.0f;
				bgColor.a(i, 1) = bgPainter->a(i).g / 255.0f;
				bgColor.a(i, 2) = bgPainter->a(i).b / 255.0f;
			}
		}
	}

	ArrayMgr<float> PixelData(4, 10, 10);
	for (int i = 0; i < 10; ++i)
		for (int j = 0; j < 10; ++j)
		{
			PixelData.a(0, i, j) = sliceColor[0];
			PixelData.a(1, i, j) = sliceColor[1];
			PixelData.a(2, i, j) = sliceColor[2];
			PixelData.a(3, i, j) = sliceAlpha;
		}
	for (int i = 0; i < 3; ++i)
	{
		glBindTexture(GL_TEXTURE_2D, pv[i].gltexture);//tell openGL this is the texture it will deal with
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 10, 10, 0, GL_RGBA, GL_FLOAT, PixelData.getP());
	}
	renders();
}

void View3D::setSliceRGBA(float r, float g, float b, float a)
{
	sliceColor[0] = r;
	sliceColor[1] = g;
	sliceColor[2] = b;
	sliceAlpha = a;
	show();
}
void View3D::drawCuboid()//first call to start drawing box in 2D and second call to draw cuboid in another dimension
{
	if (drawingType == DRAW_CUBOID)
	{
		if (drawingStatus == 3) drawingStatus = 4; // start drawing in another direction
		else
		{
			drawingType = DRAW_NOTHING;
			drawingStatus = 0;
		}
	}
	else
	{
		drawingType = DRAW_CUBOID;
		drawingStatus = 0;
	}
	renders();
}
void View3D::profileH() //first call to enable horizontal line sampling, and second call to disable
{
	if (drawingType == DRAW_LINE_HORIZONTAL) drawingType = DRAW_NOTHING;
	else drawingType = DRAW_LINE_HORIZONTAL;
	drawingStatus = 0;
	renders();
}
void View3D::profileV() //first call to enable vertical line sampling, and second call to disable
{
	if (drawingType == DRAW_LINE_VERTICAL) drawingType = DRAW_NOTHING;
	else drawingType = DRAW_LINE_VERTICAL;
	drawingStatus = 0;
	renders();
}
void View3D::profileT() //first call to enable arbitrary line sampling, and second call to disable
{
	if (drawingType == DRAW_LINE_ANY) drawingType = DRAW_NOTHING;
	else drawingType = DRAW_LINE_ANY;
	drawingStatus = 0;
	renders();
}

void View3D::nextSlice()
{
	if (activeID < VIEW3D)
	{
		cross[activeID] += 1;
		int NSlice = getNSlice(activeID);
		if (cross[activeID] >= NSlice) cross[activeID] = NSlice - 1;

		wxCommandEvent event_notify(VIEW3D_NOTIFY);
		event_notify.SetInt(VIEW3D_CROSS_MOVED);
		wxPostEvent(this, event_notify); // to View3D
		show();
	}
}
void View3D::previousSlice()
{
	if (activeID < VIEW3D)
	{
		cross[activeID] -= 1;
		int NSlice = getNSlice(activeID);
		if (cross[activeID] < 0) cross[activeID] = 0;

		// Notify to the parent window that the crossing point has been updated
		wxCommandEvent event_notify(VIEW3D_NOTIFY);
		event_notify.SetInt(VIEW3D_CROSS_MOVED);
		wxPostEvent(this, event_notify); // to View3D
		show();
	}
}

/****************************************************************************************************************************/
void View3D::setMesh(ArrayMgr<SFloat>* node, ArrayMgr<int>* elem, ArrayMgr<int>* f2tm)
{
	if (f2tm->getInnerLength() == 0) f2tm = NULL;
	wxBusyInfo wait("Please wait, establishing face map...");
	this->node = node;
	this->elem = elem;
	//need to calculate the range and center 
	NTET = elem->getWidth(1);
	int nNode = node->getWidth(1);
	xmin = xmax = node->a(0, 0);
	ymin = ymax = node->a(0, 1);
	zmin = zmax = node->a(0, 2);
	for (int i = 1; i < nNode; ++i)
	{
		xmin = min(xmin, node->a(i, 0));
		xmax = max(xmax, node->a(i, 0));
		ymin = min(ymin, node->a(i, 1));
		ymax = max(ymax, node->a(i, 1));
		zmin = min(zmin, node->a(i, 2));
		zmax = max(zmax, node->a(i, 2));
	}
	xcenter = (xmin + xmax) / 2;
	ycenter = (ymin + ymax) / 2;
	zcenter = (zmin + zmax) / 2;
	// determine DX = DY = DZ and NX, NY, NZ
	SFloat xL = xmax - xmin;
	SFloat yL = ymax - ymin;
	SFloat zL = zmax - zmin;
	const int NMAX = 100;
	SFloat maxSide = max(xL, max(yL, zL));
	DX = DY = DZ = maxSide / NMAX;
	NX = xL / DX;
	if (xL > NX*DX) ++NX;
	NY = yL / DY;
	if (yL > NY*DY) ++NY;
	NZ = zL / DZ;
	if (zL > NZ*DZ) ++NZ;

	//std::unordered_map<string, FaceValue> fmap; //face map
	auto pmap = new std::unordered_map<string, FaceValue>(); //face map
	std::unordered_map<string, FaceValue>& fmap = *pmap;
	//scan the tetrahedral to get the ranges
	tetLim.resize(NTET, 6);
	fNorm.resize(NTET, 4);
	for (int i = 0; i < NTET; ++i)
	{
		int jn = elem->a(i, 0);
		tetLim.a(i, 0) = tetLim.a(i, 1) = node->a(jn, 0);
		tetLim.a(i, 2) = tetLim.a(i, 3) = node->a(jn, 1);
		tetLim.a(i, 4) = tetLim.a(i, 5) = node->a(jn, 2);
		for (int j = 1; j < 4; ++j)
		{
			jn = elem->a(i, j);
			tetLim.a(i, 0) = min(tetLim.a(i, 0), node->a(jn, 0));
			tetLim.a(i, 1) = min(tetLim.a(i, 1), node->a(jn, 0));
			tetLim.a(i, 2) = min(tetLim.a(i, 2), node->a(jn, 1));
			tetLim.a(i, 3) = min(tetLim.a(i, 3), node->a(jn, 1));
			tetLim.a(i, 4) = min(tetLim.a(i, 4), node->a(jn, 2));
			tetLim.a(i, 5) = min(tetLim.a(i, 5), node->a(jn, 2));
		}

		//calculate the 4 normal vectors
		int v0 = elem->a(i, 0);
		int v1 = elem->a(i, 1);
		int v2 = elem->a(i, 2);
		int v3 = elem->a(i, 3);
		Vector v12(node->a(v2, 0) - node->a(v1, 0), node->a(v2, 1) - node->a(v1, 1), node->a(v2, 2) - node->a(v1, 2));
		Vector v13(node->a(v3, 0) - node->a(v1, 0), node->a(v3, 1) - node->a(v1, 1), node->a(v3, 2) - node->a(v1, 2));
		fNorm.a(i, 0) = v12&v13;
		fNorm.a(i, 0).normalizeToUnitLength();

		Vector v03(node->a(v3, 0) - node->a(v0, 0), node->a(v3, 1) - node->a(v0, 1), node->a(v3, 2) - node->a(v0, 2));
		Vector v02(node->a(v2, 0) - node->a(v0, 0), node->a(v2, 1) - node->a(v0, 1), node->a(v2, 2) - node->a(v0, 2));
		fNorm.a(i, 1) = v03&v02;
		fNorm.a(i, 1).normalizeToUnitLength();

		Vector v01(node->a(v1, 0) - node->a(v0, 0), node->a(v1, 1) - node->a(v0, 1), node->a(v1, 2) - node->a(v0, 2));
		fNorm.a(i, 2) = v01&v03;
		fNorm.a(i, 2).normalizeToUnitLength();

		fNorm.a(i, 3) = v02&v01;
		fNorm.a(i, 3).normalizeToUnitLength();

		if (f2tm == NULL)
		{
			FaceKey fkey;
			fkey.set(v1, v2, v3);
			fmap[fkey()].set(i, 0);

			fkey.set(v0, v2, v3);
			fmap[fkey()].set(i, 1);

			fkey.set(v1, v0, v3);
			fmap[fkey()].set(i, 2);

			fkey.set(v1, v2, v0);
			fmap[fkey()].set(i, 3);
		}
			
	}

	//iterate the face map to set the adjacent tetrahedral index for each face
	f2t.resize(NTET, 4);
	f2t = -1; //initially no adjacent face

	if (f2tm == NULL)
	{
		for (auto it = fmap.begin(); it != fmap.end(); ++it)
		{
			FaceValue& fv = it->second;
			if (fv.idt0 == -1)
			{
				printf("Error: a face doesn't belong to any tetrahedral.\n");
				getchar();
				return;
			}
			f2t.a(fv.idt0, fv.idf0) = fv.idt1;
			if (fv.idt1 != -1) f2t.a(fv.idt1, fv.idf1) = fv.idt0;
		}
	}
	else
	{
		int nFace = f2tm->getWidth(1);
		for (int i = 0; i < nFace; ++i)
		{
			if (f2tm->a(i,0) == -1)
			{
				printf("Error: a face doesn't belong to any tetrahedral.\n");
				getchar();
				return;
			}
			f2t.a(f2tm->a(i, 0), f2tm->a(i, 1)) = f2tm->a(i, 2);
			if (f2tm->a(i, 2) != -1) f2t.a(f2tm->a(i, 2), f2tm->a(i, 3)) = f2tm->a(i, 0);
		}
	}

#ifdef NDEBUG
		delete pmap; //it's very slow monitored by debugger
#endif
}
void View3D::setMatrix(ArrayMgr<SFloat>* mat)
{
	this->mat = mat;	
}
void View3D::setBGImage(ArrayMgr<SFloat>* bgImage)
{
	this->bgImage = bgImage;
}
void View3D::setBGPainter(ArrayMgr<RGBPixel>* bgPainter) // (0,0,0) means not to replace
{
	 this->bgPainter = bgPainter;
}
void View3D::setVoxelSize(int NX, int NY, int NZ, double DX, double DY, double DZ)
{
	this->NX = NX;
	this->NY = NY;
	this->NZ = NZ;
	this->DX = DX;
	this->DY = DY;
	this->DZ = DZ;
}
bool View3D::initShow(bool b_reinit) // init to execute a showing job
{
	g_BGInit = g_MATInit = false; // indicate to recalculate global max and min value
	if (b_windowCreated&& !b_reinit)// only check if they're empty
	{
		if (mat != NULL && mat->empty()) mat = NULL;
		if (bgImage != NULL && bgImage->empty()) bgImage = NULL;
		if (bgPainter != NULL && bgPainter->empty()) bgPainter = NULL;
		show();
		return true;
	}

	if (mat != NULL && !mat->empty())
	{
		stdValue = mat->Max();
		int maxi, mini;
		mat->getMaxMinIndex(maxi, mini);
		
		//To do: cross point based on the max value
		cross[XAXIS] = NX / 2; //crossing point index, default: the index of max value of the matrix
		cross[YAXIS] = NY / 2;
		cross[ZAXIS] = NZ / 2;
	}
	else mat = NULL;

	if (mat != NULL) //check if the matrix dimensions match
	{
		if (bgImage != NULL && !bgImage->empty())
		{
			if (!bgImage->equalDims(*mat))
			{
				bgImage = NULL;
				return false; // unequal dimensions
			}
		}
		else bgImage = NULL;
	}
	else // no mat
	{
		if (bgImage != NULL) // only have background image
		{
			cross[XAXIS] = NX / 2; //crossing point index, default: the index of max value of the matrix
			cross[YAXIS] = NY / 2;
			cross[ZAXIS] = NZ / 2;
		}
		else return false; // nothing to show
	}

	if (mat != NULL&&mat->isUniform()) wxMessageBox("Warning: the matrix to show is uniform. Make sure that's what you want to see!");

	if (bgImage != NULL && bgPainter != NULL && !bgPainter->equalDimsT<SFloat>(*bgImage))
	{
		bgPainter = NULL;
		return false; // dimension doesn't match
	}

	b_windowCreated = true; // enable the resizing
	for (int i = 0; i < 4; ++i) subw[i].picRect.setEmpty();
	//calculate the window size and picture position
	calcSubWinRect(GetSize());

	show();
	return true;
}

int View3D::getSeq(int istart, int mod) //auxilary function for rotate()
{
	int nmod = (istart + mod) % 4;
	int seq[4] = { 0, 1, 0, -1 };
	return seq[nmod];
}
int View3D::getNSlice(int wid)
{
	if (ZAXIS == wid) return NZ;
	else if (YAXIS == wid) return NY;
	else if (XAXIS == wid) return NX;
	else return 1;
}
void View3D::rotate(int dir) // dir + means counter-clockwise, - means clockwise
{
	if (dir == 0) return;
	int wid = activeID;
	if (wid < 0 || wid >= VIEW3D) return;
	if (dir > 0) ++pv[wid].rotateCounter;
	else if (dir < 0) --pv[wid].rotateCounter;
	subw[wid].picRect.transpose();
	
	//change it to equivalent positive integer
	int mod = pv[wid].rotateCounter % 4;
	if (mod < 0) mod += 4;
	//change the view angle
	if (wid == ZAXIS)
	{
		pv[wid].upv[0] = getSeq(2, mod);
		pv[wid].upv[1] = getSeq(1, mod);
	}
	else if (wid == YAXIS)
	{
		pv[wid].upv[0] = getSeq(2, mod);
		pv[wid].upv[2] = getSeq(3, mod);
	}
	else if (wid == XAXIS)
	{
		pv[wid].upv[1] = getSeq(1, mod);
		pv[wid].upv[2] = getSeq(2, mod);
	}
}
float View3D::getScale(int wid) // get the scale from 3D coordinates (cm) to pixel coordinates. used in render(int wid)
{
	if (ZAXIS == wid)
	{
		float lx = NX*DX;
		if (pv[wid].upv[0] == 0) return subw[wid].picRect.width / lx;
		else return subw[wid].picRect.height / lx;
	}
	else if (YAXIS == wid)
	{
		float lx = NX*DX;
		if (pv[wid].upv[0] == 0) return subw[wid].picRect.width / lx;
		else return subw[wid].picRect.height / lx;
	}
	else if (XAXIS == wid)
	{
		float lz = NZ*DZ;
		if (pv[wid].upv[2] == 0) return subw[wid].picRect.width / lz;
		else return subw[wid].picRect.height / lz;
	}
	else
	{
		double lx = NX*DX;
		double ly = NY*DY;
		double lz = NZ*DZ;
		double maxL = sqrt(lx*lx + ly*ly + lz*lz);
		return subw[wid].picRect.width / maxL;
	}
}
void View3D::getDWDH(int wid, double& DW, double& DH) // I may not need to exchange DW and DH ??? used when exporting the dose in DICOM ?
{
	if (ZAXIS == wid)
	{
		if (pv[wid].upv[0] == 0)
		{
			DW = DY;
			DH = DX;
		}
		else
		{
			DW = DX;
			DH = DY;
		}
	}
	else if (YAXIS == wid)
	{
		if (pv[wid].upv[0] == 0)
		{
			DW = DZ;
			DH = DX;
		}
		else
		{
			DW = DX;
			DH = DZ;
		}
	}
	else if (XAXIS == wid)
	{
		if (pv[wid].upv[2] == 0)
		{
			DW = DY;
			DH = DZ;
		}
		else
		{
			DW = DZ;
			DH = DY;
		}
	}
}
void View3D::prepare(int wid) // prepare contour and glTexture
{
	
}

void View3D::render(int wid) // render one sub-window
{
	if (NULL == mat && NULL == bgImage) return;

	FloatRect& glRect = subw[wid].glRect;
	FloatRect& picRect = subw[wid].picRect;
	glViewport(glRect.x, glRect.y, glRect.width, glRect.height);
	glScissor(glRect.x, glRect.y, glRect.width, glRect.height);
	
	//---------------------------------------Begin 3D drawing----------------------------------------
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	//apply the picking matrix if it's in selecting mode
	int gl_render_mode = GL_RENDER;
	glGetIntegerv(GL_RENDER_MODE, &gl_render_mode);
	if (GL_SELECT == gl_render_mode)
	{
		gluPickMatrix(cMouse.x, cMouse.y, 2, 2, subw[wid].viewport);
	}

	glOrtho(0, glRect.width, 0, glRect.height, -10000.0f, 10000.0f); // define the 3D space in pixels
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glColor3f(1, 1, 1);
	// translate the origin to the center of picRect
	glTranslatef(picRect.centerX(), picRect.centerY(), 0);

	if (wid == VIEW3D) return draw3D();

	// change the look at angle
	double R = 1000;
	gluLookAt(R*pv[wid].eyev[0], R*pv[wid].eyev[1], R*pv[wid].eyev[2], 0.0f, 0.0f, 0.0f, pv[wid].upv[0], pv[wid].upv[1], pv[wid].upv[2]);

	// need to calculate the ratio to scale from cm to pixel
	float scale = getScale(wid);
	glScalef(scale, scale, scale);

	// --------Note the rest plot use 3D world coordinates now-------- //
	double lx = xmax - xmin;
	double ly = ymax - ymin;
	double lz = zmax - zmin;
	glTranslatef(-lx / 2, -ly / 2, -lz / 2);

	glLineWidth(lineWidth);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_POLYGON_OFFSET_LINE);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glPolygonOffset(-1.0f, -1.0f);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA);
// 	glEnable(GL_LIGHTING);
// 	glEnable(GL_LIGHT0);
// 	glEnable(0x803A);
// 	GLfloat lightAmbient[] = { 0.3, 0.3, 0.3, 1.0 };
// 	GLfloat lightDiffuse[] = { 0.5, 0.5, 0.5, 1.0 };
// 	GLfloat lightSpecular[] = { 0.5, 0.5, 0.5, 1.0 };
// 	glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
// 	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
// 	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
// 	GLfloat lightPosition[] = { 10.0, 0, 0.0, 0.0 };   // w=0.0
// 	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
// 	//glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lightAmbient);
// 	glEnable(GL_COLOR_MATERIAL);
// 	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
// 	glColorMaterial(GL_FRONT, GL_SPECULAR);
	glPushMatrix();
	//glTranslatef(xcenter, ycenter, zcenter);
	SFloat zlim = cross[ZAXIS] * DZ + zmin;
	SFloat ylim = cross[YAXIS] * DY + ymin;
	SFloat xlim = cross[XAXIS] * DX + xmin;
#define canRender(i) (wid == ZAXIS&&tetLim.a(i, 4) > zlim || wid == YAXIS&&tetLim.a(i, 2) > ylim || wid == XAXIS&&tetLim.a(i, 1) < xlim)
	for (int i = 0; i < NTET; ++i)
	{
		if (canRender(i))
		{
			int v0 = elem->a(i, 0);
			int v1 = elem->a(i, 1);
			int v2 = elem->a(i, 2);
			int v3 = elem->a(i, 3);
			int mc = elem->a(i, 4);
			if (mc >= 7)
			{
				int im = mc % 7;
				mc = mc - 7 * im;
			}
			//judge if this face is adjacent to any face
			int iad = f2t.a(i, 0);
			if (iad == -1 || !canRender(iad)) //surface triangle
			{
				if (b_showTriangleEdges)
				{
					glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
					glColor4f(MatColor[mc][0], MatColor[mc][1], MatColor[mc][2], 0.5f);
					glBegin(GL_POLYGON);
					glVertex3f(node->a(v1, 0), node->a(v1, 1), node->a(v1, 2));
					glVertex3f(node->a(v2, 0), node->a(v2, 1), node->a(v2, 2));
					glVertex3f(node->a(v3, 0), node->a(v3, 1), node->a(v3, 2));
					glEnd();
				}
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				glColor3f(bgColor.a(i, 0), bgColor.a(i, 1), bgColor.a(i, 2));
				//glBegin(GL_TRIANGLES);
				glBegin(GL_POLYGON);
				//glNormal3fv((GLfloat*)&fNorm.a(i, 0));
				glVertex3f(node->a(v1, 0), node->a(v1, 1), node->a(v1, 2));
				glVertex3f(node->a(v2, 0), node->a(v2, 1), node->a(v2, 2));
				glVertex3f(node->a(v3, 0), node->a(v3, 1), node->a(v3, 2));
				glEnd();
			}

			iad = f2t.a(i, 1);
			if (iad == -1 || !canRender(iad)) //surface triangle
			{
				if (b_showTriangleEdges)
				{
					glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
					glColor4f(MatColor[mc][0], MatColor[mc][1], MatColor[mc][2], 0.5f);
					glBegin(GL_POLYGON);
					glVertex3f(node->a(v0, 0), node->a(v0, 1), node->a(v0, 2));
					glVertex3f(node->a(v3, 0), node->a(v3, 1), node->a(v3, 2));
					glVertex3f(node->a(v2, 0), node->a(v2, 1), node->a(v2, 2));
					glEnd();
				}
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				glColor3f(bgColor.a(i, 0), bgColor.a(i, 1), bgColor.a(i, 2));
				//glBegin(GL_TRIANGLES);
				glBegin(GL_POLYGON);
				//glNormal3fv((GLfloat*)&fNorm.a(i, 1));
				glVertex3f(node->a(v0, 0), node->a(v0, 1), node->a(v0, 2));
				glVertex3f(node->a(v3, 0), node->a(v3, 1), node->a(v3, 2));
				glVertex3f(node->a(v2, 0), node->a(v2, 1), node->a(v2, 2));
				glEnd();
			}

			iad = f2t.a(i, 2);
			if (iad == -1 || !canRender(iad)) //surface triangle
			{
				if (b_showTriangleEdges)
				{
					glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
					glColor4f(MatColor[mc][0], MatColor[mc][1], MatColor[mc][2], 0.5f);
					glBegin(GL_POLYGON);
					glVertex3f(node->a(v0, 0), node->a(v0, 1), node->a(v0, 2));
					glVertex3f(node->a(v1, 0), node->a(v1, 1), node->a(v1, 2));
					glVertex3f(node->a(v3, 0), node->a(v3, 1), node->a(v3, 2));
					glEnd();
				}
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				glColor3f(bgColor.a(i, 0), bgColor.a(i, 1), bgColor.a(i, 2));
				//glBegin(GL_TRIANGLES);
				glBegin(GL_POLYGON);
				//glNormal3fv((GLfloat*)&fNorm.a(i, 2));
				glVertex3f(node->a(v0, 0), node->a(v0, 1), node->a(v0, 2));
				glVertex3f(node->a(v1, 0), node->a(v1, 1), node->a(v1, 2));
				glVertex3f(node->a(v3, 0), node->a(v3, 1), node->a(v3, 2));
				glEnd();
			}

			iad = f2t.a(i, 3);
			if (iad == -1 || !canRender(iad)) //surface triangle
			{
				if (b_showTriangleEdges)
				{
					glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
					glColor4f(MatColor[mc][0], MatColor[mc][1], MatColor[mc][2], 0.5f);
					glBegin(GL_POLYGON);
					glVertex3f(node->a(v0, 0), node->a(v0, 1), node->a(v0, 2));
					glVertex3f(node->a(v2, 0), node->a(v2, 1), node->a(v2, 2));
					glVertex3f(node->a(v1, 0), node->a(v1, 1), node->a(v1, 2));
					glEnd();
				}
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				glColor3f(bgColor.a(i, 0), bgColor.a(i, 1), bgColor.a(i, 2));
				//glBegin(GL_TRIANGLES);
				glBegin(GL_POLYGON);
				//glNormal3fv((GLfloat*)&fNorm.a(i, 3));
				glVertex3f(node->a(v0, 0), node->a(v0, 1), node->a(v0, 2));
				glVertex3f(node->a(v2, 0), node->a(v2, 1), node->a(v2, 2));
				glVertex3f(node->a(v1, 0), node->a(v1, 1), node->a(v1, 2));
				glEnd();
			}
		}
	}
	glPopMatrix();
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_POLYGON_OFFSET_LINE);
	glDisable(GL_CULL_FACE);
	glDisable(GL_BLEND);
	glClear(GL_DEPTH_BUFFER_BIT);

	if (gl_render_mode == GL_SELECT || b_showBGImg)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		//prepare for picking
		glInitNames();
		glPushName(VIEW3D);

		//draw each phantom and dose wash
		glEnable(GL_TEXTURE_2D);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glBindTexture(GL_TEXTURE_2D, pv[wid].gltexture);
		if (ZAXIS == wid)
		{
			glLoadName(ZAXIS);
			float z = cross[ZAXIS] * DZ;
			glBegin(GL_QUADS);
			glTexCoord2f(1.0, 0.0);
			glVertex3f(lx, 0, z);
			glTexCoord2f(0.0, 0.0);
			glVertex3f(0, 0, z);
			glTexCoord2f(0.0, 1.0);
			glVertex3f(0, ly, z);
			glTexCoord2f(1.0, 1.0);
			glVertex3f(lx, ly, z);
			glEnd();
		}
		else if (YAXIS == wid)
		{
			glLoadName(YAXIS);
			float y = cross[YAXIS] * DY;
			glBegin(GL_QUADS);
			glTexCoord2f(1.0, 1.0);
			glVertex3f(lx, y, lz);
			glTexCoord2f(0.0, 1.0);
			glVertex3f(0, y, lz);
			glTexCoord2f(0.0, 0.0);
			glVertex3f(0, y, 0);
			glTexCoord2f(1.0, 0.0);
			glVertex3f(lx, y, 0);
			glEnd();
		}
		else if (XAXIS == wid)
		{
			glLoadName(XAXIS);
			float x = cross[XAXIS] * DX;
			glBegin(GL_QUADS);
			glTexCoord2f(1.0, 0.0);
			glVertex3f(x, 0, lz);
			glTexCoord2f(1.0, 1.0);
			glVertex3f(x, ly, lz);
			glTexCoord2f(0.0, 1.0);
			glVertex3f(x, ly, 0);
			glTexCoord2f(0.0, 0.0);
			glVertex3f(x, 0, 0);
			glEnd();
		}
		glLoadName(VIEW3D);//give other irrelevant object name
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_BLEND);
	}
	

	//store matrix for 2D to 3D coordinates reverse projection
	if (GL_RENDER == gl_render_mode)
	{
		glGetIntegerv(GL_VIEWPORT, subw[wid].viewport);
		glGetDoublev(GL_MODELVIEW_MATRIX, subw[wid].modelview);
		glGetDoublev(GL_PROJECTION_MATRIX, subw[wid].projection);
	}
	else return; // the rest drawing is not relevant to object picking, so we can return earlier

	

	// draw crossing lines
	if (b_showCrossLine) 
	{
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glColor3f(1.0f, 1.0f, 1.0f);
		glLineWidth(1.0f);
		glPushMatrix();

		if (ZAXIS == wid)
		{
			glTranslatef(0, cross[YAXIS] * DY, cross[ZAXIS] * DZ);
			glBegin(GL_LINES); //draw horizontal line
			glVertex2f(0, 0);
			glVertex2f(lx, 0);
			glEnd();
			glPopMatrix();

			glPushMatrix();
			glTranslatef(cross[XAXIS] * DX, 0, cross[ZAXIS] * DZ);
			glBegin(GL_LINES); //draw vertical line
			glVertex2f(0, 0);
			glVertex2f(0, ly);
			glEnd();
		}
		else if (YAXIS == wid) // I may need to change this ???
		{
			glTranslatef(0, cross[YAXIS] * DY, cross[ZAXIS] * DZ);
			glBegin(GL_LINES);
			glVertex3f(0, 0, 0);
			glVertex3f(lx, 0, 0);
			glEnd();
			glPopMatrix();

			glPushMatrix();
			glTranslatef(cross[XAXIS] * DX, cross[YAXIS] * DY, 0);
			glBegin(GL_LINES);
			glVertex3f(0, 0, 0);
			glVertex3f(0, 0, lz);
			glEnd();
		}
		else if (XAXIS == wid)
		{
			glTranslatef(cross[XAXIS] * DX, 0, cross[ZAXIS] * DZ);
			glBegin(GL_LINES);
			glVertex3f(0, 0, 0);
			glVertex3f(0, ly, 0);
			glEnd();
			glPopMatrix();

			glPushMatrix();
			glTranslatef(cross[XAXIS] * DX, cross[YAXIS] * DY, 0);
			glBegin(GL_LINES);
			glVertex3f(0, 0, 0);
			glVertex3f(0, 0, lz);
			glEnd();
		}
		glPopMatrix();
		glPopAttrib();
	}

	// draw sampling lines
	if (drawingType == DRAW_LINE_ANY && wid == activeID && (drawingStatus == 2 || drawingStatus == 3))
	{
		INT3D rp = pPos - cPos;
		if (rp.distance() > 10) // only draw the line when the line isn't too short
		{
			glPushMatrix();
			glColor3f(1, 1, 0);
			GLdouble len = 1;
			if (ZAXIS == wid)
			{
				glTranslated(pObjX, pObjY, cross[ZAXIS]*DZ);
				GLdouble dx = (cObjX - pObjX);
				GLdouble dy = (cObjY - pObjY);
				len = sqrt(dx*dx + dy*dy);
				// calculate the rotation angle
				double ra = 0;
				if (0 == dx)
				{
					if (dy > 0) ra = 90;
					else ra = -90;
				}
				else
				{
					ra = atan(dy / dx) * 180 / PI;
					if (dx < 0) ra += 180;
				}
				glRotatef(ra, 0, 0, 1);
				glRotatef(90, 0, 1, 0); // turn it to x axis first
			}
			else if (YAXIS == wid)
			{
				glTranslated(pObjX, cross[YAXIS] * DY, pObjZ);
				GLdouble dx = (cObjX - pObjX);
				GLdouble dz = (cObjZ - pObjZ);
				len = sqrt(dx*dx + dz*dz);
				// calculate the rotation angle
				double ra = 0;
				if (0 == dx)
				{
					if (dz > 0) ra = 90;
					else ra = -90;
				}
				else
				{
					ra = atan(dz / dx) * 180 / PI;
					if (dx < 0) ra += 180;
				}
				ra = 90 - ra;
				glRotatef(ra, 0, 1, 0);
			}
			else // XAXIS
			{
				glTranslated(cross[XAXIS]*DX, pObjY, pObjZ);
				GLdouble dy = (cObjY - pObjY);
				GLdouble dz = (cObjZ - pObjZ);
				len = sqrt(dy*dy + dz*dz);
				// calculate the rotation angle
				double ra = 0;
				if (0 == dz)
				{
					if (dy > 0) ra = 90;
					else ra = -90;
				}
				else
				{
					ra = atan(dy / dz) * 180 / PI;
					if (dz < 0) ra += 180;
				}
				glRotatef(ra, -1, 0, 0);
			}
			glPushMatrix();
			glScalef(len, len, len);
			glBegin(GL_LINES);
			glVertex3f(0, 0, 0);
			glVertex3f(0, 0, 1);
			glEnd();
			glPopMatrix();
			double la = 5 / pixelSize() / scale;
			glTranslated(0, 0, len-la);
			
			gluCylinder(gluQuad, 0.15*la, 0, la, 4, 4);
			glPopMatrix();
		}
	}

	if (drawingType == DRAW_LINE_VERTICAL && wid == activeID && (drawingStatus == 2 || drawingStatus == 3))
	{
		INT3D rp = pPos - cPos;
		if (rp.distance() > 10) // only draw the line when the line isn't too short
		{
			glPushMatrix();
			glColor3f(1, 1, 0);
			GLdouble len = 1;
			if (ZAXIS == wid)
			{
				if (pv[wid].transposed()) cObjY = pObjY;
				else cObjX = pObjX;
				glTranslated(pObjX, pObjY, cross[ZAXIS] * DZ);
				GLdouble dx = (cObjX - pObjX);
				GLdouble dy = (cObjY - pObjY);
				len = sqrt(dx*dx + dy*dy);
				// calculate the rotation angle
				double ra = 0;
				if (0 == dx)
				{
					if (dy > 0) ra = 90;
					else ra = -90;
				}
				else
				{
					ra = atan(dy / dx) * 180 / PI;
					if (dx < 0) ra += 180;
				}
				glRotatef(ra, 0, 0, 1);
				glRotatef(90, 0, 1, 0); // turn it to x axis first
			}
			else if (YAXIS == wid)
			{
				if (pv[wid].transposed()) cObjZ = pObjZ;
				else cObjX = pObjX;
				glTranslated(pObjX, cross[YAXIS] * DY, pObjZ);
				GLdouble dx = (cObjX - pObjX);
				GLdouble dz = (cObjZ - pObjZ);
				len = sqrt(dx*dx + dz*dz);
				// calculate the rotation angle
				double ra = 0;
				if (0 == dx)
				{
					if (dz > 0) ra = 90;
					else ra = -90;
				}
				else
				{
					ra = atan(dz / dx) * 180 / PI;
					if (dx < 0) ra += 180;
				}
				ra = 90 - ra;
				glRotatef(ra, 0, 1, 0);
			}
			else // XAXIS
			{
				if (pv[wid].transposed()) cObjY = pObjY;
				else cObjZ = pObjZ;
				glTranslated(cross[XAXIS] * DX, pObjY, pObjZ);
				GLdouble dy = (cObjY - pObjY);
				GLdouble dz = (cObjZ - pObjZ);
				len = sqrt(dy*dy + dz*dz);
				// calculate the rotation angle
				double ra = 0;
				if (0 == dz)
				{
					if (dy > 0) ra = 90;
					else ra = -90;
				}
				else
				{
					ra = atan(dy / dz) * 180 / PI;
					if (dz < 0) ra += 180;
				}
				glRotatef(ra, -1, 0, 0);
			}
			glPushMatrix();
			glScalef(len, len, len);
			glBegin(GL_LINES);
			glVertex3f(0, 0, 0);
			glVertex3f(0, 0, 1);
			glEnd();
			glPopMatrix();
			double la = 5 / pixelSize() / scale;
			glTranslated(0, 0, len - la);

			gluCylinder(gluQuad, 0.15*la, 0, la, 4, 4);
			glPopMatrix();
		}
	}

	if (drawingType == DRAW_LINE_HORIZONTAL && wid == activeID && (drawingStatus == 2 || drawingStatus == 3))
	{
		INT3D rp = pPos - cPos;
		if (rp.distance() > 10) // only draw the line when the line isn't too short
		{
			glPushMatrix();
			glColor3f(1, 1, 0);
			GLdouble len = 1;
			if (ZAXIS == wid)
			{
				if (!pv[wid].transposed()) cObjY = pObjY;
				else cObjX = pObjX;
				glTranslated(pObjX, pObjY, cross[ZAXIS] * DZ);
				GLdouble dx = (cObjX - pObjX);
				GLdouble dy = (cObjY - pObjY);
				len = sqrt(dx*dx + dy*dy);
				// calculate the rotation angle
				double ra = 0;
				if (0 == dx)
				{
					if (dy > 0) ra = 90;
					else ra = -90;
				}
				else
				{
					ra = atan(dy / dx) * 180 / PI;
					if (dx < 0) ra += 180;
				}
				glRotatef(ra, 0, 0, 1);
				glRotatef(90, 0, 1, 0); // turn it to x axis first
			}
			else if (YAXIS == wid)
			{
				if (!pv[wid].transposed()) cObjZ = pObjZ;
				else cObjX = pObjX;
				glTranslated(pObjX, cross[YAXIS] * DY, pObjZ);
				GLdouble dx = (cObjX - pObjX);
				GLdouble dz = (cObjZ - pObjZ);
				len = sqrt(dx*dx + dz*dz);
				// calculate the rotation angle
				double ra = 0;
				if (0 == dx)
				{
					if (dz > 0) ra = 90;
					else ra = -90;
				}
				else
				{
					ra = atan(dz / dx) * 180 / PI;
					if (dx < 0) ra += 180;
				}
				ra = 90 - ra;
				glRotatef(ra, 0, 1, 0);
			}
			else // XAXIS
			{
				if (!pv[wid].transposed()) cObjY = pObjY;
				else cObjZ = pObjZ;
				glTranslated(cross[XAXIS] * DX, pObjY, pObjZ);
				GLdouble dy = (cObjY - pObjY);
				GLdouble dz = (cObjZ - pObjZ);
				len = sqrt(dy*dy + dz*dz);
				// calculate the rotation angle
				double ra = 0;
				if (0 == dz)
				{
					if (dy > 0) ra = 90;
					else ra = -90;
				}
				else
				{
					ra = atan(dy / dz) * 180 / PI;
					if (dz < 0) ra += 180;
				}
				glRotatef(ra, -1, 0, 0);
			}
			glPushMatrix();
			glScalef(len, len, len);
			glBegin(GL_LINES);
			glVertex3f(0, 0, 0);
			glVertex3f(0, 0, 1);
			glEnd();
			glPopMatrix();
			double la = 5 / pixelSize() / scale;
			glTranslated(0, 0, len - la);

			gluCylinder(gluQuad, 0.15*la, 0, la, 4, 4);
			glPopMatrix();
		}
	}

	// drawing the cuboid to select a sub-matrix
	if (drawingType == DRAW_CUBOID)
	{
		if (drawingStatus == 0 || drawingStatus == 4) //need to indicate where's the mouse
		{
			if (!(drawingStatus == 4 && wid == boxDirection))
			{
				wxPoint ov(glRect.x, glRect.y);
				wxPoint r = cMouse - ov;
				ov.x = picRect.x;
				ov.y = picRect.y;
				r -= ov;
				if (0 <= r.x && r.x <= picRect.width && 0 <= r.y && r.y <= picRect.height)
				{
					glPushAttrib(GL_ALL_ATTRIB_BITS);
					const wxColour& color = *wxCYAN;
					glColor3ub(color.Red(), color.Green(), color.Blue());
					glLineWidth(1.0f);
					glPushMatrix();
					glLoadIdentity();
					glTranslatef(picRect.x, picRect.y, 1.0f);
					glPushMatrix();
					glTranslatef(0.0f, r.y, 0.0f);
					glBegin(GL_LINES);
					glVertex2f(0.0f, 0.0f);
					glVertex2f(picRect.width, 0.0f);
					glEnd();
					glPopMatrix();
					glTranslatef(r.x, 0.0f, 0.0f);
					glBegin(GL_LINES);
					glVertex2f(0.0f, 0.0f);
					glVertex2f(0.0f, picRect.height);
					glEnd();
					glPopMatrix();
					glPopAttrib();
				}
			}
		}
		if (drawingStatus == 2 || drawingStatus == 3 || drawingStatus == 4 || drawingStatus == 6 || drawingStatus == 7)
		{
			BOX3D b = box;
			b.sort();
			glPushMatrix();
			glScalef(DX, DY, DZ);
			glTranslatef(0.5f, 0.5f, 0.5f);
			glPushAttrib(GL_ALL_ATTRIB_BITS);
			glColor3f(0.0f, 1.0f, 0.0f);
			glEnable(GL_LINE_STIPPLE);
			glLineStipple(1, 0x5555);
			glBegin(GL_LINE_LOOP);
			glVertex3f(b.x1, b.y1, b.z1);
			glVertex3f(b.x2, b.y1, b.z1);
			glVertex3f(b.x2, b.y2, b.z1);
			glVertex3f(b.x1, b.y2, b.z1);
			glEnd();
			glBegin(GL_LINE_LOOP);
			glVertex3f(b.x1, b.y1, b.z2);
			glVertex3f(b.x2, b.y1, b.z2);
			glVertex3f(b.x2, b.y2, b.z2);
			glVertex3f(b.x1, b.y2, b.z2);
			glEnd();
			glBegin(GL_LINES);
			glVertex3f(b.x1, b.y1, b.z1);
			glVertex3f(b.x1, b.y1, b.z2);
			glVertex3f(b.x2, b.y1, b.z1);
			glVertex3f(b.x2, b.y1, b.z2);
			glVertex3f(b.x2, b.y2, b.z1);
			glVertex3f(b.x2, b.y2, b.z2);
			glVertex3f(b.x1, b.y2, b.z1);
			glVertex3f(b.x1, b.y2, b.z2);
			glEnd();
			glDisable(GL_LINE_STIPPLE);
			glPopAttrib();
			glPopMatrix();
		}
		
		
	}
	//------------------------------------End 3D drawing-----------------------------------

	//------------------------------------Start 2D drawing---------------------------------
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glColor3f(1, 1, 1);

	// show the 3d matrix boundary
	if (b_showBox) 
	{
		glPushMatrix();
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glTranslatef(picRect.x, picRect.y, 0.0f);
		glColor3f(1.0f, 1.0f, 1.0f);
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
	}

	//draw the legend
	{
		FloatRect bf = getBestFitRect(wid, glRect);
		float len = 0.08f*glRect.width;
		float spacek = 0.2f;
		char ch[10];
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glLineWidth(2.0f);
		glPushMatrix();
		//pos it to the right side
		glTranslatef(glRect.width*0.8f, glRect.height*0.9f, 0.0f);
		//let's amplify it 
		glScalef(len, len, len);
		int nl = NISOLine;
		for (int i = 0; i < nl; ++i)
		{
			//float kc = float(nl - 1 - i) / float(gp->NISOLevel);
			float kc = float(nl - i) / float(NISOLevel + 1)*stdColor;
			float r, g, b;
			colorGradient(r, g, b, kc, colorScheme);
			glColor3f(r, g, b);
			glBegin(GL_LINES);
			glVertex2f(0.0f, 0.0f);
			glVertex2f(1.0f, 0.0f);
			glEnd();
			float percent = float(nl - i) / float(NISOLevel + 1);
			sprintf(ch, "%.4g", percent*stdValue);

			glColor3f(1.0f, 1.0f, 0.0f);
			glText.glText2f(1.1f, 0.0f, ch);
			glTranslatef(0.0f, -0.3f, 0.0f);
		}
		glPopMatrix();
		glPopAttrib();
	}
	
	// draw yellow box indicating which window is currently active
	if (wid == activeID) 
	{
		glColor3f(1.0f, 1.0f, 0.0f);
		glBegin(GL_LINE_LOOP);
		glVertex2f(0, 0);
		glVertex2f(0, glRect.height);
		glVertex2f(glRect.width, glRect.height);
		glVertex2f(glRect.width, 0);
		glEnd();
	}
	//-------------------------------------End 2D drawing----------------------------------

	//drawing the axis
	{
		glDisable(GL_DEPTH_TEST);
		glLoadIdentity();
		glColor3f(1, 1, 1);
		// translate the origin to the upper left corner
		glTranslatef(0.1*glRect.width, 0.8*glRect.height, 0);

		// change the look at angle
		gluLookAt(R*pv[wid].eyev[0], R*pv[wid].eyev[1], R*pv[wid].eyev[2], 0.0f, 0.0f, 0.0f, pv[wid].upv[0], pv[wid].upv[1], pv[wid].upv[2]);
		//draw axis
		float axisScale = min(glRect.width, glRect.height)*0.1f;
		glScalef(axisScale, axisScale, axisScale);
		for (int i = 0; i < 3; ++i)
		{
			glPushMatrix();
			if (XAXIS == i)
			{
				glColor3f(1, 1, 0);
				if (wid != i) glText.glText3f(1.35f, 0, 0, "x");
				glRotatef(90, 0, 1, 0);
			}
			else if (YAXIS == i)
			{
				glColor3f(0, 1, 0);
				if (wid != i) glText.glText3f(0, 1.35f, 0, "y");
				glRotatef(-90, 1, 0, 0);
			}
			else
			{
				glColor3f(1, 0, 0);
				if (wid != i) glText.glText3f(0, 0, 1.35f, "z");
			}
			glBegin(GL_LINES);
			glVertex3f(0, 0, 0);
			glVertex3f(0, 0, 1);
			glEnd();
			glTranslatef(0, 0, 1);
			gluCylinder(gluQuad, 0.02, 0, 0.2, 4, 4);
			glPopMatrix();
		}
	}
}
void View3D::draw3D() // function to draw the 4th sub-window
{
	// change the look at angle
	double R = 1000;
	double ks1 = sin(theta*PI / 180.0), ks2 = sin(phi*PI / 180.0), kc1 = cos(theta*PI / 180.0), kc2 = cos(phi*PI / 180.0);
	gluLookAt(R*ks1*kc2, R*ks1*ks2, R*kc1, 0.0f, 0.0f, 0.0f, -kc1*kc2, -kc1*ks2, ks1);
	glRotatef(-90, 1, 0, 0);
	glRotatef(-90, 0, 1, 0);

	// need to calculate the ratio to scale from cm to pixel
	float scale = getScale(VIEW3D);
	glScalef(scale, scale, scale);

	// --------Note the rest plot use 3D world coordinates now-------- //
	double lx = NX*DX;
	double ly = NY*DY;
	double lz = NZ*DZ;
	glTranslatef(-lx / 2, -ly / 2, -lz / 2);

	glEnable(GL_POLYGON_OFFSET_LINE);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glPolygonOffset(-1.0f, -1.0f);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA);

	//glEnable(GL_LIGHTING);
	//	glEnable(GL_LIGHT0);
	GLfloat lightAmbient[] = { 1.0, 1.0, 1.0, 1.0 };
	// 	GLfloat lightDiffuse[] = { 1.0, 1.0, 1.0, 1.0 };
	// 	GLfloat lightSpecular[] = { 0.5, 0.5, 0.5, 1.0 };
	// 	glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
	// 	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
	// 	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
	// 	GLfloat lightPosition[] = { 0, 0, 10.0, 0.0 };   // w=0.0
	// 	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
	//glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lightAmbient);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glPushMatrix();
	//glTranslatef(xcenter, ycenter, zcenter);
	SFloat zlim = cross[ZAXIS] * DZ + zmin;
	SFloat ylim = cross[YAXIS] * DY + ymin;
	SFloat xlim = cross[XAXIS] * DX + xmin;
#define canRender3D(i) (tetLim.a(i, 4) > zlim &&tetLim.a(i, 2) > ylim &&tetLim.a(i, 1) < xlim)
	for (int i = 0; i < NTET; ++i)
	{
		if (canRender3D(i))
		{
			int v0 = elem->a(i, 0);
			int v1 = elem->a(i, 1);
			int v2 = elem->a(i, 2);
			int v3 = elem->a(i, 3);
			int mc = elem->a(i, 4);
			if (mc >= 7)
			{
				int im = mc % 7;
				mc = mc - 7 * im;
			}
			//judge if this face is adjacent to any face
			int iad = f2t.a(i, 0);
			if (iad == -1 || !canRender3D(iad)) //surface triangle
			{
				if (b_showTriangleEdges)
				{
					glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
					glColor4f(MatColor[mc][0], MatColor[mc][1], MatColor[mc][2], 0.5f);
					glBegin(GL_POLYGON);
					glVertex3f(node->a(v1, 0), node->a(v1, 1), node->a(v1, 2));
					glVertex3f(node->a(v2, 0), node->a(v2, 1), node->a(v2, 2));
					glVertex3f(node->a(v3, 0), node->a(v3, 1), node->a(v3, 2));
					glEnd();
				}
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				glColor3f(bgColor.a(i, 0), bgColor.a(i, 1), bgColor.a(i, 2));
				//glBegin(GL_TRIANGLES);
				glBegin(GL_POLYGON);
				//glNormal3fv((GLfloat*)&fNorm.a(i, 0));
				glVertex3f(node->a(v1, 0), node->a(v1, 1), node->a(v1, 2));
				glVertex3f(node->a(v2, 0), node->a(v2, 1), node->a(v2, 2));
				glVertex3f(node->a(v3, 0), node->a(v3, 1), node->a(v3, 2));
				glEnd();
			}

			iad = f2t.a(i, 1);
			if (iad == -1 || !canRender3D(iad)) //surface triangle
			{
				if (b_showTriangleEdges)
				{
					glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
					glColor4f(MatColor[mc][0], MatColor[mc][1], MatColor[mc][2], 0.5f);
					glBegin(GL_POLYGON);
					glVertex3f(node->a(v0, 0), node->a(v0, 1), node->a(v0, 2));
					glVertex3f(node->a(v3, 0), node->a(v3, 1), node->a(v3, 2));
					glVertex3f(node->a(v2, 0), node->a(v2, 1), node->a(v2, 2));
					glEnd();
				}
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				glColor3f(bgColor.a(i, 0), bgColor.a(i, 1), bgColor.a(i, 2));
				//glBegin(GL_TRIANGLES);
				glBegin(GL_POLYGON);
				//glNormal3fv((GLfloat*)&fNorm.a(i, 1));
				glVertex3f(node->a(v0, 0), node->a(v0, 1), node->a(v0, 2));
				glVertex3f(node->a(v3, 0), node->a(v3, 1), node->a(v3, 2));
				glVertex3f(node->a(v2, 0), node->a(v2, 1), node->a(v2, 2));
				glEnd();
			}

			iad = f2t.a(i, 2);
			if (iad == -1 || !canRender3D(iad)) //surface triangle
			{
				if (b_showTriangleEdges)
				{
					glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
					glColor4f(MatColor[mc][0], MatColor[mc][1], MatColor[mc][2], 0.5f);
					glBegin(GL_POLYGON);
					glVertex3f(node->a(v0, 0), node->a(v0, 1), node->a(v0, 2));
					glVertex3f(node->a(v1, 0), node->a(v1, 1), node->a(v1, 2));
					glVertex3f(node->a(v3, 0), node->a(v3, 1), node->a(v3, 2));
					glEnd();
				}
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				glColor3f(bgColor.a(i, 0), bgColor.a(i, 1), bgColor.a(i, 2));
				//glBegin(GL_TRIANGLES);
				glBegin(GL_POLYGON);
				//glNormal3fv((GLfloat*)&fNorm.a(i, 2));
				glVertex3f(node->a(v0, 0), node->a(v0, 1), node->a(v0, 2));
				glVertex3f(node->a(v1, 0), node->a(v1, 1), node->a(v1, 2));
				glVertex3f(node->a(v3, 0), node->a(v3, 1), node->a(v3, 2));
				glEnd();
			}

			iad = f2t.a(i, 3);
			if (iad == -1 || !canRender3D(iad)) //surface triangle
			{
				if (b_showTriangleEdges)
				{
					glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
					glColor4f(MatColor[mc][0], MatColor[mc][1], MatColor[mc][2], 0.5f);
					glBegin(GL_POLYGON);
					glVertex3f(node->a(v0, 0), node->a(v0, 1), node->a(v0, 2));
					glVertex3f(node->a(v2, 0), node->a(v2, 1), node->a(v2, 2));
					glVertex3f(node->a(v1, 0), node->a(v1, 1), node->a(v1, 2));
					glEnd();
				}
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				glColor3f(bgColor.a(i, 0), bgColor.a(i, 1), bgColor.a(i, 2));
				//glBegin(GL_TRIANGLES);
				glBegin(GL_POLYGON);
				//glNormal3fv((GLfloat*)&fNorm.a(i, 3));
				glVertex3f(node->a(v0, 0), node->a(v0, 1), node->a(v0, 2));
				glVertex3f(node->a(v2, 0), node->a(v2, 1), node->a(v2, 2));
				glVertex3f(node->a(v1, 0), node->a(v1, 1), node->a(v1, 2));
				glEnd();
			}
		}
	}
	glPopMatrix();
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_POLYGON_OFFSET_LINE);
	glDisable(GL_CULL_FACE);
	glDisable(GL_BLEND);
	glClear(GL_DEPTH_BUFFER_BIT);

	int gl_render_mode;
	glGetIntegerv(GL_RENDER_MODE, &gl_render_mode);
	if (gl_render_mode == GL_SELECT || b_showBGImg)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		//glBlendFunc(GL_ZERO, GL_ONE);
		//prepare for picking
		glInitNames();
		glPushName(VIEW3D);

		//draw each phantom and dose wash
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_TEXTURE_2D);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		for (int id = 0; id < 3; ++id)
		{
			glBindTexture(GL_TEXTURE_2D, pv[id].gltexture);
			if (ZAXIS == id)
			{
				glLoadName(ZAXIS);
				float z = cross[ZAXIS] * DZ;
				glBegin(GL_QUADS);
				glTexCoord2f(1.0, 0.0);
				glVertex3f(lx, 0, z);
				glTexCoord2f(0.0, 0.0);
				glVertex3f(0, 0, z);
				glTexCoord2f(0.0, 1.0);
				glVertex3f(0, ly, z);
				glTexCoord2f(1.0, 1.0);
				glVertex3f(lx, ly, z);
				glEnd();
			}
			else if (YAXIS == id)
			{
				glLoadName(YAXIS);
				float y = cross[YAXIS] * DY;
				glBegin(GL_QUADS);
				glTexCoord2f(1.0, 1.0);
				glVertex3f(lx, y, lz);
				glTexCoord2f(0.0, 1.0);
				glVertex3f(0, y, lz);
				glTexCoord2f(0.0, 0.0);
				glVertex3f(0, y, 0);
				glTexCoord2f(1.0, 0.0);
				glVertex3f(lx, y, 0);
				glEnd();
			}
			else if (XAXIS == id)
			{
				glLoadName(XAXIS);
				float x = cross[XAXIS] * DX;
				glBegin(GL_QUADS);
				glTexCoord2f(1.0, 0.0);
				glVertex3f(x, 0, lz);
				glTexCoord2f(1.0, 1.0);
				glVertex3f(x, ly, lz);
				glTexCoord2f(0.0, 1.0);
				glVertex3f(x, ly, 0);
				glTexCoord2f(0.0, 0.0);
				glVertex3f(x, 0, 0);
				glEnd();
			}
		}
		glLoadName(VIEW3D);//give other irrelevant object name
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_BLEND);
	}
	
	//store matrix for 2D to 3D coordinates reverse projection
	
	if (GL_RENDER == gl_render_mode)
	{
		glGetIntegerv(GL_VIEWPORT, subw[VIEW3D].viewport);
		glGetDoublev(GL_MODELVIEW_MATRIX, subw[VIEW3D].modelview);
		glGetDoublev(GL_PROJECTION_MATRIX, subw[VIEW3D].projection);
	}
	else return;

	

	// drawing the cuboid to select a sub-matrix
	if (drawingType == DRAW_CUBOID && (drawingStatus == 2 || drawingStatus == 3 || drawingStatus == 4 || drawingStatus == 6 || drawingStatus == 7))
	{
		BOX3D b = box;
		b.sort();
		glPushMatrix();
		glScalef(DX, DY, DZ);
		glTranslatef(0.5f, 0.5f, 0.5f);
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glColor3f(0.0f, 1.0f, 0.0f);
		glEnable(GL_LINE_STIPPLE);
		glLineStipple(1, 0x5555);
		glBegin(GL_LINE_LOOP);
		glVertex3f(b.x1, b.y1, b.z1);
		glVertex3f(b.x2, b.y1, b.z1);
		glVertex3f(b.x2, b.y2, b.z1);
		glVertex3f(b.x1, b.y2, b.z1);
		glEnd();
		glBegin(GL_LINE_LOOP);
		glVertex3f(b.x1, b.y1, b.z2);
		glVertex3f(b.x2, b.y1, b.z2);
		glVertex3f(b.x2, b.y2, b.z2);
		glVertex3f(b.x1, b.y2, b.z2);
		glEnd();
		glBegin(GL_LINES);
		glVertex3f(b.x1, b.y1, b.z1);
		glVertex3f(b.x1, b.y1, b.z2);
		glVertex3f(b.x2, b.y1, b.z1);
		glVertex3f(b.x2, b.y1, b.z2);
		glVertex3f(b.x2, b.y2, b.z1);
		glVertex3f(b.x2, b.y2, b.z2);
		glVertex3f(b.x1, b.y2, b.z1);
		glVertex3f(b.x1, b.y2, b.z2);
		glEnd();
		glDisable(GL_LINE_STIPPLE);
		glPopAttrib();
		glPopMatrix();
	}

	//draw axis
	glPushMatrix();
	glScalef(lx / 2, lx / 2, lx / 2);
	for (int i = 0; i < 3; ++i)
	{
		glPushMatrix();
		if (XAXIS == i)
		{
			glColor3f(1, 1, 0);
			glText.glText3f(1.25f, 0, 0, "x");
			glRotatef(90, 0, 1, 0);
		}
		else if (YAXIS == i)
		{
			glColor3f(0, 1, 0);
			glText.glText3f(0, 1.25f, 0, "y");
			glRotatef(-90, 1, 0, 0);
		}
		else
		{
			glColor3f(1, 0, 0);
			glText.glText3f(0, 0, 1.25f, "z");
		}
		glBegin(GL_LINES);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, 1);
		glEnd();
		glTranslatef(0, 0, 1);
		gluCylinder(gluQuad, 0.02, 0, 0.2, 4, 4);
		glPopMatrix();
	}
	glPopMatrix();

	//draw the boundary of these slices
	for (int i = 0; i < 3; ++i) //draw each phantom and dose
	{
		glBegin(GL_LINE_LOOP);
		if (ZAXIS == i)
		{
			glColor3f(1, 0, 0);
			float z = cross[ZAXIS] * DZ;
			glVertex3f(lx, 0, z);
			glVertex3f(0, 0, z);
			glVertex3f(0, ly, z);
			glVertex3f(lx, ly, z);
		}
		else if (YAXIS == i)
		{
			glColor3f(0, 1, 0);
			float y = cross[YAXIS] * DY;
			glVertex3f(lx, y, 0);
			glVertex3f(0, y, 0);
			glVertex3f(0, y, lz);
			glVertex3f(lx, y, lz);
		}
		else if (XAXIS == i)
		{
			glColor3f(0, 0, 1);
			float x = cross[XAXIS] * DX;
			glVertex3f(x, ly, 0);
			glVertex3f(x, 0, 0);
			glVertex3f(x, 0, lz);
			glVertex3f(x, ly, lz);
		}
		glEnd();
	}

	glDisable(GL_DEPTH_TEST);

	glLoadIdentity();
	glColor3f(1, 1, 1);
	if (b_showBox) // show the 3d matrix boundary
	{
		FloatRect& picRect = subw[VIEW3D].picRect;
		glPushMatrix();
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glTranslatef(picRect.x, picRect.y, 0.0f);
		glColor3f(1.0f, 1.0f, 1.0f);
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
	}

	// draw yellow box indicating which window is currently active
	if (VIEW3D == activeID)
	{
		glColor3f(1.0f, 1.0f, 0.0f);
		glBegin(GL_LINE_LOOP);
		glVertex2f(0, 0);
		glVertex2f(0, subw[VIEW3D].glRect.height);
		glVertex2f(subw[VIEW3D].glRect.width, subw[VIEW3D].glRect.height);
		glVertex2f(subw[VIEW3D].glRect.width, 0);
		glEnd();
	}
// 	if (b_showCrossLine) // draw crossing lines
// 	{
// 		glPushAttrib(GL_ALL_ATTRIB_BITS);
// 		glColor3f(1.0f, 1.0f, 1.0f);
// 		glLineWidth(1.0f);
// 		glPushMatrix();
// 
// 		if (ZAXIS == wid)
// 		{
// 			glTranslatef(0, cross[YAXIS] * DY, cross[ZAXIS] * DZ);
// 			glBegin(GL_LINES); //draw horizontal line
// 			glVertex2f(0, 0);
// 			glVertex2f(lx, 0);
// 			glEnd();
// 			glPopMatrix();
// 
// 			glPushMatrix();
// 			glTranslatef(cross[XAXIS] * DX, 0, cross[ZAXIS] * DZ);
// 			glBegin(GL_LINES); //draw vertical line
// 			glVertex2f(0, 0);
// 			glVertex2f(0, ly);
// 			glEnd();
// 		}
// 		else if (YAXIS == wid) // I may need to change this ???
// 		{
// 			glTranslatef(0, cross[YAXIS] * DY, cross[ZAXIS] * DZ);
// 			glBegin(GL_LINES);
// 			glVertex3f(0, 0, 0);
// 			glVertex3f(lx, 0, 0);
// 			glEnd();
// 			glPopMatrix();
// 
// 			glPushMatrix();
// 			glTranslatef(cross[XAXIS] * DX, cross[YAXIS] * DY, 0);
// 			glBegin(GL_LINES);
// 			glVertex3f(0, 0, 0);
// 			glVertex3f(0, 0, lz);
// 			glEnd();
// 		}
// 		else if (XAXIS == wid)
// 		{
// 			glTranslatef(cross[XAXIS] * DX, 0, cross[ZAXIS] * DZ);
// 			glBegin(GL_LINES);
// 			glVertex3f(0, 0, 0);
// 			glVertex3f(0, ly, 0);
// 			glEnd();
// 			glPopMatrix();
// 
// 			glPushMatrix();
// 			glTranslatef(cross[XAXIS] * DX, cross[YAXIS] * DY, 0);
// 			glBegin(GL_LINES);
// 			glVertex3f(0, 0, 0);
// 			glVertex3f(0, 0, lz);
// 			glEnd();
// 		}
// 		glPopMatrix();
// 		glPopAttrib();
// 	}
	//------------------------------------End 3D drawing-----------------------------------
}

void View3D::calcSubWinRect(wxSize wxsz) // this can been override in its derived class. It will invoke the resized method for sub-windows
{
	int NW = wxsz.x;
	int NH = wxsz.y;
	if (b_maxGLWindow)
	{
		resized(activeID, FloatRect(0, 0, NW - 1, NH - 1));
		return;
	}
	int npx = (NW - 3) / 2;
	int npy = (NH - 3) / 2;
	int shiftx = NW - npx - 1;
	int shifty = NH - npy - 1;
	FloatRect rc(1, 1, npx - 1, npy - 1);
	FloatRect rct = rc;
	//it will case resize of each sub-window
	resized(ZAXIS, rct.shift(0, shifty));
	rct = rc;
	resized(YAXIS, rct.shift(shiftx, shifty));
	resized(XAXIS, rc);
	resized(VIEW3D, rc.shift(shiftx, 0));
}
int View3D::getSubCoordinates(wxPoint& p) //return which window it locates in. This can been override in its derived class
{
	wxRect rc = GetClientRect();
	p.y = rc.height - 1 - p.y; // change to left-bottom corner origin
	if (b_maxGLWindow) return activeID;
	else
	{
		for (int i = 0; i < 4; ++i)
		{
			if (subw[i].glRect.x <= p.x && p.x <= subw[i].glRect.Right() &&
				subw[i].glRect.y <= p.y && p.y <= subw[i].glRect.Up()) return i;
		}
		return -1;
	}
}
FloatRect View3D::getBestFitRect(int wid, FloatRect& rect) // get the best fitted region for drawing
{
	FloatRect bf;
	//calculate the image ratio height/width
	double imageRatio = 1;
	if (VIEW3D != wid)
	{
		if (ZAXIS == wid)
		{
			imageRatio = double(NY*DY) / double(NX*DX);
			if (pv[wid].upv[0] != 0) imageRatio = 1.0 / imageRatio;
		}
		else if (YAXIS == wid)
		{
			imageRatio = double(NZ*DZ) / double(NX*DX);
			if (pv[wid].upv[0] != 0) imageRatio = 1.0 / imageRatio;
		}
		else // XAXIS == wid
		{
			imageRatio = double(NY*DY) / double(NZ*DZ);
			if (pv[wid].upv[2] != 0) imageRatio = 1.0 / imageRatio;
		}
	}
	// get window ratio
	double windowRatio = double(rect.height) / double(rect.width);
	// get the best fit region
	if (windowRatio < imageRatio)
	{
		bf.x = (rect.width - rect.height / imageRatio) / 2;
		bf.width = rect.height / imageRatio;
		bf.y = 0;
		bf.height = rect.height;
	}
	else
	{
		bf.x = 0;
		bf.width = rect.width;
		bf.y = (rect.height - rect.width*imageRatio) / 2;
		bf.height = rect.width*imageRatio;
	}
	return bf;
}
SFloat View3D::getPointValue(SFloat x, SFloat y, SFloat z)
{
	return 0;
}
bool View3D::getLineData(ArrayMgr<double>& xin, ArrayMgr<double>& yin)
{
	if (mat == NULL&&bgImage == NULL) return false;
	SFloat vx = cObjX - pObjX;
	SFloat vy = cObjY - pObjY;
	SFloat vz = cObjZ - pObjZ;
	SFloat norm = sqrt(vx*vx + vy*vy + vz*vz);
	//get normalized direction vector
	vx /= norm;
	vy /= norm;
	vz /= norm;
	SFloat step = min(DX, min(DY, DZ));
	int NSample = int(norm / step);
	xin.resize(NSample);
	yin.resize(NSample);
	for (int i = 0; i < NSample; ++i)
	{
		xin.a(i) = i*step;
		yin.a(i) = getPointValue(pObjX + i*step*vx, pObjY + i*step*vy, pObjZ + i*step*vz);
	}
	return true;
}
//mouse and keyboard events for each window
void View3D::mouseDown(int wid, wxMouseEvent& event)
{
	b_mouseLeftDown = true;
	
	pMouse = event.GetPosition();
	cMouse = pMouse;
	pRect = subw[wid].picRect;
	// find out the coordinates we clicked in the 3D world
	const unsigned int BUFFER_LENGTH = 32;
	static GLuint selectBuffer[BUFFER_LENGTH];
	glSelectBuffer(BUFFER_LENGTH, selectBuffer); //set selection buffer
	glRenderMode(GL_SELECT); // switch the selecting mode
	render(wid);
	int hits = glRenderMode(GL_RENDER);
	GLuint zmin = selectBuffer[1];
	selectedSlice = VIEW3D;
	for (int i = 0; i < hits; ++i) //get the selected object
	{
		if (selectBuffer[4 * i + 3] != VIEW3D && selectBuffer[4 * i + 1] <= zmin)//need to compare depth
		{
			zmin = selectBuffer[4 * i + 1];
			selectedSlice = selectBuffer[4 * i + 3];
		}
	}
	if (selectedSlice != VIEW3D)
	{
		glReadPixels(cMouse.x, cMouse.y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &pWinZ);
		gluUnProject(cMouse.x, cMouse.y, pWinZ, subw[wid].modelview, subw[wid].projection, subw[wid].viewport, &pObjX, &pObjY, &pObjZ);
		getMatIndex(pObjX, pObjY, pObjZ, pPos);
	}
	else pObjX = pObjY = pObjZ = -1; // didn't clicked the object

	if (isObjectClicked() && drawingType == DRAW_NOTHING)
	{
		wxCommandEvent event_notify(VIEW3D_NOTIFY);
		event_notify.SetInt(VIEW3D_VOXEL_CLICKED);
		wxPostEvent(this, event_notify); // to View3D
		if (VIEW3D == activeID) // Module drag the 3D slice to adjust the crossing point
		{
			drawingType = DRAW_DRAG_3D_SLICE;
			//backup the original cross point
			pcross[0] = cross[0];
			pcross[1] = cross[1];
			pcross[2] = cross[2];
			renders();
			return;
		}
		else // in the first 3 windows,  Module drag the crossing point
		{
			int dx = pPos.ix - cross[XAXIS];
			int dy = pPos.iy - cross[YAXIS];
			int dz = pPos.iz - cross[ZAXIS];
			if (dx*dx + dy*dy + dz*dz <= 4 * 4 * 4)
			{
				drawingType = DRAW_ADJUST_CROSSPOINT;
				pcross[0] = cross[0];
				pcross[1] = cross[1];
				pcross[2] = cross[2];
				renders();
				return;
			}
		}
		
	}
	//Module drag the picture: just need above position
	//Module zoom the picture: just need above position

	//Module DRAW_LINE_ANY turned on externally
	if ((drawingType == DRAW_LINE_ANY || drawingType == DRAW_LINE_VERTICAL || drawingType == DRAW_LINE_HORIZONTAL) && activeID != VIEW3D)
	{
		if(isObjectClicked()&&activeID!=VIEW3D) drawingStatus = 1; // clicked the object
		return;
	}

	if (drawingType == DRAW_CUBOID && activeID != VIEW3D && isObjectClicked())
	{
		// first click
		if (drawingStatus == 0 || drawingStatus == 3 )
		{
			drawingStatus = 1; // clicked the object
			boxDirection = wid;
			box.x1 = box.x2 = pPos.ix;
			box.y1 = box.y2 = pPos.iy;
			box.z1 = box.z2 = pPos.iz;
			return;
		}
		
		// one plan has finished drawing, start drawing in the perpendicular direction
		if ((drawingStatus == 4 || drawingStatus == 7) && boxDirection != wid)
		{
			drawingStatus = 5; // clicked the object
			return;
		}
	}

	renders();
}
void View3D::mouseMoved(int wid, wxMouseEvent& event)
{
	cMouse = event.GetPosition();

	//Module drag the picture: ctrl + drag the picture + no other drawing
	if (b_mouseLeftDown && event.ControlDown() && isPictureClicked(wid) && drawingType == DRAW_NOTHING)
	{
		wxPoint shiftv = cMouse - pMouse;
		FloatRect temp = pRect;
		subw[wid].picRect = temp.shift(shiftv.x, shiftv.y);
		renders();
		return;
	}

	//Module adjusting the crossing point: drag at the crossing point
	if (b_mouseLeftDown && drawingType == DRAW_ADJUST_CROSSPOINT)
	{
		//get current object position
		gluUnProject(cMouse.x, cMouse.y, pWinZ, subw[wid].modelview, subw[wid].projection, subw[wid].viewport, &cObjX, &cObjY, &cObjZ);
		//INT3D cPos;
		getMatIndex(cObjX, cObjY, cObjZ, cPos);
		if (!cPos.withinMatrix(NX, NY, NZ)) return mouseReleased(wid, event);
		cPos -= pPos;
		cross[XAXIS] = pcross[XAXIS] + cPos.ix;
		cross[YAXIS] = pcross[YAXIS] + cPos.iy;
		cross[ZAXIS] = pcross[ZAXIS] + cPos.iz;
		show();
		wxCommandEvent event_notify(VIEW3D_NOTIFY);
		event_notify.SetInt(VIEW3D_CROSS_MOVED);
		wxPostEvent(this, event_notify); // to View3D
		return;
	}

	if (b_mouseLeftDown && (drawingType == DRAW_LINE_ANY || drawingType == DRAW_LINE_VERTICAL
		|| drawingType == DRAW_LINE_HORIZONTAL) && (drawingStatus == 1 || drawingStatus == 2))
	{
		drawingStatus = 2; // indicating mouse moving now
		//get current object position
		gluUnProject(cMouse.x, cMouse.y, pWinZ, subw[wid].modelview, subw[wid].projection, subw[wid].viewport, &cObjX, &cObjY, &cObjZ);
		getMatIndex(cObjX, cObjY, cObjZ, cPos);
		if (cPos.trim(NX, NY, NZ)) return mouseReleased(wid, event);
		renders();
		return;
	}

	if (b_mouseLeftDown && drawingType == DRAW_CUBOID)
	{
		if (drawingStatus == 1 || drawingStatus == 2)
		{
			drawingStatus = 2; // moving for the first plan drawing
			//get current object position
			gluUnProject(cMouse.x, cMouse.y, pWinZ, subw[wid].modelview, subw[wid].projection, subw[wid].viewport, &cObjX, &cObjY, &cObjZ);
			getMatIndex(cObjX, cObjY, cObjZ, cPos);
			if (cPos.trim(NX, NY, NZ)) return mouseReleased(wid, event);
			box.x2 = cPos.ix;
			box.y2 = cPos.iy;
			box.z2 = cPos.iz;
			if (XAXIS == wid)
			{
				box.x1 = 0;
				box.x2 = NX-1;
			}
			else if (YAXIS == wid)
			{
				box.y1 = 0;
				box.y2 = NY-1;
			}
			else if (ZAXIS == wid)
			{
				box.z1 = 0;
				box.z2 = NZ-1;
			}

			renders();
			wxCommandEvent event_notify(VIEW3D_NOTIFY);
			event_notify.SetInt(VIEW3D_CUBOID_SAMPLING_CHANGED);
			wxPostEvent(this, event_notify); // to View3D
			return;
		}

		if (drawingStatus == 5 || drawingStatus == 6)
		{
			drawingStatus = 6; // moving for the perpendicular plan drawing
			//get current object position
			gluUnProject(cMouse.x, cMouse.y, pWinZ, subw[wid].modelview, subw[wid].projection, subw[wid].viewport, &cObjX, &cObjY, &cObjZ);
			getMatIndex(cObjX, cObjY, cObjZ, cPos);
			if (cPos.trim(NX, NY, NZ)) return mouseReleased(wid, event);
			if (XAXIS == boxDirection)
			{
				box.x1 = pPos.ix;
				box.x2 = cPos.ix;
			}
			else if (YAXIS == boxDirection)
			{
				box.y1 = pPos.iy;
				box.y2 = cPos.iy;
			}
			else if (ZAXIS == boxDirection)
			{
				box.z1 = pPos.iz;
				box.z2 = cPos.iz;
			}

			renders();
			wxCommandEvent event_notify(VIEW3D_NOTIFY);
			event_notify.SetInt(VIEW3D_CUBOID_SAMPLING_CHANGED);
			wxPostEvent(this, event_notify); // to View3D
			return;
		}
	}
	if (!b_mouseLeftDown && drawingType == DRAW_CUBOID && (drawingStatus == 0 || drawingStatus == 4 )) //need to indicate where's the mouse
	{
		FloatRect& glRect = subw[wid].glRect;
		FloatRect& picRect = subw[wid].picRect;
		wxPoint ov(glRect.x, glRect.y);
		wxPoint r = cMouse - ov;
		if (picRect.within(r.x, r.y))
		{
			if (drawingStatus == 4 && wid == boxDirection) return;
			renders();

			return;
		}
	}


	if (b_mouseLeftDown && drawingType == DRAW_DRAG_3D_SLICE)
	{
		//get current object position
		gluUnProject(cMouse.x, cMouse.y, pWinZ, subw[wid].modelview, subw[wid].projection, subw[wid].viewport, &cObjX, &cObjY, &cObjZ);
		getMatIndex(cObjX, cObjY, cObjZ, cPos);
		int dnx = 0, dny = 0, dnz = 0;
		if (ZAXIS == selectedSlice) dnz = (int)round((cObjZ - pObjZ) / DZ);
		else if (YAXIS == selectedSlice) dny = (int)round((cObjY - pObjY) / DY);
		else if (XAXIS == selectedSlice) dnx = (int)round((cObjX - pObjX) / DX);
		
		if (abs(dnx) > 0 || abs(dny) > 0 || abs(dnz) > 0)
		{
			dnx += pcross[XAXIS];
			dny += pcross[YAXIS];
			dnz += pcross[ZAXIS];
			if (0 <= dnx&&dnx < NX && 0 <= dny&&dny < NY && 0 <= dnz&&dnz < NZ)
			{
				cross[XAXIS] = dnx;
				cross[YAXIS] = dny;
				cross[ZAXIS] = dnz;
				//show();

				wxCommandEvent event_notify(VIEW3D_NOTIFY);
				event_notify.SetInt(VIEW3D_CROSS_MOVED);
				wxPostEvent(this, event_notify); // to View3D
				show();
			}
		}
		return;
	}
}
void View3D::mouseReleased(int wid, wxMouseEvent& event)
{
	if (b_mouseLeftDown == false) return;
	b_mouseLeftDown = false;

	//Module adjusting the crossing point:
	if (drawingType == DRAW_ADJUST_CROSSPOINT)
	{
		drawingType = DRAW_NOTHING;
		return;
	}

	//Module drag the slice in the 3D view:
	if (drawingType == DRAW_DRAG_3D_SLICE)
	{
		drawingType = DRAW_NOTHING;
		return;
	}

	if (drawingType == DRAW_LINE_ANY || drawingType == DRAW_LINE_VERTICAL || drawingType == DRAW_LINE_HORIZONTAL)
	{
		if (drawingStatus == 2)
		{
			drawingStatus = 3; // finished drawing the line. The line should be drew when drawingStatus == 2||3
			// notify the window that line sampling is finished. This message may be used to draw a sample graph
			wxCommandEvent event_notify(VIEW3D_NOTIFY);
			event_notify.SetInt(VIEW3D_LINE_SAMPLING_FINISHED);
			wxPostEvent(this, event_notify); // to View3D
			if (b_showSampledLine)
			{
				ArrayMgr<double> xin, yin;
				getLineData(xin, yin);
				FIGURE& fig = FIGURE::Figure(this, false);
				
				fig.Plot(MG(xin), MG(yin));
				fig.g->Label('x', "x", 0);
				fig.g->Label('y', "y", 0);
				fig.g->Title("line sampling");
				fig.renders();
			}
		}
		else if (drawingStatus == 3) // do nothing
		{
			renders();
		}
		else
		{
			drawingStatus = 0;
			renders();
		}
		return;
	}

	if (drawingType == DRAW_CUBOID)
	{
		if (drawingStatus == 2 || drawingStatus == 3) drawingStatus = 3; // finished drawing in the first plan
		else if (drawingStatus == 6 || drawingStatus == 7)
		{
			drawingStatus = 7;
			wxCommandEvent event_notify(VIEW3D_NOTIFY);
			event_notify.SetInt(VIEW3D_CUBOID_SAMPLING_FINISHED);
			wxPostEvent(this, event_notify); // to View3D
		}
		else
		{
			drawingStatus = 0;
			renders();
		}
		return;
	}
}
void View3D::mouseLeftDouble(int wid, wxMouseEvent& event)
{
	if (b_maxGLWindow) set4WinViews();
	else setMaxView(wid);
}
void View3D::keyPressed(int wid, wxKeyEvent& event)
{
	if (VIEW3D == wid)
	{
		switch (event.GetKeyCode())
		{
		case 'A':
			phi += 5;
			break;
		case 'D':
			phi -= 5;
			break;
		case 'W':
			theta += 5;
			break;
		case 'S':
			theta -= 5;
			break;
		default:
			break;
		}
		renders();
	}
	else
	{
		switch (event.GetKeyCode())
		{
		case 'A':
			moveLeft();
			break;
		case 'D':
			moveRight();
			break;
		case 'W':
			moveUp();
			break;
		case 'S':
			moveDown();
			break;
		default:
			break;
		}
	}
}
void View3D::rightClick(int wid, wxMouseEvent& event)
{

}
void View3D::mouseLeftWindow(int wid, wxMouseEvent& event)
{
	mouseReleased(wid, event); // simulate releasing the left mouse button
}
void View3D::mouseEnterWindow(int wid, wxMouseEvent& event)
{

}
void View3D::mouseWheelMoved(int wid, wxMouseEvent& event)
{
	//Module zoom the picture: ctrl + scrolling + no other drawing
	if (event.ControlDown() && drawingType == DRAW_NOTHING) //enlarge or shrink the image
	{
		if (event.GetWheelRotation() > 0) zoomIn();
		else zoomOut();
		return;
	}

	//Module change the slice view
	if (wid < VIEW3D && drawingType == DRAW_NOTHING)//scrolled in the first 3 windows
	{
		//move the slice in id direction
		int add = 1;
		if (event.GetWheelRotation() < 0) previousSlice();
		else nextSlice();
		return;
	}
}
void View3D::resized(int wid, FloatRect newGLRect) // the size of the sub-window has been changed, need to recalculate picRect
{
	if (subw[wid].picRect.empty()) // do the first time calculation
	{
		subw[wid].picRect = getBestFitRect(wid, newGLRect);
		subw[wid].picRect.zoom(0.95f); // make it a little smaller
	}
	else // has been calculated. Need to zoom it
	{
		FloatRect& glRect = subw[wid].glRect;
		FloatRect& picRect = subw[wid].picRect;
		//resize picRect if it exist
		float cx = glRect.width*0.5;
		float cy = glRect.height*0.5;
		float new_cx = newGLRect.width*0.5;
		float new_cy = newGLRect.height*0.5;

		FloatRect bf = getBestFitRect(wid, glRect);
		FloatRect newbf = getBestFitRect(wid, newGLRect);
		float sf = newbf.width / bf.width; // scaling factor
		picRect.x = new_cx + sf*(picRect.x - cx);
		picRect.y = new_cy + sf*(picRect.y - cy);
		picRect.width *= sf;
		picRect.height *= sf;
	}
	subw[wid].glRect = newGLRect; // update glRect anyway
}