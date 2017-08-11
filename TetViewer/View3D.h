#ifndef _VIEW3D_H_
#define _VIEW3D_H_
#pragma warning(disable: 4800)

#define THETA_INITIAL 60
#define PHI_INITIAL 50

wxDECLARE_EVENT(VIEW3D_NOTIFY, wxCommandEvent); //give the reference of the customized event const

void colorGradient(float& r, float& g, float&b, float kc, int scheme = 0);

class View3D;
class FigureFrame;

struct MarkPoint
{
	MarkPoint(double xi, double yi, double zi, bool shown = false, float radiusi = 0.5, uint8_t r = 255, uint8_t g = 0, uint8_t b = 0)
	{
		x = xi;
		y = yi;
		z = zi;
		b_show = shown;
		radius = radiusi;
		color[0] = r;
		color[1] = g;
		color[2] = b;
	}
	MarkPoint()
	{
		x = 0;
		y = 0;
		z = 0;
		b_show = false;
		radius = 0.5;
		color[0] = 255;
		color[1] = 0;
		color[2] = 0;
	}
	double x, y, z; // it's relative to the cuboid origin
	bool b_show;
	float radius; // default 0.5 cm
	uint8_t color[3]; // default red
};

class GLText // must call init() before using any text function
{
public:
	GLText() :gllist_id(0){}
	void init()
	{
		gllist_id = glGenLists(256);
#ifdef WIN32
		//warning this is windows specific only; need to find linux implementation
		wglUseFontBitmaps(wglGetCurrentDC(), 0, 256, gllist_id);
#else
		Display *dpy;          /* Our current X display */
		XFontStruct *fontInfo; /* Our font info */
		dpy = XOpenDisplay( NULL );

		/* Get the font information */

		fontInfo = XLoadQueryFont( dpy, "fixed" );
		/* If that font doesn't exist, something is wrong */
		if ( fontInfo == NULL )
		{
			fprintf( stderr, "no X font available?\n" );
			exit(2);
		}

		/* generate the list */
		glXUseXFont( fontInfo->fid, 0, 256, gllist_id );

		/* Recover some memory */
		XFreeFont( dpy, fontInfo );

		/* close the display now that we're done with it */
		XCloseDisplay( dpy );
#endif
	}
	void glText2f(float x, float y, const char* textstring)
	{
		glRasterPos2f(x, y);
		glListBase(gllist_id);
		glCallLists(strlen(textstring), GL_UNSIGNED_BYTE, (const GLvoid*)textstring);
	}
	void glText3f(float x, float y, float z, const char* textstring)
	{
		glRasterPos3f(x, y, z);
		glListBase(gllist_id);
		glCallLists(strlen(textstring), GL_UNSIGNED_BYTE, (const GLvoid*)textstring);
	}
	int gllist_id;
};

#define GlueF(a) a##f
class Vector {
public:

	SFloat x;  //!< x-coordinate
	SFloat y;  //!< y-coordinate
	SFloat z;  //!< z-coordinate

	/*! \brief Default constructor */
	Vector() : x(0), y(0), z(0) {}

	/*! \brief Construct from input values \a aX, \a aY and \a aZ */
	Vector(SFloat aX, SFloat aY, SFloat aZ) : x(aX), y(aY), z(aZ) {}

	/*! \brief Assignment operator */
	Vector &operator=(const Vector &v) { x = v.x; y = v.y; z = v.z; return *this; };

	/*! \brief Add \a v to this vector and return the result */
	Vector operator+(const Vector &v) const { return Vector(x + v.x, y + v.y, z + v.z); };

	/*! \brief Add vector \a v and assign the result to this vector */
	Vector &operator+=(const Vector &v) { x += v.x; y += v.y; z += v.z; return *this; };

	/*! \brief Returns the scalar product of vector \a v and this vector */
	SFloat operator*(const Vector &v) const { return x*v.x + y*v.y + z*v.z; };

	/*! \brief Returns the cross product of vector \a v and this vector */
	Vector operator&(const Vector &v) const { return Vector(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x); };

	/*! \brief Multiply all co-ordinates with \a f and assign the result to this vector */
	Vector &operator*=(SFloat f) { x *= f; y *= f; z *= f; return *this; };

	/*! \brief Multiply \a v with \a f and return the result */
	friend Vector operator*(SFloat f, Vector &v) { return v*f; };

	/*! \brief Subtract the vector \a v from this vector and return the result */
	Vector operator-(const Vector &v) const { return Vector(x - v.x, y - v.y, z - v.z); };

	/*! \brief Subtract \a v and assign the result to this vector */
	Vector &operator-=(const Vector &v) { x -= v.x; y -= v.y; z -= v.z; return *this; };

	/*! \brief Return the multiplication of this vector with \a f */
	Vector operator*(SFloat f) const { return Vector(x*f, y*f, z*f); };

	/*! \brief Get the squared length of the vector */
	SFloat getLengthSquared() const { return x*x + y*y + z*z; };

	/*! \brief Get the length of the vector */
	SFloat getLength() const { return GlueF(sqrt)(getLengthSquared()); };

	/*! \brief Normalize the vector to unit length */
	void normalizeToUnitLength() {
		SFloat length2 = getLengthSquared();
		if (length2) {
			SFloat invLength = 1 / GlueF(sqrt)(length2);
			x *= invLength; y *= invLength; z *= invLength;
		}
	};

	/*! \brief Rotate the vector by the polar angle Theta and the azimutha angle Phi */
	void rotate(SFloat cosTheta, SFloat sinTheta, SFloat cosPhi, SFloat sinPhi) {
		SFloat sinz = x*x + y*y;
		if (sinz > 1e-20) {
			sinz = GlueF(sqrt)(sinz);
			SFloat c1 = sinTheta / sinz;
			SFloat cphi = z*cosPhi;
			SFloat cx = x*cosTheta, cy = y*cosTheta;
			SFloat cx1 = cphi*x - y*sinPhi;
			SFloat cy1 = cphi*y + x*sinPhi;
			x = c1*cx1 + cx; y = c1*cy1 + cy;
			z = z*cosTheta - sinz*sinTheta*cosPhi;
		}
		else { x = sinTheta*cosPhi; y = sinTheta*sinPhi; z *= cosTheta; }
	};
};

//this is a subclass of wxGLCanvas, and is designed to provide 3D view of matrix mat
//and it may contain a background image, which will be rendered by pixel
class View3D: public wxGLCanvas
{
	// in this class, I will define some drawing process, starting from 1. In derived class, you can define other drawing processes
	// with number greater than DRAW_NEW_PROCESS = 10
	enum DRAWTYPE { DRAW_NOTHING = 0, DRAW_ADJUST_CROSSPOINT, DRAW_DRAG_3D_SLICE, DRAW_CUBOID, DRAW_LINE_ANY, DRAW_LINE_VERTICAL, DRAW_LINE_HORIZONTAL, DRAW_NEW_PROCESS = 10 };
	/*method*/
public:

	//Interface functions
	View3D(wxWindow* parent); //only need the parent windows pointer to create a window
	virtual ~View3D();
	void setMesh(ArrayMgr<SFloat>* node, ArrayMgr<int>* elem, ArrayMgr<int>* f2tm = NULL);
	void setMatrix(ArrayMgr<SFloat>* mat); // no dimension checking!
	void setBGImage(ArrayMgr<SFloat>* bgImage);
	void setBGPainter(ArrayMgr<RGBPixel>* bgPainter);// (0,0,0) means not to replace
	void setVoxelSize(int NX, int NY, int NZ, double DX, double DY, double DZ);
	bool initShow(bool b_reinit = false);// only call after all parameters has been set
	
	
	void renders(); // only do the rendering job without doing any necessary adjustment. Used in OnPaint() function
	void show();   // do necessary recalculation and then call renders()

	// events
	void OnPaint(wxPaintEvent& evt);
	void mouseMoved(wxMouseEvent& event);
	void mouseDown(wxMouseEvent& event);
	void mouseLeftDouble(wxMouseEvent& event);
	void mouseWheelMoved(wxMouseEvent& event);
	void mouseReleased(wxMouseEvent& event);
	void rightClick(wxMouseEvent& event);
	void mouseEnterWindow(wxMouseEvent& event);
	void mouseLeftWindow(wxMouseEvent& event);
	void keyPressed(wxKeyEvent& event);
	void keyReleased(wxKeyEvent& event);
	void resized(wxSizeEvent& evt);

	void OnView3DNotify(wxCommandEvent& event);

	//interface
	void showBox(unsigned int b = 2){ if (b > 1) b_showBox = !b_showBox; else b_showBox = b; renders(); }//change whether to draw the box
	void showBGImage(unsigned int b = 2){ if (b > 1) b_showBGImg = !b_showBGImg; else b_showBGImg = b; renders(); }//change whether to draw the box
	void showMat(unsigned int b = 2){ if (b > 1) b_showMat = !b_showMat; else b_showMat = b; renders(); }//change whether to show the mat
	void showGlobalContrast(unsigned int b = 2){ if (b > 1) b_globalContrast = !b_globalContrast; else b_globalContrast = b; show(); }
	void showTriangleEdges(unsigned int b = 2){ if (b > 1) b_showTriangleEdges = !b_showTriangleEdges; else b_showTriangleEdges = b; renders(); }//change whether to show the mat

	bool isMaxView(){ return b_maxGLWindow; }
	void setMaxView(int wid)
	{
		if (wid < 0 || wid>3) return;
		if (b_maxGLWindow && activeID == wid) return;
		b_maxGLWindow = true;
		activeID = wid;
		calcSubWinRect(GetSize());
		renders();
	}
	void set4WinViews()
	{
		if (b_maxGLWindow)
		{
			b_maxGLWindow = false;
			calcSubWinRect(GetSize());
			renders();
		}
	}
	//void setActiveID(int aid){ activeID = aid; }
	void nextSlice();
	void previousSlice();
	void setISOLevel(int NISOLevel) { this->NISOLevel = NISOLevel; show(); }
	void jump2Slice(int index) //this index starts from 0
	{	
		if (activeID < VIEW3D)
		{
			int NSlice = getNSlice(activeID);
			if (index < 0) index = 0;
			if (index >= NSlice) index = NSlice - 1;
			cross[activeID] = index;
			show();
		}
	}
	void addMarkPoint(const string& name, MarkPoint point, bool add = true) //if add == false, try to delete this mark
	{
		if (add) //should validate the position first
		{
			if (point.x <= 0 || point.x >= DX*NX || point.y <= 0 || point.y >= DY*NY || point.z <= 0 || point.z >= DZ*NZ) return;
			markPoints[name] = point;
		}
		else markPoints.erase(name);
	}
	void showMarkPoint(unsigned int b = 2){ if (b > 1) b_showMarkPoint = !b_showMarkPoint; else b_showMarkPoint = b; renders(); }
	
	void zoomIn();
	void zoomOut();
	void changeView(int wid) // change the active view to wid
	{
		if (!b_maxGLWindow)
		{
			activeID = wid;
			renders();
		}
		else setMaxView(wid);
	}
	void resetView();
	void moveLeft();
	void moveRight();
	void moveUp();
	void moveDown();
	void activeLeft();
	void activeRight();
	void rotateLeft();
	void rotateRight();
	void setSliceRGBA(float r = 0.1, float g = 0.1, float b = 0.1, float a = 0.1); //don't call show() again after calling this

	void drawCuboid();//first call to start drawing box in 2D and second call to draw cuboid in another dimension
	void profileH(); //first call to enable horizontal line sampling, and second call to disable
	void profileV(); //first call to enable vertical line sampling, and second call to disable
	void profileT(); //first call to enable arbitrary line sampling, and second call to disable

	bool getLineData(ArrayMgr<double>& xin, ArrayMgr<double>& yin);
	SFloat getPointValue(SFloat x, SFloat y, SFloat z);
	void bindGLContext(){ SetCurrent(*m_context); } // there are multiple openGL windows, so we need to switch context frequently

	/*data*/
protected:
	int NTET; // the number of tetrahedral
	int NX, NY, NZ; // must be specified from external matrix
	SFloat DX, DY, DZ; // default 0.1 cm
	SFloat xmin, ymin, zmin, xmax, ymax, zmax;
	SFloat xcenter, ycenter, zcenter;

	// mat and bgImage must have one presented, and bgPainter is optional
	ArrayMgr<SFloat>* mat; // a 3D matrix that will be shown
	ArrayMgr<SFloat>* bgImage; // a mono-color background. Usually from CT or MRI
	ArrayMgr<RGBPixel>* bgPainter; // a colorful background that may override the color generated from bgImag
	ArrayMgr<SFloat>* node; // the coordinates is relative to iso center (nNode*3 matrix)
	ArrayMgr<int>* elem; // (NTET * 5 matrix, the last colum is the material id)


	ArrayMgr<SFloat> tetLim; // NTET * 6 matrix, xmin,xmax, ymin,ymax, zmin,zmax, indicating the range of each tetrehedral
	ArrayMgr<SFloat> bgColor; // the grey level of the phantom (NTET*3 matrix)
	ArrayMgr<Vector> fNorm; // the normal vector of each face (NTET*4 matrix)
	ArrayMgr<int> f2t;

	map<string, MarkPoint> markPoints;

	// related to the matrix
	bool b_showMat;   // default true
	bool b_globalISOLines; // default true, show the isolines in different views with the same level standard
	int  NISOLevel; // default 9, how many level should we divide the matrix range
	int NISOLine; // the acutal iso lines to be rendered
	double stdColor; //default = 0.9, portion of the color gradient that defines the standard color
	double stdValue; // default = max of matrix.This isoline will be rendered with standard color
	float lineWidth; // default 2, line width of the isolines, range 1-32 in openGL
	double matAlpha; // default 0.2, transparent level of the matrix wash, range 0.0-1.0
	int colorScheme; // default 0, choose one color scheme to render the matrix
	float sliceAlpha; // default 0.1
	float sliceColor[3]; // default (0.1,0.1,0.1)
	bool b_showTriangleEdges; //default true

	// related to the background image
	bool b_showBGImg; // default true
	bool b_globalContrast; // default true. All sub-views have the same brightness standard
	double brightness; // default 1.0, brightness level of the background image, range 0.0-1.0
	double contrast_L; // default 0.0, low cut-off portion when showing background image, range 0.0-1.0
	double contrast_R; // default 1.0, high cut-off portion when showing background image, range 0.0-1.0

	bool b_showBGPainter; // default true
	bool b_showBox;  // default true, whether to show the matrix boundary
	bool b_showCrossLine; // default true, whether to show crossing line. It should be false when maximizing one of the sub-windows
	bool b_showMarkPoint; // default false
	bool b_maxGLWindow; //default false, if the view is maximized
	bool b_windowCreated; // default false. Used to deal with size event correctly. set it to true when calling initShow()
	//data related to the sub-windows
	int activeID; //current active window's id, default 0 (the first sub window), range 0-3
	int cross[3]; //crossing point index, default: the index of max value of the matrix
	struct GLWindow // data that each sub-windows will have 
	{
		FloatRect glRect;   // drawing region relative to wxGLCanvas' origin, used to dispatch mouse message and translate openGL drawing
		FloatRect picRect; // The region to show the picture, relative to the sub-window's upper-left corner, not to wxGLCanvas' origin
		//used for position picking
		GLint    viewport[4];
		GLdouble modelview[16];
		GLdouble projection[16];
	};
	GLWindow subw[4];
	struct PlanView
	{
		PlanView()
		{
			gltexture = 0;
		}
		GLuint gltexture;
		// the view angle for openGL
		int eyev[3];
		int upv[3];
		int rotateCounter; // record how the sub-window is rotated
		void initViewAngle(int wid) // can also be called to reset PlanView
		{
			if (wid == ZAXIS)
			{
				eyev[0] = 0; eyev[1] = 0, eyev[2] = -1;
				upv[0] = 0; upv[1] = -1; upv[2] = 0;
			}
			else if (wid == YAXIS)
			{
				eyev[0] = 0; eyev[1] = -1, eyev[2] = 0;
				upv[0] = 0; upv[1] = 0; upv[2] = 1;
			}
			else if (wid == XAXIS)
			{
				eyev[0] = 1; eyev[1] = 0, eyev[2] = 0;
				upv[0] = 0; upv[1] = -1; upv[2] = 0;
			}
			rotateCounter = 0;
		}
		bool transposed()
		{
			if (rotateCounter % 2 == 0) return false;
			else return true;
		}
	};
	PlanView pv[3];
	float theta, phi; // default theta = 60, phi = 50; used in the 3D rendering

	void calcSubWinRect(wxSize wxsz);
	int getSubCoordinates(wxPoint& p);
	FloatRect getBestFitRect(int wid, FloatRect& rect);
	//mouse and keyboard events for each window
	void mouseMoved(int wid, wxMouseEvent& event);
	void mouseDown(int wid, wxMouseEvent& event);
	void mouseLeftDouble(int wid, wxMouseEvent& event);
	void mouseReleased(int wid, wxMouseEvent& event);
	void keyPressed(int wid, wxKeyEvent& event);
	void rightClick(int wid, wxMouseEvent& event);
	void mouseLeftWindow(int wid, wxMouseEvent& event);
	void mouseEnterWindow(int wid, wxMouseEvent& event);
	void mouseWheelMoved(int wid, wxMouseEvent& event);
	void resized(int wid, FloatRect newGLRect);

	void prepare(int wid);

	int getSeq(int istart, int mod);
	int getNSlice(int wid);
	void rotate(int dir);
	float getScale(int wid);
	void getDWDH(int wid, double& DW, double& DH);
	virtual void render(int wid);
	void draw3D();
	void InitGL();

	// data related to drawing something. Note at a time, there's only one drawing process
	int drawingType; // default 0 (DRAW_NOTHING = 0), drawing nothing
	int drawingStatus; // default 0, before drawing anything. Different drawing types will have different meanings.
	wxPoint pMouse; // mouse position when clicking the left mouse button
	GLfloat pWinZ; // the projected depth in screen corresponding to pMouse
	GLdouble pObjX, pObjY, pObjZ; // the object coordinates when clicking the left mouse button. Note they're always >= 0
	GLdouble cObjX, cObjY, cObjZ; // the current object coordinates
	bool isPictureClicked(int wid)
	{
		FloatRect absRC = subw[wid].picRect;
		absRC.shift(subw[wid].glRect.x, subw[wid].glRect.y);
		if (absRC.x <= cMouse.x && cMouse.x <= absRC.Right() && absRC.y <= cMouse.y && cMouse.y <= absRC.Up()) return true;
		else return false;
	}
	inline bool isObjectClicked()
	{
		if (pObjX == -1 && pObjY == -1 && pObjZ == -1) return false; 
		else return true;
	}
	bool b_mouseLeftDown; // default false. whether the left mouse button is pressed down.
	bool b_showSampledLine; // default true. Draw the sampled line after mouse is released in corresponding mode
	int selectedSlice;// default VIEW3D, the slice index when clicking the left mouse button. No initialization needed.
	wxPoint cMouse; // current mouse position.
	FloatRect pRect; // the picRect when clicking the left button

	INT3D cPos; //current 3D position of the mouse ????????
	INT3D pPos; //previous 3D position of the mouse
	int pcross[3]; // previous crossing point index
	void getMatIndex(GLdouble& x, GLdouble& y, GLdouble& z, INT3D& Pos)
	{
		if (ZAXIS == activeID) z = (cross[ZAXIS] + 0.5)*DZ;
		else if (YAXIS == activeID) y = (cross[YAXIS] + 0.5)*DY;
		else if (XAXIS == activeID) x = (cross[XAXIS] + 0.5)*DX;
		Pos.ix = int(x / DX);
		Pos.iy = int(y / DY);
		Pos.iz = int(z / DZ);
	}
	BOX3D box; //it contains the 3D box vertex
	int boxDirection; //which direction is fully filled; it's auxilary varible

	GLText glText; // to print openGL text. call init() only once before printing openGL text
	GLUquadric* gluQuad; // openGL quadric object pointer

	wxGLContext* m_context; // default 0
	//>>
private: // private auxiliary data, very useful when rendering big 3D matrix
	SFloat g_BGMax, g_BGMin;
	bool g_BGInit; // default false. Reset it to false if new background image is loaded
	SFloat g_MatMax, g_MatMin;
	bool g_MATInit; // default false. Reset it to false if new matrix is loaded
	DECLARE_EVENT_TABLE()
};

#endif