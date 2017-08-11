#ifndef _PRECOMPILE_H_
#define _PRECOMPILE_H_

#ifdef _DEBUG   
//#pragma comment( linker, "/subsystem:console /entry:WinMainCRTStartup" )
#endif   

#define WXUSINGDLL

#include "wx/wx.h"
#include <wx/xrc/xmlres.h>
#include <wx/glcanvas.h>
#include <wx/sizer.h>
#include <wx/splitter.h>
#include <wx/bookctrl.h>
#include "wx/filedlg.h"
#include <wx/dir.h>
#include "wx/grid.h"
#include <wx/aboutdlg.h>
#include <wx/stdpaths.h>
#include <wx/busyinfo.h>
#include <GL/glu.h>


#include "../Tools/Tools.h"
#undef at

const wxString MODULENAME = "3DMatViewer";

class wxBitmapScale :public wxBitmap // this class is designed to show bitmaps in consistent physical size regardless of the DPI
{
public:
	wxBitmapScale(const wxString& name, int dw = 32, int dh = 32, double scale = 0)
	{
		if (scale == 0)
		{
			if (DPI_Scale == 1.25) scale = 1;
			else if (DPI_Scale == 2.25) scale = 2;
			else if (DPI_Scale == 2.5) scale = 2;
			else if (DPI_Scale == 2.75) scale = 3;
			else scale = DPI_Scale;
		}
		wxImage img(name);
		int w = int(dw * scale);
		int h = int(dh * scale);
		img.Rescale(w, h);
		*this = img;
	}
	wxBitmapScale(const wxImage& img) :wxBitmap(img){}
	static void setDPIScale(wxWindow* pw) { DPI_Scale = pw->GetContentScaleFactor(); }
	static double getDPIScale() { return DPI_Scale; }
private:
	static double DPI_Scale;
};

float pixelSize();

struct LineVector
{
	LineVector(double x1, double y1, double x2, double y2)
	{
		this->x1 = x1;
		this->y1 = y1;
		this->x2 = x2;
		this->y2 = y2;
	}
	double x1, y1;
	double x2, y2;
};

class Contour
{
public:
	Contour() { startLevel = 0; }
	void clear() //clear the contours manually
	{
		contourLine.resize(0);
	}
	void generate(ArrayMgr<SFloat> &d, int nc)
	{
		float maxv, minv;
		d.getMaxMin(maxv, minv);
		double dz = (maxv - minv) / (nc+1);
		z.resize(nc);
		for (int i = 0; i < nc; ++i)
		{
			z[i] = minv + (i + 1)*dz;
		}
		generate(d, z);
	}
	void generate(ArrayMgr<SFloat> &d, vector<double> &z)
	{
#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])

		//first check if z is in the right range

		int iub = d.getWidth(1) - 1;
		int jub = d.getWidth(2) - 1;
		int nc = z.size();
		contourLine.clear();
		contourLine.resize(nc);
		this->z =  z;//copy the Z array in case

		int m1, m2, m3, case_value;
		double dmin, dmax, x1 = 0, x2 = 0, y1 = 0, y2 = 0;
		int i, j, k, m;
		double h[5];
		int sh[5];
		double xh[5], yh[5];
		int im[4] = { 0, 1, 1, 0 }, jm[4] = { 0, 0, 1, 1 };
		int castab[3][3][3] = {
				{ { 0, 0, 8 }, { 0, 2, 5 }, { 7, 6, 9 } },
				{ { 0, 3, 4 }, { 1, 3, 1 }, { 4, 3, 0 } },
				{ { 9, 6, 7 }, { 5, 2, 0 }, { 8, 0, 0 } }
		};
		double temp1, temp2;

		for (j = (jub - 1); j >= 0; j--) {
			for (i = 0; i <= iub - 1; i++) {
				temp1 = min(d.a(i, j), d.a(i, j + 1));
				temp2 = min(d.a(i + 1, j), d.a(i + 1, j + 1));
				dmin = min(temp1, temp2);
				temp1 = max(d.a(i, j), d.a(i, j + 1));
				temp2 = max(d.a(i + 1, j), d.a(i + 1, j + 1));
				dmax = max(temp1, temp2);
				if (dmax < z[0] || dmin > z[nc - 1])
					continue;
				for (k = 0; k<nc; k++) {
					if (z[k] < dmin || z[k] > dmax)
						continue;
					for (m = 4; m >= 0; m--) {
						if (m > 0) {
							h[m] = d.a(i + im[m - 1], j + jm[m - 1]) - z[k];
							xh[m] = i + im[m - 1];
							yh[m] = j + jm[m - 1];
						}
						else {
							h[0] = 0.25 * (h[1] + h[2] + h[3] + h[4]);
							xh[0] = 0.50 * (i + i + 1);
							yh[0] = 0.50 * (j + j + 1);
						}
						if (h[m] > 0.0)
							sh[m] = 1;
						else if (h[m] < 0.0)
							sh[m] = -1;
						else
							sh[m] = 0;
					}

					/*
					Note: at this stage the relative heights of the corners and the
					centre are in the h array, and the corresponding coordinates are
					in the xh and yh arrays. The centre of the box is indexed by 0
					and the 4 corners by 1 to 4 as shown below.
					Each triangle is then indexed by the parameter m, and the 3
					vertices of each triangle are indexed by parameters m1,m2,and m3.
					It is assumed that the centre of the box is always vertex 2
					though this isimportant only when all 3 vertices lie exactly on
					the same contour level, in which case only the side of the box
					is drawn.
					vertex 4 +-------------------+ vertex 3
					| \               / |
					|   \    m-3    /   |
					|     \       /     |
					|       \   /       |
					|  m=2    X   m=2   |       the centre is vertex 0
					|       /   \       |
					|     /       \     |
					|   /    m=1    \   |
					| /               \ |
					vertex 1 +-------------------+ vertex 2
					*/
					/* Scan each triangle in the box */
					for (m = 1; m <= 4; m++) {
						m1 = m;
						m2 = 0;
						if (m != 4)
							m3 = m + 1;
						else
							m3 = 1;
						if ((case_value = castab[sh[m1] + 1][sh[m2] + 1][sh[m3] + 1]) == 0)
							continue;
						switch (case_value) {
						case 1: /* Line between vertices 1 and 2 */
							x1 = xh[m1];
							y1 = yh[m1];
							x2 = xh[m2];
							y2 = yh[m2];
							break;
						case 2: /* Line between vertices 2 and 3 */
							x1 = xh[m2];
							y1 = yh[m2];
							x2 = xh[m3];
							y2 = yh[m3];
							break;
						case 3: /* Line between vertices 3 and 1 */
							x1 = xh[m3];
							y1 = yh[m3];
							x2 = xh[m1];
							y2 = yh[m1];
							break;
						case 4: /* Line between vertex 1 and side 2-3 */
							x1 = xh[m1];
							y1 = yh[m1];
							x2 = xsect(m2, m3);
							y2 = ysect(m2, m3);
							break;
						case 5: /* Line between vertex 2 and side 3-1 */
							x1 = xh[m2];
							y1 = yh[m2];
							x2 = xsect(m3, m1);
							y2 = ysect(m3, m1);
							break;
						case 6: /* Line between vertex 3 and side 1-2 */
							x1 = xh[m3];
							y1 = yh[m3];
							x2 = xsect(m1, m2);
							y2 = ysect(m1, m2);
							break;
						case 7: /* Line between sides 1-2 and 2-3 */
							x1 = xsect(m1, m2);
							y1 = ysect(m1, m2);
							x2 = xsect(m2, m3);
							y2 = ysect(m2, m3);
							break;
						case 8: /* Line between sides 2-3 and 3-1 */
							x1 = xsect(m2, m3);
							y1 = ysect(m2, m3);
							x2 = xsect(m3, m1);
							y2 = ysect(m3, m1);
							break;
						case 9: /* Line between sides 3-1 and 1-2 */
							x1 = xsect(m3, m1);
							y1 = ysect(m3, m1);
							x2 = xsect(m1, m2);
							y2 = ysect(m1, m2);
							break;
						default:
							break;
						}

						/* Finally draw the line */
						//ConrecLine(x1, y1, x2, y2, z[k]);
						contourLine[k].push_back(LineVector(x1, y1, x2, y2));
						if (contourLine[k].size() > 100000) //get ride of uniform large matrix
						{
							contourLine.clear();
							return;
						}
					} /* m */
				} /* k - contour */
			} /* i */
		} /* j */
	}
	vector<double>& getZArray() { return z; }
	vector< vector<LineVector> > & getContourLine(){ return contourLine; }
	size_t size(){ return contourLine.size(); }
	size_t depthSize(){ return z.size(); }
private:

	vector< vector<LineVector> >  contourLine;
	vector<double> z;
	int startLevel;
};

struct FloatRect //aimed to resolve continuous resizing problem with wxRect
{
	float x, y, width, height;
	FloatRect() :x(0), y(0), width(0), height(0){}
	FloatRect(float inx, float iny, float inw, float inh) :x(inx), y(iny), width(inw), height(inh){}
	FloatRect& operator=(const wxRect& r)
	{ 
		x = r.x;
		y = r.y; 
		width = r.width;
		height=r.height; 
		return *this;
	}
	FloatRect& shift(float dx, float dy)
	{
		x += dx;
		y += dy;
		return *this;
	}
	FloatRect& transpose()
	{
		float cx = centerX();
		float cy = centerY();
		float temp = height;
		height = width;
		width = temp;
		x = cx - width / 2;
		y = cy - height / 2;
		return *this;
	}
	float centerX(){ return x + width / 2; }
	float centerY(){ return y + height / 2; }
	float Right(){ return x + width; }
	float Up(){ return y + height; } // openGL coordinates
	bool empty()
	{ 
		if (width == 0 && height == 0) return true;
		else return false;
	}
	void setEmpty()
	{
		x = y = 0;
		width = height = 0;
	}
	FloatRect& zoom(float scale)
	{
		if (scale <= 0) return *this;
		float cx = centerX();
		float cy = centerY();
		width *= scale;
		height *= scale;
		
		x = cx - width / 2;
		y = cy - height / 2;
		return *this;
	}
};

struct VERTEX3D
{
	VERTEX3D(double ex, double ey, double ez) { x = ex; y = ey; z = ez; }
	VERTEX3D() { x = y = z = 0; }
	double length() { return sqrt(x*x + y*y + z*z); }
	double x, y, z;
};

struct RGBPixel
{
	RGBPixel(uint8_t r, uint8_t g, uint8_t b)
	{
		this->r = r;
		this->g = g;
		this->b = b;
	}
	RGBPixel(){ r = 0; g = 0; b = 0; }
	uint8_t r, g, b;
	bool isBlack(){ if (r == 0 && g == 0 && b == 0) return true; else return false; }
	RGBPixel operator*(double mul)
	{
		RGBPixel p;
		p.r = r*mul;
		p.g = g*mul;
		p.b = b*mul;
		return p;
	}
};
inline bool operator <(const VERTEX3D vl, const VERTEX3D vr) { return vl.x < vr.x; } //this is only designed for std::set

struct INT3D
{
	INT3D() :ix(0), iy(0), iz(0){}
	int ix, iy, iz;
	INT3D& operator-=(INT3D& r)
	{
		ix -= r.ix;
		iy -= r.iy;
		iz -= r.iz;
		return *this;
	}
	INT3D operator - (INT3D& r)
	{
		INT3D temp = *this;
		temp.ix -= r.ix;
		temp.iy -= r.iy;
		temp.iz -= r.iz;
		return temp;
	}
	bool withinMatrix(int NX, int NY, int NZ)
	{
		if (ix < 0 || ix >= NX || iy < 0 || iy >= NY || iz < 0 || iz >= NZ) return false;
		else return true;
	}
	bool trim(int NX, int NY, int NZ)
	{
		if (ix < 0 || ix >= NX || iy < 0 || iy >= NY || iz < 0 || iz >= NZ)
		{
			if (ix < 0) ix = 0;
			else if (ix >= NX) ix = NX - 1;
			if (iy < 0) iy = 0;
			else if (iy >= NY) iy = NY - 1;
			if (iz < 0) iz = 0;
			else if (iz >= NZ) iz = NZ - 1;
			return true; // trimming has been applied
		}
		else return false; // no need trim
	}
	double distance()
	{
		return sqrt(ix*ix + iy*iy + iz*iz);
	}
};

struct BOX3D
{
	BOX3D() :x1(-1), x2(-1), y1(-1), y2(-1), z1(-1), z2(-1){}

	void sort()
	{
		if (x2 < x1) swapData<int>(x1, x2);
		if (y2 < y1) swapData<int>(y1, y2);
		if (z2 < z1) swapData<int>(z1, z2);
	}
	int x1, x2, y1, y2, z1, z2;
};

enum ViewDirection
{
	ZAXIS, YAXIS, XAXIS, VIEW3D //this must have the same order shown in the choice
};

const int wxGLCanvas_args[] = { WX_GL_RGBA, WX_GL_DOUBLEBUFFER, WX_GL_DEPTH_SIZE, 16, 0 };

struct ROICONTOUR
{
	wxString name;
	unsigned short r, g, b; //color
	bool show;
	vector<vector<VERTEX3D>> contour;
	struct SortAlgrithm
	{
		bool operator() (vector<VERTEX3D> left, vector<VERTEX3D> right) { return (left[0].z < right[0].z); }
	} sa;
	void sortbyz()
	{
		sort(contour.begin(), contour.end(), sa);
	}
};

class SliceInfo
{
public:
	double * density;
	short NRow, NCol;
	double thickness;
	double x, y, z; //coordinates of the upper left corner in patient coordinate system, unit mm;
	float pixSizeX, pixSizeY;
	bool operator<(const SliceInfo& r) const
	{
		return z < r.z;
	}
	void release(){ delete[] density; }
};

bool readFromDicom(VIEWRAY_FORMAT& red, const char* fname);

//pop up menu IDs
enum
{
	//the following is for class  GraphCurve
	ID_NOMALIZE_DOSE = 2001,
	ID_CHANGE_LEGEND,
	ID_ADJUST_DOSE,

	//the following is for class Gamma3D
	ID_TRIM_GAMMA_ABOVE_ONE,
	ID_TRIM_GAMMA_BELOW_ONE,
	ID_SHOW_GAMMA_BITMAP,
	ID_RESTORE_GAMMA,
	ID_SHOW_GAMMA_ISO_LINE,
	ID_SHOW_GAMMA_INFO
};

enum VIEW3D_EVENT //possible event View3D class may emitt.
{
	VIEW3D_CROSS_MOVED, //crossing point has moved by manual dragging
	VIEW3D_VOXEL_CLICKED, //one valid voxel was clicked
	VIEW3D_ACTIVE_VIEW_CHANGED,
	VIEW3D_GRAPH_DESTROYED,
	VIEW3D_LINE_SAMPLING_FINISHED, // the line sampling has been finished
	VIEW3D_CUBOID_SAMPLING_FINISHED // the cuboid samplign has been finished
};

#endif