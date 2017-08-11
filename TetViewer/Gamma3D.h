#ifndef _GAMMA3D_H_
#define _GAMMA3D_H_
#include "View3D.h"
#include "Graph.h"

class Gamma3D:public View3D
{
public:
	static Gamma3D& newGamma3D(wxWindow* parent, bool b_create_new)
	{
		FigureFrame* ff = FigureFrame::Frame(parent, b_create_new, FIGURE_3D);
		Gamma3D* gw = (Gamma3D*)(ff->getSubWindow());
		if (gw == NULL)
		{
			gw = new Gamma3D(ff);
			ff->setSubWindow(gw);
		}
		gw->bindGLContext(); // it may not create a new window, so we need to rebind the context for existed window
		ff->Show(true);//show the figure at this moment will reduce the blank window flashing
		return *gw;
	}
	Gamma3D(FigureFrame* parent);
	virtual ~Gamma3D();

	void setGamma(ArrayMgr<SFloat>* in_mat)//will copy the data into inner structure
	{
		inner_mat = *in_mat;
		src_mat = in_mat; //remember the source pointer
		setMatrix(&inner_mat); //we show showing the inner_mat, which could be modified
	}

	void setDoseRef(ArrayMgr<SFloat>* in_doseRef)
	{
		doseRef = in_doseRef; //remember the pointer of the reference dose
	}
	void setFilter(SFloat threshold, SFloat minDensity)
	{
		this->threshold = threshold;
		this->minDensity = minDensity;
	}
	void getPassingRate(int& tNum, int& tPass);

	void OnSave(wxCommandEvent& event);
	void OnExport(wxCommandEvent& event);
	void OnZoomIn(wxCommandEvent& event);
	void OnZoomOut(wxCommandEvent& event);
	void OnResetSize(wxCommandEvent& event);
	void OnMoveLeft(wxCommandEvent& event);
	void OnMoveRight(wxCommandEvent& event);
	void OnMoveUp(wxCommandEvent& event);
	void OnMoveDown(wxCommandEvent& event);
	void OnRotateLeft(wxCommandEvent& event);
	void OnRotateRight(wxCommandEvent& event);
	void OnProfileH(wxCommandEvent& event);
	void OnProfileV(wxCommandEvent& event);
	void OnProfileT(wxCommandEvent& event);

	void rightClick(wxMouseEvent& event);
	void OnPopupClick(wxCommandEvent& event);
	
protected:

	ArrayMgr<SFloat> inner_mat; // the gamma matrix may be modified, so we need an inner matrix
	ArrayMgr<SFloat>* src_mat; //source pointer of the gamma matrix, as backup for future restoration
	ArrayMgr<SFloat>* doseRef; // pointer of the external doseRef. We use it to calculate the passing rate
	ArrayMgr<RGBPixel> failVoxel;

	// parameters to calculate the passing rate
	SFloat threshold;//just a percentage
	SFloat minDensity;

	DECLARE_EVENT_TABLE()
};


#endif