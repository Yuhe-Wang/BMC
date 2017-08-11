#ifndef _GLVIEW_H_
#define _GLVIEW_H_

#include "View3D.h"

class GLView : public View3D
{
public:
	enum{PANEL_CONTROLS, PANEL_DOSE, PANEL_CONTOURS, PANEL_EXTERNAL};
	GLView(wxWindow* parent, wxFrame * frame);
	virtual ~GLView();

	void trimDose(); //called by background thread function
	void processCmdLine(int argc, wxChar ** argv);

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

	void OnNext(wxCommandEvent& event);
	void OnPrevious(wxCommandEvent& event);
	void OnShowBox(wxCommandEvent& event);
	void OnShowDose(wxCommandEvent& event);
	void OnShowPhantom(wxCommandEvent& event);
	void OnShowISOCenter(wxCommandEvent& event);
	void OnDoseLevel(wxCommandEvent& event);
	void OnViewDirection(wxCommandEvent& event);
	void OnJumpSlice(wxCommandEvent& event);
	void OnDropFiles(wxDropFilesEvent& event);
	void OnOpen(wxCommandEvent& event);
	void OnSaveSlice(wxCommandEvent& event);
	void OnAdjustDose(wxCommandEvent& event);
	void OnPDD(wxCommandEvent& event);
	void OnDoseProfile(wxCommandEvent& event);
	void OnZoomIn(wxCommandEvent& event);
	void OnZoomOut(wxCommandEvent& event);
	void OnResetSize(wxCommandEvent& event);
	void OnMoveLeft(wxCommandEvent& event);
	void OnActiveLeft(wxCommandEvent& event);
	void OnMoveRight(wxCommandEvent& event);
	void OnActiveRight(wxCommandEvent& event);
	void OnMoveUp(wxCommandEvent& event);
	void OnMoveDown(wxCommandEvent& event);
	void OnRotateLeft(wxCommandEvent& event);
	void OnRotateRight(wxCommandEvent& event);
	void OnSliderIndex(wxCommandEvent& event);
	void OnSliderLineWidth(wxCommandEvent& event);
	void OnSlidermatAlpha(wxCommandEvent& event);
	void OnSliderstdColor(wxCommandEvent& event);
	void OnSliderBrightness(wxCommandEvent& event);
	void OnSliderCenter(wxCommandEvent& event);
	void OnSliderWidth(wxCommandEvent& event);
	void OnSliderLeft(wxCommandEvent& event);
	void OnSliderRight(wxCommandEvent& event);
	void OnChangeC(wxCommandEvent& event);
	void OnChangeW(wxCommandEvent& event);
	void OnChangeL(wxCommandEvent& event);
	void OnChangeR(wxCommandEvent& event);
	void OnGlobalContrast(wxCommandEvent& event);
	void OnAutoContrast(wxCommandEvent& event);
	void OnShowContours(wxCommandEvent& event);
	void OnTrimDose(wxCommandEvent& event);
	void OnTrimLVOL(wxMouseEvent& event);
	void OnView3DNotify(wxCommandEvent& event);
	void OnViewDoseDifference(wxCommandEvent& event);
	void OnSwapDose(wxCommandEvent& event);
	void OnAverageDose(wxCommandEvent& event);
	void OnAveragePhantom(wxCommandEvent& event);
	void OnDrawBox(wxCommandEvent& event);
	void OnProfileH(wxCommandEvent& event);
	void OnProfileV(wxCommandEvent& event);
	void OnProfileT(wxCommandEvent& event);

	void OnSetPrescriptionDose(wxCommandEvent& event);
	void OnSamplePointDose(wxCommandEvent& event);
	void OnGetExpDTA(wxCommandEvent& event);

	void OnTrimDoseByDensity(wxCommandEvent& event);

	void OnSampleLineDose(wxCommandEvent& event);
	void OnSamplePlanDose(wxCommandEvent& event);
	void OnXPlanDistanceChanged(wxCommandEvent& event);
	void OnYPlanDistanceChanged(wxCommandEvent& event);
	void OnZPlanDistanceChanged(wxCommandEvent& event);

	void OnFlipDose(wxCommandEvent& event);
	void OnColorScheme(wxCommandEvent& event);
	void OnExportDiff(wxCommandEvent& event);

	void OnShowRelativeISOLines(wxCommandEvent& event);

	DECLARE_EVENT_TABLE()
	
/*data*/
private: 
	wxFrame * frame;

	//phantom info
	int NX, NY, NZ;
	double DX, DY, DZ;
	int uniform;
	double Hist;
	double COPX, COPY, COPZ; //position of phantom cuboid origin in patient coordinates
	//note that NX, DX, COPX, etc has identical variables in View3D.h, but here they're independent
	double xo, yo, zo; // position of the cuboid origin in ISO coordinates

	double Bs; //magnetic field strength
	double Bx, By, Bz; //magnetic field direction

	//3D matrix storage
	ArrayMgr<SFloat> phantom, dose, uncertainty;
	ArrayMgr<SFloat> doseRef, phantomRef; //may be used to show external dose
	//because the matrix of gamma may be a sub matrix of dose, so we need to store the following separately
	ArrayMgr<SFloat> gamma;
	ArrayMgr<SFloat> gammaBGImage;
	ArrayMgr<SFloat> gammaDoseRef;
	ArrayMgr<SFloat> gammaUncertainty;
	vector<SFloat> gammaDoseDiff; //used to export
	VIEWRAY_FORMAT RED; //relative electron density from ViewRay

	//contour set info
	
	//controls in the panel
	//wxStaticText *imcInfo;

	//ViewDirection viewAxis;
	//int imc; //current shown slice index, starting from 0
	wxString currentFile; //record current file name

//<< private Method
	void copyData(vector<SliceInfo> &slices3D); //called by loadCT
	bool loadDose(const wxString& name);

	void extraSettings();
	bool PNPoly(const vector<VERTEX3D>& vert, int ix, int iy);
	SFloat samplePointDose(double x, double y, double z);
	VERTEX3D searchExpDose(double x, double y, double z, double expDose, int NScan, int NAround);

	int getWidth() { return GetSize().x; }
	int getHeight(){ return GetSize().y; }
//>>
};

#endif 