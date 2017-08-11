#ifndef _GLVIEW_H_
#define _GLVIEW_H_

#include "View3D.h"
#include "ViewRayBeams.h"

class GLView : public View3D
{
public:
	wxFrame * frame;

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
	void OnShowCT(wxCommandEvent& event);
	void OnShowBox(wxCommandEvent& event);
	void OnShowDose(wxCommandEvent& event);
	void OnShowPhantom(wxCommandEvent& event);
	void OnShowAbsEnergy(wxCommandEvent& event);
	void OnShowISOCenter(wxCommandEvent& event);
	void OnDoseLevel(wxCommandEvent& event);
	void OnViewDirection(wxCommandEvent& event);
	void OnJumpSlice(wxCommandEvent& event);
	void OnDropFiles(wxDropFilesEvent& event);
	void OnOpen(wxCommandEvent& event);
	void OnOpenBeams(wxCommandEvent& event);
	void OnSave(wxCommandEvent& event);
	void OnSaveAs(wxCommandEvent& event);
	void OnSaveSlice(wxCommandEvent& event);
	void OnAdjustDose(wxCommandEvent& event);
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
	void OnEditData(wxCommandEvent& event);
	void OnSliderIndex(wxCommandEvent& event);
	void OnSliderLineWidth(wxCommandEvent& event);
	void OnSlidermatAlpha(wxCommandEvent& event);
	void OnSliderstdColor(wxCommandEvent& event);
	void OnSliderBrightness(wxCommandEvent& event);
	void OnSliderSliceAlpha(wxCommandEvent& event);
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
	void OnView3DNotify(wxCommandEvent& event);
	void OnShowExtDose(wxCommandEvent& event);
	void OnDelExtDose(wxCommandEvent& event);
	void OnViewDoseDifference(wxCommandEvent& event);
	void OnSwapDose(wxCommandEvent& event);
	void OnGammaAnalysis(wxCommandEvent& event);
	void OnShowPreviousGamma(wxCommandEvent& event);
	void OnAveragePhantom(wxCommandEvent& event);
	void OnDrawBox(wxCommandEvent& event);
	void OnProfileH(wxCommandEvent& event);
	void OnProfileV(wxCommandEvent& event);
	void OnProfileT(wxCommandEvent& event);

	void OnSetPrescriptionDose(wxCommandEvent& event);
	void OnSamplePointDose(wxCommandEvent& event);
	void OnGetExpDTA(wxCommandEvent& event);

	void OnTrimDoseByDensity(wxCommandEvent& event);
	void OnShowRelativeISOLines(wxCommandEvent& event);
	void OnShowTriEdges(wxCommandEvent& event);
	void OnShowBeams(wxCommandEvent& event);

	void OnSampleLineDose(wxCommandEvent& event);
	void OnSamplePlanDose(wxCommandEvent& event);
	void OnXPlanDistanceChanged(wxCommandEvent& event);
	void OnYPlanDistanceChanged(wxCommandEvent& event);
	void OnZPlanDistanceChanged(wxCommandEvent& event);

	void OnColorScheme(wxCommandEvent& event);
	void OnExportDiff(wxCommandEvent& event);
	void OnUpdateSelected(wxCommandEvent& event);
	void OnCalcUncertainty(wxCommandEvent& event);
	void OnSetSliceColor(wxCommandEvent& event);
	
	void showBF(BinaryFile& BF);
	void showCustomizedPhantom(const wxString& path);
/*data*/
private: 
	

	//phantom info
	int NX, NY, NZ;
	double DX, DY, DZ;
	double CT_kVp;
	int uniform;
	double Hist;
	double COPX, COPY, COPZ; //position of phantom cuboid origin in patient coordinates
	//note that NX, DX, COPX, etc has identical variables in View3D.h, but here they're independent
	double xo, yo, zo; // position of the cuboid origin in ISO coordinates (-xo, -yo, -zo) is the isocenter position in cuboid coordinate system

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

	ArrayMgr<int> tet; //should be nTet*5 matrix
	ArrayMgr<SFloat> vert; // should be nNode*3 matrix
	ArrayMgr<int> f2t; //optional nFace*4 matrix

	ViewRayBeams VRBeams;

	//ViewDirection viewAxis;
	bool b_showAbsE; // default false. show absolute energy instead of dose

	bool b_showBeams; // default true. Show beams in XY sub-window if any beam is loaded
	bool b_showDoseRef; //default false. Show reference dose instead of the primary dose
	int treatmentFraction; // default 1
	wxString currentFile; //record current file name
	long id_gammaWindow; // default -1, meaningless window identifier

//<< private Method
	virtual void render(int wid);
	void showDoseRef(unsigned int b = 2)
	{
		if (b > 1) b_showDoseRef = !b_showDoseRef; 
		else b_showDoseRef = b;
		if (b_showDoseRef) setMatrix(&doseRef);
		else setMatrix(&dose);
		show();
	}
	void save(wxString path, int saveDose);
	bool loadDose(const wxString& name);
	bool loadExtDose(const wxString& path, int ftype = 0);
	bool loadViewRayBeams(const wxString& path);
	
	bool hasExtDose();
	void extraSettings();
	SFloat samplePointDose(double x, double y, double z);
	VERTEX3D searchExpDose(double x, double y, double z, double expDose, int NScan, int NAround);

	int getWidth() { return GetSize().x; }
	int getHeight(){ return GetSize().y; }
//>>
	DECLARE_EVENT_TABLE()
};

extern GLView* ggl; //global pointer of GLView
#endif 