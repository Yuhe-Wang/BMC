#include "PreCompile.h"
#include <atomic>
#include "GLView.h"
#include "Graph.h"
#include "Gamma3D.h"
#include "CTOptionDlg.h"
#include "DataEditDlg.h"
#include "../GPUGamma/GPUGamma.h"

LogProvider Log;
GLView* ggl = NULL; //global pointer of GLView

#define CHECKCODE 3.14 //used to check the file format

enum DICOM_TYPE
{
	DCM_CT,
	DCM_RTSTRUCT,
	DCM_RTDOSE,
	DCM_INVALID
};

DICOM_TYPE getDicomType(wxString& path)
{
	//first identify the content. It should only open contours or external dose dicom files
	DcmFileFormat fileformat;
	fileformat.loadFile(OFFilename(path.c_str()));
	DcmDataset * dcmData = fileformat.getDataset();
	OFString ofs;
	dcmData->findAndGetOFString(DCM_Modality, ofs);
	DICOM_TYPE type;
	if (ofs == "CT")  type = DCM_CT;
	else if (ofs == "RTSTRUCT") type = DCM_RTSTRUCT;
	else if (ofs == "RTDOSE") type = DCM_RTDOSE;
	else type = DCM_INVALID;
	fileformat.clear();
	return type;
}

bool readFromDicom(VIEWRAY_FORMAT& red, const char* fname)
{
	double ex, ey, ez;
	double sliceThickness;
	double pixelSpacingX, pixelSpacingY;
	DicomFile df;
	if (!df.open(fname))
	{
		wxMessageBox("Cannot open the dicom file in function readFromDicom()");
	}
	try
	{
		MatlabType info = df.getAllItems();
		ex = info["ImagePositionPatient"].a(0);
		ey = info["ImagePositionPatient"].a(1);
		ez = info["ImagePositionPatient"].a(2);
		sliceThickness = info["SliceThickness"].a(0);
		pixelSpacingY = info["PixelSpacing"].a(0);
		pixelSpacingX = info["PixelSpacing"].a(1);
		//change the unit from mm to cm
		ex *= 0.1, ey *= 0.1, ez *= 0.1;
		sliceThickness *= 0.1;
		pixelSpacingX *= 0.1;
		pixelSpacingY *= 0.1;

		red.nx = info["Columns"].scast<int>();
		red.ny = info["Rows"].scast<int>();
		red.nz = info["NumberOfFrames"].scast<int>();

		//get the dose
		int pixpresent = info["PixelRepresentation"].scast<int>();
		if (pixpresent != 0)
		{
			wxMessageBox("Unsupported dicom file pixel representation!");
			return false;
		}
		ArrayMgr<uint16_t> pixels(red.nx, red.ny, red.nz); //rely on the MATLAB matrix memory layout

		pixels.copy((uint16_t*)df.getPixelData());
		ArrayMgr<SFloat> farray;
		pixels.cast<SFloat>(farray);
		double scale = info["DoseGridScaling"].scast<double>();
		farray *= scale;

		///////////////////////////////////////////
		red.Minor_Header = 1;
		red.Minor_Header = 0;
		red.offset_x = ex;
		red.offset_y = ey;
		red.offset_z = ez;
		red.dx = pixelSpacingX;
		red.dy = pixelSpacingY;
		red.dz = sliceThickness;
		//get the main data
		red.m = farray;
		red.uniform = red.m.isUniform();
	}
	catch (std::runtime_error& e)
	{
		wxMessageBox(e.what());
		return false;
	}
	
	return true;
}
int searchViewRayCouch(ArrayMgr<SFloat>& den) //this may change with ViewRay's definition
{
	int nx = den.getWidth(1);
	int ny = den.getWidth(2);
	int nz = den.getWidth(3);
	int ty = ny - 1;//scan for the position of ViewRa's nominal couch
	for (; ty >= 0; --ty) //the couch's height takes about 48 voxels
	{
		if (fabs(den.a(nx / 2, ty, nz / 2) / den.a(nx / 2, ty + 1, nz / 2) - 23.38) < 1e-3)
		{
			//do further test
			bool pass = true;
			for (int ix = 0; ix < nx / 2; ++ix)
				for (int iz = 0; iz < nz / 2; ++iz)
				{
					if (fabs(den.a(nx / 4 + ix, ty, nz / 4 + iz) - den.a(nx / 2, ty, nz / 2)) > 1e-4)
					{
						pass = false;
						break;
					}
				}
			if (pass) break;
		}
	}
	return ty; //ty==-1 means cannot find the ViewRay's couch
}

//C style function>>

wxDEFINE_EVENT(MY_NEW_TYPE, wxCommandEvent);

BEGIN_EVENT_TABLE(GLView, View3D)
//command event from frame
EVT_BUTTON(XRCID("m_button_previous"), GLView::OnPrevious)
EVT_BUTTON(XRCID("m_button_next"), GLView::OnNext)
EVT_BUTTON(XRCID("m_button_jump_slice"), GLView::OnJumpSlice)
EVT_BUTTON(XRCID("m_button_C"), GLView::OnChangeC)
EVT_BUTTON(XRCID("m_button_W"), GLView::OnChangeW)
EVT_BUTTON(XRCID("m_button_L"), GLView::OnChangeL)
EVT_BUTTON(XRCID("m_button_R"), GLView::OnChangeR)
EVT_BUTTON(XRCID("m_button_showExtDose"), GLView::OnShowExtDose)
EVT_BUTTON(XRCID("m_button_delExtDose"), GLView::OnDelExtDose)
EVT_BUTTON(XRCID("m_button_viewDoseDifference"), GLView::OnViewDoseDifference)
EVT_BUTTON(XRCID("m_button_swapDose"), GLView::OnSwapDose)
EVT_BUTTON(XRCID("m_button_gammaAnalysis"), GLView::OnGammaAnalysis)
EVT_BUTTON(XRCID("m_button_showPreviousGamma"), GLView::OnShowPreviousGamma)
EVT_BUTTON(XRCID("m_button_setPrescriptionDose"), GLView::OnSetPrescriptionDose)
EVT_BUTTON(XRCID("m_button_samplePointDose"), GLView::OnSamplePointDose)
EVT_BUTTON(XRCID("m_button_getExpDTA"), GLView::OnGetExpDTA)
EVT_BUTTON(XRCID("m_button_sampleLineDose"), GLView::OnSampleLineDose)
EVT_BUTTON(XRCID("m_button_samplePlanDose"), GLView::OnSamplePlanDose)
EVT_BUTTON(XRCID("m_button_trimDoseByDensity"), GLView::OnTrimDoseByDensity)
EVT_BUTTON(XRCID("m_button_exportDiff"), GLView::OnExportDiff)
EVT_BUTTON(XRCID("m_button_updateSelected"), GLView::OnUpdateSelected)
EVT_BUTTON(XRCID("m_button_calc_uncertainty"), GLView::OnCalcUncertainty)
EVT_BUTTON(XRCID("m_button_sliceColor"), GLView::OnSetSliceColor)
EVT_MENU(XRCID("m_menuItem_open"), GLView::OnOpen)
EVT_MENU(XRCID("m_menuItem_save"), GLView::OnSave)
EVT_MENU(XRCID("m_menuItem_saveSlice"), GLView::OnSaveSlice)
EVT_MENU(XRCID("m_menuItem_zoom_in"), GLView::OnZoomIn)
EVT_MENU(XRCID("m_menuItem_zoom_out"), GLView::OnZoomOut)
EVT_MENU(XRCID("m_menuItem_rotate_left"), GLView::OnRotateLeft)
EVT_MENU(XRCID("m_menuItem_rotate_right"), GLView::OnRotateRight)
EVT_TOOL(XRCID("m_tool_open"), GLView::OnOpen)
EVT_TOOL(XRCID("m_tool_beams"), GLView::OnOpenBeams)
EVT_TOOL(XRCID("m_tool_save"), GLView::OnSaveAs)
EVT_TOOL_RCLICKED(XRCID("m_tool_save"), GLView::OnSave)
EVT_TOOL(XRCID("m_tool_saveSlice"), GLView::OnSaveSlice)
EVT_TOOL(XRCID("m_tool_adjustDose"), GLView::OnAdjustDose)
EVT_TOOL(XRCID("m_tool_zoom_in"), GLView::OnZoomIn)
EVT_TOOL(XRCID("m_tool_zoom_out"), GLView::OnZoomOut)
EVT_TOOL(XRCID("m_tool_resetSize"), GLView::OnResetSize)
EVT_TOOL(XRCID("m_tool_left"), GLView::OnMoveLeft)
EVT_TOOL_RCLICKED(XRCID("m_tool_left"), GLView::OnActiveLeft)
EVT_TOOL(XRCID("m_tool_right"), GLView::OnMoveRight)
EVT_TOOL_RCLICKED(XRCID("m_tool_right"), GLView::OnActiveRight)
EVT_TOOL(XRCID("m_tool_up"), GLView::OnMoveUp)
EVT_TOOL(XRCID("m_tool_down"), GLView::OnMoveDown)
EVT_TOOL(XRCID("m_tool_rotate_left"), GLView::OnRotateLeft)
EVT_TOOL(XRCID("m_tool_rotate_right"), GLView::OnRotateRight)
EVT_TOOL(XRCID("m_tool_edit_data"), GLView::OnEditData)
EVT_TOOL(XRCID("m_tool_drawBox"), GLView::OnDrawBox)
EVT_TOOL(XRCID("m_tool_profile_horizontal"), GLView::OnProfileH)
EVT_TOOL(XRCID("m_tool_profile_vertical"), GLView::OnProfileV)
EVT_TOOL(XRCID("m_tool_profile_tilted"), GLView::OnProfileT)

EVT_CHECKBOX(XRCID("m_checkBox_showDose"), GLView::OnShowDose)
EVT_CHECKBOX(XRCID("m_checkBox_showPhantom"), GLView::OnShowPhantom)
EVT_CHECKBOX(XRCID("m_checkBox_showAbsE"), GLView::OnShowAbsEnergy)
EVT_CHECKBOX(XRCID("m_checkBox_showBox"), GLView::OnShowBox)
EVT_CHECKBOX(XRCID("m_checkBox_showISO"), GLView::OnShowISOCenter)
EVT_CHECKBOX(XRCID("m_checkBox_showBeams"), GLView::OnShowBeams)
EVT_CHECKBOX(XRCID("m_checkBox_globalContrast"), GLView::OnGlobalContrast)
EVT_CHECKBOX(XRCID("m_checkBox_autoContrast"), GLView::OnAutoContrast)
EVT_CHECKBOX(XRCID("m_checkBox_relativeISOLines"), GLView::OnShowRelativeISOLines)
EVT_CHECKBOX(XRCID("m_checkBox_showTriEdges"), GLView::OnShowTriEdges)
EVT_CHOICE(XRCID("m_choice_view_direction"), GLView::OnViewDirection)
EVT_CHOICE(XRCID("m_choice_doseLevel"), GLView::OnDoseLevel)
EVT_CHOICE(XRCID("m_choice_colorScheme"), GLView::OnColorScheme)
EVT_SLIDER(XRCID("m_slider_slice_index"), GLView::OnSliderIndex)
EVT_SLIDER(XRCID("m_slider_lineWidth"), GLView::OnSliderLineWidth)
EVT_SLIDER(XRCID("m_slider_matAlpha"), GLView::OnSlidermatAlpha)
EVT_SLIDER(XRCID("m_slider_stdColor"), GLView::OnSliderstdColor)
EVT_SLIDER(XRCID("m_slider_brightness"), GLView::OnSliderBrightness)
EVT_SLIDER(XRCID("m_slider_sliceAlpha"), GLView::OnSliderSliceAlpha)
EVT_SLIDER(XRCID("m_slider_C"), GLView::OnSliderCenter)
EVT_SLIDER(XRCID("m_slider_W"), GLView::OnSliderWidth)
EVT_SLIDER(XRCID("m_slider_L"), GLView::OnSliderLeft)
EVT_SLIDER(XRCID("m_slider_R"), GLView::OnSliderRight)
EVT_TEXT(XRCID("m_textCtrl_pointX3"), GLView::OnXPlanDistanceChanged)
EVT_TEXT(XRCID("m_textCtrl_pointY3"), GLView::OnYPlanDistanceChanged)
EVT_TEXT(XRCID("m_textCtrl_pointZ3"), GLView::OnZPlanDistanceChanged)

EVT_COMMAND(wxID_ANY, VIEW3D_NOTIFY, GLView::OnView3DNotify)
END_EVENT_TABLE()

GLView::GLView(wxWindow* parent, wxFrame* pframe) :View3D(parent)
{
	frame = pframe;

	b_showAbsE = false;
	b_showBeams = true;
	b_showDoseRef = false;
	b_showSampledLine = false; // turn off the internal function to draw the sampling line
	CT_kVp = 120;
	Hist = 0;
	xo = yo = zo = 0;
	id_gammaWindow = -1;

	//enable dragging file
	DragAcceptFiles(true);
	Connect(wxEVT_DROP_FILES, wxDropFilesEventHandler(GLView::OnDropFiles), NULL, this);

	ggl = this;
}

GLView::~GLView()
{
	//closeMCR(); //should release here, but it affects the graphics so I let it done automatically by system
}

//<<--------------------------------Event functions---------------------------------->>//

//command events
void GLView::OnDropFiles(wxDropFilesEvent& event)
{
	if (event.GetNumberOfFiles() > 1)
	{
		wxMessageBox("I can only load one file each time!");
		return;
	}
	//((wxNotebook*)FindWindowById(XRCID("m_notebook_panels")))->SetSelection(PANEL_CONTROLS);
	bool success = false;
	wxString path = event.GetFiles()[0];

	wxFileName fname(path);
	if (fname.IsDir())
	{

	}
	else //it's a file
	{
		if (wxFileExists(path))
		{
			wxString ext = fname.GetExt();
			if (ext == "mdose") success = loadDose(path);
			else if (ext == "beams") success = loadViewRayBeams(path);
			else if (ext == "txt") success = loadViewRayBeams(path);
		}
	}
	
	if(!success) wxMessageBox("Cannot open this file due to unsupported format or other reasons.");
	else extraSettings();
}

void GLView::OnOpen(wxCommandEvent& event)
{
	enum{ DOSE,BEAMS_VIEWRAY, PLANOVERVIEW_VIEWRAY }; //must have same order as the wildcard
	wxString caption = wxT("Choose a file");
	wxString wildcard = wxT("Dose files(*.mdose)|*.mdose|Beams file(*.beams)|*.beams|Plan Overview text(.txt)|*.txt");
	wxString defaultDir = wxEmptyString;
	wxString defaultFilename = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_MULTIPLE);
	if (dialog.ShowModal() == wxID_OK)
	{
		int ftype = dialog.GetFilterIndex();
		
		if (DOSE == ftype)
		{
			if (!loadDose(dialog.GetPath()))
			{
				wxMessageBox("The dose file doesn't have expected format");
				return;
			}
		}
		else if (BEAMS_VIEWRAY == ftype)
		{
			if (!loadViewRayBeams(dialog.GetPath()))
			{
				wxMessageBox("The ViewRay beams file doesn't have expected format");
				return;
			}
		}
		else if (PLANOVERVIEW_VIEWRAY == ftype)
		{
			if (!loadViewRayBeams(dialog.GetPath()))
			{
				wxMessageBox("The ViewRay beams file doesn't have expected format");
				return;
			}
		}
		extraSettings();
	}
}

void GLView::OnOpenBeams(wxCommandEvent& event)
{
	wxString caption = wxT("Choose beams configuration file");
	wxString wildcard = wxT("ViewRay beams file(*.beams)|*.beams|PlanOverview text file(*.txt)|*.txt");
	wxString defaultDir = wxEmptyString;
	wxString defaultFilename = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_OPEN);
	if (dialog.ShowModal() == wxID_OK)
	{
		if (!loadViewRayBeams(dialog.GetPath()))
		{
			wxMessageBox("The ViewRay beams file doesn't have expected format");
			return;
		}
	}
}

void GLView::OnSave(wxCommandEvent& event) //try to write the content to the original file
{
	wxFileName fname(currentFile);
	if (fname.GetExt() != "phtm" && fname.GetExt() != "dose"&& fname.GetExt() != "red") //need to call save as
	{
		wxCommandEvent evt;
		OnSaveAs(evt);
		return;
	}

	int ctype = -1;
	if (fname.GetExt() == "phtm") ctype = 1;
	else if (fname.GetExt() == "red") ctype = 3;
	else if (fname.GetExt() == "dose")
	{
		//it may be in my format or viewray's format
		BinaryFile BF;
		bool isMyFormat = BF.read(currentFile, true);
		if (isMyFormat) ctype = 0;
		else ctype = 2;
	}
	int action = wxMessageBox("Going to overwrite the original file, are you sure?", "overwrite confirmation",
		wxYES_NO | wxCANCEL);
	if (wxYES == action) save(currentFile, ctype);
}

void GLView::OnSaveAs(wxCommandEvent& event)
{
	if (dose.empty() && phantom.empty())
	{
		wxMessageBox("Nothing to save!");
		return;
	}

	wxString caption = wxT("Input a file name");
	wxString defaultFilename = wxT("out.dose");
	wxString wildcard = wxT("Dose files(*.dose)|*.dose|Phantom files(*.phtm)|*.phtm");
	wildcard += wxT("|ViewRay dose files(*.dose)|*.dose|RED files(*.red)|*.red"); //support for ViewRay
	wxString defaultDir = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_SAVE);
	if (dialog.ShowModal() == wxID_OK)
	{
		int ctype = dialog.GetFilterIndex();
		save(dialog.GetPath(), ctype);
	}
}

void GLView::OnSaveSlice(wxCommandEvent& event)
{
	if (phantom.empty() && dose.empty()) return;
	//save current slice to file that can be read by matlab
	wxString caption = wxT("Input a file name");
	wxString defaultFilename = wxT("export for matlab");
	wxString wildcard = wxT("Phantom slice(*.txt)|*.txt|Dose slice(*.txt)|*.txt|Dicom Dose(*.dcm)|*.dcm");

	wxString defaultDir = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_SAVE);
	if (dialog.ShowModal() == wxID_OK)
	{
		
	}
}

void GLView::OnAdjustDose(wxCommandEvent& event) //only adjust the current shown one
{
	static double lastFactor = 1.0; //using static variable is not safe if reloading dose
	double adjustK = 0;
	wxString str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_manualAjustment")))->GetLineText(0);
	if (!str.ToDouble(&adjustK))
	{
		wxMessageBox("incorrect adjustment!");
		return;
	};
	adjustK = 1 + adjustK*0.01;

	if (b_showDoseRef) doseRef *= adjustK / lastFactor;
	else dose *= adjustK / lastFactor;
	lastFactor = adjustK;
	wxMessageBox("overall dose adjustment done!");
}

void GLView::OnZoomIn(wxCommandEvent& event)
{
	zoomIn();
}
void GLView::OnZoomOut(wxCommandEvent& event)
{
	zoomOut();
}
void GLView::OnResetSize(wxCommandEvent& event)
{
	resetView();
}
void GLView::OnMoveLeft(wxCommandEvent& event)
{
	moveLeft();
}
void GLView::OnActiveLeft(wxCommandEvent& event)
{
	activeLeft();
}
void GLView::OnMoveRight(wxCommandEvent& event)
{
	moveRight();
}
void GLView::OnActiveRight(wxCommandEvent& event)
{
	activeLeft();
}
void GLView::OnMoveUp(wxCommandEvent& event)
{
	moveUp();
}
void GLView::OnMoveDown(wxCommandEvent& event)
{
	moveDown();
}
void GLView::OnNext(wxCommandEvent& event)
{
	nextSlice();
}
void GLView::OnPrevious(wxCommandEvent& event)
{
	previousSlice();
}
void GLView::OnRotateLeft(wxCommandEvent& event)
{
	rotateLeft();
}
void GLView::OnRotateRight(wxCommandEvent& event)
{
	rotateRight();
}

void GLView::OnEditData(wxCommandEvent& event)
{
	DataEditDlg DE(currentFile);
	DE.ShowModal();
}

void GLView::OnShowBox(wxCommandEvent& event)
{
	showBox();
}

void GLView::OnShowISOCenter(wxCommandEvent& event)
{
	showMarkPoint();
}

void GLView::OnShowDose(wxCommandEvent& event)
{
	showMat();
}

void GLView::OnShowPhantom(wxCommandEvent& event)
{
	b_showBGImg= !b_showBGImg;
	show();
}

void GLView::OnShowAbsEnergy(wxCommandEvent& event)
{
	b_showAbsE = !b_showAbsE;
	int NDose = dose.getInnerLength();
	int NDoseRef = doseRef.getInnerLength();
	if (b_showAbsE) //change dose to absolute energy
	{
		for (int i = 0; i < NDose; ++i) dose.a(i) *= phantom.a(i);
		for (int i = 0; i < NDoseRef; ++i) doseRef.a(i) *= phantom.a(i);
	}
	else //change absolution energy to dose
	{
		for (int i = 0; i < NDose; ++i)
		{
			if (phantom.a(i)>0) dose.a(i) /= phantom.a(i);
			else dose.a(i) = 0;
		}
		for (int i = 0; i < NDoseRef; ++i)
		{
			if (phantom.a(i) > 0) doseRef.a(i) /= phantom.a(i);
			else doseRef.a(i) = 0;
		}
	}
	show();
}

void GLView::OnDoseLevel(wxCommandEvent& event)
{
	int ichoice= event.GetSelection();
	setISOLevel(ichoice + 5);
}

void GLView::OnColorScheme(wxCommandEvent& event)
{
	colorScheme = event.GetSelection();
	show();
}

void GLView::OnViewDirection(wxCommandEvent& event)
{
	int wid = event.GetSelection();
	changeView(wid);
}

void GLView::OnJumpSlice(wxCommandEvent& event)
{
	wxString is=((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_slice_index")))->GetLineText(0);
	long index;
	is.ToLong(&index);
	jump2Slice(index - 1);
}

void GLView::OnSliderIndex(wxCommandEvent& event)
{
	static wxSlider* p_sliceIndex = (wxSlider*)(FindWindowByName("m_slider_slice_index"));
	static wxStaticText* p_crossText = (wxStaticText*)(FindWindowByName("m_staticText_cross"));

	int val = p_sliceIndex->GetValue();
	if (activeID < VIEW3D)
	{
		int NSlice = 0;
		if (XAXIS == activeID) NSlice = NX;
		else if (YAXIS == activeID) NSlice = NY;
		else if (ZAXIS == activeID) NSlice = NZ;
		int si = (NSlice - 1)*val / 100;
		cross[activeID] = si;
		show();
		wxString crossInfo;
		crossInfo.Printf("cross position =(%3d,%3d,%3d)", cross[XAXIS], cross[YAXIS], cross[ZAXIS]);
		p_crossText->SetLabelText(crossInfo);
	}
}
void GLView::OnSliderLineWidth(wxCommandEvent& event)
{
	int i = ((wxSlider*)FindWindowById(XRCID("m_slider_lineWidth")))->GetValue();
	lineWidth = i * 0.1f;
	renders();
}
void GLView::OnSlidermatAlpha(wxCommandEvent& event)
{
	int i = ((wxSlider*)FindWindowById(XRCID("m_slider_matAlpha")))->GetValue();
	matAlpha = i * 0.01f;
	show();
}
void GLView::OnSliderstdColor(wxCommandEvent& event)
{
	int i = ((wxSlider*)FindWindowById(XRCID("m_slider_stdColor")))->GetValue();
	stdColor = i*0.01;
	renders();
}
void GLView::OnSliderBrightness(wxCommandEvent& event)
{
	int i = ((wxSlider*)FindWindowById(XRCID("m_slider_brightness")))->GetValue();
	brightness = i * 0.01;
	show();
}
void GLView::OnSliderSliceAlpha(wxCommandEvent& event)
{
	int sv = ((wxSlider*)FindWindowById(XRCID("m_slider_sliceAlpha")))->GetValue();
	wxColour c = ((wxButton*)FindWindowById(XRCID("m_button_sliceColor")))->GetBackgroundColour();
	setSliceRGBA(c.Red() / 255.0f, c.Green() / 255.0f, c.Blue() / 255.0f, sv / 100.0f);
}
void GLView::OnSliderCenter(wxCommandEvent& event)
{
	int i=((wxSlider*)FindWindowById(XRCID("m_slider_C")))->GetValue();
	//assume width doesn't change
	double width = contrast_R - contrast_L;
	double center = i*0.01;
	if (center - width / 2 < 0 || center + width / 2 > 1) wxMessageBox("bad center position, please adjust width first!");
	else
	{
		contrast_L = center - width / 2;
		contrast_R = center + width / 2;

		wxString is;
		is.Printf("%d", i);
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_C")))->SetLabelText(is);

		is.Printf("%d", int(contrast_L*100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_L")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_L")))->SetValue(int(contrast_L * 100));
		is.Printf("%d", int(contrast_R * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_R")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_R")))->SetValue(int(contrast_R * 100));
		show();
	}
}
void GLView::OnSliderWidth(wxCommandEvent& event)
{
	int i = ((wxSlider*)FindWindowById(XRCID("m_slider_W")))->GetValue();
	//assume center doesn't change
	
	double center = (contrast_R + contrast_L) / 2;
	double width = i*0.01;
	if (center - width / 2 < 0 || center + width / 2 > 1) wxMessageBox("bad width, please adjust center first!");
	else
	{
		contrast_L = center - width / 2;
		contrast_R = center + width / 2;

		wxString is;
		is.Printf("%d", i);
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_W")))->SetLabelText(is);

		is.Printf("%d", int(contrast_L * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_L")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_L")))->SetValue(int(contrast_L * 100));
		is.Printf("%d", int(contrast_R * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_R")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_R")))->SetValue(int(contrast_R * 100));
		show();
	}
}
void GLView::OnSliderLeft(wxCommandEvent& event)
{
	int i = ((wxSlider*)FindWindowById(XRCID("m_slider_L")))->GetValue();
	if (i>contrast_R * 100) wxMessageBox("bad L value!");
	else
	{
		//assume R doesn't change
		contrast_L = i*0.01;
		double width = contrast_R - contrast_L;
		double center = (contrast_L + contrast_R) / 2;

		wxString is;
		is.Printf("%d", i);
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_L")))->SetLabelText(is);

		is.Printf("%d", int(center * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_C")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_C")))->SetValue(int(center * 100));
		is.Printf("%d", int(width * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_W")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_W")))->SetValue(int(width * 100));
		show();
	}
}
void GLView::OnSliderRight(wxCommandEvent& event)
{
	int i = ((wxSlider*)FindWindowById(XRCID("m_slider_R")))->GetValue();
	if (i<contrast_L * 100) wxMessageBox("bad R value!");
	else
	{
		//assume L doesn't change
		contrast_R = i*0.01;
		double width = contrast_R - contrast_L;
		double center = (contrast_L + contrast_R) / 2;

		wxString is;
		is.Printf("%d", i);
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_R")))->SetLabelText(is);

		is.Printf("%d", int(center * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_C")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_C")))->SetValue(int(center * 100));
		is.Printf("%d", int(width * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_W")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_W")))->SetValue(int(width * 100));
		show();
	}
}

void GLView::OnChangeC(wxCommandEvent& event)
{
	wxString is = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_C")))->GetLineText(0);
	long i;
	is.ToLong(&i);
	//assume width doesn't change
	double width = contrast_R - contrast_L;
	double center = i*0.01;
	if (center - width / 2 < 0 || center + width / 2 > 1) wxMessageBox("bad center position, please adjust width first!");
	else
	{
		contrast_L = center - width / 2;
		contrast_R = center + width / 2;

		((wxSlider*)FindWindowById(XRCID("m_slider_C")))->SetValue(i);

		is.Printf("%d", int(contrast_L * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_L")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_L")))->SetValue(int(contrast_L * 100));
		is.Printf("%d", int(contrast_R * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_R")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_R")))->SetValue(int(contrast_R * 100));
		show();
	}
}
void GLView::OnChangeW(wxCommandEvent& event)
{
	wxString is = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_W")))->GetLineText(0);
	long i;
	is.ToLong(&i);
	//assume center doesn't change
	double center = (contrast_R + contrast_L) / 2;
	double width = i*0.01;
	if (center - width / 2 < 0 || center + width / 2 > 1) wxMessageBox("bad center position, please adjust width first!");
	else
	{
		contrast_L = center - width / 2;
		contrast_R = center + width / 2;

		((wxSlider*)FindWindowById(XRCID("m_slider_W")))->SetValue(i);

		is.Printf("%d", int(contrast_L * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_L")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_L")))->SetValue(int(contrast_L * 100));
		is.Printf("%d", int(contrast_R * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_R")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_R")))->SetValue(int(contrast_R * 100));
		show();
	}
}
void GLView::OnChangeL(wxCommandEvent& event)
{
	wxString is = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_L")))->GetLineText(0);
	long i;
	is.ToLong(&i);
	if (i < 0||i>contrast_R*100) wxMessageBox("bad L value!");
	else
	{
		//assume R doesn't change
		contrast_L = i*0.01;
		double width = contrast_R - contrast_L;
		double center = (contrast_L + contrast_R) / 2;

		((wxSlider*)FindWindowById(XRCID("m_slider_L")))->SetValue(i);

		is.Printf("%d", int(center * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_C")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_C")))->SetValue(int(center * 100));
		is.Printf("%d", int(width * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_W")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_W")))->SetValue(int(width * 100));
		show();
	}
	
}
void GLView::OnChangeR(wxCommandEvent& event)
{
	wxString is = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_R")))->GetLineText(0);
	long i;
	is.ToLong(&i);
	if (i >100 || i<contrast_L * 100) wxMessageBox("bad R value!");
	else
	{
		//assume L doesn't change
		contrast_R = i*0.01;
		double width = contrast_R - contrast_L;
		double center = (contrast_L + contrast_R) / 2;

		((wxSlider*)FindWindowById(XRCID("m_slider_R")))->SetValue(i);

		is.Printf("%d", int(center * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_C")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_C")))->SetValue(int(center * 100));
		is.Printf("%d", int(width * 100));
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_W")))->SetLabelText(is);
		((wxSlider*)FindWindowById(XRCID("m_slider_W")))->SetValue(int(width * 100));
		show();
	}
}

void GLView::OnGlobalContrast(wxCommandEvent& event)
{
	showGlobalContrast();
	show();
}
void GLView::OnAutoContrast(wxCommandEvent& event)
{
	
}

void GLView::OnShowExtDose(wxCommandEvent& event)
{
	wxStaticText* edLoad = (wxStaticText *)FindWindowById(XRCID("m_staticText_extDoseLoaded"));
	wxStaticText* edOn = (wxStaticText *)FindWindowById(XRCID("m_staticText_extDoseOn"));
	wxButton* edShow = (wxButton *)FindWindowById(XRCID("m_button_showExtDose"));

	if (b_showDoseRef)
	{
		showDoseRef(false);
		edOn->SetBackgroundColour(*wxYELLOW);
		edOn->SetForegroundColour(*wxBLACK);
		edOn->SetLabelText("  NO  ");
		edShow->SetLabelText("Show External");
		show();
	}
	else //currently we're showing internal dose
	{
		if (edLoad->GetLabelText() == "  YES  ") // doseRef stores the external dose
		{
			showDoseRef(true);
			edOn->SetBackgroundColour(*wxGREEN);
			edOn->SetForegroundColour(*wxWHITE);
			edOn->SetLabelText("  Yes  ");
			edShow->SetLabelText("Show Internal");
			show();
		}
		else wxMessageBox("must load external dose through tool button first");	
	}
}
void GLView::OnDelExtDose(wxCommandEvent& event)
{
	static wxStaticText* edLoad = (wxStaticText *)FindWindowById(XRCID("m_staticText_extDoseLoaded"));
	static wxStaticText* edOn = (wxStaticText *)FindWindowById(XRCID("m_staticText_extDoseOn"));
	static wxButton* edShow = (wxButton *)FindWindowById(XRCID("m_button_showExtDose"));
	if (edLoad->GetLabelText() == "  YES  ") // we have external dose
	{
		if (b_showDoseRef) //currently show the external dose;
		{
			edOn->SetBackgroundColour(*wxYELLOW);
			edOn->SetForegroundColour(*wxBLACK);
			edOn->SetLabelText("  NO  ");
		}
		//need to reset the external loaded label
		edLoad->SetBackgroundColour(*wxYELLOW);
		edLoad->SetForegroundColour(*wxBLACK);
		edLoad->SetLabelText("  NO  ");
		showDoseRef(false);
		doseRef.release();
	}
	else wxMessageBox("Don't have external dose yet!");
}

void GLView::OnGammaAnalysis(wxCommandEvent& event)
{
	if (doseRef.empty()) return;

	wxStaticText* edLoad = (wxStaticText *)FindWindowById(XRCID("m_staticText_extDoseLoaded"));
	if (!(dose.getInnerLength() > 0 && edLoad->GetLabelText() == "  YES  "))
	{
		wxMessageBox("Must have internal and external dose loaded first!");
		return;
	}

	//some basic check
	RunTimeCounter tm;
	int nPass = 0;
	double percent = 0.03, DTA = 0.9;
	double threshold = 0;
	long NInterp = 1;
	wxString str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_deltaDose")))->GetLineText(0);
	if (!str.ToDouble(&percent) || percent < 0)
	{
		wxMessageBox("incorrect double format of the delta dose percentage");
		return;
	};
	percent *= 0.01;
	str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_DTA")))->GetLineText(0);
	if (!str.ToDouble(&DTA) || DTA <= 0)
	{
		wxMessageBox("incorrect double format of the DTA");
		return;
	};
	if (((wxCheckBox*)FindWindowById(XRCID("m_checkBox_gammaThreshold")))->IsChecked())
	{
		str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_gammaThreshold")))->GetLineText(0);
		if (!str.ToDouble(&threshold) || threshold <= 0 || threshold >= 100)
		{
			wxMessageBox("incorrect gammaThreshold");
			return;
		}; 
		threshold *= 0.01;
	}
	if (((wxCheckBox*)FindWindowById(XRCID("m_checkBox_interpolateDose")))->IsChecked())
	{
		str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_interpN")))->GetLineText(0);
		if (!str.ToLong(&NInterp) || NInterp < 1 || NInterp > 10)
		{
			wxMessageBox("incorrect interpolation parameter");
			return;
		};
	}
	
	bool b_onlyBox = ((wxCheckBox*)FindWindowById(XRCID("m_checkBox_selectedRegion")))->IsChecked();
	if (b_onlyBox && drawingStatus != 7)
	{
		wxMessageBox("Must have selected a small region first!");
		return;
	}
	bool b_local = ((wxCheckBox*)FindWindowById(XRCID("m_checkBox_localPercent")))->IsChecked();

	/*******************Start: get the two dose matrix to compare*******************/
	
	
	//do2 always points to the external dose
	ArrayMgr<SFloat>* d1 = NULL;
	ArrayMgr<SFloat>* d2 = NULL;
	//get gammaDose, gammaDoseRef, gammaBGImage
	ArrayMgr<SFloat> gammaDose;
	if (b_onlyBox)//need to resize new memory
	{
		//check the sampled cubic size
		BOX3D b = box;
		if (b.x2 < b.x1) swapData<int>(b.x1, b.x2);
		if (b.y2 < b.y1) swapData<int>(b.y1, b.y2);
		if (b.z2 < b.z1) swapData<int>(b.z1, b.z2);
		int mx = b.x2 - b.x1 + 1;
		int my = b.y2 - b.y1 + 1;
		int mz = b.z2 - b.z1 + 1;
		if (mx < 3 || mx < 3 || mz < 3)
		{
			wxMessageBox("warning: too few selected voxel!");
			return;
		}
		gammaDose.resize(mx, my, mz);
		gammaDoseRef.resize(mx, my, mz);
		gammaBGImage.resize(mx, my, mz);
		gammaUncertainty.resize(mx, my, mz);
		for (int iz = 0; iz < mz; ++iz)
			for (int iy = 0; iy < my; ++iy)
				for (int ix = 0; ix < mx; ++ix)
				{
					gammaDose.a(ix, iy, iz) = dose.a(ix + b.x1, iy + b.y1, iz + b.z1);
					gammaDoseRef.a(ix, iy, iz) = doseRef.a(ix + b.x1, iy + b.y1, iz + b.z1);
					gammaBGImage.a(ix, iy, iz) = phantom.a(ix + b.x1, iy + b.y1, iz + b.z1);
					gammaUncertainty.a(ix, iy, iz) = uncertainty.a(ix + b.x1, iy + b.y1, iz + b.z1);
				}
	}
	else
	{
		gammaDose = dose;
		gammaDoseRef = doseRef;
		gammaBGImage = phantom;
		gammaUncertainty = uncertainty;
	}

	//we may adjust the gammaDose manually
	if (((wxCheckBox*)FindWindowById(XRCID("m_checkBox_manualAdjust")))->IsChecked())
	{
		double adjustK = 0;
		str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_manualAjustment")))->GetLineText(0);
		if (!str.ToDouble(&adjustK))
		{
			wxMessageBox("incorrect adjustment!");
			return;
		};
		adjustK = 1 + adjustK*0.01;
		gammaDose *= adjustK;
	}

	//get the min density to consider
	double minDensity = 0.02;
	if (((wxCheckBox*)FindWindowById(XRCID("m_checkBox_gammaMinDensity")))->IsChecked())
	{
		wxString str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_gammaMinDensity")))->GetLineText(0);
		if (!str.ToDouble(&minDensity))
		{
			wxMessageBox("incorrect minimum density value!");
			return;
		}
	}

	//get the gamma dose difference for distribution analysis
	ArrayMgr<SFloat> doseDiff = gammaDose;
	doseDiff -= gammaDoseRef;
	//consider the low dose threshold, and the min phantom density
	int len = gammaDose.getInnerLength();
	SFloat maxv, minv;
	gammaDoseRef.getMaxMin(maxv, minv);
	SFloat minDose = maxv*threshold;
	gammaDoseDiff.resize(0);
	for (int i = 0; i < len; ++i)
	{
		if (gammaDose.a(i) > minDose && gammaBGImage.a(i) > minDensity && 0 <= gammaUncertainty.a(i) && gammaUncertainty.a(i)<0.01)
		{
			gammaDoseDiff.push_back(-doseDiff.a(i) / gammaDose.a(i));
		}
	}

	/*******************End: get the two dose matrix to compare*******************/
	gamma = gammaDose;//resize memory for the output matrix

	wxBusyCursor busyCursor;
	wxWindowDisabler disabler;
	//check whether to use GPU acceleration
	bool useGPU = true;
	if (((wxCheckBox*)FindWindowById(XRCID("m_checkBox_GPUGamma")))->IsChecked() && GPUAvailable())
	{
		if (GPUAvailable()) useGPU = true;
		else
		{
			wxMessageBox("Cannot find any compatible GPU to accelerate.\nWill use CPU to do calculation.");
			useGPU = false;
		}
	}
	else useGPU = false;

	if (useGPU)
	{
		wxBusyInfo busyInfo(wxT("doing gamma analysis with GPU, please wait..."));
		nPass = GPUGammaAnalysis(gammaDose, gammaDoseRef, gamma, DX, DY, DZ, percent, DTA, b_local, NInterp);
	}
	else
	{
		wxBusyInfo busyInfo(wxT("doing gamma analysis with CPU, please wait..."));
		nPass = gammaAnalysis(gammaDose, gammaDoseRef, gamma, DX, DY, DZ, percent, DTA, b_local, NInterp);
	}
	if (nPass < 0)
	{
		wxMessageBox("Failed to calculate gamma matrix. It may be caused by large memory allocation. Try to restart this program.");
		return;
	}

	//create the 3D gamma window
	Gamma3D& vg = Gamma3D::newGamma3D(this, false);
	vg.setGamma(&gamma);
	vg.setBGImage(&gammaBGImage);
	vg.setDoseRef(&gammaDoseRef);
	vg.setFilter(threshold, minDensity);
	vg.setVoxelSize(NX, NY, NZ, DX, DY, DZ);
	vg.initShow();
	id_gammaWindow = vg.GetParent()->GetId();
	if (!((wxCheckBox*)FindWindowById(XRCID("m_checkBox_show3DGamma")))->IsChecked()) vg.GetParent()->Hide();

	//count the passing rate, and print it out
	int tNum = 0; //the voxel number to be considered
	int tPass = 0; //the passed voxel number in considered voxels
	vg.getPassingRate(tNum, tPass);
	string strGPU = useGPU ? "yes" : "no";
	str.Printf("fitted voxel percent = %f%%, passed number = %d, passing rate = %f%%\n useGPU = %s, time cost = %.2f s", 100 * double(tNum) / gamma.getInnerLength(), tPass, 100 * double(tPass) / tNum, strGPU.c_str(), tm.stop());
	wxMessageBox(str);//print out the result
}
void GLView::OnExportDiff(wxCommandEvent& event)
{
	
	int len = gammaDoseDiff.size();
	if (len <= 0)
	{
		wxMessageBox("You should dose gamma analysis first");
		return;
	}
	wxString caption = wxT("Export dose difference");
	wxString wildcard = wxT("text file(*.txt)|*.txt");
	wxString defaultDir = wxEmptyString;
	wxString defaultFilename = wxT("Diff.txt");
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_SAVE);
	if (dialog.ShowModal() == wxID_OK)
	{
		wxString path = dialog.GetPath();
		FILE* fp = fopen(path.c_str(), "w");
		if (NULL == fp)
		{
			wxMessageBox("cannot create the file!");
			return;
		}
		for (int i = 0; i < len; ++i)
			fprintf(fp, "%f\n", gammaDoseDiff[i]);
		fclose(fp);
	}
}

void GLView::OnShowPreviousGamma(wxCommandEvent& event)
{
	if (!gamma.empty() && gamma.equalDims(gammaBGImage))
	{
		wxWindow* gw = FindWindowById(id_gammaWindow);
		if (gw != NULL)
		{
			gw->Show();
			return;
		}
		Gamma3D& vg = Gamma3D::newGamma3D(this, false);
		vg.setGamma(&gamma);
		vg.setBGImage(&gammaBGImage);
		vg.setDoseRef(&gammaDoseRef);
		vg.setVoxelSize(NX, NY, NZ, DX, DY, DZ);

		//get the filter
		wxString str;
		double threshold = 0;
		if (((wxCheckBox*)FindWindowById(XRCID("m_checkBox_gammaThreshold")))->IsChecked())
		{
			str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_gammaThreshold")))->GetLineText(0);
			if (!str.ToDouble(&threshold) || threshold <= 0 || threshold >= 100)
			{
				wxMessageBox("incorrect gammaThreshold");
				return;
			};
			threshold *= 0.01;
		}

		//get the min density to consider
		double minDensity = 0.04;
		if (((wxCheckBox*)FindWindowById(XRCID("m_checkBox_gammaMinDensity")))->IsChecked())
		{
			wxString str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_gammaMinDensity")))->GetLineText(0);
			if (!str.ToDouble(&minDensity))
			{
				wxMessageBox("incorrect minimum density value!");
				return;
			}
		}
		vg.setFilter(threshold, minDensity);
		vg.initShow();
	}
	else wxMessageBox("Haven't generate gamma matrix yet or it has been corrupted");
}

void GLView::OnDrawBox(wxCommandEvent& event)
{
	drawCuboid();
}
void GLView::OnProfileH(wxCommandEvent& event)
{
	profileH();
}
void GLView::OnProfileV(wxCommandEvent& event)
{
	profileV();
}
void GLView::OnProfileT(wxCommandEvent& event)
{
	profileT();
}

void GLView::OnView3DNotify(wxCommandEvent& event)
{
	static wxStaticText* p_crossText = (wxStaticText*)(FindWindowByName("m_staticText_cross"));
	static wxChoice* p_viewDir = (wxChoice*)(FindWindowByName("m_choice_view_direction"));
	static wxSlider* p_sliceIndex = (wxSlider*)(FindWindowByName("m_slider_slice_index"));

	int NSlice = 0;
	if (XAXIS == activeID) NSlice = NX;
	else if (YAXIS == activeID) NSlice = NY;
	else if (ZAXIS == activeID) NSlice = NZ;

	int eventType = event.GetInt();
	if (VIEW3D_CROSS_MOVED == eventType)
	{
		wxString crossInfo;
		crossInfo.Printf("cross position =(%3d,%3d,%3d)", cross[XAXIS], cross[YAXIS], cross[ZAXIS]);
		p_crossText->SetLabelText(crossInfo);
		if (activeID < VIEW3D) p_sliceIndex->SetValue(cross[activeID] * 100 / (NSlice - 1));
	}
	else if (VIEW3D_VOXEL_CLICKED == eventType)
	{
		//set the status bar text
		int kx = pPos.ix;
		int ky = pPos.iy;
		int kz = pPos.iz;
		wxString sbtext;
		sbtext.Printf("the clicked voxel index is (%3d, %3d, %3d)", kx, ky, kz);
// 		if (phantom.getInnerLength() > 0)
// 		{
// 			wxString den;
// 			den.Printf(", density = %4.2fg / cm ^ 3", phantom.a(kx, ky, kz));
// 			sbtext += den;
// 		}
// 		
// 		if (dose.getInnerLength() > 0)
// 		{
// 			wxString ds;
// 			if (b_showDoseRef) ds.Printf(", dose = %5.3f Gray", doseRef.a(kx, ky, kz));
// 			else ds.Printf(", dose = %5.3f Gray", dose.a(kx, ky, kz));
// 			sbtext += ds;
// 		}
// 		if (uncertainty.getInnerLength() > 0)
// 		{
// 			wxString err;
// 			err.Printf(", uncertainty = %4.2f%%", 100 * uncertainty.a(kx, ky, kz));
// 			sbtext += err;
// 		}
		wxStatusBar* sb = (wxStatusBar*)(FindWindowByName("m_statusBar"));
		sb->SetStatusText(sbtext);
	}
	else if (VIEW3D_ACTIVE_VIEW_CHANGED == eventType)
	{
		p_viewDir->SetSelection(activeID);

		if (activeID < VIEW3D) p_sliceIndex->SetValue(cross[activeID] * 100 / (NSlice - 1));
	}
	else if (VIEW3D_CUBOID_SAMPLING_CHANGED == eventType)
	{
		wxString sbox;
		wxTextCtrl* tc = (wxTextCtrl*)(FindWindowByName("m_textCtrl_selectedMinX"));
		sbox.Printf("%d", box.x1);
		tc->SetLabelText(sbox);	
		tc = (wxTextCtrl*)(FindWindowByName("m_textCtrl_selectedMinY"));
		sbox.Printf("%d", box.y1);
		tc->SetLabelText(sbox);
		tc = (wxTextCtrl*)(FindWindowByName("m_textCtrl_selectedMinZ"));
		sbox.Printf("%d", box.z1);
		tc->SetLabelText(sbox);
		tc = (wxTextCtrl*)(FindWindowByName("m_textCtrl_selectedMaxX"));
		sbox.Printf("%d", box.x2);
		tc->SetLabelText(sbox);
		tc = (wxTextCtrl*)(FindWindowByName("m_textCtrl_selectedMaxY"));
		sbox.Printf("%d", box.y2);
		tc->SetLabelText(sbox);
		tc = (wxTextCtrl*)(FindWindowByName("m_textCtrl_selectedMaxZ"));
		sbox.Printf("%d", box.z2);
		tc->SetLabelText(sbox);
	}
	else if (VIEW3D_LINE_SAMPLING_FINISHED == eventType)
	{
		ArrayMgr<double> xin, yin;
		FIGURE& fig = FIGURE::Figure(this, false);
		if (!doseRef.empty())
		{
			ArrayMgr<SFloat>* cmat = mat;
			setMatrix(&dose);
			getLineData(xin, yin);
			ArrayMgr<double> yin2;
			setMatrix(&doseRef);
			getLineData(xin, yin2);
			fig.Plot(MG(xin), MG(yin), MG(yin2), "b", "g");
			fig.g->AddLegend(L"PENELOPE", "b*");
			fig.g->AddLegend(L"KMC", "gd*");
			fig.g->Legend();
			setMatrix(cmat);
		}
		else //only one image
		{
			getLineData(xin, yin);
			fig.Plot(MG(xin), MG(yin), "b");
		}
		fig.g->Box();
		fig.g->Axis();
		fig.g->Grid("xy", "Y;");
		fig.g->Label('x', "distance(cm)", 0);
		if(!dose.empty()) fig.g->Label('y', "dose(Gy)", 0);
		else fig.g->Label('y', "density(g/cm^3)", 0);
		fig.g->Title("line sampling");
		if (fig.IsShown()) fig.renders();
	}
}

void GLView::OnSetPrescriptionDose(wxCommandEvent& event)
{
	wxString str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_prescriptionDose")))->GetLineText(0);
	double newPrescription = 0;
	if (!str.ToDouble(&newPrescription) || newPrescription <= 0)
	{
		wxMessageBox("Wrong new prescription dose!");
		return;
	}
	stdValue = newPrescription;
	b_globalISOLines = true;
	show();
}
void GLView::OnSamplePointDose(wxCommandEvent& event)
{
	if (dose.empty()) return;

	double DIC = 0.49, LIC = 0.64;
	double x = 0, y = 0, z = 0;
	int NScan = 100;
	wxString str;
	str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_IC_D")))->GetLineText(0);
	if (!str.ToDouble(&DIC))
	{
		wxMessageBox("Wrong diameter of ion chamber!");
		return;
	}
	str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_IC_L")))->GetLineText(0);
	if (!str.ToDouble(&LIC))
	{
		wxMessageBox("Wrong length of ion chamber!");
		return;
	}

	str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointX")))->GetLineText(0);
	if (!str.ToDouble(&x))
	{
		wxMessageBox("Wrong x coordinates!");
		return;
	}
	str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointY")))->GetLineText(0);
	if (!str.ToDouble(&y))
	{
		wxMessageBox("Wrong y coordinates!");
		return;
	}
	str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointZ")))->GetLineText(0);
	if (!str.ToDouble(&z))
	{
		wxMessageBox("Wrong z coordinates!");
		return;
	}
	long LongTemp;
	str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_NScan")))->GetLineText(0);
	if (!str.ToLong(&LongTemp) || LongTemp < 0)
	{
		wxMessageBox("Wrong NScan range!");
		return;
	}
	NScan = LongTemp;
	double dz = LIC / NScan;
	double sigma_weight = 0;
	double sigma_dose = 0;
	double max_dose = 0;
	double min_dose = 1e99;
	for (int iz = 0; iz <= NScan; ++iz)
		for (int ix = 0; ix*dz <= DIC; ++ix)
			for (int iy = 0; iy*dz <= DIC; ++iy)
			{
				//coordinates relative to the center of ion chamber
				double xc = dz*ix - DIC / 2;
				double yc = dz*iy - DIC / 2;
				double zc = dz*iz - LIC / 2;
				double ratio = (xc*xc + yc*yc) * 4 / (DIC*DIC);
				if (ratio<1) //inside the ion chamber
				{
					//coordinates relative to isocenter
					double d = samplePointDose(xc + x, yc + y, zc + z);
					double weight = 1; //exp(-ratio / 2);  //this should be a function of (xc,yc,zc)
					sigma_dose += d*weight;
					sigma_weight += weight;
					if (d > max_dose) max_dose = d;
					if (d < min_dose) min_dose = d;
				}
			}
	double average_dose = sigma_dose / sigma_weight;
	double center_dose = samplePointDose(x, y, z);
	wxString info;
	info.Printf("The center dose = %.3f Gy. In the ion chamber, max dose = %.3f Gy, min dose = %.3f Gy, average dose = %.3f Gy", center_dose, max_dose, min_dose, average_dose);
	wxMessageBox(info);
}
void GLView::OnGetExpDTA(wxCommandEvent& event)
{
	if (dose.empty()) return;
	double x = 0, y = 0, z = 0;
	double expDose = 0;
	int NScan = 100, NAround = 3;
	wxString str;
	str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointX")))->GetLineText(0);
	if (!str.ToDouble(&x))
	{
		wxMessageBox("Wrong coordinates!");
		return;
	}
	str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointY")))->GetLineText(0);
	if (!str.ToDouble(&y))
	{
		wxMessageBox("Wrong coordinates!");
		return;
	}
	str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointZ")))->GetLineText(0);
	if (!str.ToDouble(&z))
	{
		wxMessageBox("Wrong coordinates!");
		return;
	}
	str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_expPointDose")))->GetLineText(0);
	if (!str.ToDouble(&expDose) || expDose < 0)
	{
		wxMessageBox("Wrong experiment dose input!");
		return;
	}
	long LongTemp;
	str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_NScan")))->GetLineText(0);
	if (!str.ToLong(&LongTemp) || LongTemp < 0)
	{
		wxMessageBox("Wrong NScan range!");
		return;
	}
	NScan = LongTemp;
	str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_NAround")))->GetLineText(0);
	if (!str.ToLong(&LongTemp) || LongTemp < 0)
	{
		wxMessageBox("Wrong NAround range!");
		return;
	}
	NAround = LongTemp;

	VERTEX3D nearest = searchExpDose(x, y, z, expDose, NScan, NAround);
	double dist = nearest.length();
	if (dist < 1e10)
	{
		wxString info;
		info.Printf("Found the nearest point %.2f mm away. The deviation is (%.3f,%.3f,%.3f) mm", dist * 10, nearest.x * 10, nearest.y * 10, nearest.z * 10);
		wxMessageBox(info);
	}
	else wxMessageBox("Cannot find one point in nearby voxel");
}
void GLView::OnSampleLineDose(wxCommandEvent& event)
{
	if (dose.empty()) return;
	double x1 = 0, y1 = 0, z1 = 0;
	double x2 = 0, y2 = 0, z2 = 0;
	double ds;
	wxString strx1,strx2,stry1,stry2,strz1,strz2,strds;
	strx1 = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointX1")))->GetLineText(0);	
	stry1 = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointY1")))->GetLineText(0);
	strz1 = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointZ1")))->GetLineText(0);
	strx2 = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointX2")))->GetLineText(0);
	stry2 = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointY2")))->GetLineText(0);
	strz2 = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointZ2")))->GetLineText(0);
	strds = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointsSpacing")))->GetLineText(0);
	if (!strx1.ToDouble(&x1) || !stry1.ToDouble(&y1) || !strz1.ToDouble(&z1) || !strx2.ToDouble(&x2) || !stry2.ToDouble(&y2) || !strz2.ToDouble(&z2) || !strds.ToDouble(&ds) || ds <= 0)
	{
		wxMessageBox("Wrong parameters!");
		return;
	}
	//get the sampling points
	double xv = x2 - x1;
	double yv = y2 - y1;
	double zv = z2 - z1;
	double v = sqrt(xv*xv + yv*yv + zv*zv);
	if (v < ds)
	{
		wxMessageBox("Inappropriate line or points spacing! You may need to sample a longer line");
		return;
	}
	wxString caption = wxT("Input a file name");
	wxString defaultFilename = wxT("line dose.txt");
	wxString wildcard = wxT("text file(*.txt)|*.txt");

	wxString defaultDir = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_SAVE);
	if (dialog.ShowModal() == wxID_OK)
	{
		FILE *fp = fopen(dialog.GetPath().c_str(), "w");
		if (NULL == fp) wxMessageBox("Cannot create output file for line dose!");
		for (int i = 0; i*ds <= v; ++i)
		{
			SFloat pdose = samplePointDose(x1 + i*ds*xv / v, y1 + i*ds*yv / v, z1 + i*ds*zv / v);
			if (pdose < 0)
			{
				fclose(fp);
				wxMessageBox("sampling point reached out of the phantom!");
				return;
			}
			fprintf(fp, "%.2f\t%.3f\n", i*ds, pdose);
		}
		fclose(fp);
		strds.Printf("Line dose has been saved in file %s", dialog.GetPath().c_str());
		wxMessageBox(strds);
	}

}

void GLView::OnXPlanDistanceChanged(wxCommandEvent& event)
{
	wxString val = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointX3")))->GetLineText(0);
	if (((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointX3")))->GetLineText(0) != "") //clear the rest
	{
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointY3")))->SetLabelText("");
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointZ3")))->SetLabelText("");
	}
}
void GLView::OnYPlanDistanceChanged(wxCommandEvent& event)
{
	if (((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointY3")))->GetLineText(0) != "") //clear the rest
	{
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointX3")))->SetLabelText("");
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointZ3")))->SetLabelText("");
	}
}
void GLView::OnZPlanDistanceChanged(wxCommandEvent& event)
{
	if (((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointZ3")))->GetLineText(0) != "") //clear the rest
	{
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointX3")))->SetLabelText("");
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointY3")))->SetLabelText("");
	}
}

void GLView::OnSamplePlanDose(wxCommandEvent& event)
{
	if (dose.empty()) return;

	double x = 0, y = 0, z = 0;
	wxString strx, stry, strz;
	strx = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointX3")))->GetLineText(0);
	stry = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointY3")))->GetLineText(0);
	strz = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_pointZ3")))->GetLineText(0);
	int dimx, dimy;
	int direction = -1;
	float* mdata = NULL;
	double dw, dh;
	double dox = COPX + DX / 2, doy = COPY + DY / 2, doz = COPZ + DZ / 2;
	if (strx != "")
	{
		if (!strx.ToDouble(&x))
		{
			wxMessageBox("Wrong x distance!");
			return;
		}
		if (x-xo<0 || x-xo>NX*DX)
		{
			wxMessageBox("x outside the phantom!");
			return;
		}
		direction = XAXIS;
		dimx = NZ;
		dimy = NY;
		dw = DZ;
		dh = DY;
		mdata = new float[dimx*dimy];
		for (int ix = 0; ix < dimx; ++ix)
			for (int iy = 0; iy < dimy; ++iy)
				mdata[ix + dimx*iy] = samplePointDose(x, DY*(iy + 0.5) + yo, DZ*(ix + 0.5) + zo);

		dox += cross[XAXIS] * DX;
	}
	else if (stry != "")
	{
		if (!stry.ToDouble(&y))
		{
			wxMessageBox("Wrong y distance!");
			return;
		}
		if (y - yo<0 || y - yo>NY*DY)
		{
			wxMessageBox("y outside the phantom!");
			return;
		}
		direction = YAXIS;
		dimx = NX;
		dimy = NZ;
		dw = DX;
		dh = DZ;
		mdata = new float[dimx*dimy];
		for (int ix = 0; ix < dimx; ++ix)
			for (int iy = 0; iy < dimy; ++iy)
				mdata[ix + dimx*(dimy-1-iy)] = samplePointDose(DX*(ix + 0.5) + xo, y, DZ*(iy + 0.5) + zo);

		doy += cross[YAXIS] * DY;
		doz += (NZ - 1)*DZ; //there's only (NZ-1)*DZ shift
	}
	else if (strz != "")
	{
		if (!strz.ToDouble(&z))
		{
			wxMessageBox("Wrong z distance!");
			return;
		}
		if (z - zo<0 || x - zo>NZ*DZ)
		{
			wxMessageBox("z outside the phantom!");
			return;
		}
		direction = ZAXIS;
		dimx = NX;
		dimy = NY;
		dw = DX;
		dh = DY;
		mdata = new float[dimx*dimy];
		for (int ix = 0; ix < dimx; ++ix)
			for (int iy = 0; iy < dimy; ++iy)
				mdata[ix + dimx*iy] = samplePointDose(DX*(ix + 0.5) + xo, DY*(iy + 0.5) + yo, z);

		doz += cross[ZAXIS] * DZ;
	}

	if (((wxCheckBox*)FindWindowById(XRCID("m_checkBox_perFraction")))->IsChecked())
	{
		for (int i = 0; i < dimx*dimy; ++i) mdata[i] /= treatmentFraction;
	}

	//save current slice to file that can be read by matlab
	wxString caption = wxT("Input a file name");
	wxString defaultFilename = wxT("plan dose");
	wxString wildcard = wxT("Dicom dose(*.dcm)|*.dcm|Dose slice(*.txt)|*.txt");

	wxString defaultDir = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_SAVE);
	if (dialog.ShowModal() == wxID_OK)
	{
		int selection = dialog.GetFilterIndex();
		wxString path = dialog.GetPath();
		if (0 == selection) //call Matlab to write the dicom dose
		{
			//need cd to module directory
			wxString cwd = wxGetCwd();
			wxSetWorkingDirectory(wxPathOnly(wxStandardPaths::Get().GetExecutablePath()));

			DicomFile df;
			if (!df.open("2D_dose_sample.dcm"))
			{
				wxMessageBox("Cannot open the 2D_dose_sample.dcm!");
				return;
			}
			df.setItem("Width", dimx);
			df.setItem("Height", dimy);
			df.setItem("SliceThickness", 3);
			MatlabType mm;
			mm.resize(3);
			mm.a(0) = dox;
			mm.a(1) = doy;
			mm.a(2) = doz;
            MatlabType mul(10.0);
			mm *= mul; // convert cm to mm
			df.setItem("ImagePositionPatient", mm);
			df.setItem("Rows", dimy);
			df.setItem("Columns", dimx);
			mm.resize(2);
			mm.a(0) = dh;
			mm.a(1) = dw;
			df.setItem("PixelSpacing", mm);
			//find the max dose, and the scale factor
			int len = dimx*dimy;
			SFloat md = mdata[0];
			for (int i = 1; i < len; ++i)
			{
				if (mdata[i] > md) md = mdata[i];
			}
			SFloat scale = 61440 / md;
			df.setItem("DoseGridScaling", scale);
			//scale the dose
			uint16_t* pin = new uint16_t[len];
			for (int i = 1; i < len; ++i)
			{
				pin[i] = mdata[i] * scale;
			}
			delete[] mdata;
			
			// get the ImageOrientationPatient
			mm.resize(3);
			if (direction == ZAXIS)
			{
				const double iop[6] = { 1, 0, 0, 0, 1, 0 };
				df.setItem("ImageOrientationPatient", MatlabType(iop,6));
			}
			else if (direction == YAXIS)
			{
				const double iop[6] = { 1, 0, 0, 0, 0, -1 };
				df.setItem("ImageOrientationPatient", MatlabType(iop, 6));
			}
			else if (direction == XAXIS)
			{
				const double iop[6] = { 0, 0, 1, 0, 1, 0 };
				df.setItem("ImageOrientationPatient", MatlabType(iop, 6));
			}
			
			df.setPixelData(pin, len*sizeof(uint16_t));
			delete[] pin;
			df.save(path.ToStdString());
			wxSetWorkingDirectory(cwd); //recover current working directory
			wxMessageBox("DICOM file has been saved!");
		}
		else
		{
			FILE* fp = fopen(path.c_str(), "w");
			if (fp == NULL)
			{
				wxMessageBox("cannot create the output file!");
				return;
			}

			for (int iy = 0; iy < dimy; ++iy)
			{
				for (int ix = 0; ix < dimx; ++ix)
				{
					if (ix != dimx - 1) fprintf(fp, "%f,", mdata[ix + dimx*iy]);
					else fprintf(fp, "%f\n", mdata[ix + dimx*iy]);
				}
			}
			delete[] mdata;
			fclose(fp);
			wxMessageBox("Matlab compatible file has been saved!");
		}
	}
}

void GLView::OnTrimDoseByDensity(wxCommandEvent& event)
{
	wxString str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_densityThreshold")))->GetLineText(0);
	double denThresh = 0;
	if (!str.ToDouble(&denThresh) || denThresh <= 0)
	{
		wxMessageBox("Wrong density threshold!");
		return;
	}
	int len = dose.getInnerLength();
	for (int i = 0; i < len; ++i)
	{
		if (phantom.a(i) < denThresh) dose.a(i) = 0;
	}
	show();
}

void GLView::OnShowRelativeISOLines(wxCommandEvent& event)
{
	if (dose.empty()) return;
	b_globalISOLines = !b_globalISOLines;
	if (0 == stdValue&&b_globalISOLines)
	{
		//Let's assign a reasonable prescription dose
		SFloat maxv, minv;
		dose.getMaxMin(maxv, minv);
		b_globalISOLines = true;
		stdValue = maxv;

		if (0 == stdValue)
		{
			b_globalISOLines = false;
			wxMessageBox("You must provide a prescription dose before showing global ISOLines");
			((wxCheckBox*)FindWindowById(XRCID("m_checkBox_relativeISOLines")))->SetValue(true);
		}
	}
	show();
}

void GLView::OnShowTriEdges(wxCommandEvent& event)
{
	showTriangleEdges();
}
void GLView::OnShowBeams(wxCommandEvent& event)
{
	b_showBeams = !b_showBeams;
	renders();
}

void GLView::OnViewDoseDifference(wxCommandEvent& event)
{
	dose -= doseRef;
	int len = dose.getInnerLength();
	SFloat cutoff = 0.3f;
	for (int i = 0; i < len; ++i)
	{
		dose.a(i) = fabsf(dose.a(i));
		if (dose.a(i) > cutoff) dose.a(i) = cutoff;
	}
	show();
}
void GLView::OnSwapDose(wxCommandEvent& event)
{
	if (!dose.empty() && !doseRef.empty() &&dose.swap(doseRef))
	{
		wxButton* swapButton = (wxButton *)FindWindowById(XRCID("m_button_swapDose"));
		if (swapButton->GetBackgroundColour() == *wxYELLOW) swapButton->SetBackgroundColour(*wxWHITE);
		else swapButton->SetBackgroundColour(*wxYELLOW);
		show();
	}
}

void GLView::OnUpdateSelected(wxCommandEvent& event)
{
	wxTextCtrl* tc = (wxTextCtrl*)(FindWindowByName("m_textCtrl_selectedMinX"));
	wxString sbox = tc->GetLineText(0);
	long ind = 0;
	bool valid = true;
	if (sbox.ToLong(&ind) && 0 <= ind&&ind < NX) { box.x1 = ind; }
	else valid = false;

	tc = (wxTextCtrl*)(FindWindowByName("m_textCtrl_selectedMinY"));
	sbox = tc->GetLineText(0);
	if (sbox.ToLong(&ind) && 0 <= ind&&ind < NY) { box.y1 = ind; }
	else valid = false;

	tc = (wxTextCtrl*)(FindWindowByName("m_textCtrl_selectedMinZ"));
	sbox = tc->GetLineText(0);
	if (sbox.ToLong(&ind) && 0 <= ind&&ind < NZ) { box.z1 = ind; }
	else valid = false;

	tc = (wxTextCtrl*)(FindWindowByName("m_textCtrl_selectedMaxX"));
	sbox = tc->GetLineText(0);
	if (sbox.ToLong(&ind) && 0 <= ind&&ind < NX) { box.x2 = ind; }
	else valid = false;
	tc = (wxTextCtrl*)(FindWindowByName("m_textCtrl_selectedMaxY"));
	sbox = tc->GetLineText(0);
	if (sbox.ToLong(&ind) && 0 <= ind&&ind < NY) { box.y2 = ind; }
	else valid = false;
	tc = (wxTextCtrl*)(FindWindowByName("m_textCtrl_selectedMaxZ"));
	sbox = tc->GetLineText(0);
	if (sbox.ToLong(&ind) && 0 <= ind&&ind < NZ) { box.z2 = ind; }
	else valid = false;
	if (!valid)
	{
		wxMessageBox("Wrong selection range!");
		return;
	}

	renders();
}

void GLView::OnCalcUncertainty(wxCommandEvent& event)
{
	if (uncertainty.empty()) return;
	wxString str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_uncertainty_threshold")))->GetLineText(0);
	double threshold = 10;
	if (!str.ToDouble(&threshold) || threshold > 100)
	{
		wxMessageBox("Invalid threshold value");
		return;
	}
	threshold *= 0.01;
	int n = 0;
	double err = 0;
	int NVoxel = dose.getInnerLength();
	double maxDose = dose.Max();
	for (int i = 0; i < NVoxel; ++i)
	{
		if (dose[i] > threshold*maxDose)
		{
			++n;
			err += uncertainty[i] * uncertainty[i];
		}
	}
	str.Printf("Average uncertainty above %.1f%% threshold = %.2f%%", threshold*100, sqrt(err / n)*100);
	wxMessageBox(str);
}

void GLView::OnSetSliceColor(wxCommandEvent& event)
{
	wxColourDialog colorPicker(this);
	if (colorPicker.ShowModal() == wxID_OK)
	{
		wxColour c = colorPicker.GetColourData().GetColour();
		//get the alpha value
		wxSlider* slider = (wxSlider *)FindWindowById(XRCID("m_slider_sliceAlpha"));
		int sv = slider->GetValue();
		setSliceRGBA(c.Red() / 255.0f, c.Green() / 255.0f, c.Blue() / 255.0f, sv / 100.0f);
		wxButton* but = (wxButton *)FindWindowById(XRCID("m_button_sliceColor"));
		but->SetBackgroundColour(c);
	}
}
//<<------------------------------End event functions-------------------------------->>//


//<<------------------------------Processing functions------------------------------->>//
void GLView::save(wxString path, int ctype)
{
// 	wxBusyCursor busyCursor;
// 	wxWindowDisabler disabler;
// 	wxBusyInfo busyInfo(wxT("saving files, wait please..."));
	if (0 == ctype || 1 == ctype)
	{
		if (phantom.empty())
		{
			wxMessageBox("should have phantom first");
			return;
		}

		BinaryFile BF;
		double checkCode = CHECKCODE;
		BF.push("checkCode", &checkCode, sizeof(double));
		BF.push("CT_kVp", &CT_kVp, sizeof(double));
		BF.push("Hist", &Hist, sizeof(int));
		BF.push("uniform", &uniform, sizeof(int));
		BF.push("NX", &NX, sizeof(int));
		BF.push("NY", &NY, sizeof(int));
		BF.push("NZ", &NZ, sizeof(int));
		BF.push("DX", &DX, sizeof(double));
		BF.push("DY", &DY, sizeof(double));
		BF.push("DZ", &DZ, sizeof(double));
		BF.push("xo", &xo, sizeof(double));
		BF.push("yo", &yo, sizeof(double));
		BF.push("zo", &zo, sizeof(double));
		BF.push("Bs", &Bs, sizeof(double));
		BF.push("Bx", &Bx, sizeof(double));
		BF.push("By", &By, sizeof(double));
		BF.push("Bz", &Bz, sizeof(double));
		BF.push("COPX", &COPX, sizeof(double));
		BF.push("COPY", &COPY, sizeof(double));
		BF.push("COPZ", &COPZ, sizeof(double));
		BF.push("stdValue", &stdValue, sizeof(stdValue));
		BF.push("treatmentFraction", &treatmentFraction, sizeof(treatmentFraction));
		BF.push("phantom", (void*)phantom.getP(), sizeof(SFloat)*NX*NY*NZ);

		if (0 == ctype)
		{
			if (dose.empty())
			{
				wxMessageBox("have no dose to save");
				return;
			}
			BF.push("dose", (void*)dose.getP(), sizeof(SFloat)*NX*NY*NZ);
			if (uncertainty.getInnerLength() > 0) BF.push("uncertainty", (void*)uncertainty.getP(), sizeof(SFloat)*NX*NY*NZ);
		}
		if (!BF.write(path.c_str()))
		{
			wxMessageBox("can not write the output file!");
			return;
		}
	}
	else //export ViewRay compatible file (it export mat instead of dose)
	{
		//output ViewRay file format
		VIEWRAY_FORMAT VR_File;
		VR_File.nx = NX;
		VR_File.ny = NY;
		VR_File.nz = NZ;
		VR_File.dx = DX;
		VR_File.dy = DY;
		VR_File.dz = DZ;
		VR_File.Marjor_Header = 1; //the version header may change in the future
		VR_File.Minor_Header = 0;
		VR_File.offset_x = COPX + DX / 2; //the position in dicom/patient coordinate system
		VR_File.offset_y = COPY + DY / 2;
		VR_File.offset_z = COPZ + DZ / 2;
		if (2 == ctype)
		{
			VR_File.m = *mat;
			if (mat->empty())
			{
				wxMessageBox("have no dose to save");
				return;
			}
		}
		else if (3 == ctype)
		{
			if (phantom.empty())
			{
				wxMessageBox("have no phantom to save");
				return;
			}
			VR_File.m = phantom;
		}
		VR_File.write(path);
	}

	wxMessageBox("File is saved successfully!");
}

bool GLView::loadDose(const wxString& path)
{
	int ans = wxYES;
	if (!dose.empty())
	{
		ans = wxMessageBox("One dose file has been loaded. Yes to overlap it, No to load to doseRef, Cancel to abort loading.",
			"how to load dose?", wxYES_NO | wxCANCEL);
		if (wxCANCEL == ans) return false;
	}
	
	//wxBusyCursor busyCursor;
	//wxWindowDisabler disabler;
	//wxBusyInfo busyInfo(wxT("loading files, wait please..."));
	
	if (wxYES == ans)
	{
		BinaryFile BF;
		if (BF.read(path.c_str(), true)) // it's my dose format
		{
			BF.read(path.c_str());
			double checkCode = CHECKCODE;
			BF.get("checkCode", &checkCode);
			if (checkCode != CHECKCODE) return false;
			if (BF.getKeyLength("Hist") == sizeof(int))
			{
				int NHist;
				BF.get("Hist", &NHist);
				Hist = NHist;
			}
			else BF.get("Hist", &Hist); //key used double variable

			BF.get("Bs", &Bs);
			BF.get("Bx", &Bx);
			BF.get("By", &By);
			BF.get("Bz", &Bz);
			BF.get("COPX", &COPX);
			BF.get("COPY", &COPY);
			BF.get("COPZ", &COPZ);
			
			if (!BF.get("treatmentFraction", &treatmentFraction)) treatmentFraction = 1;

			int nTet = 0;
			BF.get("nTet", &nTet);
			tet.resize(nTet, 5);
			BF.get("elem", tet.getP());
			int nNode = 0;
			BF.get("nNode", &nNode);
			vert.resize(nNode, 3);
			BF.get("node", vert.getP());
			int nFace = 0;
			BF.get("nFace", &nFace);
			if (nFace > 0)
			{
				f2t.resize(nFace, 4);
				BF.get("f2t", f2t.getP());
			}
			else f2t.release();

			dose.resize(nTet);
			BF.get("dose", (void*)dose.getP());

			if (BF.exist("uncertainty"))
			{
				uncertainty.resize(nTet);
				BF.get("uncertainty", (void*)uncertainty.getP());
			}
			phantom.resize(nTet);
			if (!BF.get("phantom", (void*)phantom.getP())) phantom = 1;

			if (!BF.get("prescriptionDose", &stdValue) || stdValue <= 0) stdValue = dose.Max();
			addMarkPoint("ISOCenter", MarkPoint(-xo, -yo, -zo));
		}
		else return false;


		//record current file path
		frame->SetTitle(path);
		currentFile = path;
	}
	else return loadExtDose(path);//load the dose to doseRef	

	return true;
}

void GLView::showBF(BinaryFile& BF)
{
	double checkCode = CHECKCODE;
	BF.get("checkCode", &checkCode);
	if (checkCode != CHECKCODE) return;
	BF.get("CT_kVp", &CT_kVp);
	if (BF.getKeyLength("Hist") == sizeof(int))
	{
		int NHist;
		BF.get("Hist", &NHist);
		Hist = NHist;
	}
	else BF.get("Hist", &Hist); //key used double variable

	BF.get("uniform", &uniform);
	BF.get("NX", &NX);
	BF.get("NY", &NY);
	BF.get("NZ", &NZ);
	BF.get("DX", &DX);
	BF.get("DY", &DY);
	BF.get("DZ", &DZ);
	BF.get("xo", &xo);
	BF.get("yo", &yo);
	BF.get("zo", &zo);
	BF.get("Bs", &Bs);
	BF.get("Bx", &Bx);
	BF.get("By", &By);
	BF.get("Bz", &Bz);
	BF.get("COPX", &COPX);
	BF.get("COPY", &COPY);
	BF.get("COPZ", &COPZ);

	if (!BF.get("treatmentFraction", &treatmentFraction)) treatmentFraction = 1;

	phantom.resize(NX, NY, NZ);
	BF.get("phantom", (void*)phantom.getP());

	if (BF.exist("dose"))
	{
		dose.resize(NX, NY, NZ);
		BF.get("dose", (void*)dose.getP());
	}

	if (BF.exist("uncertainty"))
	{
		uncertainty.resize(NX, NY, NZ);
		BF.get("uncertainty", (void*)uncertainty.getP());
	}

	if (!dose.empty() && (!BF.get("prescriptionDose", &stdValue) || stdValue <= 0)) stdValue = dose.Max();
	addMarkPoint("ISOCenter", MarkPoint(-xo, -yo, -zo));
	extraSettings();
	wxPostEvent(this, wxPaintEvent()); //usually this function will called by other thread, but openGL is attached to the 
	//main thread, so the screen wouldn't be updated. Post this message to main thread to force it redrawing in openGL context
}

bool GLView::loadExtDose(const wxString& path, int ftype)
{
	if (dose.empty())
	{
		wxMessageBox("Should load one dose matrix first");
		return false;
	}

	BinaryFile BF;
	if (BF.read(path.c_str(), true)) //it's my dose format, requiring that the lattice is exactly the same
	{
		BF.read(path.c_str());

		int nTet = 0;
		BF.get("nTet", &nTet);
		if (nTet != tet.getWidth(1)) return false;
		int nNode = 0;
		BF.get("nNode", &nNode);
		if (nNode != vert.getWidth(1)) return false;

		doseRef.resize(nTet);
		BF.get("dose", (void*)doseRef.getP());
	}
	
	//change the label on buttons
	wxStaticText* edLoad =  (wxStaticText *)FindWindowById(XRCID("m_staticText_extDoseLoaded"));
	wxStaticText* edOn = (wxStaticText *)FindWindowById(XRCID("m_staticText_extDoseOn"));
	wxButton* edShow = (wxButton *)FindWindowById(XRCID("m_button_showExtDose"));
	edLoad->SetBackgroundColour(*wxGREEN);
	edLoad->SetForegroundColour(*wxWHITE);
	edLoad->SetLabelText("  YES  ");
	edOn->SetBackgroundColour(*wxGREEN);
	edOn->SetForegroundColour(*wxWHITE);
	edOn->SetLabelText("  YES  ");
	edShow->SetLabelText("Show Internal");
	
	wxNotebook* tabs = (wxNotebook*)FindWindowById(XRCID("m_notebook_panels"));
	tabs->SetSelection(PANEL_EXTERNAL);
	b_showDoseRef = true;
	return true;
}

bool GLView::loadViewRayBeams(const wxString& path)
{
	if (VRBeams.getViewRayBeams(path.ToStdString()))
	{
		double x, y, z;
		VRBeams.getISO(x, y, z);
		xo = COPX - x;
		yo = COPY - y;
		zo = COPZ - z;
		addMarkPoint("ISOCenter", MarkPoint(-xo, -yo, -zo));
		renders();
		return true;
	}
	
	return false;
}

void GLView::showCustomizedPhantom(const wxString& path)
{
	BinaryFile BF;
	//must have I/O data
	int NX, NY, NZ; //voxel number
	double DX, DY, DZ; // voxel size, unit cm
	ArrayMgr<SFloat> ph; //relative density matrix; input

	//optional I/O data
	ArrayMgr<short> matid; //material id matrix; default 0
	int nBatch; // default 0
	double Hist; // default 0
	double Bs; //magnetic field strength; default 0
	double Bx, By, Bz; //unit magnetic field direction; default (0, 0, -1)
	double xo, yo, zo; // position of the cuboid corner in isocenter coordinate system; determine where's the isocenter
	double COPX, COPY, COPZ;//DICOM coordinates info; determine where's the patient coordinate system's origin
	double prescriptionDose; // used in DoseViewer; default maxDose
	int treatmentFraction; // used in DoseViewer; default 1
	double trimDensityThreshold; // default 0
	int rngSeed;

	//auxiliary data, which can be regenerated by the dose file
	double rf; //radius factor= 100/C/B
	int uniform;
	double LX, LY, LZ; // side length Lx=DX*NX, not need to output

	//read NX,NY,NZ; DX,DY,DZ(LX, LY, LZ); density matrix
	ConfigFile gcf(path.c_str());
	ConfigFile* cf = gcf.getBlock("PHANTOM");

	bool b_fromFile = true;
	string fname;
	if (!cf->getValue("phantomName", fname)) b_fromFile = false; //no phantom file specified

	cf->getValue("trim dose threshold", trimDensityThreshold);
	if (trimDensityThreshold < 0) trimDensityThreshold = 0;

	double isoPX = 0, isoPY = 0, isoPZ = 0;//iso center position in patient coordinate
	//first try to load from the file
	if (b_fromFile)
	{
		BinaryFile BF;
		if (BF.read(fname.c_str())) return;//try to load infomation from *.phtm file
		else //try to load RED file then
		{
			VIEWRAY_FORMAT VF;
			if (VF.read(fname.c_str()))
			{
				NX = VF.nx;
				NY = VF.ny;
				NZ = VF.nz;
				DX = VF.dx;
				DY = VF.dy;
				DZ = VF.dz;
				COPX = VF.offset_x - DX / 2;
				COPY = VF.offset_y - DY / 2;
				COPZ = VF.offset_z - DZ / 2;
				ph = VF.m;
				matid.resize(NX, NY, NZ);
				matid = 278; //default water
			}
			else
			{
				Log("Warning: Cannot load phantom from file. Will try to use customized phantom in config file");
				b_fromFile = false;
			}
		}
	}

	if (!b_fromFile)//need to get regular phantom from user's definition. Patient's coordinate center is the cuboid center
	{
		vector<double> dc;
		cf->getValue("DX DY DZ", dc);
		DX = dc[0]; DY = dc[1]; DZ = dc[2];

		cf->resetSearchIndex();
		cf->getValue("matrix", dc, true);
		NX = int(dc[0]); NY = int(dc[1]); NZ = int(dc[2]);
		if (NX <= 0 || NY <= 0 || NZ <= 0)
		{
			wxMessageBox("Incorrect voxel dimensions!");
			return;
		}
		//the origin of the patient coordinates lies at the center of the cuboid,
		COPX = -NX*DX / 2;
		COPY = -NY*DY / 2;
		COPZ = -NZ*DZ / 2;
		isoPX = -dc[3];
		isoPY = -dc[4];
		isoPZ = -dc[5];
		xo = COPX - isoPX;
		yo = COPY - isoPY;
		zo = COPZ - isoPZ;
		ph.resize(NX, NY, NZ);
		ph = (dc.size() >= 7 ? SFloat(dc[6]) : SFloat(1)); //1 means use the default density of the corresponding material
		matid.resize(NX, NY, NZ);
		matid = (dc.size() >= 8 ? short(dc[7]) : short(278));

		while (cf->getValue("matrix", dc, true)) //fill inside with different density
		{
			int snx = int(dc[0]), sny = int(dc[1]), snz = int(dc[2]);
			if (snx > NX || sny > NY || snz > NZ)
			{
				wxMessageBox("Incorrect voxel dimensions!");
				return;
			}
			//need to find the origin index of the new cuboid
			int cix = lround((dc[3] - snx*DX*0.5 - xo) / DX);
			int ciy = lround((dc[4] - sny*DY*0.5 - yo) / DY);
			int ciz = lround((dc[5] - snz*DZ*0.5 - zo) / DZ);

			int pix, piy, piz;
			for (int i = 0; i < snx; ++i)
				for (int j = 0; j < sny; ++j)
					for (int k = 0; k < snz; ++k)
					{
						pix = i + cix;
						piy = j + ciy;
						piz = k + ciz;
						if (0 <= pix && pix < NX && 0 <= piy && piy < NY && 0 <= piz && piz < NZ)
						{
							ph.a(pix, piy, piz) = (dc.size() >= 7 ? SFloat(dc[6]) : SFloat(1)); //1 means use the default density of the corresponding material
							matid.a(pix, piy, piz) = (dc.size() >= 8 ? short(dc[7]) : short(278)); //default water
						}
					}
		}

		//may adjust the isoPX, isoPY, isoPZ
		double SSD = 0;
		if (cf->getValue("SSD", SSD)) isoPY = SAD - SSD - NY*DY / 2;
		vector<double> centerXZ;
		if (cf->getValue("centerXZ", centerXZ))
		{
			isoPX = -centerXZ[0];
			isoPZ = -centerXZ[1];
		}

		//different from the definition in phantom.h, we need to convert from relative density to absolute density
		PENELOPE_MAT penMat;
		wxString cwd = wxGetCwd();
		wxSetWorkingDirectory(wxPathOnly(wxStandardPaths::Get().GetExecutablePath()));
		if (!penMat.load())
		{
			wxMessageBox("Cannot load the file named PENELOPE_densities.txt in the program directory");
			wxSetWorkingDirectory(cwd);
			return;
		}
		wxSetWorkingDirectory(cwd);
		int LEN = NX*NY*NZ;
		for (int i = 0; i < LEN; ++i)
		{
			ph.a(i) = ph.a(i)*penMat.density(matid.a(i));
		}
	}

	//it's correct regardless b_fromFile
	xo = COPX - isoPX;
	yo = COPY - isoPY;
	zo = COPZ - isoPZ;

	LX = NX*DX; LY = NY*DY; LZ = NZ*DZ;

	//overwrite couch density into above phantom
	string couchFile;
	if (cf->getValue("couchConfig", couchFile)) //has the density information
	{
		FILE* fp = fopen(couchFile.c_str(), "rb");
		if (NULL == fp){
			wxMessageBox("Cannot open the couch density file you assigned"); return;
		}
		int CNX, CNY, CNZ;
		float CDX, CDY, CDZ;
		float offset_x, offset_y, offset_z;
		fread(&CNX, sizeof(int), 1, fp);
		fread(&CNY, sizeof(int), 1, fp);
		fread(&CNZ, sizeof(int), 1, fp);
		fread(&CDX, sizeof(float), 1, fp);
		fread(&CDY, sizeof(float), 1, fp);
		fread(&CDZ, sizeof(float), 1, fp);
		if (fabs(CDX - DX) > 1e-5 || fabs(CDY - DY) > 1e-5 || fabs(CDZ - DZ) > 1e-5) {
			wxMessageBox("Inconsistent couch voxel size to the phantom"); return;
		}
		fread(&offset_x, sizeof(float), 1, fp);
		fread(&offset_y, sizeof(float), 1, fp);
		fread(&offset_z, sizeof(float), 1, fp);
		ArrayMgr<SFloat> cMat(CNX, CNY, CNZ);
		float elem = 0;
		for (int k = 0; k < CNZ; ++k)
			for (int j = 0; j < CNY; ++j)
				for (int i = 0; i < CNX; ++i)
				{
					fread(&elem, sizeof(float), 1, fp);
					cMat.a(i, j, k) = SFloat(elem);
				}
		fclose(fp);

		//find the couch's location
		vector<double> pos;
		cf->getValue("couchXYZ", pos);
		if (pos.size() != 3) {
			wxMessageBox("Incorrect couch position!"); return;
		}
		//fit the couch density matrix to above phantom

		//first locate the index of couch's origin
		int cix = lround((pos[0] - CNX*CDX / 2 - xo) / DX);
		int ciy = lround((pos[1] - yo) / DY);
		int ciz = lround((pos[2] - CNZ*CDZ / 2 - zo) / DZ);
		//now the indx(couch)+cix == indx(phantom), and so forth
		//we may need to extend z direction of couch to fit above phantom
		int pix, piy, piz;
		if (ciz > 0) //need a face copy
		{
			for (int piz = 0; piz < ciz; ++piz)
			{
				for (int ix = 0; ix < CNX; ++ix)
					for (int iy = 0; iy < CNY; ++iy)
					{
						pix = ix + cix;
						piy = iy + ciy;
						if (pix >= 0 && pix < NX &&piy >= 0 && piy < NY)
							ph.a(pix, piy, piz) = cMat.a(ix, iy, 0);
					}
			}
		}
		if (ciz + CNZ < NZ) //need another face copy
		{
			for (int piz = ciz + CNZ; piz < NZ; ++piz)
			{
				for (int ix = 0; ix < CNX; ++ix)
					for (int iy = 0; iy < CNY; ++iy)
					{
						pix = ix + cix;
						piy = iy + ciy;
						if (pix >= 0 && pix < NX &&piy >= 0 && piy < NY)
							ph.a(pix, piy, piz) = cMat.a(ix, iy, CNZ - 1);
					}
			}
		}
		//do the major replacement
		for (int ix = 0; ix < CNX; ++ix)
			for (int iy = 0; iy < CNY; ++iy)
				for (int iz = 0; iz < CNZ; ++iz)
				{
					pix = ix + cix;
					piy = iy + ciy;
					piz = iz + ciz;
					if (pix >= 0 && pix < NX &&piy >= 0 && piy < NY && piz >= 0 && piz < NZ)
						ph.a(pix, piy, piz) = cMat.a(ix, iy, iz);
				}
	}

	//may change the couch density
	double couchDensityFactor = 1.0;
	if (cf->getValue("couch density factor", couchDensityFactor))
	{
		int ty = NY - 45;//scan for the position of ViewRa's nominal couch
		for (; ty >= 0; --ty) //the couch's height takes about 48 voxels
		{
			if (fabs(ph.a(NX / 2, ty, NZ / 2) - 1.819) < 1e-4 && fabs(ph.a(NX / 2, ty + 1, NZ / 2) - 0.0778) < 1e-4) break;
		}
		if (ty != -1) //we got the couch top position
		{
			for (int iy = ty; iy <= ty + 48 && iy < NY; ++iy) //the couch's height takes about 48 voxels
			{
				for (int ix = 0; ix < NX; ++ix)
					for (int iz = 0; iz < NZ; ++iz)
					{
						ph.a(ix, iy, iz) *= SFloat(couchDensityFactor);
					}
			}
		}
		else Log("Warning: failed to modify the couch's density. Make sure this couch exist.");
	}
	else Log("Warning: cannot get couch density factor in config file; will not modify the couch density!");

	//read configuration of the magnetic field
	cf->getValue("magnetic field strength", Bs);
	if (Bs <= 0) rf = 0;
	else rf = 100 / (lightSpeed*Bs);

	vector<double> dir;
	cf->getValue("magnetic field direction", dir);
	if (dir.size() == 2) //give two polar angles: theta and phi
	{
		dir[0] *= PI / 180;
		dir[1] *= PI / 180;
		Bx = sin(dir[0])*cos(dir[1]);
		if (fabs(Bx) < 1e-8) Bx = 0;
		By = sin(dir[0])*sin(dir[1]);
		if (fabs(By) < 1e-8) By = 0;
		Bz = cos(dir[0]);
		if (fabs(Bz) < 1e-8) Bz = 0;
	}
	else if (dir.size() == 3) //give unit vector
	{
		if (fabs(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2] - 1.0) > 1e-8) { wxMessageBox("Wrong unit vector!"); return; }
		Bx = dir[0];
		By = dir[1];
		Bz = dir[2];
	}
	else { wxMessageBox("Wrong magnetic field direction parameters!"); return; }

	int LEN = NX*NY*NZ;
	BF.push("NX", &NX, sizeof(int));
	BF.push("NY", &NY, sizeof(int));
	BF.push("NZ", &NZ, sizeof(int));
	BF.push("DX", &DX, sizeof(double));
	BF.push("DY", &DY, sizeof(double));
	BF.push("DZ", &DZ, sizeof(double));
	BF.push("phantom", ph.getP(), sizeof(SFloat)*LEN);
	BF.push("matid", matid.getP(), sizeof(short)*LEN);

	BF.push("Hist", &Hist, sizeof(double));
	BF.push("nBatch", &nBatch, sizeof(int));
	BF.push("Bs", &Bs, sizeof(double));
	BF.push("Bx", &Bx, sizeof(double));
	BF.push("By", &By, sizeof(double));
	BF.push("Bz", &Bz, sizeof(double));
	BF.push("xo", &xo, sizeof(double)); //tell DoseViewer where is the ISO center
	BF.push("yo", &yo, sizeof(double));
	BF.push("zo", &zo, sizeof(double));
	BF.push("COPX", &COPX, sizeof(double));//keep original DICOM coordinate(may load extra DICOM data later)
	BF.push("COPY", &COPY, sizeof(double));
	BF.push("COPZ", &COPZ, sizeof(double));

	BF.push("prescriptionDose", &prescriptionDose, sizeof(double));
	BF.push("treatmentFraction", &treatmentFraction, sizeof(int));
	BF.push("trimDensityThreshold", &trimDensityThreshold, sizeof(double));
	BF.push("rngSeed", &rngSeed, sizeof(int)); // used for continuous simulation

	//optional output; just for viewer, don't read it in function readDose()
	BF.push("uniform", &uniform, sizeof(int));

	if (!dose.empty())
	{
		if (!dose.equalDims(ph)) return;
	}
	showBF(BF);
}

void GLView::extraSettings() //extra setting work after loading file
{
	if(!b_showDoseRef) setMatrix(&dose);
	else setMatrix(&doseRef);
	setBGImage(&phantom);
	setMesh(&vert, &tet, &f2t);
	initShow(true);
	wxString info;
	info.Printf("Image Information:\nNZ=%d , DZ=%.2f cm\nNX=%d , DX=%.2f cm\nNY=%d, DY=%.2f cm\nVolume =%.2f*%.2f*%.2f cm^3",
		NZ, DZ, NX, DX, NY, DY, NX*DX, NY*DY, NZ*DZ);
	if (Hist > 0)
	{
		wxString histStr;
		//histStr.Printf("\nhist num=%.0f", Hist);
		histStr.Printf("\nhist num=%s", printWithComma(Hist).c_str());
		info += histStr;
	}
	if (!dose.empty()&&!phantom.empty())
	{
		SFloat tE = 0;
		int len = NX*NY*NZ;
		for (int i = 0; i < len; ++i)
		{
			tE += dose.a(i)*phantom.a(i);
		}
		tE *= SFloat(DX*DY*DZ*1e-3);
		wxString totE;
		totE.Printf("\ntotal Energy =%g J", tE);
		info += totE;
	}
	
	((wxStaticText*)FindWindowById(XRCID("m_staticText_info")))->SetLabelText(info);
	((wxChoice*)FindWindowById(XRCID("m_choice_view_direction")))->SetSelection(ZAXIS);

	if (stdValue != 0)
	{
		info.Printf("%.2f", stdValue);
		((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_prescriptionDose")))->SetLabelText(info);
	}
}

void GLView::processCmdLine(int argc, wxChar ** argv)
{
	if (2 == argc) //only process one file
	{
		wxString path(argv[1]);
		wxFileName fname(path);
		bool success = true;
		if (wxFileExists(path))
		{
			wxString ext = fname.GetExt();
			if (ext == "mdose") success = loadDose(path);
			else success = false;
		}
		else success = false;
		if (success) extraSettings();
		else wxMessageBox("failed to open this file due to wrong format or other reason");
	}
}

bool GLView::hasExtDose()
{
	return !doseRef.empty();
}

SFloat GLView::samplePointDose(double x, double y, double z)
{
	//check if the points is in the phantom
	x -= xo;
	y -= yo;
	z -= zo;
	if (x<0 || x>NX*DX || y<0 || y>NY*DY || z<0 || z>NZ*DZ)
	{
		return -1; //it's an absurd wrong value
	}
	//choose the voxel center as the origin now
	x -= DX / 2;
	y -= DY / 2;
	z -= DZ / 2;
	int ix = int(floor(x / DX));
	int iy = int(floor(y / DY));
	int iz = int(floor(z / DZ));
	x = x/DX - ix; //decimal part
	y = y/DY - iy;
	z = z/DZ - iz;
	//fix the boundary access problem
	if (ix < 0) { ix = 0; x = 0; }
	if (ix >= NX - 1) { ix = NX - 2; x = 1; }
	if (iy < 0) { iy = 0; y = 0; }
	if (iy >= NY - 1) { iy = NY - 2; y = 1; }
	if (iz < 0) { iz = 0; z = 0; }
	if (iz >= NZ - 1) { iz = NZ - 2; z = 1; }

// 	//only for test
// 	wxString info;
// 	info.Printf("%f, %f, %f, %f, %f, %f, %f, %f", dose.a(ix, iy, iz), dose.a(ix + 1, iy, iz), dose.a(ix, iy + 1, iz), dose.a(ix + 1, iy + 1, iz),
// 		dose.a(ix, iy, iz + 1), dose.a(ix + 1, iy, iz + 1), dose.a(ix, iy + 1, iz + 1), dose.a(ix + 1, iy + 1, iz + 1));
// 	wxMessageBox(info);

	//do the interpolation
	SFloat y1 = 0, y2 = 0;
	y1 = mat->a(ix, iy, iz)*(1 - x) + mat->a(ix + 1, iy, iz)*x;
	y2 = mat->a(ix, iy + 1, iz)*(1 - x) + mat->a(ix + 1, iy + 1, iz)*x;
	SFloat z1 = y1*(1 - y) + y2*y;

	y1 = mat->a(ix, iy, iz + 1)*(1 - x) + mat->a(ix + 1, iy, iz + 1)*x;
	y2 = mat->a(ix, iy + 1, iz + 1)*(1 - x) + mat->a(ix + 1, iy + 1, iz + 1)*x;
	SFloat z2 = y1*(1 - y) + y2*y;

	return z1*(1 - z) + z2*z;
}

VERTEX3D GLView::searchExpDose(double x, double y, double z, double expDose, int NScan, int NAround)
{
	const double dsx = 1.0 / double(NScan);
	const double dsy = 1.0 / double(NScan);
	const double DX2 = DX*DX;
	const double DY2 = DY*DY;
	const double DZ2 = DZ*DZ;
	VERTEX3D nearest(1e10,1e10,1e10); //a really far initial distance

	x -= xo;
	y -= yo;
	z -= zo;
	if (x<0 || x>NX*DX || y<0 || y>NY*DY || z<0 || z>NZ*DZ)
	{
		return nearest;
	}
	//choose the voxel center as the origin now
	x -= DX / 2;
	y -= DY / 2;
	z -= DZ / 2;
	int ix = int(floor(x / DX));
	int iy = int(floor(y / DY));
	int iz = int(floor(z / DZ));
	x = x / DX - ix; //decimal part
	y = y / DY - iy;
	z = z / DZ - iz;
	//fix the boundary access problem
	if (ix < 0) { ix = 0; x = 0; }
	if (ix >= NX - 1) { ix = NX - 2; x = 1; }
	if (iy < 0) { iy = 0; y = 0; }
	if (iy >= NY - 1) { iy = NY - 2; y = 1; }
	if (iz < 0) { iz = 0; z = 0; }
	if (iz >= NZ - 1) { iz = NZ - 2; z = 1; }

	double minDist2 = 1e20;
	//search around voxel
	for (int irx = ix - NAround; irx <= ix + NAround; ++irx)
		for (int iry = iy - NAround; iry <= iy + NAround; ++iry)
			for (int irz = iz - NAround; irz <= iz + NAround; ++irz)
			{
				if (irx < 0 || irx >= NX - 1 || iy < 0 || iy >= NY - 1 || iz < 0 || iz >= NZ - 1) continue;

				//find the max and min value
				double maxv = 0, minv = 0;
				maxv = minv = mat->a(irx, iry, irz);
				for (int i = 0; i < 2; ++i)
					for (int j = 0; j < 2; ++j)
						for (int k = 0; k < 2; ++k)
						{
							if (mat->a(irx + i, iry + j, irz + k) > maxv) maxv = mat->a(irx + i, iry + j, irz + k);
							if (mat->a(irx + i, iry + j, irz + k) < minv) minv = mat->a(irx + i, iry + j, irz + k);
						}
				
				if (minv <= expDose && expDose <= maxv) //it lays within a cube
				{
					for (double sx = 0; sx <= 1; sx += dsx)
						for (double sy = 0; sy <= 1; sy += dsy)
						{
							SFloat y1 = 0, y2 = 0;
							y1 = mat->a(irx, iry, irz)*(1 - sx) + mat->a(irx + 1, iry, irz)*sx;
							y2 = mat->a(irx, iry + 1, irz)*(1 - sx) + mat->a(irx + 1, iry + 1, irz)*sx;
							SFloat z1 = y1*(1 - sy) + y2*sy;

							y1 = mat->a(irx, iry, irz + 1)*(1 - sx) + mat->a(irx + 1, iry, irz + 1)*sx;
							y2 = mat->a(irx, iry + 1, irz + 1)*(1 - sx) + mat->a(irx + 1, iry + 1, irz + 1)*sx;
							SFloat z2 = y1*(1 - sy) + y2*sy;

							int isx = irx - ix;
							int isy = iry - iy;
							int isz = irz - iz;
							if (z1 == z2)
							{
								if (z1 != expDose) continue;
								//else
								double dist2 = (sx + isx - x)*(sx + isx - x)*DX2 + (sy + isy - y)*(sy + isy - y)*DY2;
								if (dist2 < minDist2)
								{
									minDist2 = dist2;
									nearest.x = (sx + isx - x)*DX;
									nearest.y = (sy + isy - y)*DY;
									nearest.z = 0;
								}
							}
							else
							{
								double sz = (expDose - z1) / (z2 - z1);
								if (sz < 0 || sz>1) continue;
								//else
								double dist2 = (sx + isx - x)*(sx + isx - x)*DX2 + (sy + isy - y)*(sy + isy - y)*DY2 + (sz + isz - z)*(sz + isz - z)*DZ2;
								if (dist2 < minDist2)
								{
									minDist2 = dist2;
									nearest.x = (sx + isx - x)*DX;
									nearest.y = (sy + isy - y)*DY;
									nearest.z = (sz + isz - z)*DZ;
								}
							}
						}
					
				}
				else continue;
			}
	//if (minDist2 > 1e10) return -1; //cannot find any points
	return nearest;
}
//<<-----------------------------End processing functions----------------------------->>// 

void GLView::render(int wid) // do additional drawing
{
	if (NULL == mat && NULL == bgImage) return;
	View3D::render(wid);
	FloatRect& glRect = subw[wid].glRect;
	FloatRect& picRect = subw[wid].picRect;
	glViewport(glRect.x, glRect.y, glRect.width, glRect.height);
	glScissor(glRect.x, glRect.y, glRect.width, glRect.height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, glRect.width, 0, glRect.height, -10000.0f, 10000.0f); // define the 3D space in pixels
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	// translate the origin to the center of picRect
	glTranslatef(picRect.centerX(), picRect.centerY(), 0);
	// change the look-at angle
	double R = 1000;
	gluLookAt(R*pv[wid].eyev[0], R*pv[wid].eyev[1], R*pv[wid].eyev[2], 0.0f, 0.0f, 0.0f, pv[wid].upv[0], pv[wid].upv[1], pv[wid].upv[2]);
	// need to calculate the ratio to scale from cm to pixel
	float scale = getScale(wid);
	glScalef(scale, scale, scale);
	double lx = NX*DX;
	double ly = NY*DY;
	double lz = NZ*DZ;
	glTranslatef(-lx / 2, -ly / 2, -lz / 2);
	if (ZAXIS == wid && b_showBeams&& VRBeams.beams.size()>0)
	{
		if (markPoints.count("ISOCenter") <= 0) return; // must have iso-center coordinates first
		glColor3f(1.0f, 1.0f, 0.0f);
		glPushMatrix();
		vector<BEAM>& beams = VRBeams.beams;
		int nb = beams.size();
		glTranslatef(markPoints["ISOCenter"].x, markPoints["ISOCenter"].y, markPoints["ISOCenter"].z);
		float BoreRadius = 41; //unit cm
		glScalef(BoreRadius, BoreRadius, BoreRadius);
		//check the line width supported
		float width_range[2];
		glGetFloatv(GL_LINE_WIDTH_RANGE, width_range);
		double maxf, minf;
		VRBeams.getMaxMinFlux(maxf, minf);
		float k = 0;
		if (maxf != minf) k = (width_range[1] - 1) / (maxf - minf);

		for (int i = 0; i < nb; ++i)
		{
			if (maxf == minf) glLineWidth(lineWidth);
			else glLineWidth((beams[i].flux - minf)*k + 1);
			glBegin(GL_LINES);
			glVertex2f(0, 0);
			double angle = beams[i].gantryAngle*PI / 180;
			glVertex2f(sin(angle), -cos(angle));
			glEnd();
		}

		//draw the other half projection
		glLineWidth(lineWidth);
		glColor3f(0.5f, 0.5f, 0.5f);
		glLineStipple(1, 0x0F0F);
		glEnable(GL_LINE_STIPPLE);
		glBegin(GL_LINES);
		for (int i = 0; i < nb; ++i)
		{
			glVertex2f(0, 0);
			double angle = beams[i].gantryAngle*PI / 180;
			glVertex2f(-sin(angle), cos(angle));
		}
		glEnd();
		glDisable(GL_LINE_STIPPLE);

		//draw the bore circle
		glColor3f(0.5, 0.5, 0.5);
		glBegin(GL_LINE_LOOP);
		int NS = 180;
		for (int i = 0; i < NS; ++i)
		{
			double angle = 2 * i*PI / NS;
			glVertex2f(cos(angle), sin(angle));
		}
		glEnd();
		glPopMatrix();
	}
}