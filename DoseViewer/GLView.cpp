#include "PreCompile.h"
#include <atomic>
//#include "opencv2/imgproc/imgproc.hpp"
//#include "opencv2/core/core.hpp"

#include "GLView.h"
#include "Graph.h"
#include "Gamma3D.h"
#include "CTOptionDlg.h"
#include "DataEditDlg.h"
#include "../GPUGamma/GPUGamma.h"
#include "bitset"

GLView* ggl = NULL; //global pointer of GLView
bool CmdMode = false;
LogProvider Log;
#define CHECKCODE 3.14 //used to check the file format

//<<C style function
double ConvertFunc(int ctno, double kVp)//convert CT # to density
{
	const int Na = 13;//arrary length
	const int x1[] = { 0, 197, 493, 928, 971, 994, 1049, 1051, 1243, 1887, 2371, 2858, 4015 }; // 120 kVp
	const int x2[] = { 0, 222, 512, 934, 970, 988, 1036, 1040, 1201, 1768, 2122, 2621, 3933 }; // 110 kVp
	const int x3[] = { 0, 185, 489, 927, 972, 994, 1045, 1052, 1279, 2052, 2553, 3205, 4068 }; // 100 kVp
	const double y[] = { 0, 0.195, 0.51, 0.96, 0.991, 1.016, 1.062, 1.072, 1.161, 1.53, 1.82, 2.15, 4.51 };

	const int * x = NULL;
	if (kVp == 120) x = x1;
	else if (kVp == 110) x = x2;
	else if (kVp == 100) x = x3;
	else exitApp("no data base for the incident KVP.");

	//binary search where ctno is
	int i = 0, j = Na - 1, m = 0;
	do
	{
		m = (i + j) / 2;
		if (ctno > x[m]) i = m;
		else j = m;
	} while (j - i > 1);

	return (y[j] - y[i]) / (double)(x[j] - x[i]) * (double)(ctno - x[i]) + y[i];
}
int InConvertFunc(double den, double kVp)
{
	const int Na = 13;//arrary length
	const int x1[] = { 0, 197, 493, 928, 971, 994, 1049, 1051, 1243, 1887, 2371, 2858, 4015 }; // 120 kVp
	const int x2[] = { 0, 222, 512, 934, 970, 988, 1036, 1040, 1201, 1768, 2122, 2621, 3933 }; // 110 kVp
	const int x3[] = { 0, 185, 489, 927, 972, 994, 1045, 1052, 1279, 2052, 2553, 3205, 4068 }; // 100 kVp
	const double y[] = { 0, 0.195, 0.51, 0.96, 0.991, 1.016, 1.062, 1.072, 1.161, 1.53, 1.82, 2.15, 4.51 };

	const int * x = NULL;
	if (kVp == 120) x = x1;
	else if (kVp == 110) x = x2;
	else if (kVp == 100) x = x3;
	else exitApp("no data base for the incident KVP.");

	int i = 0, j = Na - 1, m = 0;
	do
	{
		m = (i + j) / 2;
		if (den > y[m]) i = m;
		else j = m;
	} while (j - i > 1);

	return int((x[j] - x[i]) / (y[j] - y[i]) * (den - y[i]) + x[i]);

}
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
EVT_BUTTON(XRCID("m_button_trimDose"), GLView::OnTrimDose)
EVT_BUTTON(XRCID("m_button_showExtDose"), GLView::OnShowExtDose)
EVT_BUTTON(XRCID("m_button_delExtDose"), GLView::OnDelExtDose)
EVT_BUTTON(XRCID("m_button_viewDoseDifference"), GLView::OnViewDoseDifference)
EVT_BUTTON(XRCID("m_button_swapDose"), GLView::OnSwapDose)
EVT_BUTTON(XRCID("m_button_gammaAnalysis"), GLView::OnGammaAnalysis)
EVT_BUTTON(XRCID("m_button_showRED"), GLView::OnShowRED)
EVT_BUTTON(XRCID("m_button_delRED"), GLView::OnDelRED)
EVT_BUTTON(XRCID("m_button_averageDose"), GLView::OnAverageDose)
EVT_BUTTON(XRCID("m_button_averagePhantom"), GLView::OnAveragePhantom)
EVT_BUTTON(XRCID("m_button_showPreviousGamma"), GLView::OnShowPreviousGamma)
EVT_BUTTON(XRCID("m_button_setPrescriptionDose"), GLView::OnSetPrescriptionDose)
EVT_BUTTON(XRCID("m_button_samplePointDose"), GLView::OnSamplePointDose)
EVT_BUTTON(XRCID("m_button_getExpDTA"), GLView::OnGetExpDTA)
EVT_BUTTON(XRCID("m_button_sampleLineDose"), GLView::OnSampleLineDose)
EVT_BUTTON(XRCID("m_button_samplePlanDose"), GLView::OnSamplePlanDose)
EVT_BUTTON(XRCID("m_button_trimDoseByDensity"), GLView::OnTrimDoseByDensity)
EVT_BUTTON(XRCID("m_button_flipDose"), GLView::OnFlipDose)
EVT_BUTTON(XRCID("m_button_exportDiff"), GLView::OnExportDiff)
EVT_BUTTON(XRCID("m_button_updateSelected"), GLView::OnUpdateSelected)
EVT_BUTTON(XRCID("m_button_calc_uncertainty"), GLView::OnCalcUncertainty)
EVT_MENU(XRCID("m_menuItem_open"), GLView::OnOpen)
EVT_MENU(XRCID("m_menuItem_openDir"), GLView::OnOpenDir)
EVT_MENU(XRCID("m_menuItem_save"), GLView::OnSave)
EVT_MENU(XRCID("m_menuItem_saveSlice"), GLView::OnSaveSlice)
EVT_MENU(XRCID("m_menuItem_zoom_in"), GLView::OnZoomIn)
EVT_MENU(XRCID("m_menuItem_zoom_out"), GLView::OnZoomOut)
EVT_MENU(XRCID("m_menuItem_rotate_left"), GLView::OnRotateLeft)
EVT_MENU(XRCID("m_menuItem_rotate_right"), GLView::OnRotateRight)
EVT_TOOL(XRCID("m_tool_open"), GLView::OnOpen)
EVT_TOOL(XRCID("m_tool_openDir"), GLView::OnOpenDir)
EVT_TOOL(XRCID("m_tool_openContours"), GLView::OnOpenContours)
EVT_TOOL(XRCID("m_tool_openRED"), GLView::OnOpenRED)
EVT_TOOL(XRCID("m_tool_beams"), GLView::OnOpenBeams)
EVT_TOOL(XRCID("m_tool_importDose"), GLView::OnImportDose)
EVT_TOOL(XRCID("m_tool_save"), GLView::OnSaveAs)
EVT_TOOL_RCLICKED(XRCID("m_tool_save"), GLView::OnSave)
EVT_TOOL(XRCID("m_tool_saveSlice"), GLView::OnSaveSlice)
EVT_TOOL(XRCID("m_tool_adjustDose"), GLView::OnAdjustDose)
EVT_TOOL(XRCID("m_tool_pdd"), GLView::OnPDD)
EVT_TOOL(XRCID("m_tool_doseProfile"), GLView::OnDoseProfile)
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
EVT_CHECKBOX(XRCID("m_checkBox_showCT"), GLView::OnShowCT)
EVT_CHECKBOX(XRCID("m_checkBox_showDose"), GLView::OnShowDose)
EVT_CHECKBOX(XRCID("m_checkBox_showPhantom"), GLView::OnShowPhantom)
EVT_CHECKBOX(XRCID("m_checkBox_showAbsE"), GLView::OnShowAbsEnergy)
EVT_CHECKBOX(XRCID("m_checkBox_showBox"), GLView::OnShowBox)
EVT_CHECKBOX(XRCID("m_checkBox_showISO"), GLView::OnShowISOCenter)
EVT_CHECKBOX(XRCID("m_checkBox_showBeams"), GLView::OnShowBeams)
EVT_CHECKBOX(XRCID("m_checkBox_globalContrast"), GLView::OnGlobalContrast)
EVT_CHECKBOX(XRCID("m_checkBox_autoContrast"), GLView::OnAutoContrast)
EVT_CHECKBOX(XRCID("m_checkBox_showContours"), GLView::OnShowStructContours)
EVT_CHECKBOX(XRCID("m_checkBox_relativeISOLines"), GLView::OnShowRelativeISOLines)
EVT_CHOICE(XRCID("m_choice_view_direction"), GLView::OnViewDirection)
EVT_CHOICE(XRCID("m_choice_doseLevel"), GLView::OnDoseLevel)
EVT_CHOICE(XRCID("m_choice_colorScheme"), GLView::OnColorScheme)
EVT_SLIDER(XRCID("m_slider_slice_index"), GLView::OnSliderIndex)
EVT_SLIDER(XRCID("m_slider_lineWidth"), GLView::OnSliderLineWidth)
EVT_SLIDER(XRCID("m_slider_matAlpha"), GLView::OnSlidermatAlpha)
EVT_SLIDER(XRCID("m_slider_stdColor"), GLView::OnSliderstdColor)
EVT_SLIDER(XRCID("m_slider_brightness"), GLView::OnSliderBrightness)
EVT_SLIDER(XRCID("m_slider_C"), GLView::OnSliderCenter)
EVT_SLIDER(XRCID("m_slider_W"), GLView::OnSliderWidth)
EVT_SLIDER(XRCID("m_slider_L"), GLView::OnSliderLeft)
EVT_SLIDER(XRCID("m_slider_R"), GLView::OnSliderRight)
EVT_CHECKLISTBOX(XRCID("m_checkList_contours"), GLView::OnCheckCountour)
EVT_TEXT(XRCID("m_textCtrl_pointX3"), GLView::OnXPlanDistanceChanged)
EVT_TEXT(XRCID("m_textCtrl_pointY3"), GLView::OnYPlanDistanceChanged)
EVT_TEXT(XRCID("m_textCtrl_pointZ3"), GLView::OnZPlanDistanceChanged)

EVT_COMMAND(wxID_ANY, VIEW3D_NOTIFY, GLView::OnView3DNotify)
END_EVENT_TABLE()

GLView::GLView(wxWindow* parent, wxFrame* pframe) :View3D(parent)
{
	frame = pframe;
	wxButton* pb_trimDose = (wxButton*)FindWindowById(XRCID("m_button_trimDose"));
	pb_trimDose->Bind(wxEVT_RIGHT_DOWN, &GLView::OnTrimLVOL, this); //bind the right click event

	imcInfo = (wxStaticText *)FindWindowById(XRCID("m_staticText_imc"));
	contourBox = (wxCheckListBox*)FindWindowById(XRCID("m_checkList_contours"));

	b_showCT = false;
	b_showAbsE = false;
	b_showExtRED = false;
	b_showBeams = true;
	b_showStructContours = true;
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
		if (wxDirExists(path)) //it contains all CT images
		{
			wxArrayString names;
			wxDir::GetAllFiles(path, &names, "*.dcm", wxDIR_FILES);
			loadCT(names);
			success = false;
		}
	}
	else //it's a file
	{
		if (wxFileExists(path))
		{
			wxString ext = fname.GetExt();
			if (ext == "dose") success = loadDose(path);
			else if (ext == "phtm") success = loadPhantom(path);
			else if (ext == "dcm")
			{
				DICOM_TYPE dcmType = getDicomType(path);
				if (dcmType == DCM_RTSTRUCT) success = loadContours(path);
				else if (dcmType == DCM_RTDOSE) success = loadDose(path);
				else success = false;
			}
			else if (ext == "red") success = loadRED(path);
			else if (ext == "beams") success = loadViewRayBeams(path);
			else if (ext == "txt") success = loadViewRayBeams(path);
			else if (ext == "py")
			{
				showCustomizedPhantom(path);
				return;
			}
			else if (ext == "lvol")
			{
				loadViewRayLVol(path);
				return;
			}

		}
	}
	
	if(!success) wxMessageBox("Cannot open this file due to unsupported format or other reasons.");
	else extraSettings();
}

void GLView::OnOpen(wxCommandEvent& event)
{
	enum{ CT, DOSE, PHANTOM, RED_VIEWRAY, BEAMS_VIEWRAY, PLANOVERVIEW_VIEWRAY }; //must have same order as the wildcard
	wxString caption = wxT("Choose a file");
	wxString wildcard = wxT("CT files(*.dcm)|*.dcm|Dose files(*.dose)|*.dose|Phantom files(*.phtm)|*.phtm|ViewRay RED file(*.red)|*.red|Plan Overview text(.txt)|*.txt");
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
		else if (PHANTOM == ftype)
		{
			if (!loadPhantom(dialog.GetPath()))
			{
				wxMessageBox("The phantom file doesn't have expected format");
				return;
			}
		}
		else if (CT == ftype)
		{
			wxArrayString names;
			dialog.GetPaths(names);

			loadCT(names);
		}
		else if (RED_VIEWRAY == ftype)
		{
			if (!loadRED(dialog.GetPath()))
			{
				wxMessageBox("The ViewRay RED phantom file doesn't have expected format");
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

void GLView::OnOpenDir(wxCommandEvent& event)
{
	wxString defaultPath = wxEmptyString;
	wxDirDialog dialog(this, wxT("Pick a directory containing all the CT dicom file you want to load"),
		defaultPath, wxDD_NEW_DIR_BUTTON);
	if (dialog.ShowModal() == wxID_OK)
	{
		wxString path = dialog.GetPath();
		wxArrayString names;
		wxDir::GetAllFiles(path, &names, "*.dcm", wxDIR_FILES);

		loadCT(names);
		extraSettings();
	}
}

void GLView::OnOpenContours(wxCommandEvent& event)
{
	wxString caption = wxT("Choose a structure set dicom file");
	wxString wildcard = wxT("structure set files(*.dcm)|*.dcm");
	wxString defaultDir = wxEmptyString;
	wxString defaultFilename = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_OPEN);
	if (dialog.ShowModal() == wxID_OK)
	{
		wxString path = dialog.GetPath();
		loadContours(path);
	}
}

void GLView::OnOpenRED(wxCommandEvent& event)
{
	static wxStaticText* edLoad = (wxStaticText *)FindWindowById(XRCID("m_staticText_REDLoaded"));
	if (edLoad->GetLabelText() == "  YES  ") // we already have external dose
	{
		wxMessageBox("The external phantom has to be deleted before loading another one");
		return;
	}

	wxString caption = wxT("Choose a RED file");
	wxString wildcard = wxT("relative electron density(*.red)|*.red");
	wxString defaultDir = wxEmptyString;
	wxString defaultFilename = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_OPEN);
	if (dialog.ShowModal() == wxID_OK)
	{
		wxString path = dialog.GetPath();
		if (!loadRED(path))
		{
			wxMessageBox("Cannot load the RED file!");
			return;
		}
	}
}

void GLView::OnOpenBeams(wxCommandEvent& event)
{
// 	wxStaticText* edLoad = (wxStaticText *)FindWindowById(XRCID("m_staticText_extDoseLoaded"));
// 	if (edLoad->GetLabelText() == "  YES  ") // we already have external dose
// 	{
// 		wxMessageBox("The external dose has to be deleted before loading another one");
// 		return;
// 	}

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

void GLView::OnImportDose(wxCommandEvent& event) //replace dose or add dose to phantom
{
	if (phantom.getInnerLength() <= 0)
	{
		wxMessageBox("Need to load phantom first!");
		return;
	}
	wxString caption = wxT("Choose a dose file");
	wxString wildcard = wxT("Dose files(*.dose)|*.dose");
	wxString defaultDir = wxEmptyString;
	wxString defaultFilename = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_OPEN);
	if (dialog.ShowModal() == wxID_OK)
	{
		BinaryFile BF;
		BF.read(dialog.GetPath().c_str());
		double checkCode = CHECKCODE;
		BF.get("checkCode", &checkCode);
		if (checkCode != CHECKCODE)
		{
			wxMessageBox("Incorrect program/data version!");
			return;
		}
		int dnx, dny, dnz;
		BF.get("NX", &dnx);
		BF.get("NY", &dny);
		BF.get("NZ", &dnz);

		ArrayMgr<SFloat> ds;
		ds.resize(dnx, dny, dnz);
		BF.get("dose", (void*)ds.getP());
		
		//now let's try to fit ds to dose
		bool match = true;
		int nshrink = 2;
		if (NZ != dnz) match = false;
		if (NX == 2 * dnx&&NY == 2 * dny) nshrink = 2;
		else if (NX == dnx&&NY == dny) nshrink = 1;
		else match = false;
		if (!match)
		{
			wxMessageBox("Unmatched dose matrix, please make sure they're from the same CT image set");
			return;
		}

		//do the replacement
		if (1 == nshrink) dose = ds;
// 		else
// 		{
// 			dose.resize(NX, NY, NZ); //make sure we have enough memory
// 			cv::Mat smallPic(dnx, dny, CV_32FC1);
// 			cv::Mat bigPic(NX, NY, CV_32FC1);
// 			for (int iz = 0; iz < NZ; ++iz)
// 			{
// 				//fetch one small picture
// 				for (int ix = 0; ix < dnx; ++ix)
// 				{
// 					for (int iy = 0; iy < dny; ++iy)
// 					{
// 						smallPic.at<SFloat>(ix, iy) = ds.a(ix, iy, iz);
// 					}
// 				}
// 				try
// 				{
// 					cv::pyrUp(smallPic, bigPic);
// 				}
// 				catch (cv::Exception& e)
// 				{
// 					const char* err_msg = e.what();
// 					wxString info;
// 					info.Printf("openCV exception: %s", err_msg);
// 					wxMessageBox(info);
// 				}
// 				
// 				for (int ix = 0; ix < NX; ++ix)
// 				{
// 					for (int iy = 0; iy < NY; ++iy)
// 					{
// 						dose.a(ix, iy, iz) = bigPic.at<SFloat>(ix, iy);
// 					}
// 				}	
// 			}
// 		}
		
		show();
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
		if (VIEW3D == activeID) return;
		int selection = dialog.GetFilterIndex();
		ArrayMgr<SFloat> data;
		if (0 == selection) getSlice(activeID, *bgImage, data);
		else getSlice(activeID, *mat, data);
		int dimx = data.getWidth(1);
		int dimy = data.getWidth(2);

		wxString path = dialog.GetPath();
		if (2 == selection) //call Matlab to write the dicom dose
		{
			//need cd to module directory
			wxString cwd = wxGetCwd();
			wxSetWorkingDirectory(wxPathOnly(wxStandardPaths::Get().GetExecutablePath()));
			double dw = 0.3, dh = 0.3;
			getDWDH(activeID, dw, dh);
			string fname = path.ToStdString();
			
			//get the position of the upper left pixel
			double dox = COPX + DX / 2, doy = COPY + DY / 2, doz = COPZ + DZ / 2;
			if (XAXIS == activeID)
			{
				dox += cross[XAXIS] * DX;
			}
			else if (YAXIS == activeID)
			{
				doy += cross[YAXIS] * DY;
				doz += (NZ-1)*DZ; //there's only (NZ-1)*DZ shift
			}
			else if (ZAXIS == activeID)
			{
				doz += cross[ZAXIS] * DZ;
			}
			
			//get Matlab column first layout data
			float* mdata = new float[dimx*dimy];
			for (int i = 0; i < dimx; ++i)
			{
				for (int j = 0; j < dimy; ++j)
				{
					mdata[i + dimx*j] = data.a(i, j);
				}
			}

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
			if (activeID == ZAXIS)
			{
				const double iop[6] = { 1, 0, 0, 0, 1, 0 };
				df.setItem("ImageOrientationPatient", MatlabType(iop, 6));
			}
			else if (activeID == YAXIS)
			{
				const double iop[6] = { 1, 0, 0, 0, 0, -1 };
				df.setItem("ImageOrientationPatient", MatlabType(iop, 6));
			}
			else if (activeID == XAXIS)
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
			if (fp == NULL) exitApp("cannot create the output file!");

			for (int iy = 0; iy < dimy; ++iy)
			{
				for (int ix = 0; ix < dimx; ++ix)
				{
					if (ix != dimx - 1) fprintf(fp, "%f,", data.a(ix, iy));
					else fprintf(fp, "%f\n", data.a(ix, iy));
				}
			}
			fclose(fp);
			wxMessageBox("Matlab compatible File is saved successfully!");
		}
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

void GLView::OnPDD(wxCommandEvent& event)
{
	if (NULL == mat || mat->empty()) return;

	wxString caption = wxT("Input a file name");
	wxString defaultFilename = wxT("PDD");
	wxString wildcard = wxT("text files(*.txt)|*.txt");
	wxString defaultDir = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_SAVE);
	if (dialog.ShowModal() == wxID_OK)
	{
		double rx = (NX - 1) / 2.0;
		double rz = (NZ - 1) / 2.0;
		int ix = (int)rx;
		int iz = (int)rz;
		rx -= ix;
		rz -= iz;
		double fiz, fiz1;
		ArrayMgr<SFloat> pdd(NY);
		for (int iy = 0; iy < NY; ++iy)
		{
			//for iz
			fiz = (1 - rx)*mat->a(ix, iy, iz) + rx*mat->a(ix + 1, iy, iz);

			//for iz+1
			fiz1 = (1 - rx)*mat->a(ix, iy, iz + 1) + rx*mat->a(ix + 1, iy, iz + 1);
			pdd.a(iy) = (1 - rz)*fiz + rz*fiz1;
		}
		SFloat maxv, minv;
		pdd.getMaxMin(maxv, minv);

		FILE *fp = fopen(dialog.GetPath().c_str(), "w");
		if (NULL == fp) exitApp("Cannot create PDD file!");
		for (int iy = 0; iy < NY; ++iy)
		{
			fprintf(fp, "%.2f\t%0.4f\t%0.4f\n", (iy + 0.5)*DY, pdd.a(iy), pdd.a(iy) / maxv*100);
		}
		fclose(fp);
		wxMessageBox("PDD file has been exported!");
	}
}

void GLView::OnDoseProfile(wxCommandEvent& event)
{
	wxString str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_profileDepth")))->GetLineText(0);
	double depth = 0;
	if (!str.ToDouble(&depth) || depth < 0)
	{
		wxMessageBox("incorrect double format of the depth");
		return;
	};

	if (NULL == mat || mat->empty()) return;

	if (activeID != ZAXIS && activeID != XAXIS)
	{
		wxMessageBox("Please check you select following window:\nactive window == ZAXIS: export x profile\nactive window == XAXIS: export z profile");
		return;
	}

	wxString caption = wxT("Input a file name");
	wxString defaultFilename = wxT("Profile.txt");
	wxString wildcard = wxT("text files(*.txt)|*.txt");
	wxString defaultDir = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_SAVE);
	if (dialog.ShowModal() == wxID_OK)
	{
		ArrayMgr<SFloat> profile;
		if (ZAXIS == activeID)
		{
			double rz = (NZ - 1) / 2.0;
			int iz = (int)rz;
			rz -= iz;
			double ry = depth / DY - 0.5;
			int iy = (int)ry;
			ry -= iy;
			profile.resize(NX);
			double f0, f1;
			for (int ix = 0; ix < NX; ++ix)
			{
				//for iy line
				f0 = (1 - rz)*mat->a(ix, iy, iz) + rz*mat->a(ix, iy, iz + 1);
				//for iy+1 line
				f1 = (1 - rz)*mat->a(ix, iy + 1, iz) + rz*mat->a(ix, iy + 1, iz + 1);

				profile.a(ix) = (1 - ry)*f0 + ry*f1;
			}
		}
		else if (XAXIS == activeID)
		{
			double rx = (NX - 1) / 2.0;
			int ix = (int)rx;
			rx -= ix;
			double ry = depth / DY - 0.5;
			int iy = (int)ry;
			ry -= iy;
			profile.resize(NZ);
			double f0, f1;
			for (int iz = 0; iz < NZ; ++iz)
			{
				//for iy line
				f0 = (1 - rx)*mat->a(ix, iy, iz) + rx*mat->a(ix + 1, iy, iz);
				//for iy+1 line
				f1 = (1 - rx)*mat->a(ix, iy + 1, iz) + rx*mat->a(ix + 1, iy + 1, iz);

				profile.a(iz) = (1 - ry)*f0 + ry*f1;
			}
		}
		SFloat maxv, minv;
		profile.getMaxMin(maxv, minv);

		FILE *fp = fopen(dialog.GetPath().c_str(), "w");
		if (NULL == fp) exitApp("Cannot create profile file!");
		if (ZAXIS == activeID)
		{
			for (int ix = 0; ix < NX; ++ix)
			{
				fprintf(fp, "%.2f\t%0.4f\t%0.4f\n", (ix - NX / 2.0 + 0.5)*DX, profile.a(ix), profile.a(ix) / maxv * 100);
			}
			wxMessageBox("Profile in x direction has been exported");
		}
		else if (XAXIS == activeID)
		{
			for (int iz = 0; iz < NZ; ++iz)
			{
				fprintf(fp, "%.2f\t%0.4f\t%0.4f\n", (iz - NZ / 2.0 + 0.5)*DZ, profile.a(iz), profile.a(iz) / maxv * 100);
			}
			wxMessageBox("Profile in z direction has been exported");
		}

		fclose(fp);
	}
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

void GLView::OnShowCT(wxCommandEvent& event)
{
	b_showCT = !b_showCT;
	int NPht = phantom.getInnerLength();
	if (b_showCT)
	{
		for (int i = 0; i < NPht; ++i) phantom.a(i) = InConvertFunc(phantom.a(i), CT_kVp);
	}
	else
	{
		for (int i = 0; i < NPht; ++i) phantom.a(i) = ConvertFunc(phantom.a(i), CT_kVp);
	}
	show();
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

void GLView::OnCheckCountour(wxCommandEvent& event) //check which contour to show
{
	int sel = event.GetSelection();
	strt_set[sel].show = !strt_set[sel].show;
	renders();
}
void GLView::OnShowStructContours(wxCommandEvent& event) //switch for structure contour 
{
	showStructContours();
	renders();
}

void GLView::OnTrimDose(wxCommandEvent& event)
{
	
}
void GLView::OnTrimLVOL(wxMouseEvent& event)
{
	wxString caption = wxT("Choose a ViewRay lvol file");
	wxString wildcard = wxT("labled volume mask file(*.lvol)|*.lvol");
	wxString defaultDir = wxEmptyString;
	wxString defaultFilename = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_OPEN);
	if (dialog.ShowModal() == wxID_OK)
	{
		wxString path = dialog.GetPath();

		FILE *fp = fopen(path.c_str(), "rb");
		if (NULL == fp)
		{
			wxMessageBox("Cannot open lvol file!");
			return;
		}

		unsigned int Major_Header, Minor_Header;
		unsigned int rnx, rny, rnz;
		double rdx, rdy, rdz;
		double off_x, off_y, off_z;
		fread(&Major_Header, sizeof(unsigned int), 1, fp);
		fread(&Minor_Header, sizeof(unsigned int), 1, fp);
		fread(&rnx, sizeof(unsigned int), 1, fp);
		fread(&rdx, sizeof(double), 1, fp);
		fread(&off_x, sizeof(double), 1, fp);
		fread(&rny, sizeof(unsigned int), 1, fp);
		fread(&rdy, sizeof(double), 1, fp);
		fread(&off_y, sizeof(double), 1, fp);
		fread(&rnz, sizeof(unsigned int), 1, fp);
		fread(&rdz, sizeof(double), 1, fp);
		fread(&off_z, sizeof(double), 1, fp);

		unsigned int* vrMask = new unsigned int[rnx*rny*rnz];
		fread(vrMask, sizeof(unsigned int)*rnx*rny*rnz, 1, fp);
		fclose(fp);

		//do the mask fitting
		if (fabs(rdx - DX) > 1e-7 || fabs(rdy - DY) > 1e-7 || fabs(rdz - DZ) > 1e-7)
		{
			wxMessageBox("Unmatched voxel size");
			return;
		}

		ArrayMgr<unsigned int> mask;
		mask.resize(NX, NY, NZ);
		mask = 0;
		int six = (off_x - COPX) / DX; //shift index number
		int siy = (off_y - COPY) / DY;
		int siz = (off_z - COPZ) / DZ;
		unsigned int rix, riy, riz; //index in the external dose matrix

		for (int iz = (siz > 0 ? siz : 0); iz < NZ; ++iz)
			for (int iy = (siy > 0 ? siy : 0); iy < NY; ++iy)
				for (int ix = (six > 0 ? six : 0); ix < NX; ++ix)
				{
					rix = (ix - six);
					riy = (iy - siy);
					riz = (iz - siz);
					if (rix < rnx && riy < rny && riz < rnz)
					{
						mask.a(ix, iy, iz) = vrMask[rix + rnx*riy + rnx*rny*riz];
					}
				}
		delete[] vrMask;

		//do the trimming work
		for (int ix = 0; ix < NX; ++ix)
			for (int iy = 0; iy < NY; ++iy)
				for (int iz = 0; iz < NZ; ++iz)
				{
					if (0 == mask.a(ix, iy, iz))
					{
						dose.a(ix, iy, iz) = 0;
						uncertainty.a(ix, iy, iz) = 0;
					}
				}
		show();
	}
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

	vector<int> IRIDs;
	bool b_IROnly = ((wxCheckBox*)FindWindowById(XRCID("m_checkBox_ImportantRegionOnly")))->IsChecked();
	if (b_IROnly)
	{
		str = ((wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_ImportantRegionIDs")))->GetLineText(0);
		string num = "";
		vector<string> nums;
		for (int i = 0; i < str.size(); ++i)
		{
			if (str[i] != ';') num += str[i];
			else
			{
				nums.push_back(num);
				num = "";
			}
		}
		if (nums.size() <= 0 || nums.size() > 32)
		{
			wxMessageBox("Please give the ids of the important region, separated by semicolons");
			return;
		}
		IRIDs.clear();
		int intTemp;
		for (int i = 0; i < nums.size(); ++i)
		{
			try
			{
				intTemp = stoi(nums[i]);
			}
			catch (std::exception)
			{
				wxMessageBox("incorrect important region id: " + nums[i]);
				return;
			}
			
			if (intTemp < 0 || intTemp >= 32)
			{
				wxMessageBox("incorrect important region id: " + nums[i]);
				return;
			}
			IRIDs.push_back(intTemp);
		}
		if (lvol.getInnerLength() == 0) //need to load the lvol file
		{
			wxMessageBox("Please load the lvol file first");
			return;
		}
	}
	/*******************Start: get the two dose matrix to compare*******************/
	
	
	//do2 always points to the external dose
	ArrayMgr<SFloat>* d1 = NULL;
	ArrayMgr<SFloat>* d2 = NULL;
	//get gammaDose, gammaDoseRef, gammaBGImage
	ArrayMgr<SFloat> gammaDose;
	BOX3D b = box;
	if (!b_onlyBox)
	{
		b.x1 = 0; b.x2 = NX - 1;
		b.y1 = 0; b.y2 = NY - 1;
		b.z1 = 0; b.z2 = NZ - 1;
	}
	else b.sort();
	//exclude the couch region
	if (b.y2 >= couchY) b.y2 = couchY - 1;

	//need to resize new memory
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
	if(b_IROnly) gammaLVol.resize(mx, my, mz);
	else gammaLVol.release();
	for (int iz = 0; iz < mz; ++iz)
		for (int iy = 0; iy < my; ++iy)
			for (int ix = 0; ix < mx; ++ix)
			{
				gammaDose.a(ix, iy, iz) = dose.a(ix + b.x1, iy + b.y1, iz + b.z1);
				gammaDoseRef.a(ix, iy, iz) = doseRef.a(ix + b.x1, iy + b.y1, iz + b.z1);
				gammaBGImage.a(ix, iy, iz) = phantom.a(ix + b.x1, iy + b.y1, iz + b.z1);
				gammaUncertainty.a(ix, iy, iz) = uncertainty.a(ix + b.x1, iy + b.y1, iz + b.z1);
				if (b_IROnly) gammaLVol.a(ix, iy, iz) = lvol.a(ix + b.x1, iy + b.y1, iz + b.z1);
			}
	//adjust gammaLVol according to IRIDs
	std::bitset<32> bitMask;
	for (int j = 0; j < IRIDs.size(); ++j) bitMask.set(IRIDs[j]);
	for (int i = 0; i < gammaLVol.getInnerLength(); ++i)
	{
		std::bitset<32> bitTemp(gammaLVol.a(i));
		bitTemp &= bitMask;
		if (bitTemp.any()) gammaLVol.a(i) = 1;
		else gammaLVol.a(i) = 0;
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
	vg.setLVol(&gammaLVol);
	vg.setFilter(threshold, minDensity);
	vg.setVoxelSize(DX, DY, DZ);
	vg.initShow();
	id_gammaWindow = vg.GetParent()->GetId();
	if (!((wxCheckBox*)FindWindowById(XRCID("m_checkBox_show3DGamma")))->IsChecked()) vg.GetParent()->Hide();

	//count the passing rate, and print it out
	int tNum = 0; //the voxel number to be considered
	int tPass = 0; //the passed voxel number in considered voxels
	vg.getPassingRate(tNum, tPass);
	const char* strGPU = useGPU ? "yes" : "no";
	str.Printf("fitted voxel percent = %f%%, passed number = %d, passing rate = %f%%\n useGPU = %s, time cost = %.2f s", 100 * double(tNum) / gamma.getInnerLength(), tPass, 100 * double(tPass) / tNum, strGPU, tm.stop());
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
		vg.setVoxelSize(DX, DY, DZ);

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
void GLView::OnShowRED(wxCommandEvent& event)
{
	wxStaticText* edLoad = (wxStaticText *)FindWindowById(XRCID("m_staticText_REDLoaded"));
	wxStaticText* edOn = (wxStaticText *)FindWindowById(XRCID("m_staticText_REDOn"));
	wxButton* edShow = (wxButton *)FindWindowById(XRCID("m_button_showRED"));

	if (b_showExtRED) //currently we're showing external phantom, need to turn off
	{
		b_showExtRED = false;
		//exchange the content with phantomRef
		setBGImage(&phantom);
		edOn->SetBackgroundColour(*wxYELLOW);
		edOn->SetForegroundColour(*wxBLACK);
		edOn->SetLabelText("  NO  ");
		edShow->SetLabelText("Show External Phantom");
		show();
	}
	else //currently we're showing internal phantom
	{
		if (edLoad->GetLabelText() == "  YES  ")
		{
			b_showExtRED = true;
			//exchange the content with phantomRef
			setBGImage(&phantomRef);
			edOn->SetBackgroundColour(*wxGREEN);
			edOn->SetForegroundColour(*wxWHITE);
			edOn->SetLabelText("  Yes  ");
			edShow->SetLabelText("Show Internal Phantom");
			show();
		}
		else wxMessageBox("must load external RED file through tool button first");
	}
}
void GLView::OnDelRED(wxCommandEvent& event)
{
	static wxStaticText* edLoad = (wxStaticText *)FindWindowById(XRCID("m_staticText_REDLoaded"));
	static wxStaticText* edOn = (wxStaticText *)FindWindowById(XRCID("m_staticText_REDOn"));
	static wxButton* edShow = (wxButton *)FindWindowById(XRCID("m_button_showRED"));
	if (edLoad->GetLabelText() == "  YES  ") // we have external dose
	{
		if (b_showExtRED) //currently show the external dose;
		{
			edOn->SetBackgroundColour(*wxYELLOW);
			edOn->SetForegroundColour(*wxBLACK);
			edOn->SetLabelText("  NO  ");
			show();
		}
		b_showExtRED = false;
		//need to reset the external loaded label
		edLoad->SetBackgroundColour(*wxYELLOW);
		edLoad->SetForegroundColour(*wxBLACK);
		edLoad->SetLabelText("  NO  ");
		phantomRef.release();
	}
	else wxMessageBox("Don't have external phantom yet!");
}
void GLView::OnAverageDose(wxCommandEvent& event)
{
	if (dose.getInnerLength() > 0)
	{
		wxChoice* cad = (wxChoice *)FindWindowById(XRCID("m_choice_averageDose"));
		int NC = cad->GetSelection() + 2;
		int nax = NX / NC;
		int nay = NY / NC;
		for (int iz = 0; iz < NZ; ++iz)
		{
			for (int ix = 0; ix < nax; ++ix)
			{
				for (int iy = 0; iy < nay; ++iy)
				{
					SFloat sum = 0;
					for (int jx = 0; jx < NC; ++jx)
						for (int jy = 0; jy < NC; ++jy)
							sum += dose.a(ix*NC + jx, iy*NC + jy, iz);

					sum /= NC*NC; //get the average

					for (int jx = 0; jx < NC; ++jx)
						for (int jy = 0; jy < NC; ++jy)
							dose.a(ix*NC + jx, iy*NC + jy, iz) = sum;
				}
			}
		}
		show(); //reload and show
	}
	else wxMessageBox("Should have dose first!");
	
}
void GLView::OnAveragePhantom(wxCommandEvent& event)
{
	if (phantom.getInnerLength() > 0)
	{
		wxChoice* cad = (wxChoice *)FindWindowById(XRCID("m_choice_averagePhantom"));
		int NC = cad->GetSelection() + 2;
		int nax = NX / NC;
		int nay = NY / NC;
		for (int iz = 0; iz < NZ; ++iz)
		{
			for (int ix = 0; ix < nax; ++ix)
			{
				for (int iy = 0; iy < nay; ++iy)
				{
					SFloat sum = 0;
					for (int jx = 0; jx < NC; ++jx)
						for (int jy = 0; jy < NC; ++jy)
							sum += phantom.a(ix*NC + jx, iy*NC + jy, iz);

					sum /= NC*NC; //get the average

					for (int jx = 0; jx < NC; ++jx)
						for (int jy = 0; jy < NC; ++jy)
							phantom.a(ix*NC + jx, iy*NC + jy, iz) = sum;
				}
			}
		}
		show(); //reload and show
	}
	else wxMessageBox("Should have phantom first!");
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
		if (phantom.getInnerLength() > 0)
		{
			wxString den;
			den.Printf(", density = %4.2fg / cm ^ 3", phantom.a(kx, ky, kz));
			sbtext += den;
		}
		
		if (dose.getInnerLength() > 0)
		{
			wxString ds;
			if (b_showDoseRef) ds.Printf(", dose = %5.3f Gray", doseRef.a(kx, ky, kz));
			else ds.Printf(", dose = %5.3f Gray", dose.a(kx, ky, kz));
			sbtext += ds;
		}
		if (uncertainty.getInnerLength() > 0)
		{
			wxString err;
			err.Printf(", uncertainty = %4.2f%%", 100 * uncertainty.a(kx, ky, kz));
			sbtext += err;
		}
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
// 			fig.g->AddLegend(L"Detailed", "b*");
// 			fig.g->AddLegend(L"DVH-constrained", "gd*");
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
				if (ratio < 1) //inside the ion chamber
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
	wxString strx1, strx2, stry1, stry2, strz1, strz2, strds;
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
		if (x - xo<0 || x - xo>NX*DX)
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
				mdata[ix + dimx*(dimy - 1 - iy)] = samplePointDose(DX*(ix + 0.5) + xo, y, DZ*(iy + 0.5) + zo);

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
				df.setItem("ImageOrientationPatient", MatlabType(iop, 6));
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

			df.setPixelData(pin, len * sizeof(uint16_t));
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
	if (!dose.empty() && !doseRef.empty() && dose.swap(doseRef))
	{
		wxButton* swapButton = (wxButton *)FindWindowById(XRCID("m_button_swapDose"));
		if (swapButton->GetBackgroundColour() == *wxYELLOW) swapButton->SetBackgroundColour(*wxWHITE);
		else swapButton->SetBackgroundColour(*wxYELLOW);
		show();
	}
}

void GLView::OnFlipDose(wxCommandEvent& event)
{
	if (dose.empty()) return;
	//check which direction to flip
	int direction = ((wxChoice*)FindWindowById(XRCID("m_choice_view_direction")))->GetSelection();
	if (ZAXIS == direction)
	{
		mat->flip3(3);
	}
	else if (YAXIS == direction)
	{
		mat->flip3(2);
	}
	else if (XAXIS == direction)
	{
		mat->flip3(1);
	}
	show();
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
	str.Printf("For %s,\n average uncertainty above %.1f%% threshold = %.2f%%", currentFile, threshold * 100, sqrt(err / n) * 100);
	// 	if (CmdMode)
	// 	{
	// 		wxPrintf("%s\t%0.2f\n", currentFile, sqrt(err / n) * 100);
	// 	}
	// 	else wxMessageBox(str);
	wxMessageBox(str);
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

void GLView::loadCT(const wxArrayString & names)
{
	frame->SetTitle(wxPathOnly(names[0]));
	currentFile = "";
	if (dose.getInnerLength() > 0)
	{
		int saveDose = wxMessageBox(wxT("Previous dose distribution has been loaded. Do you want to keep it?\n\
										If it's the phantom used in simulation, choose Yes. Otherwise No."),
			wxT("Keep previous dose?"), wxYES_NO);
		if (saveDose != wxYES) dose.release();
	}

	DcmFileFormat fileformat;
	int NSlice = names.size();
	if (NSlice < 3) exitApp("too few slices of image! Cannot construct a reasonable phantom.");
	vector<SliceInfo> slices;
	slices.resize(NSlice); //resize memory in advance
	vector<bool> needRelease(NSlice, true); //used for CT collapsing
	vector<SliceInfo> slices3D(1); //it stores the final 3D matrix

	double slope, intercept;
	int rows, colums;
	double sliceThickness;
	double kVp;

	//above variables may be overwritten by external file.
	bool extCTS = false;
	ConfigFile ext_config;
	if (ext_config.parse("DoseViewer/Settings.py")) extCTS = true;

	{
		wxBusyCursor busyCursor;
		wxWindowDisabler disabler;
		wxBusyInfo busyInfo(wxT("loading files, wait please..."));

		//load CT to slices[]
		for (int is = 0; is < NSlice; ++is)
		{
			DicomFile df;
			if (!df.open(names[is].ToStdString()))
			{
				wxMessageBox("Cannot open the CT dicom fiel");
				return;
			}
			MatlabType info = df.getAllItems();
			slope = info["RescaleSlope"].a(0);
			intercept = info["RescaleIntercept"].a(0);


			rows = info["Rows"].scast<int>();
			colums = info["Columns"].scast<int>();
			sliceThickness = info["SliceThickness"].a(0);

			if (info["KVP"].empty())
			{
				if (!ext_config.getValue("DCM_KVP", kVp)) kVp = 120;
			}
			else kVp = info["KVP"].a(0);
			if (kVp != 100 && kVp != 110 && kVp != 120)
			{
				wxMessageBox("Warning: kVp don't have corresponding converting data.Please Check the CT data or give a default value in external_CT_setting.py");
				return;
			}

			CT_kVp = kVp;//used in InConvertFunc
			double pixelSpacingX = info["PixelSpacing"].a(1);
			double pixelSpacingY = info["PixelSpacing"].a(0);
			slices[is].x = info["ImagePositionPatient"].a(0);
			slices[is].y = info["ImagePositionPatient"].a(1);
			slices[is].z = info["ImagePositionPatient"].a(2);

			int pixpresent = info["PixelRepresentation"].scast<int>();
			if (pixpresent != 0)
			{
				wxMessageBox("Unsupported pixel format. Expected Uint16 as storage");
			}
			Uint16 *pixelData = (Uint16 *)(df.getPixelData());
			long length = rows*colums;
			double *density = new double[length];

			for (long i = 0; i < length; ++i)
			{
				density[i] = ConvertFunc(pixelData[i], kVp);
			}

			slices[is].density = density;
			slices[is].NRow = rows;
			slices[is].NCol = colums;
			slices[is].thickness = sliceThickness;
			slices[is].pixSizeX = pixelSpacingX;
			slices[is].pixSizeY = pixelSpacingY;
		}

		//check the consistency of the parameters
		double minThickness = slices[0].thickness;
		bool even = true;
		for (int is = 1; is < NSlice; ++is)
		{
			if (slices[is].NRow != slices[is - 1].NRow ||
				slices[is].NCol != slices[is - 1].NCol ||
				slices[is].pixSizeX != slices[is - 1].pixSizeX ||
				slices[is].pixSizeY != slices[is - 1].pixSizeY
				) exitApp("inconsistent CT dicom files! must have some common parameters!");
			//find the min thickness if they're not uniform
			if (slices[is].thickness < minThickness)
			{
				even = false;
				minThickness = slices[is].thickness;
			}
		}

		//we need to rearrange the sequence of slices according to z coordinate
		sort(slices.begin(), slices.end());

		if (!even) //uneven images, need interpolation to make the phantom evenly divided.
		{
			exitApp("Uneven CT thickness! Will implement it later...\nPlease contact the author for newer version.");
			// 			slices3D[0] = slices[0]; //copy the content, but didn't copy what the density pointer pointing to
			// 			needRelease[0] = false;
			// 			int iu = 0;
			// 			for (int i = 1; iu < NSlice - 1; ++i)
			// 			{
			// 				double pos = slices[0].z - i*minThickness;
			// 				//find which interval of pos located, search from is
			// 				while (pos < slices[iu + 1].z)
			// 				{
			// 					++iu;
			// 					if (iu + 1 == NSlice) break;
			// 				}
			// 				if (iu + 1 == NSlice) break;
			// 				if (pos == slices[iu + 1].z)
			// 				{
			// 					slices3D.push_back(slices[iu + 1]);
			// 					needRelease[iu + 1] = false;
			// 					++iu;
			// 				}
			// 				else if (pos > slices[iu + 1].z) //do linear interpolation between slices[is] and slices[is+1]
			// 				{
			// 					SliceInfo ins = slices[iu];
			// 					ins.density = new double[rows*colums];
			// 					for (long j = 0; j < rows*colums; ++j)
			// 					{
			// 						ins.density[j] = (slices[iu + 1].density[j] - slices[iu].density[j]) / (slices[iu + 1].z - slices[iu].z)*(pos - slices[iu].z) + slices[iu].density[j];
			// 					}
			// 					ins.z = pos;
			// 					slices3D.push_back(ins);
			// 				}
			// 			}
			// 			for (int is = 0; is < NSlice; ++is) if (needRelease[is]) slices[is].release();
		}
		else slices3D = slices;//already uniform
	}

	// this is the new slice number. It's the same as NSlice if the original CT sets are evenly spaced
	int NL = slices3D.size();
	// pop up a CT compressing option dialog
	CTOptionDlg CTOption(rows, colums, NL);
	CTOption.ShowModal();
	//nShrink: how many pixel compress to one in each CT image
	//NC: how many CT images compressed to one
	int nShrink = 1, NC = 1;
	CTOption.getParameters(nShrink, NC);

	wxBusyCursor busyCursor;
	wxWindowDisabler disabler;
	wxBusyInfo busyInfo(_("processing CT data, wait please..."));

	if (nShrink > 1) //shrink the size of each image
	{
		rows /= nShrink;
		colums /= nShrink;

		for (int i = 0; i < NL; ++i)
		{
			slices3D[i].NRow /= nShrink;
			slices3D[i].NCol /= nShrink;
			slices3D[i].pixSizeX *= nShrink;
			slices3D[i].pixSizeY *= nShrink;

			double *den = new double[rows*colums];//will abandon these uncombined pixels 
			for (int n = 0; n < rows*colums; ++n)
			{
				int ir = n / colums;
				int ic = n - ir*colums;
				int index = 0;
				double sum = 0;
				for (int j = 0; j < nShrink; ++j)
				{
					for (int k = 0; k < nShrink; ++k)
					{
						index = (nShrink*ir + j)*colums*nShrink + (nShrink*ic + k);
						sum += slices3D[i].density[index];
					}
				}
				den[n] = sum / (nShrink*nShrink);
			}
			delete[] slices3D[i].density;
			slices3D[i].density = den;
		}
	}

	vector<SliceInfo> slices3Dp;
	if (1 < NC && NC <= 10) //collapse multiple images into 1
	{
		slices3Dp.resize(NL / NC);
		int length = rows*colums;
		int is = 0;
		for (is = 0; NC*(is + 1) <= NL; ++is)
		{
			slices3Dp[is] = slices3D[NC*is]; //copy the first one
			for (int i = NC*is + 1; i < NC*(is + 1); ++i) //sum the rest
			{
				for (int j = 0; j < length; ++j)	slices3Dp[is].density[j] += slices3D[i].density[j];
				slices3D[i].release();
			}
			for (int j = 0; j < length; ++j) slices3Dp[is].density[j] /= NC;//get average
			slices3Dp[is].thickness *= NC;
		}
		if (NL%NC != 0) //deal the rest undivided slices
		{
			int NRest = NL - NC*is;
			if ((double)NRest / NC >= 0.5) // the rest is bigger than half of the joined-slice's thickness
			{
				slices3Dp.push_back(slices3D[NC*is]);
				for (int i = NC*is + 1; i < NL; ++i) //sum the rest
				{
					for (int j = 0; j < length; ++j) slices3Dp[is].density[j] += slices3D[i].density[j];
					slices3D[i].release();
				}
				for (int j = 0; j < length; ++j) slices3Dp[is].density[j] /= NRest;//get average
				slices3Dp[is].thickness *= NC;
			}
			else //need to release the memory of the rest slices
			{
				for (int i = NC*is; i < NL; ++i) slices3D[i].release();
			}
		}
	}

	//copy data to other data structure
	if (slices3Dp.size() > 0) copyData(slices3Dp);
	else copyData(slices3D);
}

bool GLView::loadPhantom(const wxString& path)
{
	frame->SetTitle(path);
	currentFile = path;
	if (!dose.empty())
	{
		int saveDose = wxMessageBox(wxT("Previous dose distribution has been loaded. Do you want to keep it?\n\
										If it's the phantom used in simulation, choose Yes. Otherwise No."),
			wxT("Keep previous dose?"), wxYES_NO);
		if (saveDose != wxYES) dose.release();
	}

	wxBusyCursor busyCursor;
	wxWindowDisabler disabler;
	wxBusyInfo busyInfo(wxT("loading files, wait please..."));


	BinaryFile BF;
	BF.read(path.c_str());
	double checkCode = CHECKCODE;
	BF.get("checkCode", &checkCode);
	if (checkCode != CHECKCODE) return false;
	BF.get("CT_kVp", &CT_kVp);
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

	phantom.resize(NX, NY, NZ);
	BF.get("phantom", (void*)phantom.getP());

	Hist = 0;
	return true;
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
			dose.resize(NX, NY, NZ);
			BF.get("dose", (void*)dose.getP());

			if (BF.exist("uncertainty"))
			{
				uncertainty.resize(NX, NY, NZ);
				BF.get("uncertainty", (void*)uncertainty.getP());
			}

			if (!BF.get("prescriptionDose", &stdValue) || stdValue <= 0) stdValue = dose.Max();
			addMarkPoint("ISOCenter", MarkPoint(-xo, -yo, -zo));
		}
		else //this should be ViewRay's dose format or dicom format
		{
			CT_kVp = 110;
			Hist = 0;

			VIEWRAY_FORMAT vdose;
			wxFileName fname(path);
			if (fname.GetExt() == "dcm") readFromDicom(vdose, path.c_str());
			else vdose.read(path.c_str());

			if (NX != vdose.nx || NY != vdose.ny || NZ != vdose.nz || fabs(DX - vdose.dx) > 1e-6 || fabs(DY - vdose.ny) > 1e-6 || fabs(DZ - vdose.nz) > 1e-6)
				phantom.release(); //not the right phantom for it

			NX = vdose.nx;
			NY = vdose.ny;
			NZ = vdose.nz;
			DX = vdose.dx;
			DY = vdose.dy;
			DZ = vdose.dz;
			COPX = vdose.offset_x - vdose.dx / 2;
			COPY = vdose.offset_y - vdose.dy / 2;
			COPZ = vdose.offset_z - vdose.dz / 2;
			dose = vdose.m;
			uncertainty = vdose.err;
			if (phantom.empty()) //try to load the phantom in the same directory
			{
				wxString name;
				wxFileName::SplitPath(path, NULL, NULL, &name, NULL);
				wxString fullName = wxPathOnly(path) + wxFileName::GetPathSeparator() + name + ".red";
				if (wxFileExists(fullName) && loadRED(fullName)) setBGImage(&phantom);
				fullName = wxPathOnly(path) + wxFileName::GetPathSeparator() + name + ".beams";
				if (wxFileExists(fullName))
				{
					loadViewRayBeams(fullName);
					((wxCheckBox*)FindWindowById(XRCID("m_checkBox_showBeams")))->SetValue(false);
					b_showBeams = false;
				}
			}
		}

		//record current file path
		frame->SetTitle(path);
		currentFile = path;
	}
	else return loadExtDose(path);//load the dose to doseRef	

	//search the position of the couch
	if (!phantom.empty())
	{
		int ty = NY / 2;
		for (; ty < NY - 1; ++ty)//the couch's height takes about 48 voxels
		{
			int nAir = 0;
			for (int tx = 0; tx < NX; ++tx)
			{
				for (int tz = 0; tz < NZ; ++tz)
				{
					if (phantom.a(tx, ty - 1, tz) < 0.01) ++nAir;
				}
			}
			float airRate = nAir / float(NX*NZ); //whether above layer is air

			SFloat sc = phantom.a(NX / 2, ty, NZ / 2);
			int nPlane = 0;
			for (int tx = 0; tx < NX; ++tx)
				for (int tz = 0; tz < NZ; ++tz)
				{
					if (fabs(phantom.a(tx, ty, tz) - sc) < 1e-4 && sc > 1.2) ++nPlane;
				}
			float planRate = nPlane / float(NX*NZ);//whether this layer is almost uniform

			sc = phantom.a(NX / 2, ty + 1, NZ / 2);
			int nCavity = 0;
			for (int tx = NX / 4; tx < NX * 3 / 4; ++tx)
				for (int tz = NZ / 4; tz < NZ * 3 / 4; ++tz)
				{
					if (fabs(phantom.a(tx, ty + 1, tz) - sc) < 1e-4 && sc < 0.2) ++nCavity;
				}
			float cavityRate = nCavity / float((NX * 3 / 4 - NX / 4)*(NZ * 3 / 4 - NZ / 4));
			if (airRate > 0.6 && planRate > 0.9&& cavityRate == 1) break;
		}
		if (ty == NY - 1) ty = NY; //means no couch found
		couchY = ty;
	}

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

bool GLView::loadContours(const wxString& path)
{
	if (phantom.getInnerLength() <= 0)
	{
		wxMessageBox("Should load phantom or CT images first");
		return false;
	}
	wxNotebook* tabs = (wxNotebook*)FindWindowById(XRCID("m_notebook_panels"));

	//use dcmtk to load the structure set
	DRTStructureSet strt;
	OFCondition err;
	{
		//the loading may take a little while
		wxBusyCursor busyCursor;
		wxWindowDisabler disabler;
		wxBusyInfo busyInfo(wxT("loading structure set infomation, please wait..."));
		err = strt.loadFile(OFFilename(path.c_str()));
	}
	if (err.bad())
	{
		wxMessageBox("Failed to load the structure set dicom file!");
		return false;
	}
	contourBox->Clear();

	//get the names
	for (int i = 1;; ++i)
	{
		DRTStructureSetROISequence::Item &rItem = strt.getROI(i);
		if (rItem.isEmpty()) break;
		OFString value;
		rItem.getROIName(value);
		contourBox->Append(wxString(value.c_str()));
	}

	//get all the contours' pointers
	OFVector<DRTROIContourSequence::Item *> cItems;
	strt.getContoursForROINumber(cItems);

	int numROI = cItems.size();

	strt_set.clear(); //clear previous data
	strt_set.resize(numROI);
	for (int i = 0; i < numROI; ++i)
	{
		strt_set[i].show = false;
		strt_set[i].name = contourBox->GetString(i);
		if (strt_set[i].name == "Skin" || strt_set[i].name == "skin")
		{
			contourBox->Check(i);
			strt_set[i].show = true;
		}
		Sint32 r, g, b;
		cItems[i]->getROIDisplayColor(r, 0);
		cItems[i]->getROIDisplayColor(g, 1);
		cItems[i]->getROIDisplayColor(b, 2);
		strt_set[i].r = r;
		strt_set[i].g = g;
		strt_set[i].b = b;

		//search the contour in sequence and store it in strt_set[i].contour
		DRTContourSequence& seq = cItems[i]->getContourSequence();
		seq.gotoFirstItem();
		do
		{
			DRTContourSequence::Item& item = seq.getCurrentItem();
			if (!item.isEmpty())
			{
				//find one line
				OFVector<double> OFLine;
				item.getContourData(OFLine);
				vector<VERTEX3D> line(OFLine.size() / 3);
				for (size_t j = 0; j < line.size(); ++j)
				{
					line[j].x = OFLine[3 * j] * 0.1;
					line[j].y = OFLine[3 * j + 1] * 0.1;
					line[j].z = OFLine[3 * j + 2] * 0.1;
				}
				//store this line
				strt_set[i].contour.push_back(line);
			}
			err = seq.gotoNextItem();
		} while (err.good());
#ifdef WIN32 // contourBox->GetItem(i) only works on Windows
		contourBox->GetItem(i)->SetBackgroundColour(wxColor(r, g, b));
		contourBox->GetItem(i)->SetTextColour(*wxWHITE);
#endif
		strt_set[i].sortbyz();
	}

	tabs->SetSelection(PANEL_CONTOURS);
	renders();
	return true;
}

bool GLView::loadExtDose(const wxString& path, int ftype)
{
	if (dose.empty())
	{
		wxMessageBox("Should load one dose matrix first");
		return false;
	}

	BinaryFile bf;
	if (bf.read(path.c_str(), true)) //it's my dose format, requiring that the lattice is exactly the same
	{
		bf.read(path.c_str());
		int nx, ny, nz;
		bf.get("NX", &nx);
		bf.get("NY", &ny);
		bf.get("NZ", &nz);
		if (nx != NX || ny != NY || nz != NZ)	return false;

		doseRef.resize(NX, NY, NZ);
		bf.get("dose", (void*)doseRef.getP());
	}
	else //it's ViewRay or DICOM format
	{
		wxWindowDisabler disableAll;
		wxBusyInfo wait("Please wait, loading file...");

		VIEWRAY_FORMAT vdose;
		wxFileName fname(path);
		if (fname.GetExt() == "dcm")
		{
			if (!readFromDicom(vdose, path.c_str())) return false;
		}
		else
		{
			if (!vdose.read(path.c_str())) return false;
		}
		//origin distance between two sets of dose
		double sx = COPX + DX / 2 - vdose.offset_x;
		double sy = COPY + DY / 2 - vdose.offset_y;
		double sz = COPZ + DZ / 2 - vdose.offset_z;

		//for customized dose calculation in zeus, we need to use following origin distance
// 		double sx = xo + DX / 2 - vdose.offset_x;
// 		double sy = yo + DY / 2 - vdose.offset_y;
// 		double sz = zo + DZ / 2 - vdose.offset_z;

		//for DICOM dose or ViewRay's dose
		int ans = wxYES;
		if (vdose.nx == NX&&vdose.ny == NY&&vdose.nz == NZ
			&&fabs(vdose.dx - DX) < 1e-6 &&fabs(vdose.dy - DY) < 1e-6&&fabs(vdose.dz - DZ) < 1e-6
			)
		{
			if (fabs(sx) < 1e-6&&fabs(sy) < 1e-6&&fabs(sz) < 1e-6) doseRef = vdose.m;
			else
			{
				ans = wxMessageBox("Same lattice but different origin. Do you want to do fitting?",
					"whether to fit?", wxYES_NO | wxCANCEL);
				if (wxYES != ans) doseRef = vdose.m;
			}
		}

		if (wxYES == ans) //we need to do voxel fitting
		{
			doseRef.resize(NX, NY, NZ);
			doseRef = 0;
			for (int ix = 0; ix < NX; ++ix)
				for (int iy = 0; iy < NY; ++iy)
					for (int iz = 0; iz < NZ; ++iz)
					{
						doseRef.a(ix, iy, iz) = vdose.interp(sx + ix*DX, sy + iy*DY, sz + iz*DZ);
					}
		}
	}

	//change the label on buttons
	wxStaticText* edLoad = (wxStaticText *)FindWindowById(XRCID("m_staticText_extDoseLoaded"));
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

bool GLView::loadRED(const wxString& path)
{
	if (!RED.read(path.c_str())) return false;

	if (phantom.empty()) //Let's treat it as the primary phantom
	{
		if (dose.empty())
		{
			Hist = 0;

			NX = RED.nx;
			NY = RED.ny;
			NZ = RED.nz;
			DX = RED.dx;
			DY = RED.dy;
			DZ = RED.dz;
			//get the cuboid origin position
			COPX = RED.offset_x - DX / 2;
			COPY = RED.offset_y - DY / 2;
			COPZ = RED.offset_z - DZ / 2;

			phantom = RED.m;
			uniform = RED.uniform;
		}
		else //we have dose loaded previously
		{
			if (NX == RED.nx&&
				NY == RED.ny&&
				NZ == RED.nz&&
				abs(DX - RED.dx) < 1e-6&&
				abs(DY - RED.dy) < 1e-6&&
				abs(DZ - RED.dz) < 1e-6)//replace only when dose and phantom have the same lattice
			{
				phantom = RED.m;
				uniform = RED.uniform;
			}
			else return false; //the grid doesn't fit. We don't use interpolation for primary phantom
		}
		//record current file path
		frame->SetTitle(path);
		currentFile = path;
		return true;
	}
	else //has phantom loaded already
	{
		if (dose.empty()) //can replace the phantom directly
		{
			phantom = RED.m;
			uniform = RED.uniform;
			return true;
		}
		else //need to consider if the dimensions fit with the dose matrix
		{
			if (NX == RED.nx&&
				NY == RED.ny&&
				NZ == RED.nz&&
				DX == RED.dx&&
				DY == RED.dy&&
				DZ == RED.dz) //this phantom can be replaced
			{
				int ans = wxMessageBox("One phantom file has been loaded. Yes to overlap it, No to load to phantomRef, Cancel to abort loading.",
					"how to load dose?", wxYES_NO | wxCANCEL);
				if (wxCANCEL == ans) return false;
				if (wxYES == ans)
				{
					phantom = RED.m;
					uniform = RED.uniform;
					return true;
				}
			}
		}
	}

	//You choosed No in previous question; fit the RED file to the phantomRef
	b_showExtRED = true;
	phantom.resize(NX, NY, NZ);
	phantomRef = 0;
	if ((RED.dz - DZ) > 1e-5) return false;
	int nshrink = 1;
	if (fabs(RED.dx / DX - 2) < 1e-5) nshrink = 2;
	else if (fabs(RED.dx / DX - 1) < 1e-5) nshrink = 1;
	else return false;
	int six = (RED.offset_x - RED.dx / 2 - COPX) / DX;
	int siy = (RED.offset_y - RED.dy / 2 - COPY) / DY;
	int siz = (RED.offset_z - RED.dz / 2 - COPZ) / DZ;
	int rix, riy, riz;
	for (int iz = (siz > 0 ? siz : 0); iz < NZ; ++iz)
		for (int iy = (siy > 0 ? siy : 0); iy < NY; ++iy)
			for (int ix = (six > 0 ? six : 0); ix < NX; ++ix)
			{
				rix = (ix - six) / nshrink;
				riy = (iy - siy) / nshrink;
				riz = (iz - siz);
				if (rix < int(RED.nx) && riy < int(RED.ny) && riz < int(RED.nz)) phantom.a(ix, iy, iz) = RED.m.a(rix, riy, riz);
			}
	wxStaticText* edLoad = (wxStaticText *)FindWindowById(XRCID("m_staticText_REDLoaded"));
	wxStaticText* edOn = (wxStaticText *)FindWindowById(XRCID("m_staticText_REDOn"));
	wxButton* edShow = (wxButton *)FindWindowById(XRCID("m_button_showRED"));
	edLoad->SetBackgroundColour(*wxGREEN);
	edLoad->SetForegroundColour(*wxWHITE);
	edLoad->SetLabelText("  YES  ");
	edOn->SetBackgroundColour(*wxGREEN);
	edOn->SetForegroundColour(*wxWHITE);
	edOn->SetLabelText("  YES  ");
	edShow->SetLabelText("Show Internal Phantom");

	wxNotebook* tabs = (wxNotebook*)FindWindowById(XRCID("m_notebook_panels"));
	tabs->SetSelection(PANEL_EXTERNAL);
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

bool GLView::loadViewRayLVol(const wxString path)
{
	wxBusyInfo busyInfo(wxT("loading lvol file, wait please..."));
	FILE* fp = fopen(path.c_str(), "rb");
	if (NULL == fp) return false;
	uint32_t temp;
	uint32_t nx, ny, nz;
	double offset_x, offset_y, offset_z;
	double dx, dy, dz;
	fread(&temp, sizeof(uint32_t), 1, fp); //Marjor_Header
	fread(&temp, sizeof(uint32_t), 1, fp); //Minor_Header
	fread(&nx, sizeof(uint32_t), 1, fp);
	fread(&dx, sizeof(double), 1, fp);
	fread(&offset_x, sizeof(double), 1, fp);
	fread(&ny, sizeof(uint32_t), 1, fp);
	fread(&dy, sizeof(double), 1, fp);
	fread(&offset_y, sizeof(double), 1, fp);
	fread(&nz, sizeof(uint32_t), 1, fp);
	fread(&dz, sizeof(double), 1, fp);
	fread(&offset_z, sizeof(double), 1, fp);
	lvol.resize(nx, ny, nz);
	uint32_t* data = new uint32_t[nx*ny*nz];
	fread(data, sizeof(float), nx*ny*nz, fp);
	for (uint32_t ix = 0; ix < nx; ++ix)
		for (uint32_t iy = 0; iy < ny; ++iy)
			for (uint32_t iz = 0; iz < nz; ++iz)
			{
				lvol.a(ix, iy, iz) = data[ix + nx*iy + nx*ny*iz];
			}
	fclose(fp);

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
		if (NULL == fp) {
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
	if (!b_showDoseRef) setMatrix(&dose);
	else setMatrix(&doseRef);
	if (!b_showExtRED) setBGImage(&phantom);
	else setBGImage(&phantomRef);
	setVoxelSize(DX, DY, DZ);
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
	if (!dose.empty() && !phantom.empty())
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

void GLView::copyData(vector<SliceInfo> &slices3D)
{
	Hist = 0;

	NX = slices3D[0].NCol;
	NY = slices3D[0].NRow;
	NZ = slices3D.size();
	DX = slices3D[0].pixSizeX / 10.0;
	DY = slices3D[0].pixSizeY / 10.0;
	DZ = slices3D[0].thickness / 10.0;
	//get the cuboid origin position
	COPX = slices3D[0].x / 10.0;
	COPY = slices3D[0].y / 10.0;
	COPZ = slices3D[0].z / 10.0;

	phantom.resize(NX, NY, NZ);

	for (int iz = 0; iz < NZ; ++iz)
	{
		for (int iy = 0; iy < NY; ++iy)
		{
			for (int ix = 0; ix < NX; ++ix)
			{
				phantom.a(ix, iy, iz) = slices3D[iz].density[ix + iy*NX];
			}
		}
		slices3D[iz].release();
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
			if (ext == "dose") success = loadDose(path);
			else if (ext == "phtm") success = loadPhantom(path);
			else if (ext == "red") success = loadRED(path);
			else success = false;
		}
		else success = false;
		if (success) extraSettings();
		else wxMessageBox("failed to open this file due to wrong format or other reason");
	}
	else
	{
		/*
#ifdef WIN32
		// maximum number of lines the output console should have
		static const WORD MAX_CONSOLE_LINES = 500;
		{
			int hConHandle;
			long lStdHandle;
			CONSOLE_SCREEN_BUFFER_INFO coninfo;
			FILE *fp;
			// allocate a console for this app
			//AllocConsole();
			if (NULL == GetStdHandle(STD_OUTPUT_HANDLE))
			{
				AllocConsole();
			}
			else
			{
				BOOL suc = AttachConsole(ATTACH_PARENT_PROCESS);
				if (suc == 0) wxMessageBox("attach failed!");
			}
			//AttachConsole(ATTACH_PARENT_PROCESS);
			// set the screen buffer to be big enough to let us scroll text
			GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &coninfo);
			coninfo.dwSize.Y = MAX_CONSOLE_LINES;
			SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE), coninfo.dwSize);
			// redirect unbuffered STDOUT to the console
			lStdHandle = (long)GetStdHandle(STD_OUTPUT_HANDLE);
			hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);
			fp = _fdopen(hConHandle, "w");
			*stdout = *fp;
			setvbuf(stdout, NULL, _IONBF, 0);
			// redirect unbuffered STDIN to the console
			lStdHandle = (long)GetStdHandle(STD_INPUT_HANDLE);
			hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);
			fp = _fdopen(hConHandle, "r");
			*stdin = *fp;
			setvbuf(stdin, NULL, _IONBF, 0);
			// redirect unbuffered STDERR to the console
			lStdHandle = (long)GetStdHandle(STD_ERROR_HANDLE);
			hConHandle = _open_osfhandle(lStdHandle, _O_TEXT);
			fp = _fdopen(hConHandle, "w");
			*stderr = *fp;
			setvbuf(stderr, NULL, _IONBF, 0);
			// make cout, wcout, cin, wcin, wcerr, cerr, wclog and clog
			// point to console as well
			std::ios::sync_with_stdio();
		}
#endif
		*/

		static const wxCmdLineEntryDesc cmdLineDesc[] =
		{
			{ wxCMD_LINE_SWITCH, "h", "help", "displays help on the command line parameters",
			wxCMD_LINE_VAL_NONE, wxCMD_LINE_OPTION_HELP },
			{ wxCMD_LINE_SWITCH, "u", "uncertainty", "calculate the uncertainty of the dose" },
			{ wxCMD_LINE_PARAM, NULL, NULL, "input file", wxCMD_LINE_VAL_STRING, wxCMD_LINE_PARAM_MULTIPLE },
			{ wxCMD_LINE_NONE }
		};
		CmdMode = true;
		wxCmdLineParser cmd(argc, argv);
		cmd.SetDesc(cmdLineDesc);
		cmd.SetSwitchChars(wxT("-"));
		cmd.Parse();
		wxString file;
		if (cmd.Found("u")) //calculate the uncertainty
		{
			if (!loadDose(cmd.GetParam())) return;
			extraSettings();
			wxCommandEvent evtx;
			OnCalcUncertainty(evtx);
		}
		exit(0);
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
	VERTEX3D nearest(1e10, 1e10, 1e10); //a really far initial distance

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
	if (ZAXIS == wid && b_showBeams&& VRBeams.beams.size() > 0)
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

	if (ZAXIS == wid && b_showStructContours)
	{
		glTranslatef(-COPX, -COPY, -COPZ);
		for (size_t i = 0; i < strt_set.size(); ++i)
		{
			if (strt_set[i].show) //ROI
			{
				glColor3ub(strt_set[i].r, strt_set[i].g, strt_set[i].b);
				for (size_t j = 0; j < strt_set[i].contour.size(); ++j) //for each line
				{
					if (int((strt_set[i].contour[j][0].z - COPZ) / DZ) == cross[ZAXIS])
					{
						glBegin(GL_LINE_LOOP);
						for (size_t k = 0; k < strt_set[i].contour[j].size(); ++k)
							glVertex3f(strt_set[i].contour[j][k].x, strt_set[i].contour[j][k].y, strt_set[i].contour[j][k].z);
						glEnd();
					}
				}
			}
		}
		glPopMatrix();
	}
}
