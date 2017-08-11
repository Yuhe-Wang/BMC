#include "PreCompile.h"
#include <atomic>
#include "GLView.h"
#include "Graph.h"

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
EVT_BUTTON(XRCID("m_button_viewDoseDifference"), GLView::OnViewDoseDifference)
EVT_BUTTON(XRCID("m_button_swapDose"), GLView::OnSwapDose)
EVT_BUTTON(XRCID("m_button_averageDose"), GLView::OnAverageDose)
EVT_BUTTON(XRCID("m_button_averagePhantom"), GLView::OnAveragePhantom)
EVT_BUTTON(XRCID("m_button_setPrescriptionDose"), GLView::OnSetPrescriptionDose)
EVT_BUTTON(XRCID("m_button_samplePointDose"), GLView::OnSamplePointDose)
EVT_BUTTON(XRCID("m_button_getExpDTA"), GLView::OnGetExpDTA)
EVT_BUTTON(XRCID("m_button_sampleLineDose"), GLView::OnSampleLineDose)
EVT_BUTTON(XRCID("m_button_samplePlanDose"), GLView::OnSamplePlanDose)
EVT_BUTTON(XRCID("m_button_trimDoseByDensity"), GLView::OnTrimDoseByDensity)
EVT_BUTTON(XRCID("m_button_flipDose"), GLView::OnFlipDose)
EVT_BUTTON(XRCID("m_button_exportDiff"), GLView::OnExportDiff)
EVT_MENU(XRCID("m_menuItem_open"), GLView::OnOpen)
EVT_MENU(XRCID("m_menuItem_saveSlice"), GLView::OnSaveSlice)
EVT_MENU(XRCID("m_menuItem_zoom_in"), GLView::OnZoomIn)
EVT_MENU(XRCID("m_menuItem_zoom_out"), GLView::OnZoomOut)
EVT_MENU(XRCID("m_menuItem_rotate_left"), GLView::OnRotateLeft)
EVT_MENU(XRCID("m_menuItem_rotate_right"), GLView::OnRotateRight)
EVT_TOOL(XRCID("m_tool_open"), GLView::OnOpen)
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
EVT_TOOL(XRCID("m_tool_drawBox"), GLView::OnDrawBox)
EVT_TOOL(XRCID("m_tool_profile_horizontal"), GLView::OnProfileH)
EVT_TOOL(XRCID("m_tool_profile_vertical"), GLView::OnProfileV)
EVT_TOOL(XRCID("m_tool_profile_tilted"), GLView::OnProfileT)
EVT_CHECKBOX(XRCID("m_checkBox_showDose"), GLView::OnShowDose)
EVT_CHECKBOX(XRCID("m_checkBox_showPhantom"), GLView::OnShowPhantom)
EVT_CHECKBOX(XRCID("m_checkBox_showBox"), GLView::OnShowBox)
EVT_CHECKBOX(XRCID("m_checkBox_globalContrast"), GLView::OnGlobalContrast)
EVT_CHECKBOX(XRCID("m_checkBox_autoContrast"), GLView::OnAutoContrast)
EVT_CHECKBOX(XRCID("m_checkBox_showContours"), GLView::OnShowContours)
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
EVT_TEXT(XRCID("m_textCtrl_pointX3"), GLView::OnXPlanDistanceChanged)
EVT_TEXT(XRCID("m_textCtrl_pointY3"), GLView::OnYPlanDistanceChanged)
EVT_TEXT(XRCID("m_textCtrl_pointZ3"), GLView::OnZPlanDistanceChanged)

EVT_COMMAND(wxID_ANY, VIEW3D_NOTIFY, GLView::OnView3DNotify)
END_EVENT_TABLE()

GLView::GLView(wxWindow* parent, wxFrame* pframe) :View3D(parent)
{
	frame = pframe;

	//imcInfo = (wxStaticText *)FindWindowById(XRCID("m_staticText_imc"));

	b_showSampledLine = false; // turn off the internal function to draw the sampling line
	Hist = 0;
	xo = yo = zo = 0;

	//enable dragging file
	DragAcceptFiles(true);
	Connect(wxEVT_DROP_FILES, wxDropFilesEventHandler(GLView::OnDropFiles), NULL, this);
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
		wxMessageBox("I can only load one file per time!");
		return;
	}
	
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
			if (ext == "mats") success = loadDose(path);
		}
	}
	
	if(!success) wxMessageBox("Cannot open this file due to unsupported format or other reasons.");
	else extraSettings();
}

void GLView::OnOpen(wxCommandEvent& event)
{
	wxString caption = wxT("Choose a file");
	wxString wildcard = wxT("mat file(*.mats)|*.mats|all files(*)|*.*");
	wxString defaultDir = wxEmptyString;
	wxString defaultFilename = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_OPEN);
	if (dialog.ShowModal() == wxID_OK)
	{
		int ftype = dialog.GetFilterIndex();
		
		if (0 == ftype)
		{
			if (!loadDose(dialog.GetPath()))
			{
				wxMessageBox("The dose file doesn't have expected format");
				return;
			}
		}
	}
}

void GLView::OnSaveSlice(wxCommandEvent& event)
{
	if (phantom.empty() && dose.empty()) return;
	//save current slice to file that can be read by matlab
	wxString caption = wxT("Input a file name");
	wxString defaultFilename = wxT("export for matlab");
	wxString wildcard = wxT("bgImage slice(*.txt)|*.txt|mat slice(*.txt)|*.txt");

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
	dose *= adjustK / lastFactor;
	lastFactor = adjustK;
	wxMessageBox("overall dose adjustment done!");
}

void GLView::OnPDD(wxCommandEvent& event)
{
	if (NULL == mat || mat->empty()) return;

	wxString caption = wxT("Input a file name");
	wxString defaultFilename = wxT("PDD.txt");
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
			fprintf(fp, "%.2f\t%0.4f\t%0.4f\n", (iy + 0.5)*DY, pdd.a(iy), pdd.a(iy) / maxv);
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
	rotate(1);
	render();
}
void GLView::OnRotateRight(wxCommandEvent& event)
{
	rotate(-1);
	render();
}


void GLView::OnShowBox(wxCommandEvent& event)
{
	showBox();
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
	render();
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
	render();
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

void GLView::OnShowContours(wxCommandEvent& event) //switch for structure contour 
{
	showContours();
	render();
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
			ds.Printf(", dose = %5.3f Gray", dose.a(kx, ky, kz));
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
	else if (VIEW3D_LINE_SAMPLING_FINISHED == eventType)
	{
		ArrayMgr<double> xin, yin;
		auto fig = FIGURE::Figure(this, false);
		if (!doseRef.empty())
		{
			ArrayMgr<SFloat>* cmat = mat;
			setMatrix(&dose);
			getLineData(xin, yin);
			fig.Plot(MG(xin), MG(yin), "b");
			fig.hold();
			setMatrix(&doseRef);
			getLineData(xin, yin);
			fig.Plot(MG(xin), MG(yin), "g");
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
		fig.g->Label('x', "distance(cm)", 0);
		fig.g->Label('y', "dose(Gy)", 0);
		fig.g->Title("line sampling");
		fig.render();
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

	//save current slice to file that can be read by matlab
	wxString caption = wxT("Input a file name");
	wxString defaultFilename = wxT("plan info");
	wxString wildcard = wxT("mat slice(*.txt)|*.txt");

	wxString defaultDir = wxEmptyString;
	wxFileDialog dialog(this, caption, defaultDir, defaultFilename, wildcard, wxFD_SAVE);
	if (dialog.ShowModal() == wxID_OK)
	{
		wxString path = dialog.GetPath();
		
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
//<<------------------------------End event functions-------------------------------->>//


//<<------------------------------Processing functions------------------------------->>//
bool GLView::loadDose(const wxString& path)
{
	DX = DY = DZ = 0.1; //default 1mm
	//record current file path
	frame->SetTitle(path);
	currentFile = path;

	FILE* fp = fopen(path.c_str(), "rb");
	if (NULL == fp) return false;
	fread(&NX, sizeof(uint32_t), 1, fp);
	fread(&DX, sizeof(double), 1, fp);
	fread(&NY, sizeof(uint32_t), 1, fp);
	fread(&DY, sizeof(double), 1, fp);
	fread(&NZ, sizeof(uint32_t), 1, fp);
	fread(&DZ, sizeof(double), 1, fp);
	int hasMat = 0, hasImg = 0;
	fread(&hasMat, sizeof(uint32_t), 1, fp);
	int len = NX*NY*NZ;
	double* data = new double[len];
	if (hasMat)
	{
		dose.resize(NX, NY, NZ);
		fread(data, sizeof(double), len, fp);
		for (int ix = 0; ix < NX; ++ix)
			for (int iy = 0; iy < NY; ++iy)
				for (int iz = 0; iz < NZ; ++iz)
				{
					dose.a(ix, iy, iz) = SFloat(data[ix + NX*iy + NX*NY*iz]);
				}
	}
	else dose.resize(0);

	fread(&hasImg, sizeof(uint32_t), 1, fp);
	if (hasImg)
	{
		phantom.resize(NX,  NY, NZ);	
		fread(data, sizeof(double), len, fp);
		for (int ix = 0; ix < NX; ++ix)
			for (int iy = 0; iy < NY; ++iy)
				for (int iz = 0; iz < NZ; ++iz)
				{
					phantom.a(ix, iy, iz) = SFloat(data[ix + NX*iy + NX*NY*iz]);
				}
	}
	else phantom.resize(0);

	delete[] data;
	fclose(fp);
	
	//we need to set DX, DY, DZ
	setMatrix(&dose);
	setBGImage(&phantom);
	setVoxelSize(DX, DY, DZ);
	initShow();
	extraSettings();
	system("cmd /c del /q temp.mats");
	return true;
}

void GLView::extraSettings() //extra setting work after loading file
{
	wxString info;
	info.Printf("Image Information:\nNZ=%d , DZ=%.2f cm\nNX=%d , DX=%.2f cm\nNY=%d, DY=%.2f cm\nVolume =%.2f*%.2f*%.2f cm^3",
		NZ, DZ, NX, DX, NY, DY, NX*DX, NY*DY, NZ*DZ);
	if (Hist > 0)
	{
		wxString histStr;
		histStr.Printf("\nhist num=%.0f", Hist);
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

void GLView::copyData(vector<SliceInfo> &slices3D)
{
	Hist = 0;

	NX = slices3D[0].NCol;
	NY = slices3D[0].NRow;
	NZ = slices3D.size();
	DX = slices3D[0].pixSizeX/10.0;
	DY = slices3D[0].pixSizeY/10.0;
	DZ = slices3D[0].thickness/10.0;
	//get the cuboid origin position
	COPX = slices3D[0].x/10.0;
	COPY = slices3D[0].y/10.0;
	COPZ = slices3D[0].z/10.0;

	phantom.resize(NX, NY, NZ);
	
	for (int iz = 0; iz < NZ; ++iz)
	{
		for (int iy = 0; iy < NY; ++iy)
		{
			for (int ix = 0; ix < NX; ++ix)
			{
				phantom.a(ix, iy, iz) = slices3D[iz].density[ix+iy*NX];
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
			if (ext == "mats") success = loadDose(path);
			else success = false;
		}
		else success = false;
		if (success) extraSettings();
		else wxMessageBox("failed to open this file due to wrong format or other reason");
	}
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