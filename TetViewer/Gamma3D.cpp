#include "PreCompile.h"
#include "Gamma3D.h"


BEGIN_EVENT_TABLE(Gamma3D, View3D)
EVT_RIGHT_DOWN(Gamma3D::rightClick)
EVT_MENU(XRCID("m_menuItem_zoom_in"), Gamma3D::OnZoomIn)
EVT_MENU(XRCID("m_menuItem_zoom_out"), Gamma3D::OnZoomOut)
EVT_TOOL(XRCID("m_tool_save"), Gamma3D::OnSave)
EVT_TOOL(XRCID("m_tool_export"), Gamma3D::OnExport)
EVT_TOOL(XRCID("m_tool_zoom_in"), Gamma3D::OnZoomIn)
EVT_TOOL(XRCID("m_tool_zoom_out"), Gamma3D::OnZoomOut)
EVT_TOOL(XRCID("m_tool_resetSize"), Gamma3D::OnResetSize)
EVT_TOOL(XRCID("m_tool_left"), Gamma3D::OnMoveLeft)
EVT_TOOL(XRCID("m_tool_right"), Gamma3D::OnMoveRight)
EVT_TOOL(XRCID("m_tool_up"), Gamma3D::OnMoveUp)
EVT_TOOL(XRCID("m_tool_down"), Gamma3D::OnMoveDown)
EVT_TOOL(XRCID("m_tool_rotate_left"), Gamma3D::OnRotateLeft)
EVT_TOOL(XRCID("m_tool_rotate_right"), Gamma3D::OnRotateRight)
EVT_TOOL(XRCID("m_tool_profile_horizontal"), Gamma3D::OnProfileH)
EVT_TOOL(XRCID("m_tool_profile_vertical"), Gamma3D::OnProfileV)
EVT_TOOL(XRCID("m_tool_profile_tilted"), Gamma3D::OnProfileT)
END_EVENT_TABLE()

Gamma3D::Gamma3D(FigureFrame* parent) :View3D(parent)
{
	threshold = 0;
	minDensity = 0.04f;
}

Gamma3D::~Gamma3D()
{

}

void Gamma3D::OnSave(wxCommandEvent& event)
{

}
void Gamma3D::OnExport(wxCommandEvent& event)
{

}
void Gamma3D::OnZoomIn(wxCommandEvent& event)
{
	zoomIn();
}
void Gamma3D::OnZoomOut(wxCommandEvent& event)
{
	zoomOut();
}
void Gamma3D::OnResetSize(wxCommandEvent& event)
{
	resetView();
}
void Gamma3D::OnMoveLeft(wxCommandEvent& event)
{
	moveLeft();
}
void Gamma3D::OnMoveRight(wxCommandEvent& event)
{
	moveRight();
}
void Gamma3D::OnMoveUp(wxCommandEvent& event)
{
	moveUp();
}
void Gamma3D::OnMoveDown(wxCommandEvent& event)
{
	moveDown();
}

void Gamma3D::OnRotateLeft(wxCommandEvent& event)
{
	rotateLeft();
}
void Gamma3D::OnRotateRight(wxCommandEvent& event)
{
	rotateRight();
}
void Gamma3D::OnProfileH(wxCommandEvent& event)
{
	profileH();
}
void Gamma3D::OnProfileV(wxCommandEvent& event)
{
	profileV();
}
void Gamma3D::OnProfileT(wxCommandEvent& event)
{
	profileT();
}

void Gamma3D::rightClick(wxMouseEvent& event)
{
	wxMenu mnu;
	if(b_showMat) mnu.Append(ID_SHOW_GAMMA_ISO_LINE, "Hide gamma iso lines");
	else mnu.Append(ID_SHOW_GAMMA_ISO_LINE, "Show gamma iso lines");
	if (b_showBGImg) mnu.Append(ID_SHOW_BACKGROUND_IMAGE, "Hide background image");
	else mnu.Append(ID_SHOW_BACKGROUND_IMAGE, "Show background image");
	if (failVoxel.empty()) mnu.Append(ID_SHOW_GAMMA_BITMAP, "Show the failing bitmap"); //change the background image to pass/failing bitmap
	else mnu.Append(ID_SHOW_GAMMA_BITMAP, "Hide the failing bitmap");
	mnu.Append(ID_TRIM_GAMMA_ABOVE_ONE, "Trim gamma above 1(make gamma <= 1)"); //the distribution of the passing voxels
	mnu.Append(ID_TRIM_GAMMA_BELOW_ONE, "Trim gamma below 1(make gamma >= 1)"); //the distribution of the failing voxels
	mnu.Append(ID_RESTORE_GAMMA, "Restore the initial gamma");
	mnu.Append(ID_SHOW_GAMMA_INFO, "View the gamma matrix info");
	mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&Gamma3D::OnPopupClick, NULL, this);
	PopupMenu(&mnu);
}

void Gamma3D::OnPopupClick(wxCommandEvent& event)
{
	int ID = event.GetId();
	int len = mat->getInnerLength();
	if (ID_SHOW_BACKGROUND_IMAGE == ID)
	{
		showBGImage();
	}
	else if (ID_SHOW_GAMMA_ISO_LINE == ID)
	{
		showMat();
	}
	else if (ID_TRIM_GAMMA_ABOVE_ONE == ID)
	{
		for (int i = 0; i < len; ++i)
		{
			if (mat->a(i) > 1) mat->a(i) = 1;
		}
		initShow();
	}
	else if (ID_TRIM_GAMMA_BELOW_ONE == ID)
	{
		for (int i = 0; i < len; ++i)
		{
			if (mat->a(i) <= 1) mat->a(i) = 1;
		}
		initShow();
	}
	else if (ID_SHOW_GAMMA_BITMAP == ID)
	{
		if (failVoxel.empty())
		{
			failVoxel.resize<SFloat>(*src_mat);
			for (int i = 0; i < len; ++i)
			{
				if (src_mat->a(i) <= 1) failVoxel.a(i) = RGBPixel(0,0,0);
				else failVoxel.a(i) = RGBPixel(255, 0, 0);
			}
			setBGPainter(&failVoxel);
			showBGImage(true);
		}
		else
		{
			failVoxel.release();
			setBGPainter(NULL);
		}
		show();
	}
	else if (ID_RESTORE_GAMMA == ID)
	{
		inner_mat = *src_mat; //restore the original gamma value
		setBGPainter(NULL); // no failing map
		showMat(true);
		showBGImage(true);
		initShow();
	}
	else if (ID_SHOW_GAMMA_INFO == ID)
	{
		int tNum, tPass;
		getPassingRate(tNum, tPass);
		
		wxString info;
		info.Printf("Image Information:\nNZ=%d , DZ=%.2f cm\nNX=%d , DX=%.2f cm\nNY=%d, DY=%.2f cm\nVolume =%.2f*%.2f*%.2f cm^3\nPassint rate = %f%% ",
			NZ, DZ, NX, DX, NY, DY, NX*DX, NY*DY, NZ*DZ, 100 * double(tPass) / tNum);
		wxMessageBox(info);
	}
}

void Gamma3D::getPassingRate(int& tNum, int& tPass)
{
	tNum = tPass = 0;//reset counter
	//int put threshold is a percentage of max dose.
	SFloat maxv, minv;
	doseRef->getMaxMin(maxv, minv);

	//update threshold, minDensity first
	SFloat threshold_val = maxv*threshold;

	//search where the couch is
	int couchY = searchViewRayCouch(*bgImage);
	int nx = bgImage->getWidth(1);
	int ny = bgImage->getWidth(2);
	int nz = bgImage->getWidth(3);
	if (!src_mat->equalDims(*bgImage))
	{
		wxMessageBox("unmatched dimension between gamma and gammaBGImage");
		tNum = -1;
		tPass = -1;
		return;
	}
	if (-1 == couchY) couchY = ny; //if there's no couch, we count all voxels
	for (int ix = 0; ix < nx; ++ix)
		for (int iy = 0; iy < ny; ++iy)
			for (int iz = 0; iz < nz; ++iz)
			{
				if (iy < couchY && doseRef->a(ix, iy, iz) > threshold_val && bgImage->a(ix, iy, iz) > minDensity)
				{
					++tNum;
					if (src_mat->a(ix, iy, iz) < 1) ++tPass;
				}
			}
}