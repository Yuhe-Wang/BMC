#include "PreCompile.h"

#include "CTOptionDlg.h"


BEGIN_EVENT_TABLE(CTOptionDlg, wxDialog)
EVT_BUTTON(XRCID("m_button_shrink"), CTOptionDlg::OnShrink)
EVT_BUTTON(XRCID("m_button_quit"), CTOptionDlg::OnCancel)
END_EVENT_TABLE()

CTOptionDlg::CTOptionDlg(int rows, int cols, int NZ)
{
	wxXmlResource::Get()->LoadDialog(this, NULL, _T("CTInputDlg"));
	
	wxString res;
	res.Printf(wxT("the original resolution is %d * %d, which may be too large for simulation."), rows, cols);
	((wxControl*)FindWindowById(XRCID("m_staticText_resolution")))->SetLabelText(res);
	wxChoice *resChoice = (wxChoice *)FindWindowById(XRCID("m_choice_resolution"));	
	for (int i = 1; i <= 8; ++i)
	{
		res.Printf("%d", i);
		resChoice->Append(res);
	}
	resChoice->SetSelection(0);

	wxString nzs;
	nzs.Printf(wxT("the slice number is %d, which may be too large for simulation."), NZ);
	((wxControl*)FindWindowById(XRCID("m_staticText_NZ")))->SetLabelText(nzs);
	wxChoice *ncChoice = (wxChoice *)FindWindowById(XRCID("m_choice_nc"));		
	for (int i = 0; i < 10; ++i)
	{
		nzs.Printf(wxT("%d"), i+1);
		ncChoice->Append(nzs);
	}
	ncChoice->SetSelection(0);
}

void CTOptionDlg::OnShrink(wxCommandEvent& event)
{

	wxChoice *ncChoice = (wxChoice *)FindWindowById(XRCID("m_choice_resolution"));
	nShrink = ncChoice->GetSelection();
	nShrink += 1;

	ncChoice = (wxChoice *)FindWindowById(XRCID("m_choice_nc"));
	NC = ncChoice->GetSelection();
	NC += 1;

	EndModal(0);
}

void CTOptionDlg::OnCancel(wxCommandEvent& event)
{
	nShrink = 1;
	NC = 1;
	EndModal(0);
}

int CTOptionDlg::ShowModal()
{
	return wxDialog::ShowModal();
}