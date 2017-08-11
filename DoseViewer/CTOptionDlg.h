#ifndef _CTOPTIONDLG_H_
#define _CTOPTIONDLG_H_

class CTOptionDlg : public wxDialog
{
public:
	CTOptionDlg(int rows, int cols, int NZ);
	//virtual ~DVFrame();
	void getParameters(int& nShrink, int& NC){ nShrink = this->nShrink; NC = this->NC; }

	void OnShrink(wxCommandEvent& event);
	void OnCancel(wxCommandEvent& event);
	virtual int ShowModal(); //override to give detailed control
public:

private:
	int nShrink, NC;
	DECLARE_EVENT_TABLE()
};

#endif