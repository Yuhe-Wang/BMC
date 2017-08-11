#ifndef _DVFORM_H_
#define _DVFORM_H_

#include "GLView.h"
#include "SimuPanel.h"

class DVFrame : public wxFrame
{
public:
	DVFrame(wxWindow* parent);
	bool ProcessEvent(wxEvent& event);
	//virtual ~DVFrame();

	void OnExit(wxCommandEvent& event);
	void OnAbout(wxCommandEvent& event);
	void OnClose(wxCloseEvent& event);
	void OnSimuFinished(wxCommandEvent& event);
public:
	GLView *gl;
private:
	SimuPanel* pSimu;
	bool bMarkClose;
	DECLARE_EVENT_TABLE()
};

#endif