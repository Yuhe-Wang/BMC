#ifndef _DVFORM_H_
#define _DVFORM_H_

#include "GLView.h"

class DVFrame : public wxFrame
{
public:
	DVFrame(wxWindow* parent);
	bool ProcessEvent(wxEvent& event);
	//virtual ~DVFrame();

	void OnExit(wxCommandEvent& event);
	void OnAbout(wxCommandEvent& event);
	
public:
	GLView *gl;
private:
	
	DECLARE_EVENT_TABLE()
};

#endif