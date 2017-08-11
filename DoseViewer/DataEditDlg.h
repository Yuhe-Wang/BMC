#ifndef _DATAEDITDLG_H_
#define _DATAEDITDLG_H_

#define ID_ASSIGN_INT		2001
#define ID_ASSIGN_DOUBLE	2002
#define ID_ASSIGN_FLOAT	2003
#define ID_ASSIGN_CHAR	2004

class DataEditDlg : public wxDialog
{
public:
	DataEditDlg(const wxString& currentFile="");
	void setCurrentFile(const wxString& currentFile){ _fname = currentFile; }

	void OnLoad(wxCommandEvent& event);
	void OnOpenCurrentFile(wxMouseEvent& event);
	void OnSave(wxCommandEvent& event);
	void OnSaveAs(wxCommandEvent& event);
	void OnChunkAll(wxCommandEvent& event);
	void OnExport(wxCommandEvent& event);
	void OnAdd(wxCommandEvent& event);
	void OnDel(wxCommandEvent& event);
	void OnImport(wxCommandEvent& event);
	void OnCellChange(wxGridEvent& event);
	void OnAssignType(wxGridEvent& event);
	void OnTypeChange(wxCommandEvent& event);
	void OnPopupClick(wxCommandEvent &evt);
	void OnHelp(wxCommandEvent &evt);
	void OnExit(wxCommandEvent &evt){ Close(); }
	

private:
	void showGrid();
	

	BinaryFile _BF;

	wxGrid* _grid;
	wxTextCtrl* _key;
	wxTextCtrl* _value;
	wxComboBox* _type;
	wxStaticText* _len;
	wxString _fname;
	wxString _chkname;
	DECLARE_EVENT_TABLE()
};

#endif