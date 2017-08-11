#include "PreCompile.h"

#include "DataEditDlg.h"

BEGIN_EVENT_TABLE(DataEditDlg, wxDialog)
EVT_BUTTON(XRCID("m_button_load"), DataEditDlg::OnLoad)

EVT_BUTTON(XRCID("m_button_save"), DataEditDlg::OnSave)
EVT_BUTTON(XRCID("m_button_saveas"), DataEditDlg::OnSaveAs)
EVT_BUTTON(XRCID("m_button_chunkAll"), DataEditDlg::OnChunkAll)
EVT_BUTTON(XRCID("m_button_export"), DataEditDlg::OnExport)
EVT_BUTTON(XRCID("m_button_add"), DataEditDlg::OnAdd)
EVT_BUTTON(XRCID("m_button_del"), DataEditDlg::OnDel)
EVT_BUTTON(XRCID("m_button_import"), DataEditDlg::OnImport)
EVT_BUTTON(XRCID("m_button_help"), DataEditDlg::OnHelp)
EVT_BUTTON(XRCID("m_button_exit"), DataEditDlg::OnExit)
END_EVENT_TABLE()

DataEditDlg::DataEditDlg(const wxString& currentFile) :_fname(currentFile)
{
	wxXmlResource::Get()->LoadDialog(this, NULL, _T("DataEditDlg"));
	_grid = (wxGrid*)FindWindowById(XRCID("m_grid_data"));
	_grid->CreateGrid(100, 4);
	double DPIScale = wxBitmapScale::getDPIScale();
	_grid->SetColLabelSize(20 * DPIScale);
	_grid->SetRowLabelSize(40*DPIScale);
	_grid->SetColSize(0, 150*DPIScale);
	_grid->SetColSize(1, 80 * DPIScale);
	_grid->SetColSize(2, 60 * DPIScale);
	_grid->SetColSize(3, 170*DPIScale);
	_grid->SetColLabelValue(0, wxT("name"));
	_grid->SetColLabelValue(1, wxT("bytes"));
	_grid->SetColLabelValue(2, wxT("type"));
	_grid->SetColLabelValue(3, wxT("value"));
	//_grid->SetDefaultCellAlignment(wxALIGN_CENTER, wxALIGN_CENTER);

	_grid->Bind(wxEVT_GRID_CELL_CHANGING, &DataEditDlg::OnCellChange, this);
	_grid->Bind(wxEVT_GRID_CELL_RIGHT_CLICK, &DataEditDlg::OnAssignType, this);
	_key = (wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_key"));
	_value = (wxTextCtrl*)FindWindowById(XRCID("m_textCtrl_value"));
	_type = (wxComboBox*)FindWindowById(XRCID("m_choice_type"));
	_type->SetSelection(0);
	_type->Bind(wxEVT_CHOICE, &DataEditDlg::OnTypeChange, this);
	_len = (wxStaticText*)FindWindowById(XRCID("m_staticText_len"));

	wxButton* pb_load = (wxButton*)FindWindowById(XRCID("m_button_load"));
	pb_load->Bind(wxEVT_RIGHT_DOWN, &DataEditDlg::OnOpenCurrentFile, this);
	wxMouseEvent evtx;
	OnOpenCurrentFile(evtx);
}

void DataEditDlg::OnLoad(wxCommandEvent& event)
{
	wxFileDialog dialog(this, wxT("Open Binary Data File"), wxEmptyString, wxEmptyString, wxT("All files(*.*)|*.*"), wxFD_OPEN);
	if (dialog.ShowModal() == wxID_OK)
	{
		_fname = dialog.GetPath();
		if (_BF.read(_fname.c_str())) showGrid();
		else wxMessageBox("The file cannot be modified with this tool!");
	}
}
void DataEditDlg::OnOpenCurrentFile(wxMouseEvent& event)
{
	if (_BF.read(_fname.c_str())) showGrid();
	else wxMessageBox("Current file cannot be modified with this tool! You may need load another BF file");
}
void DataEditDlg::showGrid()
{
	_grid->ClearGrid();
	_BF.beginSearch();
	
	int row = 0;
	map<string, Chunk>& bm = _BF.getBM();
	for (map<string, Chunk>::iterator ic = bm.begin(); ic != bm.end(); ++ic)
	{
		const string& name = ic->first;
		Chunk& chk = ic->second;
		wxString slen;
		slen.Printf("%d", chk.len);
		_grid->SetCellValue(row, 0, name);
		_grid->SetCellValue(row, 1, slen);
		_grid->SetCellAlignment(row, 1, wxALIGN_CENTRE, wxALIGN_BOTTOM);
		if (chk.type == "")
		{
			//let's guess the type, may not right
			if (8 == chk.len)
			{
				double val = 0;
				memcpy(&val, chk.p, sizeof(double));
				wxString sv;
				sv.Printf("%g", val);
				chk.type = "double";
				chk.val = sv.c_str();
			}
			else if (4 == chk.len)
			{
				int val = 0;
				memcpy(&val, chk.p, sizeof(int));
				wxString sv;
				sv.Printf("%d", val);
				chk.type = "int";
				chk.val = sv.c_str();
			}
			else if (chk.len > 8)//guess it as a float array
			{
				float val = 0;
				memcpy(&val, chk.p, sizeof(float));
				wxString sv;
				sv.Printf("%g", val);
				if (chk.len > sizeof(float)) sv += ",  ...";
				chk.type = "float";
				chk.val = sv.c_str();
			}
		}
		_grid->SetCellValue(row, 2, chk.type);
		_grid->SetCellAlignment(row, 2, wxALIGN_CENTRE, wxALIGN_BOTTOM);
		_grid->SetCellValue(row, 3, chk.val);
		_grid->SetCellAlignment(row, 3, wxALIGN_CENTRE, wxALIGN_BOTTOM);
		_grid->SetReadOnly(row, 0, false);
		_grid->SetReadOnly(row, 1);
		_grid->SetReadOnly(row, 2);
		_grid->SetReadOnly(row, 3);
		if (8 == chk.len&&chk.type == "double" || 4 == chk.len && (chk.type == "int" || chk.type == "float")) _grid->SetReadOnly(row, 3, false);
		else _grid->SetReadOnly(row, 3);
		++row;
	}
	while (row < 100)
	{
		_grid->SetReadOnly(row, 0);
		_grid->SetReadOnly(row, 1);
		_grid->SetReadOnly(row, 2);
		_grid->SetReadOnly(row, 3);
		++row;
	}
}
void DataEditDlg::OnCellChange(wxGridEvent& event)
{
	if (_fname != wxEmptyString)
	{
		bool success = true;
		
		if (0 == event.GetCol()) //change key name
		{
			//get new value
			//string oldKey(event.GetString().c_str());
			string oldKey(_grid->GetCellValue(event.GetRow(), event.GetCol()).c_str());
			string newKey(event.GetString().c_str());
			
			if (_BF.changeKeyName(oldKey, newKey)) success = true;
			else success = false;
		}
		else if (3 == event.GetCol())//change the value
		{	
			string Key(_grid->GetCellValue(event.GetRow(), 0).c_str());
			wxString newVal(event.GetString().c_str());
			Chunk& chk =_BF.getChunk(Key);
			if (chk.type == "double")
			{
				double dval = 0;
				if (newVal.ToDouble(&dval))
				{
					memcpy(chk.p, &dval, sizeof(double));
					newVal.Printf("%f", dval);
					chk.val = newVal.c_str();
				}
				else success = false;
			}
			else if( chk.type == "float")
			{
				double dval = 0;
				if (newVal.ToDouble(&dval))
				{
					float fval = float(dval);
					memcpy(chk.p, &fval, sizeof(float));
					newVal.Printf("%f", fval);
					chk.val = newVal.c_str();
				}
				else success = false;
			}
			else if (chk.type == "int")
			{
				long lval = 0;
				if (newVal.ToLong(&lval))
				{
					int ival = int(lval);
					memcpy(chk.p, &ival, sizeof(int));
					newVal.Printf("%d", ival);
					chk.val = newVal.c_str();
				}
				else success = false;
			}
			else success = false;
		}
		if (success) showGrid();
		else _grid->SetCellValue(event.GetRow(), event.GetCol(), event.GetString().c_str());
	}
}
void DataEditDlg::OnSave(wxCommandEvent& event)
{
	if (_fname != wxEmptyString)
	{
		_BF.write(_fname.c_str());
		wxMessageBox("Saved the file successfully");
	}
	else
	{
		wxMessageBox("Must open a file first!");
	}
	
}

void DataEditDlg::OnSaveAs(wxCommandEvent& event)
{
	if (_BF.size())
	{
		wxFileDialog dialog(this, wxT("Save Binary Data File"), wxEmptyString, _fname, wxT("All files(*.*)|*.*"), wxFD_SAVE);
		if (dialog.ShowModal() == wxID_OK)
		{
			wxString path = dialog.GetPath();
			_BF.write(path.c_str());
			wxMessageBox("Saved as another file successfully");
		}
	}
	else
	{
		wxMessageBox("Must have data first");
	}

}

void DataEditDlg::OnChunkAll(wxCommandEvent& event)
{
	if (!_BF.size())
	{
		wxMessageBox("Must have data first");
		return;
	}
	wxString key = _key->GetLineText(0).c_str();
	if (key!=wxEmptyString)
	{
		char* cname[64];
		strcpy((char*)cname, key.c_str());
		wxString defaultFile = key + ".chunk";
		wxFileDialog dialog(this, wxT("Save the whole data as a chunk"), wxEmptyString,defaultFile, wxT("chunk files(*.chunk)|*.chunk"), wxFD_SAVE);
		if (dialog.ShowModal() == wxID_OK)
		{
			wxString path = dialog.GetPath();
			FILE *fp = fopen(path.c_str(), "wb");
			if (NULL == fp)
			{
				wxMessageBox("cannot open output file");
				return;
			}
			fwrite(cname, 64, 1, fp);
			int len = _BF.getTotBytes();
			fwrite(&len, sizeof(int), 1, fp);
			_BF.write(fp);
			fclose(fp);
			wxMessageBox("Saved as a whole chunk file successfully");
		}
	}
	else
	{
		wxMessageBox("Should give a name for this whole chunk");
	}

}

void DataEditDlg::OnExport(wxCommandEvent& event)
{
	if (_grid->GetGridCursorRow()<_BF.size())
	{
		wxString defaultName = _grid->GetCellValue(_grid->GetGridCursorRow(), 0);
		string name(defaultName.c_str());
		//defaultName += wxT(".chunk");
		wxFileDialog dialog(this, wxT("Save Chunk Data File"), wxEmptyString,defaultName, wxT("Chunk files(*.chunk)|*.chunk"), wxFD_SAVE);
		if (dialog.ShowModal() == wxID_OK)
		{
			wxString path = dialog.GetPath();
			_BF.writeChunk(path.c_str(), name);
		}
	}
	else
	{
		wxMessageBox("Invalid data to export");
	}

}

void DataEditDlg::OnTypeChange(wxCommandEvent& event)
{
	int type = _type->GetSelection();
	int len = 0;
	if (type == INT32_T) len = sizeof(int);
	else if (type == DOUBLE_T) len = sizeof(double);
	else if (type == FLOAT_T) len = sizeof(float);
	wxString slen;
	if (len > 0)slen.Printf("len=%d", len);
	else slen.Printf("len = ?");
	_len->SetLabel(slen);
}

void  DataEditDlg::OnPopupClick(wxCommandEvent &evt)
{
	wxPoint *data = (wxPoint*)(static_cast<wxMenu *>(evt.GetEventObject())->GetClientData());
	int row = data->x;
	string key(_grid->GetCellValue(row, 0).c_str());
	Chunk& chk = _BF.getChunk(key);
	switch (evt.GetId()) 
	{
	case ID_ASSIGN_INT:
		if (chk.len >= sizeof(int))
		{
			int val = 0;
			memcpy(&val, chk.p, sizeof(int));
			wxString sv;
			sv.Printf("%d", val);
			if (chk.len > sizeof(int)) sv += ",...";
			chk.type = "int";
			chk.val = sv.c_str();
		}
		break;
	case ID_ASSIGN_DOUBLE:
		if (chk.len >= sizeof(double))
		{
			double val = 0;
			memcpy(&val, chk.p, sizeof(double));
			wxString sv;
			sv.Printf("%f", val);
			if (chk.len > sizeof(double)) sv += ",...";
			chk.type = "double";
			chk.val = sv.c_str();
		}
		break;
	case ID_ASSIGN_FLOAT:
		if (chk.len >= sizeof(float))
		{
			float val = 0;
			memcpy(&val, chk.p, sizeof(float));
			wxString sv;
			sv.Printf("%f", val);
			if (chk.len > sizeof(float)) sv += ",...";
			chk.type = "float";
			chk.val = sv.c_str();
		}
		break;
	case ID_ASSIGN_CHAR:
		char cv[33];
		int len = chk.len;
		if (len > 32) 	len = 32;
		memcpy(cv, chk.p, len);
		cv[32] = 0; //make the string null terminated
		cv[len] = 0;
		wxString sv;
		sv.Printf("%s", cv);
		if (chk.len > 32) sv += "...";
		chk.type = "char*";
		chk.val = sv.c_str();
		break;
	}
	delete data;
	showGrid();
}

void DataEditDlg::OnAssignType(wxGridEvent& event)
{
	int row = event.GetRow();
	int col = event.GetCol();
	_grid->GoToCell(row, col);
	if (row < _BF.size() && col == 2)
	{
		//void *data = reinterpret_cast<void *>(evt.GetItem().GetData());
		wxPoint *data = new wxPoint;
		data->x = row;
		data->y = col;
		wxMenu mnu;
		mnu.SetClientData((void*)data);
		mnu.Append(ID_ASSIGN_INT, "Assign int type");
		mnu.Append(ID_ASSIGN_DOUBLE, "Assign double type");
		mnu.Append(ID_ASSIGN_FLOAT, "Assign float type");
		mnu.Append(ID_ASSIGN_CHAR, "Assign char* type");
		mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, (wxObjectEventFunction)&DataEditDlg::OnPopupClick, NULL, this);
		PopupMenu(&mnu);
	}
}
void DataEditDlg::OnAdd(wxCommandEvent& event)
{
	string key(_key->GetLineText(0).c_str());
	wxString value = _value->GetLineText(0);
	if (key == wxEmptyString)
	{
		wxMessageBox("key name should not be empty");
		return;
	}
	else if (_BF.exist(key))
	{
		wxMessageBox("the key name has already existed");
		return;
	}
	else
	{
		//get the chunk value
		int type = _type->GetSelection();
		Chunk chk;
		if (INT32_T == type)
		{
			chk.len = sizeof(int);
			chk.p = (char*)new int;
			long lv;
			if (!value.ToLong(&lv))
			{ 
				wxMessageBox("the value isn't int type");
				return;
			}
			*((long*)chk.p) = lv;
			chk.type = "int";
			chk.val = value.c_str();
		}
		else if (DOUBLE_T == type)
		{
			chk.len = sizeof(double);
			chk.p = (char*)new double;
			double dv;
			if (!value.ToDouble(&dv))
			{
				wxMessageBox("the value isn't double float type");
				return;
			}
			*(double*)chk.p = dv;
			chk.type = "double";
			chk.val = value.c_str();
		}
		else if (FLOAT_T == type)
		{
			chk.len = sizeof(float);
			chk.p = (char*)new float;
			double dv;
			if (!value.ToDouble(&dv))
			{
				wxMessageBox("the value isn't float type");
				return;
			}
			*(double*)chk.p = float(dv);
			chk.type = "float";
			chk.val = value.c_str();
		}
		else if (STRING_T == type)
		{
			if (value != wxEmptyString)
			{
				chk.len = value.size()+1;
				chk.type = "string";
				chk.val = value;
				chk.p = new char[value.size() + 1];
				memcpy(chk.p, value.c_str(), value.size() + 1);
			}
			else
			{
				wxMessageBox("The value cannot be empty");
				return;
			}
		}
		else if (VOID_T == type)
		{
			FILE* fp = fopen(_chkname.c_str(), "rb");
			if (NULL == fp)
			{
				wxMessageBox("cannot open the chunk file");
				return;
			}
			int len;
			char cname[64];
			fread(cname, 1, 64, fp);
			fread(&len, sizeof(int), 1, fp);
			chk.len = len;
			chk.p = new char[len];
			fread(chk.p, len, 1, fp);
			fclose(fp);

			_chkname = wxEmptyString;
		}

		if (!_BF.addChunk(key, chk))
		{
			wxMessageBox("failed to add this chunk");
			return;
		}
		else showGrid();

		_key->SetLabel("");
		_value->SetLabel(wxEmptyString);
		_type->SetSelection(INT32_T);
		_len->SetLabel("len = 4");
	}
}
void DataEditDlg::OnDel(wxCommandEvent& event)
{
	if (_grid->GetGridCursorRow() < _BF.size())
	{
		string name(_grid->GetCellValue(_grid->GetGridCursorRow(), 0).c_str());
		if (_BF.erase(name))
		{
			showGrid();
		}
		else wxMessageBox("Cannot delete the selected data");
	}
	else
	{
		wxMessageBox("Invalid data to delete");
	}
}
void DataEditDlg::OnImport(wxCommandEvent& event)
{
	wxFileDialog dialog(this, wxT("Open Chunk File"), wxEmptyString, wxEmptyString, wxT("Chunk files(*.chunk)|*.chunk"), wxFD_OPEN);
	if (dialog.ShowModal() == wxID_OK)
	{
		_chkname = dialog.GetPath();
		FILE* fp = fopen(_chkname.c_str(),"rb");
		if (NULL == fp)
		{
			wxMessageBox("cannot open the chunk file");
			return;
		}
		int len;
		char cname[64];
		fread(cname, sizeof(char), 64, fp);
		wxString name((char*)cname);
		fread(&len, sizeof(int), 1, fp);
		fclose(fp);
		_key->SetLabel(name);
		_value->SetLabel(wxEmptyString);
		_type->SetSelection(VOID_T);
		wxString slen;
		slen.Printf("len = %d",len);
		_len->SetLabel(slen);
	}
}
void DataEditDlg::OnHelp(wxCommandEvent &evt)
{
	wxMessageBox("Will write in detail later!");
}
