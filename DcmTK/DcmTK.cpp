
#include "DcmTK.h"
#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dctk.h"

#define DFF(p) ((DcmFileFormat*)p)

DicomFile::DicomFile()
{
	dcmFF = new DcmFileFormat;
}

DicomFile::DicomFile(const string& fname)
{
	dcmFF = new DcmFileFormat;
	open(fname);
}

DicomFile::~DicomFile()
{
	delete DFF(dcmFF);
}

bool DicomFile::open(const string& fname)
{
	if(DFF(dcmFF)->loadFile(fname.c_str()).bad()) return false;
	DcmDataset *dataset = DFF(dcmFF)->getDataset();
	// decompress data set if compressed
	dataset->chooseRepresentation(EXS_LittleEndianExplicit, NULL);
	return true;
}
bool DicomFile::save(const string& fname)
{
	if (DFF(dcmFF)->saveFile(fname.c_str(), EXS_LittleEndianExplicit).bad()) return false;
	if (DFF(dcmFF)->clear().bad()) return false;
	return true;
}

// get one item. If you know their data layout structure
MatlabType DicomFile::getItem(const string& TagName)
{
	DcmTag tag;
	if ((DcmTag::findTagFromName(TagName.c_str(), tag)).bad()) return MatlabType();
	DcmElement* elem;
	if (DFF(dcmFF)->getDataset()->findAndGetElement(tag, elem).bad()) return MatlabType();
	return getItem(elem);
}

MatlabType DicomFile::getAllItems()
{
	MatlabType val;
	DcmDataset* dataset = DFF(dcmFF)->getDataset();
	int num = dataset->card();
	for (int i = 0; i < num; ++i)
	{
		DcmElement* subElem = dataset->getElement(i);
		DcmTag tag  = subElem->getTag();
		val[tag.getTagName()] = getItem(subElem);
	}
	return val;
}

MatlabType DicomFile::getItem(void* e)
{
	DcmElement* elem = (DcmElement*)e;
	MatlabType val;
	if (elem->isaString())
	{
		DcmByteString* str = (DcmByteString*)elem;
		int vm = str->getVM();
		OFString ofs;
		// decide the type
		if (elem->getVR() == EVR_DS) // decimal string
		{
			val.setType(MATLAB_DOUBLE);
			T_Basic& b = val.castBasic();
			b.resize(vm);
			for (int i = 0; i < vm; ++i)
			{
				str->getOFString(ofs, i);
				b.a(i) = std::stod(ofs.c_str()); //convert to double
			}
		}
		else if (elem->getVR() == EVR_IS) // integer string
		{
			val.setType(MATLAB_INT32);
			T_Basic& b = val.castBasic();
			b.resize(vm);
			for (int i = 0; i < vm; ++i)
			{
				str->getOFString(ofs, i);
				b.a(i) = double(std::stoi(ofs.c_str())); //convert to int
			}
		}
		else // string type
		{
			if (vm == 0) val.setType(MATLAB_STRING);
			else if (vm == 1)
			{
				str->getOFString(ofs, 0);
				val = ofs.c_str();
			}
			else if (vm > 1)
			{
				val.setType(MATLAB_CELL);
				T_Cell& c = val.castCell();
				c.resize(vm);
				for (int i = 0; i < vm; ++i)
				{
					str->getOFString(ofs, i);
					c.a(i) = ofs.c_str();
				}
			}
		}
		
	}
	else if (elem->getVR() == EVR_FD) // double 
	{
		int vm = elem->getVM();
		val.setType(MATLAB_DOUBLE);
		T_Basic& b = val.castBasic();
		b.resize(vm);
		Float64 fd;
		for (int i = 0; i < vm; ++i)
		{
			elem->getFloat64(fd, i);
			b.a(i) = fd;
		}
	}
	else if (elem->getVR() == EVR_FL) // float
	{
		int vm = elem->getVM();
		val.setType(MATLAB_FLOAT);
		T_Basic& b = val.castBasic();
		b.resize(vm);
		Float32 f;
		for (int i = 0; i < vm; ++i)
		{
			elem->getFloat32(f, i);
			b.a(i) = f;
		}
	}
	else if (elem->getVR() == EVR_SL)
	{
		int vm = elem->getVM();
		val.setType(MATLAB_INT32);
		T_Basic& b = val.castBasic();
		b.resize(vm);
		Sint32 si;
		for (int i = 0; i < vm; ++i)
		{
			elem->getSint32(si, i);
			b.a(i) = si;
		}
	}
	else if (elem->getVR() == EVR_SS)
	{
		int vm = elem->getVM();
		val.setType(MATLAB_INT16);
		T_Basic& b = val.castBasic();
		b.resize(vm);
		Sint16 si;
		for (int i = 0; i < vm; ++i)
		{
			elem->getSint16(si, i);
			b.a(i) = si;
		}
	}
	else if (elem->getVR() == EVR_UL)
	{
		int vm = elem->getVM();
		val.setType(MATLAB_UINT32);
		T_Basic& b = val.castBasic();
		b.resize(vm);
		Uint32 ui;
		for (int i = 0; i < vm; ++i)
		{
			elem->getUint32(ui, i);
			b.a(i) = ui;
		}
	}
	else if (elem->getVR() == EVR_US || elem->getVR() == EVR_AT)
	{
		int vm = elem->getVM();
		val.setType(MATLAB_UINT16);
		T_Basic& b = val.castBasic();
		b.resize(vm);
		Uint16 ui;
		for (int i = 0; i < vm; ++i)
		{
			elem->getUint16(ui, i);
			b.a(i) = ui;
		}
	}
	else if (elem->getVR() == EVR_SQ) // sequence data
	{
		elem->getTag().getVRName();
		val.setType(MATLAB_CELL);
		DcmSequenceOfItems* sqItem = (DcmSequenceOfItems*)elem;
		int num = sqItem->card();
		T_Cell& c = val.castCell();
		c.resize(num);
		//get each item
		for (int i = 0; i < num; ++i) 
		{
			DcmItem* subItem = sqItem->getItem(i);
			DcmTag tag = subItem->getTag();
			string itemName = tag.getTagName();
			MatlabType mitem;
			// search all elements for each item
			int nc = subItem->card();
			for (int j = 0; j < nc; ++j)
			{
				DcmElement* subElem = subItem->getElement(j);
				DcmTag subtag = subElem->getTag();
				mitem[subtag.getTagName()] = getItem(subElem);
			}
			c.a(i) = mitem;
		}
	}
	return val;
}

unsigned char* DicomFile::getPixelData()
{
	unsigned char* pImage = NULL;
	DcmElement* element = NULL;
	if (DFF(dcmFF)->getDataset()->findAndGetElement(DCM_PixelData, element).bad()) return NULL;
	element->getUint8Array(pImage);
	return pImage;
}

bool DicomFile::setItem(const string& TagName, const MatlabType& val)
{
	DcmTag tag;
	if ((DcmTag::findTagFromName(TagName.c_str(), tag)).bad())
	{
		printf("Cannot find the tag name\n");
		return false;
	}
	DcmDataset *dataset = DFF(dcmFF)->getDataset();
	DcmVR vr = tag.getVR();
	MDataType type = val.getType();
	DcmEVR evr = vr.getEVR();
	if (vr.isaString())
	{
		if (type == MATLAB_STRING) dataset->putAndInsertOFStringArray(tag, val.cString().c_str());
		else if (type == MATLAB_CELL)
		{
			const T_Cell& c = val.cCell();
			int ns = c.getInnerLength();
			string s = c[0].cString();
			for (int i = 1; i < ns; ++i)
			{
				s += "\\";
				s += c[i].cString();
			}
			dataset->putAndInsertOFStringArray(tag, s.c_str());
		}
		else if (type < MATLAB_BASIC)
		{
			printf("Expecting string, but now get Matlab_Basic type. Will perform casting!\n");
			int ns = val.cBasic().getInnerLength();
			string s;
			char wd[60];
			if (type <= MATLAB_FLOAT) sprintf(wd, "%f", val.cBasic()[0]);
			else if (type < MATLAB_BASIC) // integer
			{
				int d = int(val.cBasic()[0]);
				sprintf(wd, "%d", d);
			}
			s += wd;
			for (int i = 1; i < ns; ++i)
			{
				if (type <= MATLAB_FLOAT) sprintf(wd, "%f", val.cBasic()[i]);
				else if (type < MATLAB_BASIC) // integer
				{
					int d = int(val.cBasic()[i]);
					sprintf(wd, "%d", d);
				}
				s += "\\";
				s += wd;
			}
			dataset->putAndInsertOFStringArray(tag, s.c_str());
		}
		else
		{
			printf("Warning: the value provided doesn't match the type expected!\n");
			return false;
		}
	}
	else if (evr == EVR_SQ) // sequence
	{
		printf("Warning: currently doesn't support SQ tag writing back!\n");
		return false;
	}
	else // basic type and sequence
	{
		if (type > MATLAB_BASIC)
		{
			printf("Warning: the value provided doesn't match the type expected!\n");
			return false;
		}
		const Float64* p = val.cBasic().getConstP();
		const unsigned long len = val.cBasic().getInnerLength();

		if (evr == EVR_FD) // double
		{
			if (type != MATLAB_DOUBLE) printf("Warning: casting other basic type to double");
			for (unsigned long i = 0; i < len; ++i)
			{
				dataset->putAndInsertFloat64(tag, Float64(p[i]), i);
			}
		}
		else if (evr == EVR_FL) // float
		{
			if (type != MATLAB_FLOAT) printf("Warning: casting other basic type to float");			
			for (unsigned long i = 0; i < len; ++i)
			{
				dataset->putAndInsertFloat32(tag, Float32(p[i]), i);
			}
		}
		else if (evr == EVR_SL)
		{
			if (type != MATLAB_INT32) printf("Warning: casting other basic type to signed long (int32)");
			for (unsigned long i = 0; i < len; ++i)
			{
				dataset->putAndInsertSint32(tag, Sint32(p[i]), i);
			}
		}
		else if (evr == EVR_SS)
		{
			if (type != MATLAB_INT16) printf("Warning: casting other basic type to signed long (int16)");
			for (unsigned long i = 0; i < len; ++i)
			{
				dataset->putAndInsertSint16(tag, Sint16(p[i]), i);
			}
		}
		else if (evr == EVR_UL)
		{
			if (type != MATLAB_UINT32) printf("Warning: casting other basic type to unsigned long (uint32)");
			for (unsigned long i = 0; i < len; ++i)
			{
				dataset->putAndInsertUint32(tag, Uint32(p[i]), i);
			}
		}
		else if (evr == EVR_US || evr == EVR_AT)
		{
			if (type != MATLAB_UINT16) printf("Warning: casting other basic type to unsigned short (uint16)");
			for (unsigned long i = 0; i < len; ++i)
			{
				dataset->putAndInsertUint16(tag, Uint16(p[i]), i);
			}
		}
	}
	
	return true;
}

bool DicomFile::setPixelData(void* pin, unsigned long bytes)
{
	Uint8* p = (Uint8*)pin;
	if(DFF(dcmFF)->getDataset()->putAndInsertUint8Array(DCM_PixelData, p, bytes).good()) return true;
	else return false;
}