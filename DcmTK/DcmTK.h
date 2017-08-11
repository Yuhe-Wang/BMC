#pragma once
#ifndef _DCMTK_WRAPPER_H_
#define _DCMTK_WRAPPER_H_
// This is wrapper of DCMTK to behave as simple as MATLAB

#ifndef WIN32 //linux platform
#define __declspec(a)
#endif

#include "../Tools/Tools.h"
#undef at

#define T_Basic ArrayMgr<double>
#define T_String string
#define T_Struct map<string, MatlabType>
#define T_Cell ArrayMgr<MatlabType>

#define M_Basic(d) (*((T_Basic*)(d)))
#define M_String(d) (*((T_String*)(d)))
#define M_Struct(d) (*((T_Struct*)(d)))
#define M_Cell(d) (*((T_Cell*)(d)))

enum MDataType
{
	// the basic c++ type
	MATLAB_DOUBLE = 0,
	MATLAB_FLOAT,
	MATLAB_INT8, MATLAB_INT16, MATLAB_INT32, MATLAB_INT64,
	MATLAB_UINT8, MATLAB_UINT16, MATLAB_UINT32, MATLAB_UINT64,
	MATLAB_BASIC, // used as a separator
	// compound type
	MATLAB_STRING,
	MATLAB_STRUCT,
	MATLAB_CELL
};

class MatlabType // It holds double float matrix, or string, or struct, or cell like MATLAB
{
	
public:
	
	MatlabType()
	{
		d = new ArrayMgr<double>;
		type = MATLAB_DOUBLE;
	}
	MatlabType(const double val)
	{
		type = MATLAB_DOUBLE;
		d = new ArrayMgr<double>;
		M_Basic(d).resize(1);
		M_Basic(d).a(0) = val;
	}
	MatlabType(const double* val, int len)
	{
		type = MATLAB_DOUBLE;
		d = new ArrayMgr<double>;
		if (len <= 0) return;
		M_Basic(d).resize(len);
		for (int i = 0; i < len; ++i) M_Basic(d).a(i) = val[i];
	}
	MatlabType(const string val)
	{
		d = new string;
		type = MATLAB_STRING;
		M_String(d) = val;
	}
	MatlabType(const MatlabType& r)
	{
		d = new ArrayMgr<double>;
		type = MATLAB_DOUBLE;
		setType(r.getType());
		if (type < MATLAB_BASIC) // basic type
		{
			M_Basic(d) = M_Basic(r.getP());
		}
		else if (MATLAB_STRING == type)
		{
			M_String(d) = M_String(r.getP());
		}
		else if (MATLAB_STRUCT == type)
		{
			M_Struct(d) = M_Struct(r.getP());
		}
		else if (MATLAB_CELL == type)
		{
			M_Cell(d) = M_Cell(r.getP());
		}
	}
	~MatlabType()
	{
		if (d != NULL)
		{
			if (type < MATLAB_BASIC) delete ((T_Basic*)d);
			else if (type == MATLAB_STRING) delete ((T_String*)d);
			else if (type == MATLAB_STRUCT) delete ((T_Struct*)d);
			else if (type == MATLAB_CELL) delete ((T_Cell*)d);
		}
	}

	MDataType getType() const { return type; }
	bool setType(MDataType newType)
	{
		if (newType < 0 || newType == MATLAB_BASIC || newType > MATLAB_CELL) return false;
		if (type < MATLAB_BASIC && newType < MATLAB_BASIC) // internal type doesn't change, still basic
		{
			type = newType;
			return true;
		}
		if (type != newType)
		{
			delete d;
			type = newType;
			if (type < MATLAB_BASIC) // basic type
			{
				d = new ArrayMgr<double>; // arrange the order by frequency
			}
			else if (MATLAB_STRING == type)
			{
				d = new string;
			}
			else if (MATLAB_STRUCT == type)
			{
				d = new map<string, MatlabType>;
			}
			else if (MATLAB_CELL == type)
			{
				d = new ArrayMgr<MatlabType>;
			}
		}
		return true;
	}
	void* getP() const { return d; }
	bool isBasic(){ if (type < MATLAB_BASIC) return true; else return false; }
	bool empty()
	{
		if (type < MATLAB_BASIC) return M_Basic(d).empty();
		else if (type == MATLAB_STRING) return M_String(d).empty();
		else if (type == MATLAB_STRUCT) return M_Struct(d).empty();
		else if (type == MATLAB_CELL) return M_Cell(d).empty();
		return true;
	}
	T_Basic& castBasic()
	{ 
		setType(MATLAB_DOUBLE);
		return M_Basic(d);
	}
	T_String& castString()
	{ 
		setType(MATLAB_STRING);
		return M_String(d); 
	}
	T_Struct& castStruct()
	{ 
		setType(MATLAB_STRUCT);
		return M_Struct(d);
	}
	T_Cell& castCell()
	{
		setType(MATLAB_CELL);
		return M_Cell(d); 
	}

	const T_Basic& cBasic() const
	{
		return M_Basic(d);
	}
	const T_String& cString() const
	{
		return M_String(d);
	}
	const T_Struct& cStruct() const
	{
		return M_Struct(d);
	}
	const T_Cell& cCell() const
	{
		return M_Cell(d);
	}

	//for single value assignment
	MatlabType& operator =(const double val) // If the original is basic, then it wouldn't change type
	{
		if (type > MATLAB_BASIC) // not basic type
		{
			delete d;
			d = new ArrayMgr<double>;
			type = MATLAB_DOUBLE;
		}
		M_Basic(d).resize(1);
		M_Basic(d).a(0) = val;
		return *this;
	}
	MatlabType& operator =(const string& val)
	{
		if (type != MATLAB_STRING)
		{
			delete d;
			d = new string;
			type = MATLAB_STRING;
		}
		M_String(d) = val;
		return *this;
	}
	MatlabType& operator =(const MatlabType& r)
	{
		setType(r.getType());
		if (type < MATLAB_BASIC) // basic type
		{
			M_Basic(d) = M_Basic(r.getP());
		}
		else if (MATLAB_STRING == type)
		{
			M_String(d) = M_String(r.getP());
		}
		else if (MATLAB_STRUCT == type)
		{
			M_Struct(d) = M_Struct(r.getP());
		}
		else if (MATLAB_CELL == type)
		{
			M_Cell(d) = M_Cell(r.getP());
		}
		return *this;
	}
	MatlabType& operator [](const string& name) // for struct access
	{
		if (type != MATLAB_STRUCT) 
		{
			setType(MATLAB_STRUCT);
		}
		return M_Struct(d)[name];
	}
	void addField(const string& name, const MatlabType& val) // for struct only
	{
		if (type != MATLAB_STRUCT) setType(MATLAB_STRUCT);
		M_Struct(d)[name] = val;
	}

	MatlabType operator +(MatlabType& r)
	{
		MatlabType ans;
		if (type < MATLAB_BASIC && r.getType() < MATLAB_BASIC)
		{
			T_Basic& b = ans.castBasic();
			b = castBasic() + r.castBasic();
		}
		else if (type == MATLAB_STRING && r.getType() == MATLAB_STRING)
		{
			T_String& s = ans.castString();
			s = castString() + r.castString();
		}
		else printf("Warning: unsupported operation + , please check the type\n");
		return ans;
	}
	MatlabType operator -(MatlabType& r)
	{
		MatlabType ans;
		if (type < MATLAB_BASIC && r.getType() < MATLAB_BASIC)
		{
			T_Basic& b = ans.castBasic();
			b = this->castBasic() - r.castBasic();
		}
		else printf("Warning: unsupported operation + , please check the type\n");
		return ans;
	}
	MatlabType operator *(MatlabType& r)
	{
		MatlabType ans;
		if (type < MATLAB_BASIC && r.getType() < MATLAB_BASIC)
		{
			T_Basic& b = ans.castBasic();
			b = this->castBasic() * r.castBasic();
		}
		else printf("Warning: unsupported operation * , please check the type\n");
		return ans;
	}
	MatlabType operator /(MatlabType& r)
	{
		MatlabType ans;
		if (type < MATLAB_BASIC && r.getType() < MATLAB_BASIC)
		{
			T_Basic& b = ans.castBasic();
			b = this->castBasic() / r.castBasic();
		}
		else printf("Warning: unsupported operation / , please check the type\n");
		return ans;
	}

	MatlabType& operator +=(MatlabType& r)
	{
		if (type < MATLAB_BASIC && r.getType() < MATLAB_BASIC)
		{
			castBasic() += r.castBasic();
		}
		else if (type == MATLAB_STRING && r.getType() == MATLAB_STRING)
		{
			castString() += r.castString();
		}
		else printf("Warning: unsupported operation + , please check the type\n");
		return *this;
	}
	MatlabType& operator -=(MatlabType& r)
	{
		if (type < MATLAB_BASIC && r.getType() < MATLAB_BASIC)
		{
			this->castBasic() -= r.castBasic();
		}
		else printf("Warning: unsupported operation + , please check the type\n");
		return *this;
	}
	MatlabType& operator *=(MatlabType& r)
	{
		if (type < MATLAB_BASIC && r.getType() < MATLAB_BASIC)
		{
			this->castBasic() *= r.castBasic();
		}
		else printf("Warning: unsupported operation * , please check the type\n");
		return *this;
	}
	MatlabType& operator /=(MatlabType& r)
	{
		if (type < MATLAB_BASIC && r.getType() < MATLAB_BASIC)
		{
			this->castBasic() /= r.castBasic();
		}
		else printf("Warning: unsupported operation / , please check the type\n");
		return *this;
	}
	template<class Type>
	Type scast() // single cast
	{
		if (type < MATLAB_BASIC)
		{
			return castBasic().scast<Type>();
		}
		else
		{
			cout << "Warning: unsupported cast type!" << endl;
			return Type();
		}
	}
	
	// basic type matrix access
	double& a(const int32_t i)
	{
		if (type >= MATLAB_BASIC) throw std::runtime_error("It's not MATLAB basic type as expected");
		T_Basic* pb = (T_Basic*)d;
		if (pb->empty()) throw std::runtime_error("Accessing empty MATLAB array");
		return M_Basic(d).a(i);
	}
	double& a(const int32_t i, const int32_t j)
	{
		if (type >= MATLAB_BASIC) throw std::runtime_error("It's not MATLAB basic type as expected");
		T_Basic* pb = (T_Basic*)d;
		if (pb->empty()) throw std::runtime_error("Accessing empty MATLAB array");
		return M_Basic(d).a(i, j);
	}
	double& a(const int32_t i, const int32_t j, const int32_t k)
	{
		if (type >= MATLAB_BASIC) throw std::runtime_error("It's not MATLAB basic type as expected");
		T_Basic* pb = (T_Basic*)d;
		if (pb->empty()) throw std::runtime_error("Accessing empty MATLAB array");
		return M_Basic(d).a(i, j, k);
	}
	double& a(const int32_t i, const int32_t j, const int32_t k, const int32_t m)
	{
		if (type >= MATLAB_BASIC) throw std::runtime_error("It's not MATLAB basic type as expected");
		T_Basic* pb = (T_Basic*)d;
		if (pb->empty()) throw std::runtime_error("Accessing empty MATLAB array");
		return M_Basic(d).a(i, j, k, m);
	}

	// cell elements access, no type checking performed for speed concern
	MatlabType& cell(const int32_t i)
	{
		return M_Cell(d).a(i);
	}
	MatlabType& cell(const int32_t i, const int32_t j)
	{
		return M_Cell(d).a(i, j);
	}
	MatlabType& cell(const int32_t i, const int32_t j, const int32_t k)
	{
		return M_Cell(d).a(i, j, k);
	}
	MatlabType& cell(const int32_t i, const int32_t j, const int32_t k, const int32_t m)
	{
		return M_Cell(d).a(i, j, k, m);
	}

	bool resize(int n)
	{
		if (type < MATLAB_BASIC) return M_Basic(d).resize(n);
		else if (type == MATLAB_STRING) 
		{
			M_String(d).resize(n); 
			return true;
		}
		else if (type == MATLAB_CELL) return M_Cell(d).resize(n);
		else return false;
	}
	template<class T>
	void cast(ArrayMgr<T>& in)
	{
		if (type > MATLAB_BASIC) return; // not basic type, so do nothing
		M_Basic(d).cast<T>(in);
	}
	
	void print(int dent=0) // print the content, if it's a long array, just print the first 3 items
	{
		if (type < MATLAB_BASIC)
		{
			int numElem = M_Basic(d).getInnerLength();
			if (numElem > 1)
			{
				printDent(dent);
				cout << "(array size " << numElem << ")" << endl;
			}
			printDent(dent);		
			for (int i = 0; i < 3 && i < numElem; ++i)
			{
				if (type <= MATLAB_FLOAT)
				{
					cout << M_Basic(d).a(i) << " ";
				}
				else if (type < MATLAB_BASIC)
				{
					int pi = int(M_Basic(d).a(i));
					cout << pi << " ";
				}
			}
			if (numElem > 3) cout << "...";
			cout << endl;
		}
		else if (type == MATLAB_STRING)
		{
			printDent(dent);
			cout << M_String(d) << endl;
		}
		else if (type == MATLAB_STRUCT)
		{
			printDent(dent);
			cout << "(struct size " << M_Struct(d).size() << ")" << endl;
			
			for (map<string, MatlabType>::iterator it = M_Struct(d).begin(); it != M_Struct(d).end(); ++it)
			{
				printDent(dent);
				cout << it->first << ":" << endl;
				it->second.print(dent + 1);
			}
		}
		else if (type == MATLAB_CELL)
		{
			printDent(dent);
			cout << "(cell ";
			M_Cell(d).printSize();
			cout << ")" << endl;
			int numElem = M_Cell(d).getInnerLength();
			for (int i = 0; i < 3 && i < numElem; ++i) M_Cell(d).a(i).print(dent + 1);
			if (numElem > 3)
			{
				printDent(dent);
				cout << "..." << endl;
			}
		}
	}
private:
	void *d; // d should be always none-null
	MDataType type;
	void printDent(int n) // each dent takes two blank spaces
	{
		for (int i = 0; i < n; ++i) cout << "  ";
	}
};

#ifdef WIN32
#define DEXPORT __declspec(dllexport)
#else
#define DEXPORT __attribute__ ((visibility ("default")))
#endif

class DEXPORT DicomFile
{
public:
	DicomFile();
	DicomFile(const string& fname);
	~DicomFile();

	bool open(const string& fname);
	bool save(const string& fname);
	
	MatlabType getItem(const string& TagName); // get one item. If you know their data layout structure
	MatlabType getAllItems();
	unsigned char* getPixelData();

	bool setItem(const string& TagName, const MatlabType& val);
	bool setPixelData(void* pin, unsigned long bytes);

private:
	void* dcmFF;
	MatlabType getItem(void *e); // it will get the element recursively
};

#endif
