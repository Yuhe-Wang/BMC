#pragma once
#ifndef _TOOLS_H_
#define _TOOLS_H_
#define _CRT_SECURE_NO_WARNINGS //close the thread warning in VC++
/*This file contains useful tools */

#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <stack>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>
#include <utility>
#include <functional>
#include <atomic>
#include <thread>
#include <mutex>
#include <chrono>
#include <random>
#include <valarray>
#include <new>
#include <regex>
#include <stdarg.h> 
#include <cstdint>

#ifdef WIN32
#define DECLSPECIFIER __declspec(dllexport)
#define EXPIMP_TEMPLATE
#pragma warning (disable : 4251)
#else // unix like
#define DECLSPECIFIER __attribute__ ((visibility ("default")))
#endif

//here are some macro definition switches
#define USE_OPENMP
#define CUDA_KERNEL_CHECK
typedef float SFloat;

#ifdef WIN32
#include <filesystem> // Will be officially available in C++ 17 standard
namespace fs = std::tr2::sys; //may be changed in the future or in other platforms
#include <windows.h> //windows specific
#define getLibAddress GetProcAddress
#define freeLib FreeLibrary
#define FileSeparator '\\'
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#include <dlfcn.h>
#include <unistd.h>
#include <signal.h>
#include <sys/sysinfo.h>
#define getLibAddress dlsym
#define freeLib dlclose
#define FileSeparator '/'
#endif

#ifdef USE_OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#ifdef _DEBUG
//#include <vld.h> //help to locate the memory leakage
#endif

extern int NProcess, pID;

using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::list;
using std::vector;
using std::stack;
using std::getline;
using std::max;
using std::min;
using std::pair;
using std::stringstream;
using std::make_pair;
using std::map;
using std::ifstream;
using std::deque;
using std::valarray;
using std::bad_alloc;
using std::regex;

typedef vector<pair<string, string>> MPS; //macro pair string

const double PI = 3.141592653589793238463;
const double Es = 5.10998902e5; //electron or positron static energy , unit eV
const double TEs = 2 * Es;
const double lightSpeed = 299792458; //unit m/s
const int PRNGBufferSize = 512;
const double SAD = 105.0; //unit cm, this is for ViewRay

void exitApp(const char *inf);
extern "C" DECLSPECIFIER bool existFile(const char* filename);
extern "C" DECLSPECIFIER long fileSize(FILE *fp);
extern "C" DECLSPECIFIER int get_thread_num();
DECLSPECIFIER void replaceAll(string& str, const string& from, const string& to);
extern "C" DECLSPECIFIER double nodeSpeedTest(int NThread);
extern "C" DECLSPECIFIER void matchFile(const string& dir, string& match);
extern "C" DECLSPECIFIER void pauseSeconds(const unsigned int t_seconds);

DECLSPECIFIER string changeExtension(const char* fname, const char* newExt);
DECLSPECIFIER string printWithComma(const double fnum);
DECLSPECIFIER string findFirstFileInDir(const string& dir, const string& match);

inline void readline(FILE *fp)
{
	while (fgetc(fp) != '\n');
}

inline int binarySearch(double *a, int n /*array length>2*/, double x /*input x to search*/, int shift = 0) //return left index <=n-1
{
	//if a is two dimensional array a[m][n], and we specify the i-th line's content a[i][?], then
	//we need to input the parameter shift= i*n
	//note that a[i][?] must in increasing order and x must in this range
	//or you can input &a[i][0] to a and set shift=0
	int iL = 0, iR = n - 1;
	int iM;
	do
	{
		iM = (iL + iR) / 2;
		if (x > a[shift + iM]) iL = iM;
		else iR = iM;
	} while (iR - iL > 1);

	return iL;
}

typedef void(*LogCallBack)(const char*);
class DECLSPECIFIER LogProvider
{
public:
	LogProvider();
	~LogProvider();
	void setLogName(int i, const string mode = "w", string cdir = "");// the log name would be Process-i.log
	void setLogName(const char* name);
	string getLogName(){ return logName; }
	void closeFile();
	void flush();
	void operator()(const char *fmt, ...);
	void setCallBack(LogCallBack pfunc);
	char* timeNow();
	void warningColor(); //change the text to alerting color
	void normalColor(); //change the text to normal color
	void shortMode(bool mode = true);// only output in one line on screen 
private:
	FILE* fp;
	string logName;
	bool b_shortMode;
	LogCallBack callBackFunc;
};
DECLSPECIFIER extern LogProvider Log;

#include "./PRNG.h"

class DECLSPECIFIER RunTimeCounter
{
	//methods
public:
	RunTimeCounter()
	{
		tstart = std::chrono::system_clock::now();
		tlast = tstart;
		tAccu = -1;
	}
	double getStoredTime()
	{
		return tAccu;
	}
	void resetStoredTime()
	{
		tAccu = 0;
	}
	void start() // It will reset tlast and tstart, but won't reset tAccu (enable time accumulation if not yet)
	{
		tstart = std::chrono::system_clock::now();
		tlast = tstart;
		if (tAccu < 0) tAccu = 0; // enable tAccu at the first time calling start()
	}
	double stop(bool sinceStart = false) // This function will always add tAccu (time from tlast to tnow)
	{
		auto tnow = std::chrono::system_clock::now();
		std::chrono::duration<double> dt;
		if (sinceStart)
		{
			if (tAccu >= 0)
			{
				dt = tnow - tlast;
				tAccu += dt.count();
			}
			tlast = tnow; // always update tlast when calling stop()
			dt = tnow - tstart;
			return dt.count();
		}
		else
		{
			dt = tnow - tlast;
			tlast = tnow; // always update tlast when calling stop()
			double tseconds = dt.count();
			if (tAccu >= 0) tAccu += tseconds;
			return tseconds;
		}
	}
	//data
private:
	std::chrono::time_point<std::chrono::system_clock> tstart, tlast;
	double tAccu; // if you want to use tAccu, you must call start() at least once
};

// this definition will be used in both CPU and GPU simulation code
#define at(x,y,z) ((x) + NX*((y)+ NY*(z))) //Column first memory layout, this should be the same as ArrayManager

const int32_t MATMGR_MAX_DIM = 4; // The max number of dimensions that ArrayMgr will directly handle
//the memory layout is column first just like MATLAB and FORTRAN
template<class Type>
class ArrayMgr
{
public:
	//shouldn't have default constructor
	ArrayMgr()
	{
		init();
	}
	ArrayMgr(const int32_t nx)
	{
		init();
		resize(nx);
	}
	bool resize(const int32_t nx)
	{
		numDimensions = 1;
		size[0] = nx;
		return resize();
	}
	ArrayMgr(const int32_t nx, const int32_t ny)
	{
		init();
		resize(nx, ny);
	}
	bool resize(const int32_t nx, const int32_t ny)
	{
		numDimensions = 2;
		size[0] = nx;
		size[1] = ny;
		return resize();
	}
	ArrayMgr(const int32_t nx, const int32_t ny, const int32_t nz)
	{
		init();
		resize(nx, ny, nz);
	}
	bool resize(const int32_t nx, const int32_t ny, const int32_t nz)
	{
		numDimensions = 3;
		size[0] = nx;
		size[1] = ny;
		size[2] = nz;
		return resize();
	}
	ArrayMgr(const int32_t nx, const int32_t ny, const int32_t nz, const int32_t ns)
	{
		init();
		resize(nx, ny, nz, ns);
	}
	bool resize(const int32_t nx, const int32_t ny, const int32_t nz, const int32_t ns)
	{
		numDimensions = 4;
		size[0] = nx;
		size[1] = ny;
		size[2] = nz;
		size[3] = ns;
		return resize();
	}
	ArrayMgr(const ArrayMgr<Type>& r) //copy constructor
	{
		init();
		canFreeData = r.canFreeData;
		//copy the dimensions
		numDimensions = r.numDimensions;
		for (int32_t i = 0; i < numDimensions; ++i) size[i] = r.size[i];
		if (canFreeData)
		{
			resize(); //request memory
			//copy the data content
			int32_t length = getInnerLength();
			for (int32_t i = 0; i < length; ++i) data[i] = r.data[i];
		}
		else data = r.data;
	}
	~ArrayMgr()
	{
		if (data&&canFreeData) delete[] data;
	}

	// Image related
	template<class T>
	bool createImage(const int32_t np, ArrayMgr<T>& img)
	{
		// create an image whose pointer can be directly accessed by C image process program such as openGL
		// np is the number of variables for each pixel
		if (numDimensions != 2) return false;
		bool result = img.resize(np, size[1], size[0]);
		for (int32_t i = 0; i < np; ++i)
			for (int32_t j = 0; j < size[1]; ++j)
				for (int32_t k = 0; k < size[0]; ++k)
					img.a(i, j, k) = a(k, j);
		return result;
	}
	inline Type& a_img(const int32_t& x, const int32_t& y, const int32_t& ip) //access the pixel corresponding to the original matrix
	{
		return a(ip, y, x);
	}

	inline Type& a(const int32_t &x)//access the element, please choose the right one only
	{
		return data[x];
	}
	inline Type& a(const int32_t &x, const int32_t &y)//access the element, please choose the right one only
	{
		return data[x + y*size[0]];
	}
	inline Type& a(const int32_t &x, const int32_t &y, const int32_t &z)//access the element, please choose the right one only
	{
		return data[x + size[0] * (y + z*size[1])];
	}
	inline Type& a(const int32_t &x, const int32_t &y, const int32_t &z, const int32_t& s)//access the element, please choose the right one only
	{
		return data[x + size[0] * (y + size[1] * (z + size[2] * s))];
	}
	inline const Type& operator[](const int32_t &x) const// for constant access, equal to a(x) for one dimension
	{
		return data[x];
	}
	
	bool deepCopy(const ArrayMgr<Type>& r) // if the old one has no inner copy, call deepCopy(*this) will change itself to deep copy.
	{
		if (&r == this) //may change the canFreeData
		{
			if (!canFreeData) // the old one has no inner copy
			{
				canFreeData = true;
				resize();  //request memory
				Type* p = data;
				int32_t length = getInnerLength();
				for (int32_t i = 0; i < length; ++i) data[i] = p[i];
			}
		}
		canFreeData = true;
		//copy the dimensions
		numDimensions = r.numDimensions;
		for (int32_t i = 0; i < numDimensions; ++i) size[i] = r.size[i];
		resize(); //request memory
		//copy the data content
		int32_t length = getInnerLength();
		for (int32_t i = 0; i < length; ++i) data[i] = r.data[i];
		return true;
	}
	bool shallowCopy(const ArrayMgr<Type>& r)
	{
		if (&r == this) //may change the canFreeData
		{
			if (canFreeData) return false;
			else return true;
		}
		canFreeData = false;
		//copy the dimensions
		numDimensions = r.numDimensions;
		for (int32_t i = 0; i < numDimensions; ++i) size[i] = r.size[i];
		if (allocatedSize > 0) // release previous allocated data
		{
			delete[] data;
			allocatedSize = 0;
		}
		data = r.data;
		return true;
	}
	ArrayMgr& operator =(const ArrayMgr<Type>& r) // note this operation will keep the attribution of canFreeData. You may instead use deepCopy, shallowCopy explicitly
	{
		if (&r == this) return *this; // copy to itself
		canFreeData = r.canFreeData;
		//copy the dimensions
		numDimensions = r.numDimensions;
		for (int32_t i = 0; i < numDimensions; ++i) size[i] = r.size[i];
		if (canFreeData)
		{
			resize(); //request memory
			//copy the data content
			int32_t length = getInnerLength();
			for (int32_t i = 0; i < length; ++i) data[i] = r.data[i];
		}
		else // now there's no inner managed array
		{
			if (allocatedSize > 0) // release previous allocated data
			{
				delete[] data;
				allocatedSize = 0;
			}
			data = r.data;
		}
		return *this;
	}
	ArrayMgr& operator +=(const ArrayMgr<Type>& r)
	{
		if (equalDims(r))
		{
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i) data[i] += r[i];
		}
		else if (r.getInnerLength() == 1)
		{
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i) data[i] += r[0];
		}
		else printf("error:dimension doesn't match for += operation between arrays!\n");
		return *this;
	}
	ArrayMgr& operator -=(const ArrayMgr<Type>& r)
	{
		if (equalDims(r))
		{
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i) data[i] -= r[i];
		}
		else if (r.getInnerLength() == 1)
		{
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i) data[i] -= r[0];
		}
		else printf("error:dimension doesn't match for -= operation between arrays!\n");
		return *this;
	}
	ArrayMgr& operator *=(const ArrayMgr<Type>& r)
	{
		if (equalDims(r))
		{
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i) data[i] *= r[i];
		}
		else if (r.getInnerLength() == 1)
		{
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i) data[i] *= r[0];
		}
		else printf("error:dimension doesn't match for *= operation between arrays!\n");
		return *this;
	}
	ArrayMgr& operator /=(const ArrayMgr<Type>& r)
	{
		if (equalDims(r))
		{
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i) data[i] /= r[i];
		}
		else if (r.getInnerLength() == 1)
		{
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i) data[i] /= r[0];
		}
		else printf("error:dimension doesn't match for *= operation between arrays!\n");
		return *this;
	}

	ArrayMgr& operator *=(const double mul) // the element must support *= double operation
	{
		int32_t len = getInnerLength();
		for (int32_t i = 0; i < len; ++i) data[i] = Type(data[i] * mul);
		return *this;
	}
	ArrayMgr& operator +=(const double addv) // the element must support *= double operation
	{
		int32_t len = getInnerLength();
		for (int32_t i = 0; i < len; ++i) data[i] = Type(data[i] + addv);
		return *this;
	}
	ArrayMgr& operator -=(const double minusv) // the element must support *= double operation
	{
		int32_t len = getInnerLength();
		for (int32_t i = 0; i < len; ++i) data[i] = Type(data[i] - minusv);
		return *this;
	}
	ArrayMgr& operator /=(const double denominator) //the element must support *= double operation
	{
		int32_t len = getInnerLength();
		for (int32_t i = 0; i < len; ++i) data[i] = Type(data[i] / denominator);
		return *this;
	}
	ArrayMgr& operator =(const Type& val) // the class must support element=val operation
	{
		int32_t length = getInnerLength();
		for (int32_t i = 0; i < length; ++i) data[i] = val;
		return *this;
	}
	ArrayMgr& setValue(const Type& val){ return *this = val; }

	ArrayMgr operator+(const ArrayMgr<Type>& r) //the return array is a deep copy
	{
		ArrayMgr ans;
		if (r.getInnerLength() == 1)
		{
			ans.deepCopy(*this);
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i)
			{
				ans.a(i) += r[0];
			}
		}
		else if (getInnerLength() == 1)
		{
			ans.deepCopy(r);
			int32_t len = r.getInnerLength();
			for (int32_t i = 0; i < len; ++i)
			{
				ans.a(i) += a(0);
			}
		}
		else if (equalDims(r))
		{
			ans.deepCopy(*this);
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i)
			{
				ans.a(i) += r[i];
			}
		}
		else printf("Warning: unsupported operator + , please check the dimensions\n");
		return ans;
	}
	ArrayMgr operator-(const ArrayMgr<Type>& r)
	{
		ArrayMgr ans;
		if (r.getInnerLength() == 1)
		{
			ans.deepCopy(*this);
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i)
			{
				ans.a(i) -= r[0];
			}
		}
		else if (getInnerLength() == 1)
		{
			ans.deepCopy(r); // get the size
			ans = a(0); // set to one value
			int32_t len = r.getInnerLength();
			for (int32_t i = 0; i < len; ++i)
			{
				ans.a(i) -= r[i];
			}
		}
		else if (equalDims(r))
		{
			ans.deepCopy(*this);
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i)
			{
				ans.a(i) *= r[i];
			}
		}
		else printf("Warning: unsupported operator - , please check the dimensions\n");
		return ans;
	}
	ArrayMgr operator*(const ArrayMgr<Type>& r)
	{
		ArrayMgr ans;
		if (r.getInnerLength() == 1)
		{
			ans.deepCopy(*this);
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i)
			{
				ans.a(i) *= r[0];
			}
		}
		else if (getInnerLength() == 1)
		{
			ans.deepCopy(r);
			int32_t len = r.getInnerLength();
			for (int32_t i = 0; i < len; ++i)
			{
				ans.a(i) *= a(0);
			}
		}
		else if (equalDims(r))
		{
			ans.deepCopy(*this);
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i)
			{
				ans.a(i) *= r[i];
			}
		}
		else printf("Warning: unsupported operator * , please check the dimensions\n");
		return ans;
	}
	ArrayMgr operator/(const ArrayMgr<Type>& r)
	{
		ArrayMgr ans;
		if (r.getInnerLength() == 1)
		{
			ans.deepCopy(*this);
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i)
			{
				ans.a(i) /= r[0];
			}
		}
		else if (getInnerLength() == 1)
		{
			ans.deepCopy(r); // get the size
			ans = a(0); // set to one value
			int32_t len = r.getInnerLength();
			for (int32_t i = 0; i < len; ++i)
			{
				ans.a(i) /= r[i];
			}
		}
		else if (equalDims(r))
		{
			ans.deepCopy(*this);
			int32_t len = getInnerLength();
			for (int32_t i = 0; i < len; ++i)
			{
				ans.a(i) /= r[i];
			}
		}
		else printf("Warning: unsupported operator / , please check the dimensions\n");
		return ans;
	}

	// Operations to exact information
	int32_t getInnerLength() const
	{
		if (numDimensions == 0) return 0;
		int32_t len = 1;
		for (int32_t i = 0; i < numDimensions; ++i) len *= size[i];
		return len;
	}
	int32_t getByteNum() const { return sizeof(Type)*getInnerLength(); }
	int32_t getDim(int32_t* s=NULL) const
	{ 
		if (s != NULL)
		{
			for (int i = 0; i < numDimensions; ++i) s[i] = size[i];
		}
		return numDimensions; 
	}
	template<class Type2>
	bool equalDimsT(const ArrayMgr<Type2>& r)
	{
		int32_t s[4];
		int r_dim = r.getDim(s);
		if (numDimensions != r_dim) return false;
		else
		{
			for (int32_t i = 0; i < numDimensions; ++i)
			{
				if (size[i] != s[i]) return false;
			}
			return true;
		}
	}
	bool equalDims(const ArrayMgr<Type>& r) // tell if the two matrix has the same dimension
	{
		int32_t s[4];
		int r_dim = r.getDim(s);
		if (numDimensions != r_dim) return false;
		else
		{
			for (int32_t i = 0; i < numDimensions; ++i)
			{
				if (size[i] != s[i]) return false;
			}
			return true;
		}
	}
	int32_t getWidth(int32_t d) //d starts from 1
	{
		if (d > 4)
		{
			printf("error: unsupported dimension in getWidth().\n");
			return 0;
		}
		return size[d - 1];
	}
	void getMaxMinIndex(int32_t &maxi, int32_t &mini) // get the internal indice of the max and min value
	{
		maxi = mini = 0;
		int32_t len = getInnerLength();
		for (int32_t i = 1; i < len; ++i)
		{
			if (data[i] > data[maxi]) maxi = i;
			if (data[i] < data[mini]) mini = i;
		}
	}
	void getMaxMin(Type *maxp, Type *minp) //get the internal pointer of the max and min value
	{
		int32_t maxi = 0, mini = 0;
		int32_t length = getInnerLength();
		for (int32_t i = 1; i < length; ++i)
		{
			if (data[i] > data[maxi]) maxi = i;
			if (data[i] < data[mini]) mini = i;
		}
		maxp = data + maxi;
		minp = data + mini;
	}
	void getMaxMin(Type &maxv, Type &minv) // get the max and min value. request Type support copy
	{
		int32_t maxi = 0, mini = 0;
		int32_t length = getInnerLength();
		for (int32_t i = 1; i < length; ++i)
		{
			if (data[i] > data[maxi]) maxi = i;
			if (data[i] < data[mini]) mini = i;
		}
		maxv = data[maxi];
		minv = data[mini];
	}
	Type& Max() // get the max value
	{
		int32_t maxi = 0;
		int32_t length = getInnerLength();
		for (int32_t i = 1; i < length; ++i)
		{
			if (data[i] > data[maxi]) maxi = i;
		}
		return data[maxi];
	}
	Type& Min() // get the min value
	{
		int32_t mini = 0;
		int32_t length = getInnerLength();
		for (int32_t i = 1; i < length; ++i)
		{
			if (data[i] < data[mini]) mini = i;
		}
		return data[mini];
	}
	bool empty()
	{
		if (getInnerLength() > 0) return false;
		else return true;
	}
	bool isUniform()
	{
		if (numDimensions == 0) return true;
		int32_t len = getInnerLength();
		Type val = data[0];
		for (int32_t i = 0; i < len; ++i)
		{
			if (data[i] != val)
			{
				return false;
			}
		}
		return true;
	}
	int32_t countAbove(const Type& Val)
	{
		int32_t len = getInnerLength();
		int32_t NC = 0;
		for (int32_t i = 0; i < len; ++i)
		{
			if (data[i] >= Val) ++NC;
		}
		return NC;
	}

	// Index transform
	void ind2sub(const int32_t in, int32_t& ox, int32_t& oy)
	{
		if (numDimensions != 2)
		{
			printf("Wrong dimensions. Expect the number of dimension to be 2\n");
			return;
		}
		oy = in / size[0];
		ox = in - size[0] * oy;
	}
	void ind2sub(const int32_t in, int32_t& ox, int32_t& oy, int32_t& oz)
	{
		if (numDimensions != 3)
		{
			printf("Wrong dimensions. Expect the number of dimension to be 3\n");
			return;
		}
		int32_t usize = size[0] * size[1];
		oz = in / usize;
		int32_t rest = in - usize*oz;
		oy = rest / size[0];
		ox = rest - size[0] * oy;
	}
	void ind2sub(const int32_t in, int32_t& ox, int32_t& oy, int32_t& oz, int32_t& os)
	{
		if (numDimensions != 4)
		{
			printf("Wrong dimensions. Expect the number of dimension to be 4\n");
			return;
		}
		int32_t usize = size[0] * size[1] * size[2];
		os = in / (usize);
		int32_t rest = in - usize*os;
		usize = size[0] * size[1];
		oz = rest / usize;
		rest -= usize*oz;
		oy = rest / size[0];
		ox = rest - size[0] * oy;
	}
	inline int32_t sub2ind(int32_t& x, int32_t& y)
	{
		return x + size[0] * y;
	}
	inline int32_t sub2ind(int32_t& x, int32_t& y, int32_t& z)
	{
		return x + size[0] * (y + size[1] * z);
	}
	inline int32_t sub2ind(int32_t& x, int32_t& y, int32_t& z, int32_t& s)
	{
		return x + size[0] * (y + size[1] * (z + size[2] * s));
	}

	// Memory management
	void release() //release the data manually, i.e. turn it into an empty matrix
	{
		if (data != NULL)
		{
			if(canFreeData) delete[] data;
			data = NULL;
			allocatedSize = 0;
		}
		numDimensions = 0;
	}
	bool setReserveMem(const int32_t nr)
	{
		if (!canFreeData) return false;
		if (nr > allocatedSize)
		{
			if (data) delete[] data;
			try
			{
				data = new Type[nr];
			}
			catch (const bad_alloc& e)
			{
				printf("bad_alloc caught:%s\n", e.what());
				return false;
			}
			allocatedSize = nr;
		}
		return true;
	}
	void toRowFirstLayout() // change the memory layout from MATLAB type column first (default) to C type row first
	{
		ArrayMgr<Type> temp;
		temp.deepCopy(*this);
		if (2 == numDimensions)
		{
			for (int ix = 0; ix < size[0]; ++ix)
				for (int iy = 0; iy < size[1]; ++iy)
				{
					data[ix*size[1] + iy] = temp.a(ix, iy);
				}
		}
		else if (3 == numDimensions)
		{

			for (int ix = 0; ix < size[0]; ++ix)
				for (int iy = 0; iy < size[1]; ++iy)
					for (int iz = 0; iz < size[2]; ++iz)
					{
						data[(ix*size[1] + iy)*size[2] + iz] = temp.a(ix, iy, iz);
					}
		}
		else if (4 == numDimensions)
		{
			for (int ix = 0; ix < size[0]; ++ix)
				for (int iy = 0; iy < size[1]; ++iy)
					for (int iz = 0; iz < size[2]; ++iz)
						for (int is = 0; is < size[3]; ++is)
						{
							data[((ix*size[1] + iy)*size[2] + iz)*size[3] + is] = temp.a(ix, iy, iz, is);
						}
		}
	}
	void toColFirstLayout() // change the memory layout from C type row first to MATLAB type column first (default)
	{
		ArrayMgr<Type> temp;
		temp.deepCopy(*this);
		Type* tdata = temp.getP();
		if (2 == numDimensions)
		{
			for (int ix = 0; ix < size[0]; ++ix)
				for (int iy = 0; iy < size[1]; ++iy)
				{
					this->a(ix, iy) = tdata[ix*size[1] + iy];
				}
		}
		else if (3 == numDimensions)
		{
			for (int ix = 0; ix < size[0]; ++ix)
				for (int iy = 0; iy < size[1]; ++iy)
					for (int iz = 0; iz < size[2]; ++iz)
					{
						this->a(ix, iy, iz) = tdata[(ix*size[1] + iy)*size[2] + iz];
					}
		}
		else if (4 == numDimensions)
		{
			for (int ix = 0; ix < size[0]; ++ix)
				for (int iy = 0; iy < size[1]; ++iy)
					for (int iz = 0; iz < size[2]; ++iz)
						for (int is = 0; is < size[3]; ++is)
						{
							this->a(ix, iy, iz, is) = tdata[((ix*size[1] + iy)*size[2] + iz)*size[3] + is];
						}
		}
	}
	Type* getP(){ return data; } // It's maybe not safe. Use with caution
	const Type* getConstP() const{ return data; }
	bool canFree(){ return canFreeData; }
	void  setP(Type* po) // It's maybe not safe. Use with caution
	{
		if (data&&canFreeData) delete[] data;
		data = po;
		canFreeData = false;
	}
	void detachP() // It's maybe not safe. Use with caution
	{
		release();
	}
	void copy(Type* p) // you should ensure the validation of p, and make sure the lengths match
	{
		if (NULL == p) return;
		int len = getInnerLength();
		for (int i = 0; i < len; ++i) data[i] = p[i];
	}
	bool resize(const int32_t nd, const int32_t* s) // only update the size array
	{
		numDimensions = nd;
		for (int32_t i = 0; i < numDimensions; ++i)
		{
			size[i] = s[i];
		}
		return resize();
	}
	template<class T>
	bool resize(ArrayMgr<T>& input_a) // resize to the same size of the external matrix
	{
		int sz[4];
		numDimensions = input_a.getDim(sz);
		for (int32_t i = 0; i < numDimensions; ++i)
		{
			size[i] = sz[i];
		}
		return resize();
	}
// 	template<class T>
// 	ArrayMgr<T>* cast() // need to delete the new array manually!!!
// 	{
// 		auto p = new ArrayMgr<T>;
// 		p->resize(numDimensions, size);
// 		int32_t len = getInnerLength();
// 		for (int i = 0; i < len; ++i) p->a(i) = (T) a(i); // the element must support this cast
// 		return p;
// 	}
	template<class T>
	void cast(ArrayMgr<T>& in) // need to delete the new array manually!!!
	{
		in.resize(numDimensions, size);
		int32_t len = getInnerLength();
		for (int i = 0; i < len; ++i) in.a(i) = (T)a(i); // the element must support this cast
	}
	template<class T>
	T scast() // single cast
	{
		return T(data[0]);
	}

	// Operations that may change the content
	bool swap(ArrayMgr<Type>& r) // swap the content if they have the same dimension
	{
		if (canFreeData&& r.canFree() && equalDims(r))
		{
			Type* temp = r.getP();
			r.setP(data);
			data = temp;
			return true;
		}
		return false;
	}
	void normalize()
	{
		if (0 == getInnerLength()) return;
		Type maxv = Max();
		(*this) /= maxv;
	}
	void flip3(const int32_t di)
	{
		if (numDimensions != 3) return;
		int32_t nx = size[0];
		int32_t ny = size[1];
		int32_t nz = size[2];

		if (1 == di)
		{
			for (int32_t ix = 0; ix <= nx / 2; ++ix)
				for (int32_t iy = 0; iy < ny; ++iy)
					for (int32_t iz = 0; iz < nz; ++iz)
					{
						//swap between m(ix, iy, iz) and m(nx-1-ix, iy, iz)
						Type temp = this->a(ix, iy, iz);
						this->a(ix, iy, iz) = this->a(nx - 1 - ix, iy, iz);
						this->a(nx - 1 - ix, iy, iz) = temp;
					}
		}
		else if (2 == di)
		{
			for (int32_t iy = 0; iy <= ny / 2; ++iy)
				for (int32_t ix = 0; ix < nx; ++ix)
					for (int32_t iz = 0; iz < nz; ++iz)
					{
						//swap between m(ix, iy, iz) and m(ix, ny - 1 - iy, iz)
						Type temp = this->a(ix, iy, iz);
						this->a(ix, iy, iz) = this->a(ix, ny - 1 - iy, iz);
						this->a(ix, ny - 1 - iy, iz) = temp;
					}
		}
		else if (3 == di)
		{
			for (int32_t iz = 0; iz <= nz / 2; ++iz)
				for (int32_t ix = 0; ix < nx; ++ix)
					for (int32_t iy = 0; iy < ny; ++iy)
					{
						//swap between m(ix, iy, iz) and m(ix, iy, nz - 1 - iz)
						Type temp = this->a(ix, iy, iz);
						this->a(ix, iy, iz) = this->a(ix, iy, nz - 1 - iz);
						this->a(ix, iy, nz - 1 - iz) = temp;
					}
		}
	}
	void flipup()
	{
		if (numDimensions != 2) return;
		for (int i = 0; i < size[0] / 2; ++i)
			for (int j = 0; j < size[1]; ++j)
			{
				Type temp = a(i, j);
				int32_t ic = size[0] - 1 - i;
				a(i, j) = a(ic, j);
				a(ic, j) = temp;
			}
	}
	void fliplr()
	{
		if (numDimensions != 2) return;
		for (int i = 0; i < size[0]; ++i)
			for (int j = 0; j < size[1]/2; ++j)
			{
				Type temp = a(i, j);
				int32_t jc = size[1] - 1 - j;
				a(i, j) = a(i, jc);
				a(i, jc) = temp;
			}
	}
	void transpose()
	{
		if (numDimensions != 2) return;
		ArrayMgr<Type> temp;
		temp.deepCopy(*this);
		resize(size[1], size[0]);
		for (int i = 0; i < size[0]; ++i)
			for (int j = 0; j < size[0]; ++j)
				a(i, j) = temp.a(j, i);
	}

	//for 3D operation
	void slice(const int id, const int index, ArrayMgr<Type>& sm) //get slice from one dimension
	{
		if (numDimensions != 3) return;
		if (1 == id)
		{
			sm.resize(size[1], size[2]);
			for (int i = 0; i < size[1]; ++i)
				for (int j = 0; j < size[2]; ++j)
					sm.a(i, j) = a(index, i, j);
		}
		else if (2 == id)
		{
			sm.resize(size[0], size[2]);
			for (int i = 0; i < size[0]; ++i)
				for (int j = 0; j < size[2]; ++j)
					sm.a(i, j) = a(i, index, j);
		}
		else // 3==id
		{
			sm.resize(size[0], size[1]);
			for (int i = 0; i < size[0]; ++i)
				for (int j = 0; j < size[1]; ++j)
					sm.a(i, j) = a(i, j, index);
		}
	}


	void printSize()
	{
		cout << "size: ";
		for (int i = 0; i < numDimensions; ++i)
		{
			cout << size[i];
			if (i != numDimensions - 1) cout << "*";
		}
	}
private:
	Type* data;
	int32_t allocatedSize;
	int32_t size[MATMGR_MAX_DIM];
	int32_t numDimensions;
	bool canFreeData; // default true

	void init()
	{
		data = NULL;
		allocatedSize = 0;
		numDimensions = 0;
		canFreeData = true;
	}
	bool resize() //we support 4 dimensional memory management
	{
		if (!canFreeData) return false;
		int32_t len = 1;
		for (int32_t i = 0; i < numDimensions; ++i) len *= size[i];
		if (len > allocatedSize) //requested size is bigger than previously allocated
		{
			if (data) delete[] data;
			try
			{
				data = new Type[len];
			}
			catch (const bad_alloc& e)
			{
				printf("bad_alloc caught:%s\n", e.what());
				return false;
			}
			allocatedSize = len;
		}
		return true;
	}
};

DECLSPECIFIER bool interp3N(ArrayMgr<SFloat>& in, ArrayMgr<SFloat>& out, int N);
DECLSPECIFIER int  gammaAnalysis(ArrayMgr<SFloat>& d1, ArrayMgr<SFloat>& d2, ArrayMgr<SFloat>& g, double dx, double dy, double dz, double percent, double DTA, bool locale = false, int NInterp = 1);
DECLSPECIFIER int  gammaAnalysis2(ArrayMgr<SFloat>& d1, ArrayMgr<SFloat>& d2, ArrayMgr<SFloat>& g, double dx, double dy, double dz, double percent, double DTA, bool locale = false, int NInterp = 1);

template<class T>
void swapData(T& d1, T& d2)
{
	T temp;
	temp = d1;
	d1 = d2;
	d2 = temp;
}

class DECLSPECIFIER ConfigFile //class used to parse the config file
{
public:
	ConfigFile() :iLast(-1){}
	ConfigFile(const string &input);
	ConfigFile(const char *filename);
	~ConfigFile();
	ConfigFile* getBlock(const string &key, bool followLast = false);
	vector<ConfigFile*> getBlockArray(const string &key);
	bool getValue(const string &key, string &value, bool followLast = false);
	bool getValue(const string &key, int &value, bool followLast = false);
	bool getValue(const string &key, double &value, bool followLast = false);
	bool getValue(const string &key, vector<string> &values, bool followLast = false);
	bool getValue(const string &key, vector<int> &values, bool followLast = false);
	bool getValue(const string &key, vector<double> &values, bool followLast = false);
	void resetSearchIndex(){ iLast = -1; }

	bool parse(const string &input);
	bool parse(const char *filename);
	void addBlock(const pair<string, ConfigFile*> &block){ blocks.push_back(block); }
	void addKey(const pair<string, string> &key){ mps.push_back(key); }
	int getBlockNum(){ return (int)blocks.size(); }

	void copyMps(MPS& cp) { cp = mps; }
	void macroReplace(const MPS &ps);
	void macroReplace(const string &from, const string& to);
private:
	vector<pair<string, string> > mps; //key-value maps
	vector<pair<string, ConfigFile*> > blocks; //children blocks
	int iLast;
	void trim(string& s);
};

class DECLSPECIFIER EmailNotify
{
public:
	EmailNotify() :notify(false){}
	void init(ConfigFile* cf);
	bool doNotify();
	void send(bool normal = true, LogProvider* pLog = NULL);
private:
	bool notify;
	string address;
	string pwd;
	bool getEncryptedAccount(string& filename);
};

extern EmailNotify emailNotify;

#ifndef __CUDACC__
#define __device__ 
#define __host__
#define __forceinline__
#define __global__
#define __constant__
#endif

namespace MonteCarlo //to avoid duplicated definition
{
	class Vector {

	public:

		double x; //!< x-component
		double y; //!< y-component
		double z; //!< z-component

		Vector(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {};
		Vector(const Vector &v) : x(v.x), y(v.y), z(v.z) {};
		Vector() : x(0), y(0), z(0) {};

		Vector &operator=(const Vector &v) {
			x = v.x;
			y = v.y;
			z = v.z;
			return *this;
		};

		// Vector additions
		//
		Vector operator+(const Vector &v) const {
			return Vector(x + v.x, y + v.y, z + v.z);
		};
		Vector &operator+=(const Vector &v) {
			x += v.x;
			y += v.y;
			z += v.z;
			return *this;
		};

		// Vector subtractions
		//
		Vector operator-(const Vector &v) const {
			return Vector(x - v.x, y - v.y, z - v.z);
		}
		Vector &operator-=(const Vector &v) {
			x -= v.x;
			y -= v.y;
			z -= v.z;
			return *this;
		};

		// Vector multiplications
		//
		Vector operator*(const double f) const {
			return Vector(x*f, y*f, z*f);
		};
		Vector &operator*=(const double f) {
			x *= f;
			y *= f;
			z *= f;
			return *this;
		};
		friend Vector operator*(double f, Vector &v) {
			return v*f;
		};
		double operator*(const Vector &v) const {
			return x*v.x + y*v.y + z*v.z;
		};

		// vector product
		Vector times(const Vector &v) const {
			return Vector(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
		};
		Vector operator%(const Vector &v) const {
			return Vector(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
		};

		// scale
		Vector getScaled(const Vector &s) const {
			return Vector(x*s.x, y*s.y, z*s.z);
		};
		void scale(const Vector &s) {
			x *= s.x;
			y *= s.y;
			z *= s.z;
		};

		// Some other useful methods.
		double length() const {
			return sqrt(x*x + y*y + z*z);
		};
		double length2() const {
			return x*x + y*y + z*z;
		};
		void normalize() {
			double tmp = 1. / length();
			x *= tmp;
			y *= tmp;
			z *= tmp;
		};
		void changeByCos(double cosTheta, double phi)
		{
			double sinPhi = sin(phi);
			double cosPhi = cos(phi);
			double dxy = x*x + y*y;
			double dxyz = dxy + z*z;
			if (fabs(dxyz - 1.0) > 1e-14)
			{
				double norm = 1.0 / sqrt(dxyz);
				x *= norm;
				y *= norm;
				z *= norm;
				dxy = x*x + y*y;
			}

			if (dxy > 1.0e-28)
			{
				double sinTheta = sqrt((1 - cosTheta*cosTheta) / dxy);
				double up = x;
				x = x*cosTheta + sinTheta*(up*z*cosPhi - y*sinPhi);
				y = y*cosTheta + sinTheta*(y*z*cosPhi + up*sinPhi);
				z = z*cosTheta - dxy*sinTheta*cosPhi;
			}
			else
			{
				double sinTheta = sqrt(1 - cosTheta*cosTheta);
				if (z > 0)
				{
					x = sinTheta*cosPhi;
					y = sinTheta*sinPhi;
					z = cosTheta;
				}
				else
				{
					x = -sinTheta*cosPhi;
					y = -sinTheta*sinPhi;
					z = -cosTheta;
				}
			}
		}
		void rotate(double cos_t, double sin_t,
			double c_phi, double s_phi) {
			double sin_z = x*x + y*y;
			if (sin_z > 1e-10) {
				sin_z = sqrt(sin_z);
				double temp = sin_t / sin_z;
				double temp_phi = z*c_phi;
				double temp_x = x*cos_t;
				double temp_y = y*cos_t;
				double temp_x1 = temp_phi*x - y*s_phi;
				double temp_y1 = temp_phi*y + x*s_phi;
				x = temp*temp_x1 + temp_x;
				y = temp*temp_y1 + temp_y;
				z = z*cos_t - sin_z*sin_t*c_phi;
			}
			else {
				x = sin_t*c_phi;
				y = sin_t*s_phi;
				z *= cos_t;
			}
		};

	};
}


enum ParticleType { electron, photon, positron };

class Particle //store status information of a particle
{

public:
	__device__ __host__ void changeByCos(double cosTheta, double phi)
	{
		double sinPhi = sin(phi);
		double cosPhi = cos(phi);
		double dxy = u*u + v*v;
		double dxyz = dxy + w*w;
		if (fabs(dxyz - 1.0) > 1e-14)
		{
			double norm = 1.0 / sqrt(dxyz);
			u *= norm;
			v *= norm;
			w *= norm;
			dxy = u*u + v*v;
		}

		if (dxy > 1.0e-28)
		{
			double sinTheta = sqrt((1 - cosTheta*cosTheta) / dxy);
			double up = u;
			u = u*cosTheta + sinTheta*(up*w*cosPhi - v*sinPhi);
			v = v*cosTheta + sinTheta*(v*w*cosPhi + up*sinPhi);
			w = w*cosTheta - dxy*sinTheta*cosPhi;
		}
		else
		{
			double sinTheta = sqrt(1 - cosTheta*cosTheta);
			if (w > 0)
			{
				u = sinTheta*cosPhi;
				v = sinTheta*sinPhi;
				w = cosTheta;
			}
			else
			{
				u = -sinTheta*cosPhi;
				v = -sinTheta*sinPhi;
				w = -cosTheta;
			}
		}
	}
	__device__ __host__ void newDirection(double cosTheta, double phi, double& ux, double& vx, double& wx)
	{
		double cosPhi = cos(phi);
		double sinPhi = sqrt(1 - cosPhi*cosPhi);
		if (phi > 3.141592653589793238463) sinPhi = -sinPhi;
		double dxy = u*u + v*v;
		double dxyz = dxy + w*w;
		if (fabs(dxyz - 1.0) > 1e-14)
		{
			double norm = 1.0 / sqrt(dxyz);
			u *= norm;
			v *= norm;
			w *= norm;
			dxy = u*u + v*v;
		}

		if (dxy > 1.0e-28)
		{
			double sinTheta = sqrt((1 - cosTheta*cosTheta) / dxy);
			ux = u*cosTheta + sinTheta*(u*w*cosPhi - v*sinPhi);
			vx = v*cosTheta + sinTheta*(v*w*cosPhi + u*sinPhi);
			wx = w*cosTheta - dxy*sinTheta*cosPhi;
		}
		else
		{
			double sinTheta = sqrt(1 - cosTheta*cosTheta);
			if (w > 0)
			{
				ux = sinTheta*cosPhi;
				vx = sinTheta*sinPhi;
				wx = cosTheta;
			}
			else
			{
				ux = -sinTheta*cosPhi;
				vx = -sinTheta*sinPhi;
				wx = -cosTheta;
			}
		}
	}
	__device__ __host__ void normDirection()
	{
		double length2 = u*u + v*v + w*w;
		if (length2) 
		{
			double invLength = 1 / sqrt(length2);
			u *= invLength; v *= invLength; w *= invLength;
		}
	}
	//data
	double x, y, z; //position, unit cm
	double u, v, w; //direction vector
	int ivx, ivy, ivz, iabsv; //voxel index for current particle

	ParticleType type;

	double E; //energy, unit eV
	double weight;
};


struct VIEWRAY_FORMAT //load ViewRay's density file
{
	VIEWRAY_FORMAT() :Marjor_Header(1), Minor_Header(0){}
	bool read(const char* fname)
	{
		FILE* fp = fopen(fname, "rb");
		if (NULL == fp) return false;
		filename = fname;
		fread(&Marjor_Header, sizeof(uint32_t), 1, fp); //Marjor_Header
		fread(&Minor_Header, sizeof(uint32_t), 1, fp); //Minor_Header
		fread(&nx, sizeof(uint32_t), 1, fp);
		fread(&dx, sizeof(double), 1, fp);
		fread(&offset_x, sizeof(double), 1, fp);
		fread(&ny, sizeof(uint32_t), 1, fp);
		fread(&dy, sizeof(double), 1, fp);
		fread(&offset_y, sizeof(double), 1, fp);
		fread(&nz, sizeof(uint32_t), 1, fp);
		fread(&dz, sizeof(double), 1, fp);
		fread(&offset_z, sizeof(double), 1, fp);
		m.resize(nx, ny, nz);
		float* data = new float[nx*ny*nz];
		fread(data, sizeof(float), nx*ny*nz, fp);
		SFloat initVal = data[0];
		uniform = true;
		for (uint32_t ix = 0; ix < nx; ++ix)
			for (uint32_t iy = 0; iy < ny; ++iy)
				for (uint32_t iz = 0; iz < nz; ++iz)
				{
					m.a(ix, iy, iz) = SFloat(data[ix + nx*iy + nx*ny*iz]);
					if (initVal != m.a(ix, iy, iz)) uniform = false;
				}

		if (fread(data, sizeof(float), nx*ny*nz, fp))
		{
			err.resize(nx, ny, nz);
			for (uint32_t ix = 0; ix < nx; ++ix)
				for (uint32_t iy = 0; iy < ny; ++iy)
					for (uint32_t iz = 0; iz < nz; ++iz)
					{
						err.a(ix, iy, iz) = SFloat(data[ix + nx*iy + nx*ny*iz] / m.a(ix, iy, iz));
					}
		}
		delete[] data;
		fclose(fp);
		return true;
	}
	bool write(const char* fname = NULL)
	{
		FILE* fp = NULL;
		if (NULL == fname) fp = fopen(filename.c_str(), "wb");
		else fp = fopen(fname, "wb");
		if (NULL == fp) return false;
		fwrite(&Marjor_Header, sizeof(uint32_t), 1, fp); //Marjor_Header
		fwrite(&Minor_Header, sizeof(uint32_t), 1, fp); //Minor_Header
		fwrite(&nx, sizeof(uint32_t), 1, fp);
		fwrite(&dx, sizeof(double), 1, fp);
		fwrite(&offset_x, sizeof(double), 1, fp);
		fwrite(&ny, sizeof(uint32_t), 1, fp);
		fwrite(&dy, sizeof(double), 1, fp);
		fwrite(&offset_y, sizeof(double), 1, fp);
		fwrite(&nz, sizeof(uint32_t), 1, fp);
		fwrite(&dz, sizeof(double), 1, fp);
		fwrite(&offset_z, sizeof(double), 1, fp);

		float* data = new float[nx*ny*nz];
		for (uint32_t ix = 0; ix < nx; ++ix)
			for (uint32_t iy = 0; iy < ny; ++iy)
				for (uint32_t iz = 0; iz < nz; ++iz)
					data[ix + nx*iy + nx*ny*iz] = m.a(ix, iy, iz);
		fwrite(data, sizeof(float), nx*ny*nz, fp);
		if (!err.empty()) //may output the err
		{
			for (uint32_t ix = 0; ix < nx; ++ix)
				for (uint32_t iy = 0; iy < ny; ++iy)
					for (uint32_t iz = 0; iz < nz; ++iz)
						data[ix + nx*iy + nx*ny*iz] = err.a(ix, iy, iz)*m.a(ix, iy, iz);
			fwrite(data, sizeof(float), nx*ny*nz, fp);
		}
		delete[] data;
		fclose(fp);
		return true;
	}

	void release(){ m.release(); err.release(); } //release memory
	SFloat interp(double x, double y, double z) //linear 3D interpolation
	{
		//x,y,z is the coordinates relative to upper left voxel center
		if (x + dx / 2 < 0 || x + dx / 2 > nx*dx || y + dy / 2 < 0 || y + dy / 2 > ny*dy || z + dz / 2 < 0 || z + dz / 2 > nz*dz)
		{
			return 0; //set to zero for points beyond the boundary
		}
		int ix = int(floor(x / dx));
		int iy = int(floor(y / dy));
		int iz = int(floor(z / dz));
		x = x / dx - ix; //decimal part
		y = y / dy - iy;
		z = z / dz - iz;
		//fix the boundary access problem
		if (ix < 0) { ix = 0; x = 0; }
		if (ix >= int(nx - 1)) { ix = nx - 2; x = 1; }
		if (iy < 0) { iy = 0; y = 0; }
		if (iy >= int(ny - 1)) { iy = ny - 2; y = 1; }
		if (iz < 0) { iz = 0; z = 0; }
		if (iz >= int(nz - 1)) { iz = nz - 2; z = 1; }
		//do the interpolation
		double y1 = 0, y2 = 0;
		y1 = m.a(ix, iy, iz)*(1 - x) + m.a(ix + 1, iy, iz)*x;
		y2 = m.a(ix, iy + 1, iz)*(1 - x) + m.a(ix + 1, iy + 1, iz)*x;
		double z1 = y1*(1 - y) + y2*y;

		y1 = m.a(ix, iy, iz + 1)*(1 - x) + m.a(ix + 1, iy, iz + 1)*x;
		y2 = m.a(ix, iy + 1, iz + 1)*(1 - x) + m.a(ix + 1, iy + 1, iz + 1)*x;
		double z2 = y1*(1 - y) + y2*y;

		return SFloat(z1*(1 - z) + z2*z);
	}

	ArrayMgr<SFloat> m; //main data matrix
	ArrayMgr<SFloat> err; //err, may be empty
	uint32_t nx, ny, nz;
	uint32_t Marjor_Header, Minor_Header;
	double offset_x, offset_y, offset_z;
	double dx, dy, dz;

	//auxiliary variable
	bool uniform;//uniform is set in read().
	string filename; //remember which file we just opened.
};


const double BinaryFile_MagicNum = 3.1415926;
const double BinaryFile_CurrentVerison = 1.1;
const int32_t HEADER_SECTION_LENGTH = 128 * sizeof(double); // the real header may have smaller length

enum DataType{ INT32_T, DOUBLE_T, FLOAT_T, STRING_T, VOID_T };
struct Chunk
{
	Chunk() :p(NULL), len(0){}
	Chunk(const Chunk& r)
	{		
		len = r.len;
		p = new char[len];
		memcpy(p, r.p, len);

		type = r.type;
		val = r.val;
	}
	Chunk& operator = (const Chunk& r)
	{
		if (len != r.len)
		{
			if (p) delete[] p; // delete old memory
			len = r.len;
			p = new char[len]; // allocate new memory
		}

		memcpy(p, r.p, len); // copy the content

		type = r.type;
		val = r.val;
		return *this;
	}
	~Chunk()
	{
		if (p) delete[] p;
	}
	void set(void* rp, int rlen)
	{
		if (len != rlen)
		{
			if (p) delete[] p; // delete old memory
			len = rlen;
			p = new char[len]; // allocate new memory
		}
		memcpy(p, rp, len); // copy the content
	}
	int32_t len;
	char* p;

	//the following is used in DoseViewer
	string type; //printed type
	string val; //printed value
};

class DECLSPECIFIER BinaryFile
{
public:
	
	struct FileHeader // The file header type may change in the future. All added variable should follow behind to ensure compatibility
	{
		double MagicNum; //as an indicator of whether this file follow the format I define
		int8_t b_compress; //whether the data is compressed
		int32_t compressMethod; // LZ4 =0 , QuickLZ = 1;
		double version;

		bool read(FILE* fp)
		{
			if (fileSize(fp) < HEADER_SECTION_LENGTH) return false;
			int size_header = sizeof(FileHeader);
			if (size_header > HEADER_SECTION_LENGTH) return false;
			fread(this, size_header, 1, fp);
			fseek(fp, HEADER_SECTION_LENGTH, SEEK_SET);
			if (MagicNum != BinaryFile_MagicNum) return false;
			if (BinaryFile_CurrentVerison == version) return true;
			else printf("Warning: expect file version %.2f, but the actual version is %.2f\n", BinaryFile_CurrentVerison, version);
			return true;
		}
		void write(FILE* fp)
		{
			fwrite(this, sizeof(FileHeader), 1, fp); //write the header
			char chz = 0; //complement the rest to HEADER_LENGTH
			for (int i = 0; i < HEADER_SECTION_LENGTH - sizeof(FileHeader); ++i) fwrite(&chz, 1, 1, fp);
		}
	};
	~BinaryFile();
	bool read(const char* fname, bool onlyTest = false);
	double getVersion();
	bool write(const char* fname, bool compress = true, int32_t compressMethod = 0);
	bool write(FILE *fp);
	bool writeChunk(const char* fname, string& key);
	bool get(const string& key, void* rp);
	bool getPointer(const string& key, void**p, int32_t &len);
	int  getKeyLength(const string& key); //get key length in byte
	bool push(const string& key, void* rp, int32_t len);
	bool exist(const string& key);
	int  size();
	void beginSearch();
	bool next(string& name, Chunk& chk);
	bool changeKeyName(string& oldKey, string& newKey);
	bool erase(const string& key);
	bool addChunk(string& key, Chunk& chk);
	Chunk& getChunk(const string& key);
	int  getTotBytes();
	void delAll();
	map<string, Chunk>& getBM();
private:
	map<string, Chunk> bm; //string must have <= 63 characters
	map<string, Chunk>::iterator ic;
	FileHeader header;
};

class DECLSPECIFIER PENELOPE_MAT
{
	struct MAT
	{
		string name;
		double density;
	};
public:
	bool load()
	{
		matlist.resize(0);
		string line;
		ifstream myfile("PENELOPE_densities.txt");
		if (myfile.is_open())
		{
			MAT m;
			while (getline(myfile, line) && !line.empty())
			{
				size_t p1 = line.find_first_of('|');
				size_t p2 = line.find_last_of('|');
				m.name = line.substr(p1 + 1, p2 - p1 - 1);
				string strden = line.substr(p2 + 1, string::npos);
				m.density = std::stod(strden);
				matlist.push_back(m);
			}
			myfile.close();
			return true;
		}
		else return false;
	}
	double density(int id) // note id start from 1
	{
		return matlist[id - 1].density;
	}
	string name(int id) // note id start from 1
	{
		return matlist[id - 1].name;
	}
private:
	vector<MAT> matlist;
};

//several callback function definitions
typedef void(*ProgressCallBack)(double /*percent:%*/, double /*speed;h/s*/, double /*time used;s*/, double /*time left;s*/, double /*uncertainty;%*/);
typedef void(*PeekDoseCallBack)(BinaryFile&);
typedef void(*JobFinishedCallBack)(bool /*0 normal, 1 abort*/, BinaryFile&);
#endif
