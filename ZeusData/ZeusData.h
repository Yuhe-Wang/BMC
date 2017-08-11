#pragma once
#ifndef _ZEUSDATA_H_
#define _ZEUSDATA_H_
#include <vector>
using std::vector;

#ifdef WIN32
#define DEXPORT __declspec(dllexport)
#else
#define DEXPORT __attribute__ ((visibility ("default")))
#endif

//dataFolder is the folder where the cross section files is stored
extern "C" DEXPORT bool ZeusData_load(const char* dataFolder);

//don't forget to call this after finishing simulation to avoid memory leakage
extern "C" DEXPORT void ZeusData_unload();

extern "C" DEXPORT void ZeusData_lambdaWoodcock(int aNumElements, const float *dens, const char *mat, int nvoxels, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b);

extern "C" DEXPORT void ZeusData_BremConst(int aNumElements, double f0[3], double& bcb);

extern "C" DEXPORT void ZeusData_lambdaPhoton(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d);

extern "C" DEXPORT void ZeusData_lambdaCo(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d);

extern "C" DEXPORT void ZeusData_lmbdapp(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d);

extern "C" DEXPORT void ZeusData_ScreeningParameter(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d);

extern "C" DEXPORT void ZeusData_QSurface(int matid, double& fq, double& ieq0, vector<vector<double>>& qs);

extern "C" DEXPORT void ZeusData_Range(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b);

extern "C" DEXPORT void ZeusData_InverseRange(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b);

#endif