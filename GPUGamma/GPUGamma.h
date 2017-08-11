#pragma once
#ifndef _GPUGAMMA_H_
#define _GPUGAMMA_H_

#include "../Tools/Tools.h"

#ifdef WIN32
#define DEXPORT __declspec(dllexport)
#else
#define DEXPORT __attribute__ ((visibility ("default")))
#endif

DEXPORT int GPUAvailable();

DEXPORT int GPUGammaAnalysis(ArrayMgr<SFloat>& d1, ArrayMgr<SFloat>& d2, ArrayMgr<SFloat>& g, double dx, double dy, double dz, double percent, double DTA, bool locale = false, int NInterp = 1);

#endif
