//
// MATLAB Compiler: 5.2 (R2014b)
// Date: Fri Aug 28 00:21:25 2015
// Arguments: "-B" "macro_default" "-v" "-W" "cpplib:MatlabTools" "-T"
// "link:lib" "wdicom.m" "writeDose2D.m" 
//

#ifndef __MatlabTools_h
#define __MatlabTools_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#include "mclcppclass.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_MatlabTools
#define PUBLIC_MatlabTools_C_API __global
#else
#define PUBLIC_MatlabTools_C_API /* No import statement needed. */
#endif

#define LIB_MatlabTools_C_API PUBLIC_MatlabTools_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_MatlabTools
#define PUBLIC_MatlabTools_C_API __declspec(dllexport)
#else
#define PUBLIC_MatlabTools_C_API __declspec(dllimport)
#endif

#define LIB_MatlabTools_C_API PUBLIC_MatlabTools_C_API


#else

#define LIB_MatlabTools_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_MatlabTools_C_API 
#define LIB_MatlabTools_C_API /* No special import/export declaration */
#endif

extern LIB_MatlabTools_C_API 
bool MW_CALL_CONV MatlabToolsInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_MatlabTools_C_API 
bool MW_CALL_CONV MatlabToolsInitialize(void);

extern LIB_MatlabTools_C_API 
void MW_CALL_CONV MatlabToolsTerminate(void);



extern LIB_MatlabTools_C_API 
void MW_CALL_CONV MatlabToolsPrintStackTrace(void);

extern LIB_MatlabTools_C_API 
bool MW_CALL_CONV mlxWdicom(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

extern LIB_MatlabTools_C_API 
bool MW_CALL_CONV mlxWriteDose2D(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);


#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__BORLANDC__)

#ifdef EXPORTING_MatlabTools
#define PUBLIC_MatlabTools_CPP_API __declspec(dllexport)
#else
#define PUBLIC_MatlabTools_CPP_API __declspec(dllimport)
#endif

#define LIB_MatlabTools_CPP_API PUBLIC_MatlabTools_CPP_API

#else

#if !defined(LIB_MatlabTools_CPP_API)
#if defined(LIB_MatlabTools_C_API)
#define LIB_MatlabTools_CPP_API LIB_MatlabTools_C_API
#else
#define LIB_MatlabTools_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_MatlabTools_CPP_API void MW_CALL_CONV wdicom(const mwArray& info);

extern LIB_MatlabTools_CPP_API void MW_CALL_CONV writeDose2D(const mwArray& info);

#endif
#endif
