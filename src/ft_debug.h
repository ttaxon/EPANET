#ifndef _FT_DEBUG_H_
#define _FT_DEBUG_H_

// these defines must appear before the include of transport.h
// the header file uses them, but I feel they shouldn't be put
// in the header.
// ENABLE_DBG_PRINT and DBG_PRINT control the inclusion ofthe logging code
//#define ENABLE_DBG_PRINT
#ifdef ENABLE_DBG_PRINT
extern int *dbgPrintNodes;
extern int *dbgPrintLinks;

// DBG_PRINT macro is to reduce the amount of #ifdef/#endif pairs in an attempt to keep the code "cleaner"
// where there are only one or two lines that need to be excluded.
#define DBG_PRINT(code) code
#else
#define DBG_PRINT(code)
#endif

#endif
