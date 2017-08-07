#ifndef _QUALITY_DBG_H_
#include "types.h"
#include "funcs.h"
#define  EXTERN  extern
#include "vars.h"
#include "transport.h"

#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)
#define VAR_NAME_VALUE(var) #var "="  VALUE(var)

#pragma message("QDBG_EXT=" VAR_NAME_VALUE(QDBG_EXT))



//#define DBG_PRINT_QUAL
#define FLOW_DBG_PRINT
#define DBG_PRINT_DETAIL

QDBG_EXT FILE *fpMassInfo;
QDBG_EXT FILE *fpLV;
QDBG_EXT double extraMass;
QDBG_EXT double lostMass;
QDBG_EXT double extraVol;
QDBG_EXT double cLostMass;
QDBG_EXT double cExtraMass;
QDBG_EXT double cExtraVol;

void dbgPrintDetail(long stime, int loc, long dt, double *x);
void dbgPrintLinkQual(long stime, char *type);
void dbgPrintTankQual(long stime, char *type);
int LinkHasNZConcentration(Pseg seg);
void dbgPrintNodeQual(long stime, char *type);
void writeLinkVolumes(FILE *fplv, long stime, char *id, long dt);

#endif