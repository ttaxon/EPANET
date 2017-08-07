#ifndef _FT_DBG_PRINT_H_
#define _FT_DBG_PRINT_H_
#include <stdio.h>

#include "transport.h"

//#define DEBUG

#ifdef ENABLE_DBG_PRINT
static long CurTime;
void ftPrintInfoHeader(FILE *fpi,int maxIn, int maxOut);
void ftPrintNodeInfo(FILE *fpi, PNodeTransport n, int maxIn, int maxOut);
void ftPrintDemandHeader(FILE *fpi);
void ftPrintFlowLinksHeader(FILE *fpi,int maxOut);
void ftPrintMoveDSHeaders(FILE *fpi);
void ftPrintInflowHeader(FILE *fpi,int maxIn);
void ftPrintSegmentsToMoveHeader(FILE *fpi,int maxIn);
void ftPrintSegmentsToMove(FILE *fpi,PSegList segs);
void printSegStats(long stime,char *str);
void ftPrintInflowOutflow(long stime);
void ftPrintWarningsFile(long stime);
void ftPrintIncompleteSegments(PIncompleteSeg *s, int n, FILE *fp);
void ftPrintIncompleteSegment(PIncompleteSeg s, int level, FILE *fp, char *indent);
void ftCloseDebugFiles();
void printMassBalance(FILE *fp, long step);
void printmassbalanceDEBUG(long t, long dt);
void ftPrintDetail(long t, char *id, PSegList segsToMove, int sequence, int loc);
void ftPrintAllIncompleteSegments(PNodeTransport n);
void ftPrintAllIncSegs(long stime);
void ftCheckForIncompleteSegments(PNodeTransport n, long stime);
void ftPrintSource(FILE *fp, long long segSrc);

char *getSegmentType(Pseg s);
char *ftGetOwnerObject(PSegOwner o);
char *ftGetOwnerData(PSegOwner o);
void ftOpenMBFiles();
void ftdAccumulateMassAdded(int i);
void ftPrintFlow(long stime, long dt);
void ftOpenDebugFiles(long stime);
void ftOpenInfoFile(long stime, long dt);
void ftPrintAllConcentrations();
void ftPrintMassBalanceData(long stime);
void ftPrintMBArrays(long stime);

void ftPrintInflowOutflowData();


void ftDocPrintNode(long stime, SNodeTransport *node, char *desc);
void ftDocPrintSegments(long stime, SNodeTransport *node, PSegList segs, char *desc);
void ftDocPrintLinks(long stime, SNodeTransport *node, char *desc);
void ftDocPrintLinkVolumes(long stime, SNodeTransport *node, char *desc);
void ftDocPrintConcentrations(long stime, char *desc);

extern FILE *fpDeltas;
extern FILE *fpAdj;
extern FILE *fpIncTrace;
extern FILE *demFP;
extern FILE *qualFP;
extern FILE *massFP;
extern FILE *fpAllAvgConc;
extern FILE *fpAllInstConc;

extern FILE *fpMassBalance;
extern FILE *fpcombine;
extern FILE *incfp;
extern FILE *fpTrace;
extern FILE *fpi;
extern FILE *warnfp;
extern FILE *fpOwners;
extern FILE *fpCompleteInc;
extern FILE *fpUpdateIncSeg;
extern FILE *fpLoss;
extern FILE *fpTank;
extern FILE *fpd;

extern int maxIn;
extern int maxOut;

extern int printTrace;

extern int printAllLinkSegs;
extern int openMassBalanceData;
extern int *nodesToWrite;
extern int *linksToWrite;
#endif

#endif
