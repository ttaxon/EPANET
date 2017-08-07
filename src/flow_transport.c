#include <stdlib.h>
#include <stdio.h>
#include <io.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <float.h>

#ifdef _WIN32
#define isnan(n) _isnan(n)
#define isinf(n) (!_finite(n))
#endif
#include "ft_dbg_print.h"
#include "transport.h"
#include "ft_checks.h"

int traceProcessIncomplete;

// enable logging (to a file) of competing incomplete segments
//#define LOG_COMPLETE_INC

// ENABLE_TIMING controls whether or not the timing code is added and an output file is created.
//#define ENABLE_TIMING
#include <time.h>
#ifdef ENABLE_TIMING
#define TIMING(code) code
char *sects[7]={"Init","Loop","WarnFile","Incomplete completion","Finalization","PrintLinkSegs","CloseFiles"};
#else
#define TIMING(code)
#endif

double FT_QZERO=1e-6;
// set FT_RTOL to 0 because sometimes (especially with long WQ time steps),
// segments were being ignored when they shouldn't have been.  This seems
// to be working for now.
//double FT_RTOL=1e-4;
double FT_RTOL=0;
// FT_ATOL_BASE is the basis for computing FT_ATOL based on the dt that is being
// used for any given step
double FT_ATOL_BASE=1e-6;
double FT_ATOL=0;
// these two are used to aid in determining if two segments can be combined
// for more detail, see the comment in ftCombineSegemnts
double FT_CTOL=1e-2;
double FT_VTOL=1.0;


char *OwnerTypes[]={"NONE","THROUGH","INFLOW","DEMAND","LINK","MERGED","INCOMPLETE"};
/*
 * since the SegOwner references in the incomplete and merged segment types
 * are not pointers, these objects are used when NULL needs to be assigned to
 * one of them.
*/
SegOwnerObject NULLOwnerObject = {0};
SegOwnerData   NULLOwnerData = {0};

double volAdded;
double volLost;
double massLost;

SNodeTransport *TransportNodes; /* Flow transport data structures */

/* Function prototypes for methonds used only in this file */
void ftFinalizeMoveSegs();
void ftInitMoveSegs(long stime, long dt, SNodeTransportList**nodesToProcess);

void ftAddLink(SNodeLinks *links, int linkIdx, double flow);
void ftComputeRatios(SNodeLinks *links);
Pseg ftMergeSegments(PNodeSegs volSegs, PNodeSegs nodeSegs, PNodeLinks inflow);
int ftFlowComplete(PNodeLinks flow);
void ftAddThroughFlow(Pseg seg, PNodeTransport tNode,Pseg throughFlow);
void ftFreeIncompleteSegments(PIncompleteSegmentData incomplete, int nlinks, PIncompleteSegmentData *head);
double ftMassAdded(Psource source, long dt, SNodeTransport *node);
void ftAddNodeToProcess(SNodeTransport *node,SNodeTransportList **nodesToProcess);
void ftAddNodeSeg(SegList *segs, Pseg seg, int loc);
void ftGetSegmentsToMove(SegList * segs,SNodeTransport *node);
void ftGetSegmentsToMoveOld(SegList * segs,SNodeTransport *node);
int ftGetLinkIndex(SNodeLinks *inflow,int linkIdx);
void ftRemoveFirstLinkSegment(int linkIdx);
void ftReleaseSeg(Pseg seg);
void ftRemoveFirstSegment(SegList *segList);
Pseg ftUnlinkFirstSegment(PSegList segList);
SNodeTransport *ftGetNextNodeToProcess(SNodeTransportList **nodesToProcess);
Pseg ftRemoveDemand(Pseg seg,PNodeTransport tNode, Pseg throughFlow);
void ftMoveToDownstreamNodes(PNodeTransport tNode,PNodeTransportList *nodesToProcess);
void ftAdjustOutflow(Pseg seg, double outflowAdj);
void ftAdjustIncompleteSegment(PIncompleteSeg incSeg,double adj);
void ftAdjustMergedSeg(PMergedSeg seg,double adj);
void ftMoveFlowToLinks(Pseg seg,PNodeTransport tNode,Pseg throughFlow);
void ftSplitIncompleteSegment(PIncompleteSeg incSeg,int nparts,SegOwner owners[],VolSplitData volSplit[]);
double ftProcessIncomplete(PIncompleteSeg incSeg,Pseg *seg, int idx, int applyAdjustment, int completing);
void ftSplitMergedSegment(PMergedSeg mergedSeg,int nparts,PSegOwner types,VolSplitData volSplit[],PMergedSeg *segs);
void ftUpdateMergedSegment(PMergedSeg ms,PIncompleteSeg incSeg,double mass, int completing);
void ftRemoveSegFromList(Pseg seg,Pseg *list);
void ftFreeAllSegs(PSegList segList);
void ftCombineSegs(Pseg *firstSeg, Pseg *lastSeg,int loc);
void ftReplacePseg(Pseg oldSeg,Pseg newSeg, Pseg *first, Pseg *last);
void ftReplaceLinkSeg(Pseg oldSeg,Pseg newSeg, int linkIdx,int loc);
void ftReplaceSeg(Pseg oldSeg,Pseg newSeg,SegList *segs, int loc);
void ftFreeIncompleteSegment(PIncompleteSeg incSeg);
int ftIsComplete(PIncompleteSeg incSeg);
void ftAddMass(double segMass,PIncompleteSeg is);
void ftReleaseSegList(PSegList segs);
int ftHaveAllInflowSegs(SNodeSegs *segs);
int ftIncompleteDone(PIncompleteSegmentData incomplete,int nlinks);
double ftSumMass(PSegList segs);
void ftMoveSegList(PSegList from, PSegList to, int loc);
Pseg ftUnlinkFirstLinkSegment(int linkIdx);
void ftAddTankOutflow(PSegList segs,int tankIdx,double vol,int dt);
void ftTankOutflow1(int tankIdx, double vol, PSegList systemInflow, int dt);
void ftTankOutflow2(int tankIdx, double vol, PSegList systemInflow, int dt);
void ftTankOutflow3(int tankIdx, double vol, PSegList systemInflow, int dt);
void ftTankOutflow4(int tankIdx, double vol, PSegList systemInflow, int dt);
void ftAddTankInflow(PSegList segs,int tankIdx);
void ftTankInflow1(int tankIdx,PSegList inflow);
void ftTankInflow2(int tankIdx,PSegList inflow);
void ftTankInflow3(int tankIdx,PSegList inflow);
void ftTankInflow4(int tankIdx,PSegList inflow);
int ftUpdatePartSegment(PIncSegPart incpart,PIncompleteSeg incpartOwner,Pseg incSeg,double mass, int completing);
int ftUpdateIncompleteSegs(PNodeTransport node,int inflowIdx,double c,double *volMoved, int completing);
void ftAddIncompleteSegment(PIncompleteSegmentData incSeg,PIncompleteSegmentData *list);
int ftValidateMergedSegment(PMergedSeg ms,char *src);
void ftPrintAllLinkSegs(long stime,char *mode);
int ftCombineSegments(double v1, double c1, double v2, double c2,int loc);
void ftFreeNodeLinks(PNodeLinks nl);
void ftResetFlows(PNodeLinks f);
void ftCompleteIncompleteSegments(PNodeTransport n, FILE *fp);

void ftCheckForInvalidSegments();

void ftCompleteAllIncompleteSegments(long stime);

void ftAdjustLinkSegmentVolumes();
void ftComputeMassBalance();
void checkConc(double c);

static long SegSequence=0;
long ftGetSegSequence() {
	long rv=SegSequence;
    if(rv==3844) {
        printf("");
    }
	SegSequence++;
	return rv;
}

#define FTInitNoneOwnerType(owner,sourceID) \
	(owner).type=SO_NONE; \
	(owner).object=NULLOwnerObject; \
	(owner).data=NULLOwnerData; \
	(owner).source=sourceID;

#define FTInitThroughOwnerType(owner,ownerObject,sourceID) \
	(owner).type=SO_THROUGH; \
	(owner).object.node=ownerObject; \
	(owner).data=NULLOwnerData; \
	(owner).source=sourceID; \
	(owner).sequence=ftGetSegSequence();

#define FTInitInflowPartOwner(owner,ownerObject,ownerData,sourceID) \
	(owner).type=SO_INFLOW; \
	(owner).object.node=ownerObject; \
	(owner).data.inflowIdx=ownerData; \
	(owner).source=sourceID; \
	(owner).sequence=ftGetSegSequence();

#define FTInitDemandPartOwner(owner,ownerObject, sourceID) \
	(owner).type=SO_DEMAND; \
	(owner).object.node=ownerObject; \
	(owner).data=NULLOwnerData; \
	(owner).source=sourceID; \
	(owner).sequence=ftGetSegSequence();

#define FTInitLinkPartOwner(owner,ownerObject,sourceID) \
	(owner).type=SO_LINK; \
	(owner).object.linkIdx=ownerObject; \
	(owner).data=NULLOwnerData; \
	(owner).source=sourceID; \
	(owner).sequence=ftGetSegSequence();

#define FTInitMergedPartOwner(owner,ownerObject,sourceID) \
	(owner).type=SO_MERGED; \
	(owner).object.mergedSeg=ownerObject; \
	(owner).data=NULLOwnerData; \
	(owner).source=(owner).source*100+sourceID; \
	(owner).sequence=ftGetSegSequence();

#define FTInitIncompletePartOwner(owner,ownerObject,ownerData,sourceID) \
	(owner).type=SO_INCOMPLETE_PART; \
	(owner).object.incSeg=ownerObject; \
	(owner).data.incSeg=ownerData; \
	(owner).source=sourceID; \
	(owner).sequence=ftGetSegSequence();

#define FTInitCopyOwner(owner,originalOwner,newSourceID) \
	memcpy(&(owner), originalOwner, sizeof(SegOwner)); \
	(owner).source=(owner).source*100+newSourceID; \
	(owner).sequence=ftGetSegSequence();


/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: Initializes data structures for a new WQ step.
**            called from initqual()
**--------------------------------------------------------------
*/
void ftInitQual() {
	int i;
	for (i=1; i<=Nnodes; i++) {
		PNodeTransport node=&TransportNodes[i];
		node->demand=0;
		node->massadded=0;
		node->massused=0;
		node->totalVolInflow=0;
		node->totalVolOutflow=0;
		MB_MassAdded[i]=0;
		MB_MassMoved[i]=0;
		MB_MassRemoved[i]=0;
		MB_MassNoOutflow[i]=0;
	}
	for (i=1; i<=Ntanks; i++) {
		MB_MassInTanks[i]=0;
		MB_MassNegTankVol[i]=0;
	}
	for (i=1; i<=Nlinks; i++) {
		MB_MassInPipes[i]=0;
	}
	DBG_PRINT(ftPrintAllConcentrations(););
}
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: Initializes data structures for the FlowTransport method.
**            called from openqual()
**--------------------------------------------------------------
*/
void ftOpenMoveSegs() {
	int i;
	int openMBFiles=0;
	int errcode = 0;
	if (transportMethod==TM_FLOW) {
		AvgNodeQual = (double *)calloc(Nnodes + 1, sizeof(double));
		Mass        = (double *)calloc(Nnodes+1, sizeof(double));
		ERRCODE(MEMCHECK(AvgNodeQual));
		ERRCODE(MEMCHECK(Mass));
	}
	MB_MassAdded = (double *)calloc(Nnodes+1, sizeof(double));
	MB_MassMoved = (double *)calloc(Nnodes+1, sizeof(double));
	MB_MassRemoved = (double *)calloc(Nnodes+1, sizeof(double));
	MB_MassNoOutflow = (double *)calloc(Nnodes+1, sizeof(double));
	MB_MassInTanks = (double *)calloc(Ntanks+1, sizeof(double));
	MB_MassNegTankVol = (double *)calloc(Ntanks+1, sizeof(double));
	MB_MassInPipes = (double *)calloc(Nlinks+1, sizeof(double));
	ERRCODE(MEMCHECK(MB_MassAdded));
	ERRCODE(MEMCHECK(MB_MassMoved));
	ERRCODE(MEMCHECK(MB_MassRemoved));
	ERRCODE(MEMCHECK(MB_MassInTanks));
	ERRCODE(MEMCHECK(MB_MassInPipes));
	ERRCODE(MEMCHECK(MB_MassNegTankVol));
	ERRCODE(MEMCHECK(MB_MassNoOutflow));
	TransportNodes=(SNodeTransport*)calloc(Nnodes+1, sizeof(SNodeTransport));
	for (i=1; i<=Nnodes; i++) {
		TransportNodes[i].node=&Node[i];
		TransportNodes[i].nodeIdx=i;
	}
#ifdef ENABLE_DBG_PRINT
	if (openMBFiles) {
		ftOpenMBFiles();
	}
#endif
}
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: Releases data structures for the FlowTransport method.
**            called from closequal()
**--------------------------------------------------------------
*/
void ftCloseMoveSegs() {
	int i;
	free(AvgNodeQual);
	free(Mass);
	for (i=1; i<=Nnodes; i++) {
		PNodeTransport node=&TransportNodes[i];
		ftFreeNodeLinks(&node->inflow);
		ftFreeNodeLinks(&node->outflow);
		free(node->inflowSegs.segs);
	}
	free(TransportNodes);
	free(MB_MassAdded);
	free(MB_MassMoved);
	free(MB_MassRemoved);
	free(MB_MassInTanks);
	free(MB_MassInPipes);
	free(MB_MassNegTankVol);
	free(MB_MassNoOutflow);
#ifdef ENABLE_DBG_PRINT
	fclose(fpAllAvgConc);
	fclose(fpAllInstConc);
#endif
}

/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: Performs necessary initialization when a new WQ run
**            is being started.
**--------------------------------------------------------------
*/
void ftInitSegs() {
	// a new iter is starting.
	ftNewHydStep();
}
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: A new set of hydraulic results was reached.  Perform
**            necessary updates
**--------------------------------------------------------------
*/
void ftNewHydStep() {
	int i;
	// a new hydraulic step has been reached. rebuild the inflow & outflow links
	for (i=1; i<=Nnodes; i++) {
		TransportNodes[i].outflow.nlinks=0;
		TransportNodes[i].inflow.nlinks=0;
	}
    // accumulate total flowrates (in & out) due to link flow
	for (i=1; i<=Nlinks; i++) {
		double q;
		int ds=DOWN_NODE(i);
		q=ABS(Q[i]);             /* Flow rate */
		if (q >= FT_QZERO) {
			int us=UP_NODE(i);
            // add the flowrate to the outflow links for the upstream node
			ftAddLink(&TransportNodes[UP_NODE(i)].outflow, i, q);
            // add the flowrate to the inflow links for the downstream node
			ftAddLink(&TransportNodes[ds].inflow, i, q);
		}
	}
	for (i=1; i<=Nnodes; i++) {
		PNodeTransport n=&TransportNodes[i];
		// find all nodes where water is entering the system
		if (NodeDemand[i] < 0.0) {
            // and add the flowrate to the inflow
			ftAddLink(&n->inflow, -i, ABS(NodeDemand[i]));
		} else {
			double inflowRate=0;
			int j;
			for (j=0; j<n->inflow.nlinks; j++) { inflowRate+=n->inflow.flowRate[j]; }
            // if the inflow rate is zero, still need to check the outflow.  there have been cases
            // where the outflowRate+demandRate are > 0 even though there is no inflow
			if (inflowRate == 0) {
				double outflowRate=0;
				double demandRate=NodeDemand[i];
				for (j=0; j<n->outflow.nlinks; j++) { outflowRate+=n->outflow.flowRate[j]; }
				if (outflowRate + demandRate > 0) {
					// in order for the transport to function properly, add an inflow
					// segment equal to demand+outflow
					inflowRate=outflowRate+demandRate;
					ftAddLink(&n->inflow, 0, inflowRate);
				}
			}
		}
	}
	for (i=1; i<=Nnodes; i++) {
		SNodeTransport *tn=&TransportNodes[i];
		tn->inflowSegs.nlinks=tn->inflow.nlinks;
        // increase the size of the data structures if necessary.  Don't want to always allocate/deallocate.
		if (tn->inflowSegs.size < tn->inflowSegs.nlinks) {
			tn->inflowSegs.segs=(SegList*)realloc(tn->inflowSegs.segs, tn->inflowSegs.nlinks*sizeof(SegList));
			memset(tn->inflowSegs.segs, 0, tn->inflowSegs.nlinks*sizeof(SegList));
			tn->inflowSegs.size = tn->inflowSegs.nlinks;
		}
        // compute inflow and outflow ratios
		ftComputeRatios(&tn->outflow);
		ftComputeRatios(&tn->inflow);
	}
}

/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: Ensure that there are no demand segs or through
**            segs left from a previous iteration
**--------------------------------------------------------------
*/
void ftInitializeTransport() {
	int i;
	int anyLeft=0;
	for (i=1; i<=Nnodes; i++) {
		PNodeTransport n=&TransportNodes[i];
		anyLeft=n->demandSegs.firstSeg != NULL || n->throughSegs.firstSeg != NULL;
		ftReleaseSegList(&n->demandSegs);
		ftReleaseSegList(&n->throughSegs);
		MB_MassAdded[i]=0;
		MB_MassMoved[i]=0;
		MB_MassRemoved[i]=0;
		MB_MassNoOutflow[i]=0;
	}
	for (i=1; i<=Ntanks; i++) {
		MB_MassInTanks[i]=0;
		MB_MassNegTankVol[i]=0;
	}
	for (i=1; i<=Nlinks; i++) {
		MB_MassInPipes[i]=0;
	}
	if (anyLeft) {
		fprintf(stdout, "Are there any demandSegs or throughSegs left here?\n");  fflush(stdout);
		fprintf(stdout, "%s\n", anyLeft?"yes":"no");  fflush(stdout);
	}
}
/*
**--------------------------------------------------------------
**   Input:   stime: current simulation time
**            dt: current time step
**            nodesToProcess: the list of nodes to process
**   Output:  none
**   Purpose: Perform necessary initialization for transport for
**            one time step
**--------------------------------------------------------------
*/
void ftInitMoveSegs(long stime, long dt, SNodeTransportList**nodesToProcess) {
	int i;

	DBG_PRINT(ftOpenDebugFiles(stime););
	volAdded=0;
	volLost=0;
	massLost=0;
	/* set up SNodeTransport data structure */
	for (i=1; i<=Nnodes; i++) {
		PNodeTransport node=&TransportNodes[i];
		if (NodeDemand[i] > 0) {
			// negative demand (inflow) is handled elsewhere
			node->demand=NodeDemand[i]*dt;
		} else {
			node->demand=0;
		}
		node->totalVolInflow=0;
		node->totalVolOutflow=0;
		node->massadded=0;
		node->massused=0;
#ifdef ENABLE_DBG_PRINT
		node->visited=0;
#endif
		ftResetFlows(&node->inflow);
		ftResetFlows(&node->outflow);

		DBG_PRINT(if (fpd!=NULL) { fprintf(fpd, "%d\t%s\t%.10f\t%.10f\n", i, Node[i].ID, NodeDemand[i]>0?TransportNodes[i].demand:NodeDemand[i]*dt, NodeDemand[i]); fflush(fpd); });
	}
	DBG_PRINT(if (fpd!=NULL) fclose(fpd););

	/* Initialize inflow/outflow link data */
	for (i=1; i<=Nlinks; i++) {
		double origv, v;
		int ds=DOWN_NODE(i), us=UP_NODE(i), l;
		origv=ABS(Q[i])*dt;             /* Flow volume */
		if (ABS(Q[i]) < FT_QZERO) {
			v=0;
		} else {
			v=origv;
		}
		// accumulate the inflow volume on the ds nodes and outflow volume on the us nodes
		// and the inflow volume on each inflow and outflow link
		if (v > 0) {
			TransportNodes[ds].totalVolInflow += v;
			l=ftGetLinkIndex(&TransportNodes[ds].inflow, i);
			if (l>=0) {
				if (TransportNodes[ds].inflow.flowVolume[l] != 0) {
					printf("");
				}
				TransportNodes[ds].inflow.flowVolume[l]+=v;
			}
			TransportNodes[us].totalVolOutflow += v;
			l=ftGetLinkIndex(&TransportNodes[us].outflow, i);
			if (l>=0) {
				if (TransportNodes[us].outflow.flowVolume[l] != 0) {
					printf("");
				}
				TransportNodes[us].outflow.flowVolume[l] += v;
			}
		}
	}
	DBG_PRINT(ftPrintFlow(stime,dt););
	DBG_PRINT(maxIn=0;);
	DBG_PRINT(maxOut=0;);
	for (i=1; i<=Nnodes; i++) {
		Psource source;
		PNodeTransport n=&TransportNodes[i];
		// find all nodes where water is entering the system
		if (NodeDemand[i] < 0.0) {
			double v=-NodeDemand[i]*dt;
			n->totalVolInflow+=v;
			if (i>Njuncs) {
				// this is a tank or reservoir...
				int l, tidx;
				l=ftGetLinkIndex(&n->inflow, -i);
				if (n->inflow.flowVolume[l] != 0) {
					printf("");
				}
				n->inflow.flowVolume[l] += v;
				// tank or reservoir
				tidx=i-Njuncs;

				DBG_PRINT(if (fpTank!=NULL) { fprintf(fpTank, "%ld\tOUTFLOW\t%s\t%s\t%.10g\t%.10g\t%.10g\t%.10g",
						stime, Node[Tank[tidx].Node].ID, Tank[tidx].A==0?"Res":"Tank", Tank[tidx].V, NodeDemand[i], v, 0); });
				if (Tank[tidx].A==0) { //reservoir
					ftAddNodeSeg(&n->inflowSegs.segs[l], createseg(v, 0), 6);
					n->inflow.volumeFlowed[l] += v;
				} else { //tank
					ftAddTankOutflow(&n->inflowSegs.segs[l], tidx, v, dt);
					n->inflow.volumeFlowed[l]+=v;
				}
				DBG_PRINT(if (fpTank!=NULL) { fprintf(fpTank, "\t%.10g\t%.10g\t%.10g\n", Tank[tidx].V, Tank[tidx].C, Tank[tidx].M); fflush(fpTank); });
				// add this node to the list of nodes to process.
				ftAddNodeToProcess(n, nodesToProcess);
			} else {
				int l;
				//printf("Inflow from a node that is not a tank or reservoir\n");
				l=ftGetLinkIndex(&n->inflow, -i);
				if (n->inflow.flowVolume[l] != 0) {
					printf("");
				}
				n->inflow.flowVolume[l] += v;
				n->inflow.volumeFlowed[l]+=v;
				ftAddNodeSeg(&n->inflowSegs.segs[l], createseg(v, 0), 95);
				ftAddNodeToProcess(n, nodesToProcess);
			}
		} else {
			if (n->inflow.nlinks == 1  && n->inflow.linkIdxs[0] == 0) {
				// there must have been an imbalance in the hydraulic solution that would cause this condition
				// in order for the transport to function properly, add an inflow segment equal to demand+outflow
				// and add the node to nodes to process.
				double vol=n->totalVolOutflow+n->demand;

				volAdded += vol;
				ftAddNodeSeg(&n->inflowSegs.segs[0], createseg(vol, 0), 66);
				n->totalVolInflow=vol;
				n->inflow.volumeFlowed[0]+=vol;
				ftAddNodeToProcess(n, nodesToProcess);
				DBG_PRINT(if (fpLoss!=NULL) { fprintf(fpLoss, "%s\t%0.9f\t\t\n", n->node->ID, vol); });
			}
		}
		// get source mass
		source = Node[i].S;
		if (source != NULL && source->C0 > 0.0) {
			// compute mass added in this timestep
			TransportNodes[i].massadded=ftMassAdded(source, dt, &TransportNodes[i]);
		}
		DBG_PRINT(maxIn=MAX(maxIn, TransportNodes[i].inflow.nlinks););
		DBG_PRINT(maxOut=MAX(maxOut, TransportNodes[i].outflow.nlinks););
	}
    // compute the outflow adjustment.  This is used to ensure that inflow=outflow+demand.
	for (i=1; i<=Nnodes; i++) {
		PNodeTransport n=&TransportNodes[i];
		n->outflowAdjustment=(n->totalVolOutflow+n->demand)/n->totalVolInflow;
	}
	DBG_PRINT(ftPrintInflowOutflow(stime););
	DBG_PRINT(ftOpenInfoFile(stime,dt););
}

/*
**--------------------------------------------------------------
**   Input:   dt: the time step in seconds
**            stime: the current simulation time
**   Output:  none
**   Purpose: Move the water segments through the network
**--------------------------------------------------------------
*/
void ftMoveSegs(long dt, long stime) {
    // The list of nodes yet to process
	SNodeTransportList *nodesToProcess=NULL;
	int i;
#ifdef ENABLE_DBG_PRINT
	int seq=0;
	CurTime=stime;
#endif
#ifdef ENABLE_TIMING
	static FILE *fpt=NULL;
	clock_t t[8];
	t[0]=clock();
#endif
    // If the output is being redirected to a file, print the current time
	if (!_isatty(_fileno(stdout))) {
		fprintf(stdout, "%ld\n", stime); fflush(stdout);
	}
	else {
//		fprintf(stdout, "%ld\n", stime); fflush(stdout);
	}
	/* the absolute tolerance is dependent on the time step. */
	FT_ATOL=FT_ATOL_BASE*dt;
//    /* handle cases where TANK_DIRECT_MASS is in effect */
//	tankdirectinput(dt);
	ftInitMoveSegs(stime, dt, &nodesToProcess);
	TIMING(t[1]=clock();)
	/* start moving the water... */
	DBG_PRINT(ftPrintDetail(stime, "Begin", NULL, seq++, 0);)
	while (nodesToProcess != NULL) {
		Pseg seg;
		SegList segs;
		SNodeTransport *tNode;


        /* get the next node to process */
		tNode=ftGetNextNodeToProcess(&nodesToProcess);
		DBG_PRINT(printTrace = dbgPrintNodes[tNode->nodeIdx]);
		DBG_PRINT(if (fpTrace != NULL && printTrace) { fprintf(fpTrace, "\n"); fflush(fpTrace); });
		if (tNode->nodeIdx == 2226 && stime == 89928) {
//		if (tNode->nodeIdx == 109 && stime == 3600) {
				printf("");
		}
		if (tNode->nodeIdx == 55119) {
			traceProcessIncomplete = 1;
		} else {
			traceProcessIncomplete = 0;
		}
		DBG_PRINT(ftDocPrintNodeInfo(stime, tNode, "begin"));
		DBG_PRINT(ftDocPrintLinkVolumes(stime, tNode, "begin"));

		DBG_PRINT(ftDocPrintNode(stime, tNode, "begin"));
		DBG_PRINT(ftPrintDetail(stime, Node[tNode->nodeIdx].ID, NULL, seq++, 1);)
		DBG_PRINT(tNode->visited++;);
		DBG_PRINT(if (fpTrace != NULL && printTrace) { fprintf(fpTrace, "Processing node %s (%d) [%d]\n", tNode->node->ID, tNode->visited, tNode->nodeIdx); fflush(fpTrace); });
		DBG_PRINT(ftPrintNodeInfo(fpi, tNode, maxIn, maxOut););
		memset(&segs, 0, sizeof(segs));
		// get the segments that we can move from this node
		ftGetSegmentsToMove(&segs, tNode);
		DBG_PRINT(ftDocPrintSegments(stime,tNode,&segs, "segmentsToMove"));
		DBG_PRINT(ftPrintDetail(stime, Node[tNode->nodeIdx].ID, &segs, seq++, 2););
		if (ftHaveAnyInflowSegs(&tNode->inflowSegs)) {
            // if there are still inflow segments, add this node back to the "to proceee" list
			ftAddNodeToProcess(tNode, &nodesToProcess);
		}
		DBG_PRINT(ftPrintSegmentsToMove(fpi, &segs););
		seg=segs.firstSeg;
		while (seg != NULL) {
			// throughFlow is the amount of water that flows through the node
			// it is generally the same as seg, except in the case of incomplete or merged segments
			// because during processing, those types get split during demand removal, they don't
			// just have their volume reduced.
			Pseg throughFlow=seg;
			double inv=seg->v;
			double inc=seg->c;
            // adjust the outflow if necessary
			ftAdjustOutflow(seg, tNode->outflowAdjustment);
			DBG_PRINT(if (fpTrace != NULL && printTrace) { fprintf(fpTrace, "Moving segment:\tvol\t%.10g\tc\t%.10g\tadjusted:\tvol\t%.10g\tc\t%.10g\n", inv, inc, seg->v, seg->c); fflush(fpTrace); })
			if (tNode->demand > 0) {
				// remove demand from the segment
				throughFlow=ftRemoveDemand(seg, tNode, throughFlow);
				DBG_PRINT(ftDocPrintNode(stime,tNode, "post-demand"));
				DBG_PRINT(ftPrintDetail(stime, Node[tNode->nodeIdx].ID, NULL, seq++, 3);)
			}
			if (seg->v > 0) {
				// if there is any volume left to move, add the volume to the through flow
				ftAddThroughFlow(seg, tNode, throughFlow);
				DBG_PRINT(ftDocPrintNode(stime,tNode, "post-throughflow"));
				DBG_PRINT(ftPrintDetail(stime, Node[tNode->nodeIdx].ID, NULL, seq++, 4);)
				// and also move it to the node's outflow links.
				ftMoveFlowToLinks(seg, tNode, throughFlow);
				DBG_PRINT(ftDocPrintNode(stime,tNode, "post-moveflow"));
				DBG_PRINT(ftDocPrintLinks(stime,tNode, "post-moveflow"));
				DBG_PRINT(ftPrintDetail(stime, Node[tNode->nodeIdx].ID, NULL, seq++, 5);)
			}
			seg=seg->prev;
		}
		DBG_PRINT(ftPrintDetail(stime, Node[tNode->nodeIdx].ID, NULL, seq++, 6););
		// Move the water segemnts on the outflow links to the downstream nodes.
		ftMoveToDownstreamNodes(tNode, &nodesToProcess);
		DBG_PRINT(ftDocPrintNode(stime,tNode, "post-ds"));
		DBG_PRINT(ftDocPrintLinks(stime,tNode, "post-ds"));
		DBG_PRINT(ftPrintDetail(stime, Node[tNode->nodeIdx].ID, NULL, seq++, 7););
	}
	// now look to see if any nodes have mass that was not moved.  This likely indicates
	// that there was no flow in any attached pipe.
	for (i=1; i<=Nnodes; i++) {
		PNodeTransport n=&TransportNodes[i];
		double unusedMass = n->massadded-n->massused;
        // if there is, track it.
		if (unusedMass > 0) {
			MB_MassNoOutflow[n->nodeIdx]+=unusedMass;
		}
	}
	DBG_PRINT(ftPrintDetail(stime, "End", NULL, seq++, 8);)
	TIMING(t[2]=clock();)
#ifdef ENABLE_DBG_PRINT
	ftPrintWarningsFile(stime);
#endif
	TIMING(t[3]=clock());
	// complete any incomplete segments.  Generally, any incomplete
	// segments left have segments that are very small, so this should be OK
    // for networks that have significant hydraulic problems, there may be segments with large volumes.
	DBG_PRINT(ftPrintInflowOutflowData());
	ftCompleteAllIncompleteSegments(stime);
	TIMING(t[4]=clock();)
	// move any demand on tank nodes into the tanks.
	for (i=1; i<=Ntanks; i++) {
		int nidx=Tank[i].Node;
		PNodeTransport n=&TransportNodes[nidx];
		if (Tank[i].A > 0) {
			// this is a tank, not a reservoir
			// move all segments from allDemandSegs into the tank.
			ftAddTankInflow(&n->demandSegs, i);
		}
	}
	ftAdjustLinkSegmentVolumes();
	ftComputeMassBalance();
	DBG_PRINT(ftDocPrintConcentrations(stime, "end"));
	DBG_PRINT(printmassbalanceDEBUG(stime, dt));
	DBG_PRINT(ftPrintMassBalanceData(stime));
	DBG_PRINT(ftPrintMBArrays(stime));
    ftFinalizeMoveSegs(stime);
	TIMING(t[5]=clock();)
	DBG_PRINT(if (printAllLinkSegs) { ftPrintAllLinkSegs(stime, "flow"); })
	TIMING(t[6]=clock();)

#ifdef ENABLE_DBG_PRINT
	ftCloseDebugFiles();
	if (fpMassBalance != NULL) {
		printMassBalance(fpMassBalance, stime);
	}
#endif
#ifdef ENABLE_TIMING
	t[7]=clock();
	{
		int ii;
		if (fpt == NULL) {
			fpt=fopen("mstiming.txt", "w");
			fprintf(fpt, "Time");
			for (ii=0; ii<7; ii++) {
				fprintf(fpt, "\t%s", sects[ii]);
			}
			fprintf(fpt, "\n");
		}
		if (fpt!=NULL) {
			fprintf(fpt, "%d", stime);
			for (ii=0; ii<7; ii++) {
				fprintf(fpt, "\t%f", ((double)(t[ii+1]-t[ii]))/CLOCKS_PER_SEC);
			}
			fprintf(fpt, "\n");
		}
	}
#endif
}

/*
**--------------------------------------------------------------
**   Input:   stime: the current simulation time
**   Output:  none
**   Purpose: Perform necessary finalization
**--------------------------------------------------------------
*/
void ftFinalizeMoveSegs(long stime) {
	// move all demand segments created during this iteration to 
	// the allDemand segments list and the same for throughSegs.
	int inf=0, inc=0;
    int i;
	for (i=1; i<=Nnodes; i++) {
		int l;
		int linf = 0;
		PNodeTransport n=&TransportNodes[i];
		if (n->nodeIdx == 2225 && stime == 89928) {
			printf("");
		}
        // accumulate demand
		n->totalDemand+=n->demand;
        // move the demand segments to the allDemandSegs list
		ftMoveSegList(&n->demandSegs, &n->allDemandSegs, 4);
        // move the through segments to the allThroughSegs list
		ftMoveSegList(&n->throughSegs, &n->allThroughSegs, 5);
        // check for any inflow segments left
        for (l=0; l<n->inflow.nlinks; l++) {
			Pseg s=n->inflowSegs.segs[l].firstSeg;
			linf |= s!=NULL;
		}
		// if there are, report the issue and free them
		if (linf) {
			fprintf(stderr, "Node %s[%d]: has inflow segments in ftFinalizeMoveSegs.\n", n->node->ID, n->nodeIdx);
			for (l = 0; l<n->inflow.nlinks; l++) {
				ftFreeAllSegs(&n->inflowSegs.segs[l]);
			}
		}
		inf |= linf;
		// if there are any incomplete segments left,
		// free them
		if (n->incompleteSegs != NULL) {
			fprintf(stderr, "Incomplete segments found in ftFinalizeMoveSegs for node: %s[%d]\n",n->node->ID,n->nodeIdx);
			inc=1;
			while(n->incompleteSegs != NULL) {
				ftFreeIncompleteSegments(n->incompleteSegs,n->inflow.nlinks,&n->incompleteSegs);
			}
		}
	}
	// if HERE HERE
//	ftAlternatingSegmentCheck(stime);
	// if there were ANY nodes with inflow or incomplete segments report them
	if (inf || inc) {
		fprintf(stdout, "Are there any inflow or incomplete segements left?\n");  fflush(stdout);
		fprintf(stdout, "inflow: %s  incomplete: %s\n", inf?"yes":"no", inc?"yes":"no");  fflush(stdout);
	}
}

/*
**--------------------------------------------------------------
**   Input:   stime - current time
**   Output:  none
**   Purpose: Complete any incomplete segments that are left.
**--------------------------------------------------------------
*/
void ftCompleteAllIncompleteSegments(long stime) {
	static FILE *fp = NULL;
	int i;
	DBG_PRINT(ftPrintAllIncSegs(stime));

#ifdef ENABLE_DBG_PRINT
	for (i = 1; i <= Nnodes; i++) {
		PNodeTransport n = &TransportNodes[i];
		ftCheckForIncompleteSegments(n, stime);
	}
#endif
#ifdef LOG_COMPLETE_INC
	if (fp==NULL) {
		fp=fopen("complete_inc.txt","w");
	}
#endif
	for (i=1; i<=Nnodes; i++) {
		PNodeTransport n=&TransportNodes[i];
		if (i == 726 && stime == 89928) {
			printf("");
		}
		if (n->toProcess) {
            FILE *ofp=stdout;
			if(fp != NULL) { ofp=fp; }
			fprintf(ofp, "Didn't process node: %s [%d]\n", n->node->ID, n->nodeIdx);  fflush(fp);
		}
		if (n->incompleteSegs != NULL) {
            FILE *ofp=stdout;
			if(fp != NULL) { ofp=fp; }
			fprintf(ofp, "Completing incomplete segments for node: %s [%d]\n", n->node->ID, n->nodeIdx);  fflush(fp);
			ftCompleteIncompleteSegments(n,fp);
		}
	}
}


/*
**--------------------------------------------------------------
**   Input:   incSeg - the incomplete segment to be adjusted
**            adj - the adjustment factor
**   Output:  none
**   Purpose: Adjust flow on the incomplete segment
**--------------------------------------------------------------
*/
void ftAdjustIncompleteSegment(PIncompleteSeg incSeg, double adj) {
	int i;
	if (incSeg->nchildren > 0) {
		fprintf(stdout, "Incomplete segment with children to get adjusted.\n");  fflush(stdout);
	}
	// adjust the volume & save this adjustment
	incSeg->v *= adj;
	incSeg->adjustment = adj;
	// for each of the parts
	for (i=0; i<incSeg->nparts; i++) {
		PIncSegPart part=&incSeg->parts[i];
		// adjust the part's volume
		part->volume*=adj;
		// if there is any water flowed already
		if (part->flowed >0) { // an actual water segment on this inflow link.
			part->flowed*=adj;  // adjust the amount flowed
			part->conc/=adj;    // and the concentration (to conserver mass)
		}
		if (part->segs != NULL) {
			Pseg ps;
			// fort any segments inside this incomplete segment
			for (ps=part->segs; ps!=NULL; ps=ps->prev) {
				// perform necessary adjustments
				if (IsIncompleteSeg(ps)) {
					ftAdjustIncompleteSegment((PIncompleteSeg)ps, adj);
				} else if (IsMergedSeg(ps)) {
					ftAdjustMergedSeg((PMergedSeg)ps, adj);
				} else if (IsNormalSeg(ps)) {
					// shouldn't be here...
					fprintf(stdout, "Have a normal segment in a part...  Shouldn't happen?\n");  fflush(stdout);
				}
			}
		}
	}
}
/*
**--------------------------------------------------------------
**   Input:   seg - the merged segment to be adjusted
**            adj - the adjustment factor
**   Output:  none
**   Purpose: Adjust flow on the merged segment
**--------------------------------------------------------------
*/
void ftAdjustMergedSeg(PMergedSeg seg, double adj) {
	PIncompleteSeg incSeg;
	// adjust the volume
	seg->v*=adj;
	// and the volume already flowed
	seg->knownVol*=adj;
	// no concentration adjustment necessary - merged segments store known mass
	// adjust all of the incomplete segments as well
	for (incSeg=seg->incSegs; incSeg != NULL; incSeg=incSeg->prev) {
		ftAdjustIncompleteSegment(incSeg, adj);
	}
}
/*
**--------------------------------------------------------------
**   Input:   seg - the segment to be adjusted
**            outflowAdj - the adjustment factor
**   Output:  none
**   Purpose: Adjust the outflow on the  segment
**--------------------------------------------------------------
*/
void ftAdjustOutflow(Pseg seg, double outflowAdj) {
	if (outflowAdj==1)
		return;
	if (IsNormalSeg(seg)) {
		// maintain mass
		seg->v*=outflowAdj;
		seg->c/=outflowAdj;
	} else if (IsIncompleteSeg(seg)) {
		ftAdjustIncompleteSegment((PIncompleteSeg)seg, outflowAdj);
	} else if (IsMergedSeg(seg)) {
		ftAdjustMergedSeg((PMergedSeg)seg, outflowAdj);
	}
}
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: Adjust the link segment volumes.  This ensures that all
**    the link segments add up to the link volume.  Typically these
**    adjustments are very small.
**--------------------------------------------------------------
*/
void ftAdjustLinkSegmentVolumes() {
	int i;
	for (i=1; i<=Nlinks; i++) {
		double lv=LINKVOL(i);  // the link's volume
		double lwv=ftComputeSegmentWaterVolume(FirstSeg[i]);  // the link's segments total
		Pseg pmax=NULL, p;
		double m, vadj;
		if (FirstSeg[i]==NULL) {
			// must be a pipe with a volume smaller than FT_ATOL
			addseg(i, lv, 0);
		} else {
			vadj=lwv-lv;  // the amount to add or remove
			if (vadj != 0) {
				// find the largest segment and adjust its volume and mass;
				for (p = FirstSeg[i]; p != NULL; p = p->prev) {
					if (pmax == NULL || p->v > pmax->v) pmax=p;
				}
				m=pmax->v*pmax->c; // mass
				pmax->v-=vadj;  // add (or subtract) the volume
// !!! Need to check if there are any situations where the adjustment actually leaves a negative volume !!!
				if (pmax->v>0) { // if there is still volume left
					pmax->c=m/pmax->v;   // compute the new concentration
					checkConc(pmax->c);
				} else {
					pmax->c=0;
				}
			}
		}
	}
}
/*
**--------------------------------------------------------------
**   Input:   segs - list of segments
**   Output:  total volume of merged segments in seg list
**   Purpose: Compute the total merged segment volume in the list of segments
**--------------------------------------------------------------
*/
double ftMergedSegmentVolume(Pseg segs) {
	double rv=0;
	for(Pseg s = segs; s != NULL; s=s->prev) {
		if(IsMergedSeg(s)) rv+=s->v;
	}
	return rv;
}
/*
**--------------------------------------------------------------
**   Input:   The node with incomplete segments
**   Output:  none
**   Purpose: Complete the node's incomplete segments.
**            The volume of any remaining incomplete segments
**            at this point are generally very small, so set
**            their concentration to zero and let the completion
**            propogate
**--------------------------------------------------------------
*/
void ftCompleteIncompleteSegments(PNodeTransport n, FILE *fp) {
	PIncompleteSegmentData prevInc=NULL;
	DBG_PRINT(ftPrintAllIncompleteSegments(n));
	while (n->incompleteSegs != NULL && prevInc != n->incompleteSegs) {
		int l, complete;
		PIncompleteSeg *incsegs=n->incompleteSegs->incomplete;
		prevInc=n->incompleteSegs;
		complete=0;
		for (l=0; l<n->inflow.nlinks && !complete; l++) {
			PIncompleteSeg is=incsegs[l];
			if (is != NULL) {
				int p=l;
				PIncSegPart part=&is->parts[p];
				// INVESTIGATE: the msVol computation and check is currently bypassed - it was discovered that
				// merged segments were causing issues sometimes.  Leaving the test in for now in case it is needed again
				double msVol = ftMergedSegmentVolume(part->segs);
				if(msVol > 0) {
					printf("");
				}
				msVol=0;
				if(part->conc == -1 || part->flowed + msVol < part->volume) {
						double v=part->volume-part->flowed-msVol;
						if(fp != NULL) { fprintf(fp, "  vol: %.10g\n", v);  fflush(fp); }
						complete=ftUpdateIncompleteSegs(n, l, 0, &v, TRUE);
// INVESTIGATE:
//	ftCheckForInvalidSegments();
				}
			}
		}
	}
	if (n->incompleteSegs != NULL) {
		fprintf(stdout, "ERROR: could not complete incomplete segments.\n");  fflush(stdout);
	}
////	ftCheckForInvalidSegments();
}
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: Perfrom necessary calculations when the transport
**            is done
**--------------------------------------------------------------
*/
void ftTransportDone() {
	int i;
	for (i=1; i<=Nnodes; i++) {
		PNodeTransport n=&TransportNodes[i];
		// if the node has demand, set the instantaneous concentration
		// to the the most recently added demand segment's concentration
		// and compute the average concetration from all the demand segments
		if (NodeDemand[i]>0) {
			if (n->allDemandSegs.lastSeg != NULL) {
				Mass[i]=ftGetMass(n->allDemandSegs.firstSeg);
				NodeQual[i]=n->allDemandSegs.lastSeg->c;
				//  not using ftGetAvgConc here because the summed volume
				// may not be exsctly the same as the computed demand due to tolerance
				// issues, so simply sum the mass and divide it by the actaul demand.
				AvgNodeQual[i] = Mass[i] / n->totalDemand;
			} else if(i<=Njuncs) {
				Mass[i]=0;
				NodeQual[i]=0;
				AvgNodeQual[i] = 0;
			}
		} else {
			// if the node does not have demand, set the instantaneous concentration
			// to the the most recently added through segment's concentration
			// and compute the average concetration from all the through segments
			Mass[i]=0;
			if (n->allThroughSegs.lastSeg != NULL) {
				NodeQual[i]=n->allThroughSegs.lastSeg->c;
				AvgNodeQual[i] = ftGetAvgConc(&n->allThroughSegs);
			} else {
				NodeQual[i]=0;
				AvgNodeQual[i] = 0;
			}
		}
		// and free all the through and demand segments
		ftFreeAllSegs(&n->allThroughSegs);
		ftFreeAllSegs(&n->allDemandSegs);
		n->totalDemand=0;
		DBG_PRINT(ftdAccumulateMassAdded(i));
	}

	DBG_PRINT(appendConcFile(fpAllAvgConc, AvgC));
	DBG_PRINT(appendConcFile(fpAllInstConc, C));
}
/*
**--------------------------------------------------------------
**   Input:   links: the node's links (either inflow or outflow)
**            linkIdx: The link's global index
**            flow: the flow on the link
**   Output:  none
**   Purpose: Add the link to the links collection
**--------------------------------------------------------------
*/
void ftAddLink(SNodeLinks *links, int linkIdx, double flow) {
	int i=links->nlinks;
	links->nlinks++;
	if (links->nlinks > links->size) {
		// use size to remember how big these arrays are.  If this was not done,
		// these arrays would be allocated and reallocated every hydraulic timestep,
		// when they likely won't change size very often.
		links->linkIdxs=(int *)realloc(links->linkIdxs, links->nlinks*sizeof(int));
		links->flowRatio=(double*)realloc(links->flowRatio, links->nlinks*sizeof(double));
		links->flowRate=(double*)realloc(links->flowRate, links->nlinks*sizeof(double));
		links->flowVolume=(double*)realloc(links->flowVolume, links->nlinks*sizeof(double));
		links->volumeFlowed=(double*)realloc(links->volumeFlowed, links->nlinks*sizeof(double));
		links->size++;
	}
	links->linkIdxs[i]=linkIdx;
	links->flowRate[i]=flow;
}
/*
**--------------------------------------------------------------
**   Input:   links: the node's links (either inflow or outflow)
**   Output:  none
**   Purpose: Compute the ratio of each link to the total of all
**            links for the node
**--------------------------------------------------------------
*/
void ftComputeRatios(SNodeLinks *links) {
	if (links->nlinks==1) {
		links->flowRatio[0]=1.0;
	} else {
		double totFlow=0;
		int i;
		for (i=0; i< links->nlinks; i++) {
			totFlow+=links->flowRate[i];
		}
		for (i=0; i< links->nlinks; i++) {
			// compute ratio of link's flow to total flow
			links->flowRatio[i]=links->flowRate[i]/totFlow;
		}
	}
}
/*
**--------------------------------------------------------------
**   Input:   segs: The list of segemnts to compute the average concentration
**   Output:  the average concentration
**   Purpose: Compute the average concentration of the list of segments
**--------------------------------------------------------------
*/
double ftGetAvgConc(PSegList segs) {
	Pseg p=segs->firstSeg;
	double m=0, v=0;
	while (p != NULL) {
		m+=p->v*p->c;  // compute mass
		v+=p->v;       // sum volume
		p=p->prev;
	}
	return v!=0?m/v:0;
}
/*
**--------------------------------------------------------------
**   Input:   segs: The list of segemnts to compute the mass
**   Output:  the average concentration
**   Purpose: Compute the mass contained of the list of segments
**--------------------------------------------------------------
*/
double ftGetMass(Pseg segs) {
	Pseg p=segs;
	double m=0, v=0;
	while (p != NULL) {
		m+=p->v*p->c;
		p=p->prev;
	}
	return m;
}
/*
**--------------------------------------------------------------
**   Input:   segs: The list of segmens to release
**   Output:  none
**   Purpose: Release the segments back to the available pool
**--------------------------------------------------------------
*/
void ftReleaseSegList(PSegList segs) {
	Pseg s;
	while (segs->firstSeg != NULL) {
		s=segs->firstSeg->prev;
		ftReleaseSeg(segs->firstSeg);
		segs->firstSeg=s;
	}
	segs->firstSeg=NULL;
	segs->lastSeg=NULL;
}
/*
**--------------------------------------------------------------
**   Input:   links: the node's inflow or outflow data
**   Output:  none
**   Purpose: Release memory used by the node's inflow and outflow data
**--------------------------------------------------------------
*/
void ftFreeNodeLinks(PNodeLinks links) {
	free(links->flowRatio);
	free(links->flowRate);
	free(links->flowVolume);
	free(links->linkIdxs);
	free(links->volumeFlowed);
}
/*
**--------------------------------------------------------------
**   Input:   links: the node's links (either inflow or outflow)
**   Output:  none
**   Purpose: Reset the total flow volume and the volume flowed
**--------------------------------------------------------------
*/
void ftResetFlows(PNodeLinks links) {
	int l, n;
	n=links->nlinks;
	for (l=0; l<n; l++) {
		links->flowVolume[l]=0;
		links->volumeFlowed[l]=0;
	}
}
/*
**--------------------------------------------------------------
**   Input:   tankIdx: the tank's index
**            masadded: the amount of mass added into the tank
**   Output:  none
**   Purpose: Called from tankdirectinput in quality.c to
**            update the MassBalance arrays appropriately
**--------------------------------------------------------------
*/
void ftTankDirectInput(int tankIdx, double massadded) {
	int nodeIdx=tankIdx+Njuncs;
	TransportNodes[nodeIdx].massused+=massadded;
	MB_MassMoved[nodeIdx]+=massadded;
	MB_MassAdded[tankIdx+Njuncs]+=massadded;
}
/*
**--------------------------------------------------------------
**   Input:   segs: the water flowing into a tank
**            tankIdx: the tank's index
**   Output:  none
**   Purpose: move water into the tank and remove segments from segs.
**            After this call, segs will have no more segments
**--------------------------------------------------------------
*/
void ftAddTankInflow(PSegList segs, int tankIdx) {
	DBG_PRINT(if (fpTank!=NULL && NodeDemand[Tank[tankIdx].Node]>=0) {
		fprintf(fpTank, "%ld\tINFLOW\t%s\t%s\t%.10g\t%.10g\t%.10g",
			CurTime, Node[Tank[tankIdx].Node].ID, Tank[tankIdx].A==0?"Res":"Tank",
			Tank[tankIdx].V, NodeDemand[Tank[tankIdx].Node],
			ftSumVolume(segs), ftSumMass(segs)); })
		switch (Tank[tankIdx].MixModel) {
		case MIX2: ftTankInflow2(tankIdx, segs); break;
		case FIFO: ftTankInflow3(tankIdx, segs); break;
		case LIFO: ftTankInflow4(tankIdx, segs); break;
		default:   ftTankInflow1(tankIdx, segs); break;
	}
	DBG_PRINT(if (fpTank!=NULL && NodeDemand[Tank[tankIdx].Node]>=0) {
		fprintf(fpTank, "\t%.10g\t%.10g\t%.10g\n",
			Tank[tankIdx].V, Tank[tankIdx].C, Tank[tankIdx].M);
		fflush(fpTank);
	})
}
/*
**--------------------------------------------------------------
**   Input:   segs: the list of segs that will receive this tank's outflow
**            tankIdx: the index of the tank
**            vol: the amount of water to remove from the tank
**            dt: the timestep
**   Output:  none
**   Purpose: Move water from the tank into segs
**--------------------------------------------------------------
*/
void ftAddTankOutflow(PSegList segs, int tankIdx, double vol, int dt) {
	switch (Tank[tankIdx].MixModel) {
	case MIX2: ftTankOutflow2(tankIdx, vol, segs, dt); break;
	case FIFO: ftTankOutflow3(tankIdx, vol, segs, dt); break;
	case LIFO: ftTankOutflow4(tankIdx, vol, segs, dt); break;
	default:   ftTankOutflow1(tankIdx, vol, segs, dt); break;
	}
}
/*
**--------------------------------------------------------------
**   Input:   tankIdx: the index of the tank
**            vol: the volume to move from the tank
**            systemInflow: the segment list to move the segments into
**            dt: the timestep
**   Output:  none
**   Purpose: Move water from a fully mixed tank
**--------------------------------------------------------------
*/
void ftTankOutflow1(int tankIdx, double vol, PSegList systemInflow, int dt) {
	// remove water from the tank, create a new inflow segment
	// and add it to systemInflow

	double c;
	Stank *tank=&Tank[tankIdx];

	c=tankreact(tank->C, tank->V, tank->Kb, dt);
	if (vol<= tank->V) {
		tank->V -= vol;
		tank->M -= vol*c;
		// if there is a bulk growth reaction specified, the mass in the tank could be < 0
		// ensure that doesn't happen.
		if (tank->M < 0) tank->M=0;
		tank->C = tank->M/tank->V;
	} else {
		// the tank has been "overdrawn".  Compute the new concentration to conserve mass
		// ands set the volume, mass and concentration to 0
		c=tank->M/vol;
		tank->V=0;
		tank->M=0;
		tank->C=0;
	}
	ftAddNodeSeg(systemInflow, createseg(vol, c), 7);
}
/*
**--------------------------------------------------------------
**   Input:   tankIdx: the index of the tank
**            inflow: the segment list containing the water to
**                    move into the tank
**   Output:  none
**   Purpose: Move water into a fully mixed tank
**--------------------------------------------------------------
*/
void ftTankInflow1(int tankIdx, PSegList inflow) {
	Pseg seg;
	Stank *tank=&Tank[tankIdx];
	for (seg=inflow->firstSeg; seg != NULL; seg=seg->prev) {
		tank->V += seg->v;
		tank->M += seg->v*seg->c;
	}
	if (tank->V < 0) {
		MB_MassNegTankVol[tankIdx]+=tank->M;
		tank->C=0;
		tank->M=0;
	} else {
		tank->C=tank->M/tank->V;
	}
	NodeQual[Tank[tankIdx].Node]=Tank[tankIdx].C;
	ftReleaseSegList(inflow);
}

/*
**--------------------------------------------------------------
**   Input:   tankIdx: the index of the tank
**            vol: the volume to move from the tank
**            systemInflow: the segment list to move the segments into
**            dt: the timestep
**   Output:  none
**   Purpose: Move water from a 2-compartment tank
**--------------------------------------------------------------
*/
void ftTankOutflow2(int tankIdx, double vol, PSegList systemInflow, int dt) {
	// 2 compartment mixing.
	// remove water from the mixing zone, possibly the main zone
	// move water from the main zone to the mixing zone,
	// create new inflow segments and add them to systemInflow

	Pseg mainZone, mixingZone;
	double vremain=vol;
	Stank *tank=&Tank[tankIdx];
	double v1max=tank->V1max;
	double v;

	mixingZone=LastSeg[Nlinks+tankIdx];
	mainZone=FirstSeg[Nlinks+tankIdx];
	if (mainZone==NULL || mixingZone==NULL) return; // why?  could this happen?

	mixingZone->c=tankreact(mixingZone->c, mixingZone->v, tank->Kb, dt);
	mainZone->c=tankreact(mainZone->c, mainZone->v, tank->Kb, dt);

	// take from the mixing zone first
	v=MIN(vremain, mixingZone->v);
	ftAddNodeSeg(systemInflow, createseg(v, mixingZone->c), 8);
	mixingZone->v-=v;
	vremain-=v;
	if (vremain > 0) {
		// take from the main zone
		v=MIN(vremain, mainZone->v);
		ftAddNodeSeg(systemInflow, createseg(v, mainZone->c), 9);
		mainZone->v-=v;
		vremain-=v;
	}
	// now move water from the main zone to the mixing zone
	if (mixingZone->v < v1max) {
		double vmove=MIN(v1max-mixingZone->v, mainZone->v); // can't move more water than we have;
		double mass=mixingZone->v*mixingZone->c+vmove+mainZone->c; // new mass in mixing zone
		mixingZone->v+=vmove;  mainZone->v-=vmove;  // move the volume
		mixingZone->c=mass/mixingZone->v;  // compute the new concentration
		if (mainZone->v==0) { mainZone->c=0; }
	}
	if (mixingZone->v==0) { mixingZone->c=0; }

	tank->C=mixingZone->c;
	tank->V=mixingZone->v+mainZone->v;
	tank->M=mixingZone->v*mixingZone->c+mainZone->v*mainZone->c;
}
/*
**--------------------------------------------------------------
**   Input:   tankIdx: the index of the tank
**            inflow: the segment list containing the water to
**                    move into the tank
**   Output:  none
**   Purpose: Move water into a 2-compartment tank
**--------------------------------------------------------------
*/
void ftTankInflow2(int tankIdx, PSegList inflow) {
	// 2 compartment mixing.
	// for each incoming segment, add volume & concentration to the mixing zone
	// and move excess to main zone

	Stank *tank=&Tank[tankIdx];
	Pseg mainZone, mixingZone, seg;
	double v1max=tank->V1max;

	mixingZone=LastSeg[Nlinks+tankIdx];
	mainZone=FirstSeg[Nlinks+tankIdx];

	for (seg=inflow->firstSeg; seg!=NULL; seg=seg->prev) {
		double volTransfer=MAX(0.0, mixingZone->v + seg->v - v1max);
		mixingZone->c=(mixingZone->c*mixingZone->v + seg->v*seg->c)/(mixingZone->v+seg->v);
		mixingZone->v+=seg->v-volTransfer;
		if (volTransfer > 0.0) {
			mainZone->c=(mainZone->v*mainZone->c + volTransfer * mixingZone->c)/(mainZone->v + volTransfer);
			mainZone->v+=volTransfer;
		}
	}
	ftReleaseSegList(inflow);
}
/*
**--------------------------------------------------------------
**   Input:   tankIdx: the index of the tank
**            vol: the volume to move from the tank
**            systemInflow: the segment list to move the segments into
**            dt: the timestep
**   Output:  none
**   Purpose: Move water from a FIFO tank
**--------------------------------------------------------------
*/
void ftTankOutflow3(int tankIdx, double vol, PSegList systemInflow, int dt) {
	// FIFO tank mixing model
	int k=Nlinks+tankIdx;
	Pseg seg, outSeg;
	Stank *tank=&Tank[tankIdx];
	double vremain=vol;
	double vsum, msum;

	if (Reactflag) {
		seg=FirstSeg[k];
		while (seg != NULL) {
			seg->c=tankreact(seg->c, seg->v, tank->Kb, dt);
			seg=seg->prev;
		}
	}
	// withdraw starting at FirstSeg and add those segments to systemInflow
	seg=FirstSeg[k];
	msum=0;
	vsum=0;
	while (vremain > 0.0 && seg != NULL) {
		if (seg->v > vremain) {
			// simply unlink from FirstSeg & add it to systemInflow.
			ftUnlinkFirstLinkSegment(k);
			outSeg=seg;
		} else {
			outSeg=createseg(vremain, seg->c);
			seg->v-=vremain;
		}
		vremain-=outSeg->v;
		vsum+=outSeg->v;
		msum+=outSeg->v*outSeg->c;
		ftAddNodeSeg(systemInflow, outSeg, 10);
		seg=FirstSeg[k];
	}
	tank->M-=msum;
	tank->V-=vsum;
	// if there is a bulk growth reaction specified, the mass in the tank could be < 0
	// ensure that doesn't happen.
	if (tank->M < 0) tank->M=0;
	if (vsum>0) {
		tank->C=msum/vsum;
	} else {
		tank->C=FirstSeg[k]->c;
	}
}
/*
**--------------------------------------------------------------
**   Input:   tankIdx: the index of the tank
**            inflow: the segment list containing the water to
**                    move into the tank
**   Output:  none
**   Purpose: Move water into a FIFO tank
**--------------------------------------------------------------
*/
void ftTankInflow3(int tankIdx, PSegList inflow) {
	// FIFO tank mixing model
	int k=Nlinks+tankIdx;
	Stank *tank=&Tank[tankIdx];
	double vsum, msum;

	msum=0;
	vsum=0;
	while (inflow->firstSeg != NULL) {
		Pseg seg=ftUnlinkFirstSegment(inflow);
		insertseg(k, seg);
		vsum+=seg->v;
		msum+=seg->v*seg->c;
	}
	tank->M+=msum;
	tank->V+=vsum;
	tank->C=FirstSeg[k]->c;
}
/*
**--------------------------------------------------------------
**   Input:   tankIdx: the index of the tank
**            vol: the volume to move from the tank
**            systemInflow: the segment list to move the segments into
**            dt: the timestep
**   Output:  none
**   Purpose: Move water from a LIFO tank
**--------------------------------------------------------------
*/
void ftTankOutflow4(int tankIdx, double vol, PSegList systemInflow, int dt) {
	// LIFO
	int k=Nlinks+tankIdx;
	Pseg seg, outSeg;
	Stank *tank=&Tank[tankIdx];
	double vsum, msum;
	double vremain=vol;

	if (Reactflag) {
		seg=FirstSeg[k];
		while (seg != NULL) {
			seg->c=tankreact(seg->c, seg->v, tank->Kb, dt);
			seg=seg->prev;
		}
	}
	// withdraw starting at FirstSeg and add those segments to systemInflow
	seg=LastSeg[k];
	msum=0;
	vsum=0;
	while (vremain > 0.0 && seg != NULL) {
		if (seg->v > vremain) {
			// simply unlink from LastSeg & add it to systemInflow.
			LastSeg[k]=seg->prev;
			seg->prev=NULL;
			outSeg=seg;
		} else {
			outSeg=createseg(vremain, seg->c);
			seg->v-=vremain;
		}
		vremain-=outSeg->v;
		vsum+=outSeg->v;
		msum+=outSeg->v*outSeg->c;
		ftAddNodeSeg(systemInflow, outSeg, 11);
		seg=LastSeg[k];
	}
	tank->M-=msum;
	tank->V-=vsum;
	// if there is a bulk growth reaction specified, the mass in the tank could be < 0
	// ensure that doesn't happen.
	if (tank->M < 0) tank->M=0;
	if (vsum>0) {
		tank->C=msum/vsum;
	} else {
		tank->C=FirstSeg[k]->c;
	}
}
/*
**--------------------------------------------------------------
**   Input:   tankIdx: the index of the tank
**            inflow: the segment list containing the water to
**                    move into the tank
**   Output:  none
**   Purpose: Move water into a LIFO tank
**--------------------------------------------------------------
*/
void ftTankInflow4(int tankIdx, PSegList inflow) {
	// LIFO tank mixing model
	int k=Nlinks+tankIdx;
	Stank *tank=&Tank[tankIdx];
	double vsum, msum;

	msum=0;
	vsum=0;
	while (inflow->firstSeg != NULL) {
		Pseg seg=ftUnlinkFirstSegment(inflow);
		Pseg tmpseg=LastSeg[k];
		LastSeg[k]=NULL;
		insertseg(k, seg);
		LastSeg[k]->prev=tmpseg;
		vsum+=seg->v;
		msum+=seg->v*seg->c;
	}
	tank->M+=msum;
	tank->V+=vsum;
	tank->C=LastSeg[k]->c;
}

/*
**--------------------------------------------------------------
**   Input:   segs: The list of segments
**   Output:  the total mass of all the segments.  Returns -1
**            if there were incomplete segments, -2 if there
**            were merged segments, or -3 if there were both.
**            Otherwise returns the actual mass if all segments
**            are regular segments
**   Purpose: Compute the volume of the segments
**--------------------------------------------------------------
*/
double ftSumMass(PSegList segs) {
	double mass=0;
	int nInc=0, nMerged=0;
	Pseg seg;
	for (seg=segs->firstSeg; seg != NULL; seg=seg->prev) {
		if (IsNormalSeg(seg)) {
			mass += seg->v * seg->c;
		} else if (IsIncompleteSeg(seg)) {
			nInc++;
		} else if (IsMergedSeg(seg)) {
			nMerged++;
		} else {
			// concentration is < 0
		}
	}
	if (nMerged==0 && nInc==0) {
		return mass;
	} else if (nMerged > 0 && nInc > 0) {
		return -3;
	} else if (nInc > 0) {
		return -1;
	} else {
		return -2;
	}
}

/*
**--------------------------------------------------------------
**   Input:   segs: The list of segments
**   Output:  the total volume of all the segments
**   Purpose: Compute the volume of the segments
**--------------------------------------------------------------
*/
double ftSumVolume(PSegList segs) {
	double vol=0;
	Pseg seg;
	for (seg=segs->firstSeg; seg != NULL; seg=seg->prev) {
		vol += seg->v;
	}
	return vol;
}

/*
**--------------------------------------------------------------
**   Input:   from: The source segment list
**            to: The destination segment list
**   Output:  none
**   Purpose: Move the segments from the from list and put them
**            at the beginning of the to list
**--------------------------------------------------------------
*/
void ftMoveSegList(PSegList from, PSegList to, int loc) {
	if (to->lastSeg != NULL) {
		to->lastSeg->prev=from->firstSeg;
	}
	to->lastSeg=from->lastSeg;
	if (to->firstSeg==NULL) {
		to->firstSeg=from->firstSeg;
	}
	from->firstSeg=NULL;
	from->lastSeg=NULL;
	ftCombineSegs(&to->firstSeg, &to->lastSeg, loc);
}

/*
**--------------------------------------------------------------
**   Input:   tNode: The node
**            nodesToProcess: The list of nodes to process.  This
**                            list  may be updated by this method
**   Output:  none
**   Purpose: Move excess water from each of the node's outflow
**            links to their downstream nodes.
**--------------------------------------------------------------
*/
void ftMoveToDownstreamNodes(PNodeTransport tNode, PNodeTransportList *nodesToProcess) {
	int l;
    static int cnt=0;
	PNodeLinks out=&tNode->outflow;
    cnt++;
	// move excess water from links to downstream nodes.
	for (l=0; l<out->nlinks; l++) {
		int linkIdx=out->linkIdxs[l];
		SNodeTransport *downNode=&TransportNodes[DOWN_NODE(linkIdx)];
		int inflowIdx=ftGetLinkIndex(&downNode->inflow, linkIdx);
		// compute the excess volume that needs to be moved
		double lwv=ftComputeSegmentWaterVolume(FirstSeg[linkIdx]);
		double lv=LINKVOL(linkIdx);
		double excessVolume=lwv-lv;
		PNodeLinks in = &downNode->inflow;

		DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "MOVEDSINFO\t\t%s\t%s\t%.8f", Link[linkIdx].ID, downNode->node->ID, excessVolume); fflush(fpi); })
		DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\n");  fflush(fpi); })
		DBG_PRINT(if (fpTrace!=NULL && printTrace) { fprintf(fpTrace, "Moving excess volume\t%.10g\t from link %s to node %s\n", excessVolume, Link[linkIdx].ID, Node[DOWN_NODE(linkIdx)].ID); fflush(fpTrace); })
		while (excessVolume>0 && FirstSeg[linkIdx] != NULL) {
			Pseg seg=FirstSeg[linkIdx];
			if (IsNormalSeg(seg)) {
				// for "normal" segments
				Pseg newSeg;
				double volMoved;
				double massMoved=0;
				double delta;
				double volumeWaiting;
				DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "MOVEDS\t\t\t%.8f\t%.8f\t%.8f", excessVolume, seg->v, seg->c);  fflush(fpi); })
				// compute the amount of excess water to move
				volMoved=MIN(seg->v, excessVolume);
				massMoved=volMoved*seg->c;
				DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f", volMoved);  fflush(fpi); })
				DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f", volMoved);  fflush(fpi); })
				// remove the volume moved from the excess volume and from the segment
				excessVolume-=volMoved;
				seg->v -= volMoved;
				DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f\t%.8f", seg->v, excessVolume);  fflush(fpi); })
				// if this node has any incomplete segments, this may be a segment it needs to help complete it
				if (downNode->incompleteSegs != NULL) {
					DBG_PRINT(if(fpIncTrace!=NULL) { fprintf(fpIncTrace, "Updating incomplete segments for node %s [%d]\n", downNode->node->ID, downNode->nodeIdx); fflush(fpIncTrace); })
					ftUpdateIncompleteSegs(downNode,inflowIdx,seg->c,&volMoved, FALSE);
				}
				// INVESTIGATE: if the volume moved + the volume flowed to this point is sufficiently close to the expected flow volume
				// set the volume moved to the expected flow - the volume flowed.
				// this will also aid in stabilizing the algorithm
				// not sure if this is needed anymore with the inflow adjustment mod
				volumeWaiting=ftSumVolume(&downNode->inflowSegs.segs[inflowIdx]);
				delta = ABS(volMoved+in->volumeFlowed[inflowIdx]+volumeWaiting-in->flowVolume[inflowIdx]);
				if (delta > 1e-8 && ftWithinTolerance(volMoved+in->volumeFlowed[inflowIdx]+volumeWaiting, in->flowVolume[inflowIdx], 1e-8, FT_RTOL)) {
					DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "MTDN1\t%.10g\n", volMoved-in->volumeFlowed[inflowIdx]+volumeWaiting-in->flowVolume[inflowIdx]); fflush(fpDeltas); })
					volMoved=in->flowVolume[inflowIdx]-(in->volumeFlowed[inflowIdx]+volumeWaiting);
				}
				if (volMoved >0) {
					// create a new segment
					newSeg=createseg(volMoved, massMoved/volMoved);
					DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f\t%.8f\n", newSeg->v, newSeg->c); })
					DBG_PRINT(if (fpTrace!=NULL && printTrace) { fprintf(fpTrace, "Moving volume\t%.10g\tc\t%.10g\n", newSeg->v, newSeg->c); fflush(fpTrace); })
					// add the new segment to thew inflowSegs of the downstream node
					ftAddNodeSeg(&downNode->inflowSegs.segs[inflowIdx], newSeg, 12);
					downNode->inflow.volumeFlowed[inflowIdx] += volMoved;
					// and add the downstream node to the list of nodes to process
					ftAddNodeToProcess(downNode, nodesToProcess);
				} else {
					DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t\t\n"); })
				}
				if (seg->v == 0) {
					if (linkIdx > Npipes && seg->prev == NULL) {
						// this is a valve or pump - just set the concentration to 0 if it is the only segment
						seg->c=0;
					} else {
						// remove this segment
						ftRemoveFirstLinkSegment(linkIdx);
					}
				}
				// update the downstream node's inflow volume flowed.
				if (downNode->inflow.volumeFlowed[inflowIdx]+volumeWaiting+volMoved > downNode->inflow.flowVolume[inflowIdx]+1e-7) {
#ifdef DEBUG
					fprintf(stdout, "");  fflush(stdout);
#endif
					DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "MTDN6\t%.10g\n", downNode->inflow.volumeFlowed[inflowIdx]+volMoved - downNode->inflow.flowVolume[inflowIdx]); fflush(fpDeltas); })
				}
				// if the excess volume is sufficiently close to zero, set it to zero.
				if (excessVolume > 0 && ftWithinTolerance(excessVolume, 0, 1e-8, FT_RTOL)) {
					DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "MTDN2\t%.10g\n", excessVolume); fflush(fpDeltas); })
					DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "MOVEDS-EXCVOL\t\t%.8e\n", excessVolume); })
					excessVolume=0;
				}
				DBG_PRINT(if (fpTrace!=NULL && printTrace) { fprintf(fpTrace, "flowvol-volumeflowed\t%.10g\n", downNode->inflow.flowVolume[inflowIdx]-downNode->inflow.volumeFlowed[inflowIdx]-volumeWaiting); fflush(fpTrace); })
			} else {
				// the segment being moved is either a merged segment or incomplete segment.
				// this segment will need to be split
				double volMoved;
				DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "MOVEDS\t\t\t%.8f\t%.8f\t%.8f", excessVolume, seg->v, seg->c); })
				// compute the volume moved
				volMoved=MIN(seg->v, excessVolume);
				DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f", volMoved); })
				if (seg->v != volMoved && ftWithinTolerance(seg->v, volMoved, FT_ATOL, FT_RTOL)) {
					// if it is sufficiently close to the segement's volume, use that amount.
					// this keeps many extremely small segments from creeping into the system
					// and causing the algorithm to not work
					DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "MTDN3\t%.10g\n", seg->v-volMoved); fflush(fpDeltas); })
				}
				DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f", volMoved); })
				// decrement the excess volume by the volume moved
				excessVolume-=volMoved;
				// if the downstream node has incomplete segments, the segMoved may be one of the segments needed
				// to complete it.
				if (downNode->incompleteSegs != NULL) {
					Pseg origSeg=seg;
					PIncompleteSegmentData incompleteData, t;
					incompleteData=downNode->incompleteSegs;
					t=NULL;
					while (incompleteData != NULL && volMoved > 0) {
						// got the segment we've been waiting for, but it's either incomplete or is merged...
						// split segment into vol needed to complete and remainder (if any)
						DBG_PRINT(if (fpUpdateIncSeg != NULL) { fprintf(fpUpdateIncSeg, "ftMoveToDownstreamNode-a\n"); })
						DBG_PRINT(if (fpUpdateIncSeg != NULL) { ftPrintIncompleteSegment(incompleteData->incomplete[inflowIdx], 0, fpUpdateIncSeg, NULL); })
						if (incompleteData->incomplete[inflowIdx] != NULL) {
							double volUsed;
							int isMerged = IsMergedSeg(seg);
							int isIncSeg = IsIncompleteSeg(seg);
#ifdef ENABLE_DBG_PRINT
							char *stype=getSegmentType(seg);
#endif
							volUsed=ftProcessIncomplete(incompleteData->incomplete[inflowIdx], &seg, inflowIdx, TRUE, FALSE);
							DBG_PRINT(if (fpTrace != NULL) { fprintf(fpTrace, "Received %s segment for an incomplete segment into node %s.  vol used: %0.10g\n", stype, downNode->node->ID, volUsed); fflush(fpTrace); })
							if (isMerged || isIncSeg) {
//								downNode->inflow.volumeFlowed[inflowIdx]+=volUsed;
							}
							volMoved-=volUsed;
							if (seg == NULL) {
								ftUnlinkFirstLinkSegment(linkIdx);
							} else if (seg != origSeg) {
								ftReplacePseg(origSeg, seg, &FirstSeg[linkIdx], &LastSeg[linkIdx]);
								origSeg=seg;
							}
						}
						DBG_PRINT(if (fpUpdateIncSeg != NULL) { fprintf(fpUpdateIncSeg, "ftMoveToDownstreamNode-b\n"); })
						DBG_PRINT(if (fpUpdateIncSeg != NULL) { ftPrintIncompleteSegment(incompleteData->incomplete[inflowIdx], 0, fpUpdateIncSeg, NULL); })
						incompleteData=incompleteData->prev;

					}
					if (volMoved != 0 && ftWithinTolerance(volMoved, 0, FT_ATOL, FT_RTOL)) {
						DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "MTDN4\t%.10g\n", volMoved); fflush(fpDeltas); })
						volMoved=0;
					}
				}

				if (volMoved > 0) {
					int nparts;
					SegOwner *owners;
					VolSplitData *volSplit;
					Pseg segMoved;
					Pseg replacementSeg=NULL;

					if (volMoved != seg->v && ftWithinTolerance(volMoved, seg->v, 1e-8, FT_RTOL)) {
						DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "MTDN7\t%.10g\n", volMoved-seg->v); fflush(fpDeltas); })
						volMoved=seg->v;
					}
					// if the volume moved is the same as the segment's volume, then there will only be one
					// part to "split", otherwise there will be 2
					if (volMoved==seg->v) {
						nparts=1;
					} else {
						nparts=2;
					}
					// set up the data structures to split the segment
					owners=(PSegOwner)calloc(nparts, sizeof(SegOwner));
					volSplit=(PVolSplitData)calloc(nparts, sizeof(VolSplitData));
					// the first part is the part that will be moved, so set it to be an INFLOW type for the downstream node
					FTInitInflowPartOwner(owners[0], downNode, inflowIdx, 1)
					volSplit[0].volume=volMoved;
					volSplit[0].adjustment=1.0;
					if (nparts==2) {
						// if there will be 2 parts, then the second part will be a LINK type for the link.
						volSplit[1].volume=seg->v-volMoved; // remaining part shouldn't be adjusted
						volSplit[1].adjustment=1;
					}
					if (IsIncompleteSeg(seg)) {
						// segment to move is an incomplete segment.
						PIncompleteSeg incSeg=(PIncompleteSeg)seg;
						if (nparts==2) {
							// set the owner info for the replacement seg
							FTInitCopyOwner(owners[1], &incSeg->owner, 2)
						}
						// split the segment
						ftSplitIncompleteSegment(incSeg, nparts, owners, volSplit);
						// this segment is now no type
						FTInitNoneOwnerType(incSeg->owner, incSeg->owner.source*100+3)
						// the segment moved is the first child created
						segMoved=(Pseg)&incSeg->children[0];
						if (nparts==2) {
							// and if there were 2 parts, the replacement is the second child created
							replacementSeg=(Pseg)&incSeg->children[1];
						}
					} else if (IsMergedSeg(seg)) {
						// segment to move is a merged segment
						PMergedSeg mergedSeg=(PMergedSeg)seg;
						// allocate the return array
						PMergedSeg *segs=(PMergedSeg*)calloc(nparts, sizeof(PMergedSeg));
						if (nparts==2) {
							// set the owner info for the replacement seg
							FTInitCopyOwner(owners[1], &mergedSeg->owner, 4)
//							FTInitLinkPartOwner(owners[1],mergedSeg->owner.object.linkIdx,4)
							// finish setting owner info
						}
						// split the segment
						ftSplitMergedSegment(mergedSeg, nparts, owners, volSplit, segs);
						// the segment moved is the first element
						segMoved=(Pseg)segs[0];
						if (nparts==2) {
							// and the replacement segment is the second
							replacementSeg=(Pseg)segs[1];
						}
						// release the array
						free(segs);
					}
					// release the data structures
					free(volSplit);
					free(owners);
					DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f\t%.8f", (nparts==1?0:replacementSeg->v), excessVolume); })

					DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f\t%.8f\n", segMoved->v, segMoved->c); })
					DBG_PRINT(if (fpTrace!=NULL && printTrace) { fprintf(fpTrace, "Moving volume\t%.10g\tc\t%.10g\n", segMoved->v, segMoved->c); fflush(fpTrace); })
					// add the new segment to thew inflowSegs of the downstream node
					ftAddNodeSeg(&downNode->inflowSegs.segs[inflowIdx], (Pseg)segMoved, 13);
					downNode->inflow.volumeFlowed[inflowIdx] += segMoved->v;
					// and add the downstream node to the list of nodes to process
					ftAddNodeToProcess(downNode, nodesToProcess);

					// nparts will be 1 if the entire segment was moved
					if (nparts==1) {
						if (linkIdx > Npipes && seg->prev == NULL) {
							// this is a valve or pump - and this is the only segment left
							// remove the incSeg from FirstSeg;
							FirstSeg[linkIdx]=NULL;
							LastSeg[linkIdx]=NULL;
							addseg(linkIdx, 0, 0);
						} else {
							// unlink this segment, but do not return it to the FreeSeg list
							ftUnlinkFirstLinkSegment(linkIdx);
						}
					} else {
						// replace the original segment with the replacement seg
						ftReplacePseg(seg, replacementSeg, &FirstSeg[linkIdx], &LastSeg[linkIdx]);
						if (IsMergedSeg(seg)) {
							// if the original segment was a merged segment,
							// setit's type to NONE
							PMergedSeg ms=(PMergedSeg)seg;
							FTInitNoneOwnerType(ms->owner, ms->owner.source*100+5)
						}
					}
//					downNode->inflow.volumeFlowed[inflowIdx] += volMoved;
					// if the excess volume is sufficiently close to zero, set it to zero.
					if (excessVolume > 0 && ftWithinTolerance(excessVolume, 0, 1e-8, FT_RTOL)) {
						DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "MTDN5\t%.10g\n", excessVolume); fflush(fpDeltas); })
						DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "MOVEDSINC-EXCVOL\t\t%.8e\n", excessVolume); })
						excessVolume=0;
					}
				} else {
					DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t\t\n"); })
				}
			}
		}
	}
}
/*
**--------------------------------------------------------------
**   Input:   ms: The merged segment to validate
**   Output:  1 if it id valid, 0 otherwise
**   Purpose: Validate the merged segment's volume.
**
**   Note: This function is not likely needed anymore, as it was
**         added to check a specific issue during develeopment.
**         Once the flow transport method is accepted, this method should be removed
**--------------------------------------------------------------
*/

int ftValidateMergedSegment(PMergedSeg ms, char *src) {
	double v1, v2;
	int rv;
	v1=ms->v;
	v2=ms->knownVol+ftComputeSegmentWaterVolume((Pseg)ms->incSegs);
	rv=ftWithinTolerance(v1, v2, FT_ATOL, FT_RTOL);
	if (rv==0) {
		fprintf(stdout, "Invalid merged segment! %s[%d] %s ms->v: %.8f  ms->knownVol+segvol: %.8f  diff:  fflush(stdout); %.8f\n", OwnerTypes[ms->owner.type], ms->owner.source, src, v1, v2, ABS(v2-v1));
	}
	return rv;
}
/*
**--------------------------------------------------------------
**   Input:   node: The node whose incomplete nodes need to be
**                  updated
**            inflowIdx: the index of the inflow link the segment
**                       arrived on
**            c: the concentration
**            volMoved: The volume moved.  This value will be
**                      updated to reflect the removal by this
**                      function call
**   Output:  1: if the any incomplete segment was completed, 0 otherwise
**   Purpose: Update a node's incomplete segments using the
**            concentration and volume moved
**--------------------------------------------------------------
*/
int ftUpdateIncompleteSegs(PNodeTransport node, int inflowIdx, double c, double *volMoved, int completing) {
	static unsigned long cnt=0;  //debug variable to help stop at a specific call
	DBG_PRINT(FILE *fpi=NULL;)
	PIncompleteSegmentData incompleteData;
	PNodeLinks in = &node->inflow;
	int complete=0;
	double adjVol=*volMoved;
    
	incompleteData=node->incompleteSegs;
    if(node->nodeIdx==2230) {
        printf("");
    }
    if(cnt==39) {
        printf("");
    }
	// INVESTIGATE:  I believe this is where we may be running into some problems
	// if an incomplete segment has very "uneven" parts, ftIsComplete may erroneously return true on one side, but not the other...
	while (incompleteData != NULL && adjVol > 0) {
		PIncompleteSegmentData prevInc=incompleteData->prev;
		if (incompleteData->incomplete[inflowIdx] != NULL) {
			Pseg tseg=createseg(adjVol, c);  // new segment
			double volUsed;
			// INVESTIGATE: keep part & msVol as they are currently used for debug checks
			PIncSegPart part;
            double msVol;

			DBG_PRINT(if (fpi!= NULL) { ftPrintIncompleteSegments(incompleteData->incomplete, in->nlinks, fpi); })
			// got the segment we've been waiting for...
			DBG_PRINT(if (fpTrace!=NULL && printTrace) { fprintf(fpTrace, "Received %s segment for an incomplete segment into node %s\n", getSegmentType(tseg), node->node->ID); fflush(fpTrace); })
			DBG_PRINT(if (incfp!=NULL) ftPrintIncompleteSegment(incompleteData->incomplete[inflowIdx], 0, incfp, NULL);)
			part = &(incompleteData->incomplete[inflowIdx]->parts[inflowIdx]);
			if(part->conc==-2) {
				printf("");
			}
            msVol = ftMergedSegmentVolume(part->segs);
            if(msVol > 0) {
                printf("");
            }
			// tell the incomplete segment to process the new segment.  volUsed is the amount of the segment that it used.
			volUsed=ftProcessIncomplete(incompleteData->incomplete[inflowIdx], &tseg, inflowIdx, completing?FALSE:TRUE, completing);
//			node->inflow.volumeFlowed[inflowIdx]+=volUsed;  // adjust the volume flowed
			adjVol-=volUsed;
			ftReleaseSeg(tseg);  //free the segment
			if (ftIncompleteDone(incompleteData, in->nlinks)) {
				// if it is complete, free the segments
				ftFreeIncompleteSegments(incompleteData, in->nlinks, &node->incompleteSegs);
				complete=1;
			}
			cnt++;
		}
		incompleteData=prevInc;
	}
#ifdef DEBUG
	if (adjVol != 0) {
		// just to check to make sure volMoved is correct in this case...
		fprintf(stdout, "");  fflush(stdout);
	}
#endif
	*volMoved = adjVol;  // update volMoved
	DBG_PRINT(if (fpi != NULL) { fclose(fpi); })
	return complete;
}
/*
**--------------------------------------------------------------
**   Input:   incompleteData: The incomplete segment data to free
**            nlinks: The number of inflow links
**            head: the head of the list of incomplete segment data
**                  for the node.  This may be modified on return
**   Output:  none
**   Purpose: Release the memory allocated for this incomplete
**            segment data
**--------------------------------------------------------------
*/
void ftFreeIncompleteSegments(PIncompleteSegmentData incompleteData, int nlinks, PIncompleteSegmentData *head) {
	int i;
    int rv=0;
	for (i=0; i<nlinks; i++) {
		if (incompleteData->incomplete[i] != NULL) {
			ftFreeIncompleteSegment(incompleteData->incomplete[i]);
			free(incompleteData->incomplete[i]);
			incompleteData->incomplete[i]=NULL;
			i=nlinks; // all incomplete->incomplete[i] point to the same memory
		}
	}
	free(incompleteData->incomplete);
	if (*head==incompleteData) { //this was the first one in the list
		*head=incompleteData->prev;
	} else {
		PIncompleteSegmentData t, p;
		t=NULL;
		p=*head;
		while (p!=incompleteData) {
			t=p;
			p=p->prev;
		}
		t->prev=incompleteData->prev;
	}
	free(incompleteData);
}
/*
**--------------------------------------------------------------
**   Input:   incompleteData: The incomplete segment data
**            nlinks: the number of inflow links
**   Output:  1 if it is now complete, 0 otherwise
**   Purpose: Determine if the incomplete segment is now complete
**--------------------------------------------------------------
*/
int ftIncompleteDone(PIncompleteSegmentData incompleteData, int nlinks) {
	int i=0;
	for (i=0; i<nlinks; i++) {
		if (!ftIsComplete(incompleteData->incomplete[i])) {
			return 0;
		}
	}
	return 1;
}
/*
**--------------------------------------------------------------
**   Input:   seg: The segment to move to the node's outflow links
**            tNode: the node
**   Output:  none
**   Purpose: Move the water segment to the outflow links based
**            on the outflow ratios
**--------------------------------------------------------------
*/
void ftMoveFlowToLinks(Pseg seg, PNodeTransport tNode, Pseg throughFlow) {
	int l;
	PNodeLinks out=&tNode->outflow;
	double vf;
	// move water segments to links

	if (throughFlow==NULL || throughFlow->v==0) {
		// there is no through flow.  This would only occur if an
		// incomplete or merged segment were being processed by ftMoveSegs
		return;
	}

	if (IsNormalSeg(throughFlow)) {
		// in this case, throughflow and seg are the same
		DBG_PRINT(if (fpi != NULL) { fprintf(fpi, "FLOW2LINKS\t"); })
		if (out->nlinks == 0) {
			volLost+=seg->v;
			massLost+=seg->c*seg->v;
			DBG_PRINT(if (fpLoss!=NULL) { fprintf(fpLoss, "%s\t\t%0.9f\t%0.9f\n", tNode->node->ID, seg->v, seg->c*seg->v); })
		}
		for (l=0; l<out->nlinks; l++) {
			// compute the volume leaving the node for the link
			double vl=seg->v * out->flowRatio[l];
			int linkIdx=out->linkIdxs[l];
			Pseg linkSeg=LastSeg[linkIdx];

			// if the current amount flowed out on this link + the new volume leaving the node
			// is sufficiently close to the expected flow, make it so.
			vf=out->volumeFlowed[l]+vl;
			if ((vf != out->flowVolume[l]) && ftWithinTolerance(vf, out->flowVolume[l], FT_ATOL, FT_RTOL)) {
				if (vl != vf-out->volumeFlowed[l]) {
					DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "MFTL1\t%.10g\n", vl-(vf-out->volumeFlowed[l])); fflush(fpDeltas); })
				}
				vl=vf-out->volumeFlowed[l];
			}
			DBG_PRINT(if (fpTrace!=NULL && printTrace) { fprintf(fpTrace, "moving vol\t%.10g\tc\t%.10g\tto link\t%s", vl, seg->c, Link[linkIdx].ID); fflush(fpTrace); })
			DBG_PRINT(if (fpi != NULL) { fprintf(fpi, "\t%.8f\t%.8f\t%.8f\t%.8f", vl, seg->c, linkSeg->v, linkSeg->c); })
			// if this segment and the first segment on the link should be combined (see ftCombinedSegments for details)
			// then combine them
			if (ftCombineSegments(vl, seg->c, linkSeg->v, linkSeg->c, 1*100000+linkIdx)) {
				linkSeg->c = (linkSeg->c*linkSeg->v + seg->c*vl) / (linkSeg->v + vl);
				checkConc(linkSeg->c);
				linkSeg->v += vl;
				DBG_PRINT(if (fpTrace!=NULL && printTrace) { fprintf(fpTrace, "\tadded vol to current segment. vol\tc\t%f\t%f\n", linkSeg->v, linkSeg->c); fflush(fpTrace); })
				DBG_PRINT(if (fpi != NULL) { fprintf(fpi, "\t%s\t%.8f\t%.8f\t", "ADD", linkSeg->v, linkSeg->c); })
			} else {
				/* Otherwise add a new seg to end of link */
				checkConc(seg->c);
				addseg(linkIdx, vl, seg->c);
				DBG_PRINT(if (fpTrace!=NULL && printTrace) { fprintf(fpTrace, "\tadded new segment.\n"); fflush(fpTrace); })
				DBG_PRINT(if (fpi != NULL) { fprintf(fpi, "\t%s\t%.8f\t%.8f\t", "NEW", vl, seg->c); })
			}
			// update volumeFlowed on the outflow link
			out->volumeFlowed[l]+=vl;
		}
		DBG_PRINT(if (fpi != NULL) { fprintf(fpi, "\n"); fflush(fpi); });
	} else {
		// if it is either an incomplete or merged segment
		int l;
		// allocate the data for splitting the segment
		SegOwner *owners=(PSegOwner)calloc(out->nlinks, sizeof(SegOwner));
		VolSplitData *volSplit=(VolSplitData *)calloc(out->nlinks, sizeof(VolSplitData));
		Pseg *segsMoved=(Pseg*)calloc(out->nlinks, sizeof(Pseg));

		DBG_PRINT(if (fpi != NULL) { fprintf(fpi, "FLOW2LINKS-%s\t", IsIncompleteSeg(throughFlow)?"INC":"MERGE"); });

		for (l=0; l<out->nlinks; l++) {
			// throughFlow contains the actual amount of water flowed through the node
			double vl=throughFlow->v * out->flowRatio[l];
			DBG_PRINT(if (fpTrace!=NULL && printTrace) { fprintf(fpTrace, "moving %s vol\t%f\tc\t%f\tto link\t%s\n", IsIncompleteSeg(throughFlow)?"incomplete":"merged", vl, throughFlow->c, Link[out->linkIdxs[l]].ID); fflush(fpTrace); });
			// if the current amount flowed out on this link + the new volume leaving the node
			// is sufficiently close to the expected flow, make it so.
			vf=out->volumeFlowed[l]+vl;
			if (vf != out->flowVolume[l] && ftWithinTolerance(vf, out->flowVolume[l], FT_ATOL, FT_RTOL)) {
				if (vl != vf-out->volumeFlowed[l]) {
					DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "MFTL2\t%.10g\n", vl-(vf-out->volumeFlowed[l])); fflush(fpDeltas); });
				}
				vl=vf-out->volumeFlowed[l];
			}
			FTInitLinkPartOwner(owners[l], out->linkIdxs[l], 6)
			volSplit[l].volume=vl;
			volSplit[l].adjustment=1.0;
			out->volumeFlowed[l] += vl;
		}
		if (IsIncompleteSeg(throughFlow)) {
			PIncompleteSeg is=(PIncompleteSeg)throughFlow;
			// slpit the incomplete segment and set the segsMoved data
			ftSplitIncompleteSegment(is, out->nlinks, owners, volSplit);
			for (l=0; l<out->nlinks; l++) {
				segsMoved[l]=(Pseg)&is->children[l];
			}
		} else {
			PMergedSeg ms=(PMergedSeg)throughFlow;
			PMergedSeg *segs=(PMergedSeg*)calloc(out->nlinks, sizeof(PMergedSeg));
			// slpit the merged segment and set the segsMoved data
			ftSplitMergedSegment(ms, out->nlinks, owners, volSplit, segs);
			for (l=0; l<out->nlinks; l++) {
				segsMoved[l]=(Pseg)segs[l];
			}
			free(segs);
		}
		for (l=0; l<out->nlinks; l++) {
			// for each split segment, insert it into the list for the outflow link
			int linkIdx=out->linkIdxs[l];
			insertseg(linkIdx, segsMoved[l]);
			DBG_PRINT(if (fpi != NULL) { fprintf(fpi, "\t%.8f\t%.8f\t", segsMoved[l]->v, segsMoved[l]->c); })
		}
		DBG_PRINT(if (fpi != NULL) { fprintf(fpi, "\n"); fflush(fpi); })
		free(owners);
		free(volSplit);
		free(segsMoved);
	}
}
/*
**--------------------------------------------------------------
**   Input:   v1: The first segment's volume
**            c1: The first segment's concentration
**            v2: The second segment's volume
**            c2: The second segment's concentration
**            loc: A debugging aid to identify where this call came from
**   Output:  1 if the two segments should be combined, 0 if not
**   Purpose: Determine if two segments should be combined based
**            on the following rules:
**            - if both concentrtations are 0, then combine
**            - if either concentration is zero (but not both) then do not combine
**            - if the absolute conc difference is less than Ctol, then combine (original rule)
**            - if the sum of the two volumes is less than FT_VTOL, then combine
**               the idea behind this is that if the volumes are sufficiently small (currently < 1 cf)
**               combining them will have no significant impact on the overall results and
**               will greatly reduce the nmber of segments in the network
**            - if the absolute difference between the potential new concentration and the concentration of the larger
**              of the two segments is less than F_CTOL, then combine them.
**               the idea behind this is that if combining two segments does not
**               significantly alter the concentration (FT_CTOL),
**               then it will have negligible outcome on the overall results
**               and reduce the number of segments in the network
**--------------------------------------------------------------
*/

// INVESTIGATE:
// need to look at alternative methods to determine if two segments can be combined
// how it is done in quality.c currently is
// compare c1 & c2 and if the difference is less than CTOL, then combine, otherwise don't
// This leads to the possibility of large segments of water  with extremely small (~e-20 or smaller) concetrations
// the current method here is to combine if both concentrations are 0,
// do not combine if one is 0 and the other is not 0 (this deals with the issue above)
// if c1 and c2 are within Ctol, then combine.
//#define _FT_USE_CTOL_
int ftCombineSegments(double v1, double c1, double v2, double c2, int loc) {
	static int seq=0;
#ifdef _FT_USE_CTOL_
	return ABS(c1 - c2) < Ctol;
#else
	int combine=0;
	if (c1==0 && c2==0) return 1;
	DBG_PRINT(if (fpcombine!= NULL) { fprintf(fpcombine, "%d\t%d\t%g\t%g\t%g\t%g", seq++, loc, (v1>v2?v1:v2), (v1>v2?c1:c2), (v1>v2?v2:v1), (v1>v2?c2:c1)); })
	if (c1==0 || c2==0) {
		int tcombine;
		double bigc, newc;
		newc=(v1*c1 + v2*c2) / (v1+v2);
		if (!v1 > v2) {
			bigc=c1;
		} else {
			bigc=c2;
		}
		tcombine=ABS(newc-bigc) < FT_CTOL;

		DBG_PRINT(if (fpcombine!= NULL) { fprintf(fpcombine, "\t1\t\t\t%g\t%g\t%d\n", newc, ABS(newc-bigc), tcombine?10:0); })

		return 0;
	}
	if (ABS(c1 - c2) < Ctol) {
		DBG_PRINT(if (fpcombine!=NULL) { fprintf(fpcombine, "\t2\t\t%g\t\t\t%d\n", ABS(c1 - c2), 1); })
		return 1;
//	} else {
//		DBG_PRINT(if (fpcombine!=NULL) { fprintf(fpcombine, "\t5\t\t\t\t\t%d\n", 0); })
//		return 0;
	}
//	if (v1+v2 < FT_VTOL) {  // if the sum of the two volumes is sufficiently small, combine them regardless of their concentrations
//		combine=1;
//		DBG_PRINT(if (fpcombine!=NULL) { fprintf(fpcombine, "\t3\t%g\t\t\t\t%d\n", v1+v2, combine); })
//	}
	if (!combine) {
		// compute what the concentration *would* be if the two were combined
		// and check it against the concentration of the bigger of the two segments
		// if it is less than Ctol, combine them
		double bigc, newc;
		newc=(v1*c1 + v2*c2) / (v1+v2);
		if (!v1 > v2) {
			bigc=c1;
		} else {
			bigc=c2;
		}
		combine=ABS(newc-bigc) < FT_CTOL;
		DBG_PRINT(if (fpcombine!= NULL) { fprintf(fpcombine, "\t4\t\t\t%g\t%g\t%d\n", newc, ABS(newc-bigc), combine); })
	}
	return combine;
#endif
}
/*
**--------------------------------------------------------------
**   Input:   seg: The segment being flowed
**            tNode: the node being flowed through
**            throughFlow: the through flow if seg is not a normal seg
**   Output:  none
**   Purpose: Add the through flow to the node's through flow segment list
**--------------------------------------------------------------
*/
void ftAddThroughFlow(Pseg seg, PNodeTransport tNode, Pseg throughFlow) {
	if (throughFlow == NULL || throughFlow->v ==0)
		return;
	if (IsNormalSeg(seg)) {
		// create and add a new segment to the node's throughflow list
		ftAddNodeSeg(&tNode->throughSegs, createseg(seg->v, seg->c), 14);
	} else {
		// change the type to THROUGH and set the owner object accordingly
		if (IsIncompleteSeg(seg)) {
			PIncompleteSeg is=(PIncompleteSeg)throughFlow;
			FTInitThroughOwnerType(is->owner, tNode, is->owner.source*100+7)
		} else if (IsMergedSeg(seg)) {
			// this segment is a merged segment
			PMergedSeg ms=(PMergedSeg)throughFlow;
			FTInitThroughOwnerType(ms->owner, tNode, ms->owner.source*100+8)
		}
		// and add the throughflow segment to the node's throughflow list
		ftAddNodeSeg(&tNode->throughSegs, (Pseg)throughFlow, 15);
	}
}
/*
**--------------------------------------------------------------
**   Input:   seg: The segment being flowed
**            tNode: the node being flowed through
**            throughFlow: the through flow seg
**   Output:  the through flow segment
**   Purpose: Remove any demand from the segment and add it to
**            the node's demand segments
**--------------------------------------------------------------
*/
Pseg ftRemoveDemand(Pseg seg, PNodeTransport tNode, Pseg throughFlow) {
	// the amount of water removed from this segment for demand is
	// computed as (seg->v/node inflowvolume) * nodal demand
	double volRatio=seg->v / (tNode->totalVolInflow*tNode->outflowAdjustment);
	double demFromThis;

	DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "DEMAND\t%.8f", volRatio); })
#ifdef ENABLE_DBG_PRINT
	if (fpi!=NULL) { fprintf(fpi, "\t%.8f", volRatio); }
	if (volRatio>1) {
		if (fpTrace!=NULL && printTrace) {
			fprintf(fpTrace, "Demand Ratio > 1:\t%.14g\n", volRatio); fflush(fpTrace);
		}
	}
#endif
	demFromThis=tNode->demand*volRatio;
	DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f\t%.8f", demFromThis, seg->v); })
	if (ABS(demFromThis-seg->v) > 0 && ABS(demFromThis-seg->v) < 1e-6) {
#ifdef DEBUG
		fprintf(stdout, "");  fflush(stdout);
#endif
		demFromThis=seg->v;
	}
	if (IsNormalSeg(seg)) {
		// if it is a normal segment, create a new segment and add it to the node's demand segs
		ftAddNodeSeg(&tNode->demandSegs, createseg(demFromThis, seg->c), 16);
		// and decrease the amount left in seg
		seg->v-=demFromThis;
		DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f", seg->v); })
	} else {
		// it is either an incomplete node or merged node and needs to be split
		// into SO_DEMAND and SO_INFLOW parts
		PSegOwner owners;
		PVolSplitData volSplit;
		int nparts;
		if (seg->v == demFromThis) {
			nparts=1;
		} else {
			nparts=2;
		}
		owners=(PSegOwner)calloc(nparts, sizeof(SegOwner));
		volSplit=(PVolSplitData)calloc(nparts, sizeof(VolSplitData));
		FTInitDemandPartOwner(owners[0], tNode, 9)
		volSplit[0].volume=demFromThis;
		volSplit[0].adjustment=1.0;
		if (nparts == 2) {
			volSplit[1].volume=seg->v-demFromThis;
			volSplit[1].adjustment=1.0;
		}
		DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f", nparts==1?0:volSplit[1].volume); })
		if (IsIncompleteSeg(seg)) {
			// this is an incomplete segment.  Split it into parts and put the demand "part" on the demandSegmens
			PIncompleteSeg srcSeg=(PIncompleteSeg)seg;
			if (nparts == 2) {
				FTInitInflowPartOwner(owners[1], srcSeg->owner.object.node, srcSeg->owner.data.inflowIdx, 10)
			}
			ftSplitIncompleteSegment(srcSeg, nparts, owners, volSplit);

			FTInitNoneOwnerType(srcSeg->owner, srcSeg->owner.source*100+11)
			// add the demand seg to the demandSegs list
			ftAddNodeSeg(&tNode->demandSegs, (Pseg)&srcSeg->children[0], 17);
			if (nparts == 2) {
				// and set throughFlow to the remainder
				throughFlow=(Pseg)&srcSeg->children[1];
			} else {
				// unless there isn't one, then set it to null
				throughFlow=NULL;
			}
		} else if (IsMergedSeg(seg)) {
			// this segment is a merged segment
			PMergedSeg ms=(PMergedSeg)seg;
			PMergedSeg *segs=(PMergedSeg*)calloc(nparts, sizeof(PMergedSeg));
			if (nparts == 2) {
				FTInitInflowPartOwner(owners[1], ms->owner.object.node, ms->owner.data.inflowIdx, 12)
			}
			ftSplitMergedSegment(ms, nparts, owners, volSplit, segs);
			// add the demand seg to demandSegs list
			ftAddNodeSeg(&tNode->demandSegs, (Pseg)segs[0], 18);
			if (nparts==2) {
				// and set throughFLow to the remainder
				throughFlow=(Pseg)segs[1];
			} else {
				// unless there isn't one, then set it to null
				throughFlow=NULL;
			}
			free(segs);
			FTInitNoneOwnerType(ms->owner, ms->owner.source*100+30)

		}
		free(owners);
		free(volSplit);
	}
	DBG_PRINT(if (fpTrace!=NULL && printTrace) {
		fprintf(fpTrace, "Total demand:\t%.10g\tDemand from this seg\t%.10g\tmass from this\t%.10g\ttotal vol out\t%.10g\ttotal mass out\t%.10g\n",
			tNode->demand, demFromThis, seg->c*demFromThis/*,tNode->volout,tNode->massout*/, 0, 0); fflush(fpTrace);
	})
	DBG_PRINT(if (fpi != NULL) { fprintf(fpi, "\n"); fflush(fpi); })
	return throughFlow;
}
/*
**--------------------------------------------------------------
**   Input:   incSeg: An incomplete segment
**   Output:  1 if the segment is now complete, 0 if not
**   Purpose: Determine if the segment is now complete
**--------------------------------------------------------------
*/
int ftIsComplete(PIncompleteSeg incSeg) {
	/* if all of incSeg's parts have no more segs and the volume and flowed are equal, then this segment is "complete" */
	int i;
	if (incSeg==NULL) return 1;
	for (i=0; i<incSeg->nparts; i++) {
		PIncSegPart part=&incSeg->parts[i];
		if (part->segs != NULL)
			return 0;
		if (part->conc < 0)
			return 0;
		if (part->volume != part->flowed) {
			if (!ftWithinTolerance(part->volume, part->flowed, FT_ATOL, FT_RTOL)) {
				return 0;
			} else {
				DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "IC\t%.10g\n", part->volume-part->flowed); fflush(fpDeltas); })
			}
		}
	}
	return 1;
}
/*
**--------------------------------------------------------------
**   Input:   incSeg: The incomplete segment to release memory for
**   Output:  none
**   Purpose: Release the memory allocated for this incomplete segment
**--------------------------------------------------------------
*/
void ftFreeIncompleteSegment(PIncompleteSeg incSeg) {
	int i;
	for (i=0; i<incSeg->nchildren; i++) {
		ftFreeIncompleteSegment(&incSeg->children[i]);
	}
	free(incSeg->children);
	free(incSeg->parts);
}

/*
**--------------------------------------------------------------
**   Input:   incSeg: The incomplete segment being processed
**            seg: The seg which triggered the processing
**            idx: the inflow index
**   Output:  The amount of seg's volume that was used.
**   Purpose: to update all the parts and children with the now
**            known volume and concentration coming in on inflow link idx
**--------------------------------------------------------------
*/
double ftProcessIncomplete(PIncompleteSeg incSeg, Pseg *segRef, int idx, int applyAdjustment, int completing) {
	int i;
	double notused=0, used=0, adjUsed=0;
	double needed;
	double mass;
	Pseg seg=*segRef;
	double vol=seg->v;
	double adjVol;
	int seq;
	PIncSegPart part=&incSeg->parts[idx];
	double partTotalVol;
	double partCompleteVol;

	seq=incSeg->owner.sequence;
	if (traceProcessIncomplete) {
		if (seq == 506) {
			printf("");
		}
		fprintf(stdout, "Processing incSeg: %d\n", incSeg->owner.sequence);
	}
	DBG_PRINT(
	if (fpAdj != NULL) {
		fprintf(fpAdj, "INC\t%lld\t%d\t%0.10g\t%0.10g\t%0.10g\t%0.10g\t%0.10g\t%d",
			incSeg->owner.source, incSeg->owner.sequence, vol, vol*incSeg->adjustment, incSeg->adjustment,
			part->volume-part->flowed-vol, part->volume-part->flowed-vol*incSeg->adjustment,
			applyAdjustment);
		ftPrintSource(fpAdj, incSeg->owner.source);
		fprintf(fpAdj, "\n");
		fflush(fpAdj);
	}
	)
	// apply the adjustment if necessary
	if (applyAdjustment) {
		adjVol=vol*incSeg->adjustment;
	} else {
		adjVol=vol;
	}

	DBG_PRINT(
	if (fpIncTrace!=NULL) {
		PSegOwner o=&incSeg->owner;
		fprintf(fpIncTrace, "PI\t%ld\t%s\t%s\t%s\t%0.14g\t%d\t%0.14g\t%0.14g\n",
			o->sequence, OwnerTypes[o->type], ftGetOwnerObject(o), ftGetOwnerData(o), incSeg->v, idx, part->volume, part->flowed);
		fflush(fpIncTrace);
	}
	)
	// note that this is being processed
	incSeg->processing++;
#ifdef DEBUG
	if (incSeg->processing > 1) {
		fprintf(stdout, "");  fflush(stdout);
	}
#endif
	// how much volume is needed to complete the part
	partTotalVol=ftComputeSegmentWaterVolume(part->segs);
	partCompleteVol=ftComputeCompleteSegmentWaterVolume(part->segs);
	if(partTotalVol != partCompleteVol) {
		printf("");
	}
	// INVESTIGATE:
	// this was apprently added at some point to fix some issues, but it may be causing other problems.
	if(!completing) {
		needed=part->volume-part->flowed-partTotalVol;
	} else {
		needed=part->volume-part->flowed-partCompleteVol;
	}
	// if we don't need anything from this segment, then decrement processing flag and return
	if (needed==0) {
		incSeg->processing--;
		return 0;
	}
    if(incSeg->owner.sequence==351) {
        printf("");
    }
	// if there is (likely) excess volume
    /*****
	INVESTIGATE:
    figure out a way to possibly pass in whether or not the parent segment will be completed.
    There is an issue if one of the ratios is very close to 1.0, the remaining children may 
    be thought os as complete even though they aren't really. 
    *****/
//	if(!completing) {
		if (adjVol > (needed-FT_ATOL)) {
			// if adjVOl and needed are sufficiently close
			if (adjVol != needed && ftWithinTolerance(adjVol, needed, FT_ATOL, FT_RTOL)) {
				// set used to be the entire adjVol and notused to 0
				DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "PI1\t%.10g\n", adjVol-needed); fflush(fpDeltas); })
				notused=0;
				used=needed;
				adjUsed=adjVol;
			} else {
				// therwise set used to be what is needed and notused to be what is left
				notused = adjVol-needed;
				if (applyAdjustment) {
					notused /=incSeg->adjustment;
				}
				if (adjVol==needed) {
					notused=0;
				} else {
					notused=vol-(applyAdjustment?needed/incSeg->adjustment:needed);
				}
				used=needed;
				adjUsed=needed;
			}
		} else {
			// there is no excess volume
			notused=0;
        	used=adjVol;
			adjUsed=adjVol;
		}
//	} else {
//		// if we are completing incomplete segments, use the exact values.
//		notused=0;
//		used=needed;
//		adjUsed=adjVol;
//	}
	if (IsNormalSeg(seg)) {
		static unsigned long cnt=0;
		unsigned long lcnt=cnt;
		// the segment that came in was a normal segment
		double c=seg->c;
		cnt++;
        if(cnt>=320) {
            printf("");
        }
		if(part->conc==-2) {
			printf("");
		}
		if (part->conc==-1) {
			// signify that some normal segment processing has begun
			part->conc=0;
		}
		// compute total new mass of part
		mass=part->flowed*part->conc + used*c;
		// increment the part's flowed amount
		part->flowed+=used;
		// and now compute the part's concentrtation
		part->conc=mass/part->flowed;
		// now propogate this segment to this segs children
		if (incSeg->processing <= 1) {
            if(incSeg->nchildren>0) {
				// INVESTIGATE: see above...
				int willBeComplete = ftIsComplete(incSeg); /////////
                int anyChildrenComplete=0,allChildrenComplete=1;
			    for (i=0; i<incSeg->nchildren; i++) {
				    double tvol;
				    // create a new segment with volume based on the ratio of the child to this seg
				    Pseg tseg=createseg(adjVol * incSeg->children[i].ratio, c);
				    // and process that incomplete segment
				    tvol=ftProcessIncomplete(&incSeg->children[i], &tseg, idx, TRUE /* ?? or TRUE ?? */, completing);
                    if(ftIsComplete(&incSeg->children[i])) {
                        anyChildrenComplete = 1;
                    } else {
                        allChildrenComplete = 0;
                    }
			    }
                if(anyChildrenComplete && !allChildrenComplete) {
                    printf("");
					if (ftIsComplete(incSeg)) {
						printf("");
					}
                }
            }
			// if this segment is now complete
			if (ftIsComplete(incSeg)) {
				// compute the total mass in this segment
				mass=incSeg->massadded;
				for (i=0; i<incSeg->nparts; i++) {
					part=&incSeg->parts[i];
					mass += part->volume*part->conc;
				}
				// and properly update the "owner" 
				switch (incSeg->owner.type) {
				case SO_NONE:
					break;
				case SO_THROUGH:
				{
				   // for through segment types, replace the incSeg (will be in throughSegs list)
					// with a new normal seg with the correct volume and concentration
					PNodeTransport node=incSeg->owner.object.node;
					ftReplaceSeg((Pseg)incSeg, createseg(incSeg->v, mass/incSeg->v), &node->throughSegs, 23);
				}
					break;
				case SO_INFLOW:
				{
					// for inflow segment types, replace the incSeg (will be in the inflow segs list)
					// with a new normal segment with the correct volume and concentration
					PNodeTransport node=incSeg->owner.object.node;
					int inflowIdx=incSeg->owner.data.inflowIdx;
					ftReplaceSeg((Pseg)incSeg, createseg(incSeg->v, mass/incSeg->v), &node->inflowSegs.segs[inflowIdx], 24);
				}
					break;
				case SO_DEMAND:
				{
					// for demand segment types, replace the incSeg (will be in the demandSegs segs list)
					// with a new normal segment with the correct volume and concentration
					PNodeTransport node=incSeg->owner.object.node;
					ftReplaceSeg((Pseg)incSeg, createseg(incSeg->v, mass/incSeg->v), &node->demandSegs, 25);
				}
					break;
				case SO_LINK:
				{
					// for link segment types, replace the incSeg (will be in the FirstSeg and LastSeg segs list)
					// with a new normal segment with the correct volume and concentration
					int linkIdx = incSeg->owner.object.linkIdx;
					ftReplaceLinkSeg((Pseg)incSeg, createseg(incSeg->v, mass/incSeg->v), linkIdx, 21);
				}
					break;
				case SO_MERGED:
				{
					PMergedSeg ms = incSeg->owner.object.mergedSeg;
					// get the merged segment and if it is null, that means that it has already been processed
					// so don't do it again
					if (ms != NULL) {
						// for merged segment types, call updateMergedSegment with the 
						// inc segment and total mass
						ftUpdateMergedSegment(ms, incSeg, mass, completing);
						// added back above function call on 11/29/2016.  It had been deleted and not sure if it was intentional
						if (ms->incSegs == NULL) {
							free(ms);
							incSeg->owner.object.mergedSeg=NULL;
						}
						// INVETIGATE: is this necessary or not...
////					free(ms);
					}
				}
					break;
				case SO_INCOMPLETE_PART:
				{
					// incSeg is owned by another incomplete segment
					PSegOwner owner=&incSeg->owner;
					PIncSegPart incpart;
					PIncompleteSeg incpartOwner;
					int completed;

#ifdef ENABLE_DBG_PRINT
#if 0
										   if (fpUpdateIncSeg != NULL) {
											   PIncompleteSeg pi=incSeg;

											   PIncompleteSegmentData p;
											   while (pi->node == NULL) pi=pi->parent;
											   p=pi->node->incompleteSegs;
											   fprintf(fpUpdateIncSeg, "in ftProcessIncomplete\n");
											   while (p != NULL) {
												   ftPrintIncompleteSegments(p->incomplete, pi->node->inflow.nlinks, fpUpdateIncSeg);
												   p=p->prev;
											   }
										   }
#endif
#endif
					incpart = owner->object.incSeg;
					incpartOwner=owner->data.incSeg;
					// call updatePartSegment to handle it.
					completed=ftUpdatePartSegment(incpart, incpartOwner, (Pseg)incSeg, mass,completing);
					if (completed) {
						FTInitNoneOwnerType(incSeg->owner, incSeg->owner.source*100+33)
					} else {
//							fprintf(stdout,"Not Completed\n");  fflush(stdout);
					}
				}
					break;
				default:
				{
				   fprintf(stdout, "Unhandled incomplete owner type: %s  fflush(stdout); (%d)\n", OwnerTypes[incSeg->owner.type], incSeg->owner.type);
				}
					break;
				}
				cnt++;
// INVESTIGATE: needed?
if(completing) {
//ftCheckForInvalidSegments();
}
			}
#ifdef DEBUG
		} else {
			// if processing > 1
			fprintf(stdout, "");  fflush(stdout);
#endif
		}
	} else {
		// adding a merged or incomplete segment...
		int nparts = notused>0 ? 2 : 1;
		SegOwner *owners=(PSegOwner)calloc(nparts, sizeof(SegOwner));
		VolSplitData *volSplit=(VolSplitData*)calloc(nparts, sizeof(VolSplitData));
		PSegOwner owner;
		Pseg replacementSeg=NULL;

		if (IsIncompleteSeg(seg)) { // "incomplete" segment
			PIncompleteSeg pIncSeg=(PIncompleteSeg)seg;
			owner=&pIncSeg->owner;
		} else if (IsMergedSeg(seg)) { // "Merged" segment
			PMergedSeg mergedSeg=(PMergedSeg)seg;
			owner=&mergedSeg->owner;
		}
		FTInitCopyOwner(owners[0], owner, 31)
			volSplit[0].volume=used;
		volSplit[0].adjustment=1.0;

		if (nparts == 2) {
			FTInitCopyOwner(owners[1], owner, 32)
				volSplit[1].volume=notused;
			volSplit[1].adjustment=1.0;
		}

		if (IsIncompleteSeg(seg)) { // "incomplete" segment
			PIncompleteSeg pIncSeg=(PIncompleteSeg)seg;
			ftSplitIncompleteSegment(pIncSeg, nparts, owners, volSplit);
			seg=(Pseg)&pIncSeg->children[0];
			//				FTInitIncompletePartOwner(pIncSeg->owner,part,incSeg,pIncSeg->owner.source*100+41)
			if (nparts==2) {
				*segRef=(Pseg)&pIncSeg->children[1];
			} else {
				*segRef=NULL;
			}
		} else if (IsMergedSeg(seg)) { // "Merged" segment
			PMergedSeg *msegs=(PMergedSeg*)calloc(nparts, sizeof(PMergedSeg));
			ftSplitMergedSegment((PMergedSeg)seg, nparts, owners, volSplit, msegs);
			seg=(Pseg)msegs[0];
			if (nparts==2) {
				*segRef=(Pseg)msegs[1];
			} else {
				*segRef=NULL;
			}
			free(msegs);
		}
		// now that it has been split, set this type to NONE because nobody owns it anymore
		FTInitNoneOwnerType(*owner, owner->source*100+21)
		free(volSplit);
		free(owners);

		// set the ownership appropriately
		if (IsIncompleteSeg(seg)) {
			PIncompleteSeg pIncSeg=(PIncompleteSeg)seg;
			FTInitIncompletePartOwner(pIncSeg->owner, part, incSeg, pIncSeg->owner.source*100+14)
		} else if (IsMergedSeg(seg)) {
			PMergedSeg mergedSeg=(PMergedSeg)seg;
			FTInitIncompletePartOwner(mergedSeg->owner, part, incSeg, mergedSeg->owner.source*100+15)
		}
		// and add it to the part's seg list
		seg->prev=part->segs;
		part->segs=seg;
	}
	incSeg->processing--;
	return vol-notused;
	//	return applyAdjustment ? adjUsed/incSeg->adjustment : adjUsed;
}
/*
**--------------------------------------------------------------
**   Input:   part: The part of an incomplete segment to update
**            partOwner: the incomplete segment to which this part belongs
**            seg: the incomplete segment that is now complete
**            mass: the mass
**   Output:  none
**   Purpose: Update the incomplete segment that owns an incomplete segmnent that was just made complete
**--------------------------------------------------------------
*/
int ftUpdatePartSegment(PIncSegPart part, PIncompleteSeg partOwner, Pseg seg, double mass, int completing) {
	int completed=0;
	int n=0, c=0;
	static unsigned long cnt=0;
	// create a new segment with the seg's volume and appropriate concetration
	Pseg newSeg=createseg(seg->v, mass/seg->v);
	PNodeTransport node=partOwner->node;
	int inflowIdx=(int)(part-&partOwner->parts[0]);
	PIncompleteSegmentData incompleteData;

	cnt++;
	DBG_PRINT(
	if (fpIncTrace!=NULL) {
		PSegOwner o=&partOwner->owner;
		fprintf(fpIncTrace, "UPS\t%ld\t%s\t%s\t%s\t%0.14g\t%d\t%0.14g\t%0.14g\n",
			o->sequence, OwnerTypes[o->type], ftGetOwnerObject(o), ftGetOwnerData(o), seg->v, inflowIdx,
			part->volume, part->flowed);
		fflush(fpIncTrace);
	}
	)

	// remove the seg from the part's seg list
	ftRemoveSegFromList(seg,(Pseg *)&part->segs);
	// if the partOwner is not complete, process it 
	// it is possible for the partOwner to be already complete at this
	// time, so it doesn't need to be processed again.
	if (!ftIsComplete(partOwner)) {
#ifdef ENABLE_DBG_PRINT
		if (fpUpdateIncSeg != NULL) {
			PIncompleteSegmentData p=partOwner->node->incompleteSegs;
			fprintf(fpUpdateIncSeg, "in ftUpdatePartSegment: cnt=%d\n",cnt);
			while (p != NULL) {
				ftPrintIncompleteSegments(p->incomplete, partOwner->node->inflow.nlinks, fpUpdateIncSeg);
				p=p->prev;
			}
		}
#endif
// INVESTIGATE:
//		if(!partOwner->processing) {
			ftProcessIncomplete(partOwner,&newSeg,inflowIdx,FALSE, completing);
//		}
	} else {
		// part owner is complete
		ftProcessIncomplete(partOwner, &newSeg, inflowIdx, FALSE, completing);
	}
	// now free any of the node's incomplete segments that are now complete
	incompleteData=node->incompleteSegs;
	while (incompleteData != NULL) {
		PIncompleteSegmentData prevInc=incompleteData->prev;
		if (incompleteData->incomplete[inflowIdx] == partOwner) {
			n++;
			// this is the segment we are updating
			if (ftIncompleteDone(incompleteData, node->inflow.nlinks)) {
				DBG_PRINT(if (fpUpdateIncSeg != NULL) { fprintf(fpUpdateIncSeg, "ftUpdatePartSegment: 2\n"); })
/////			DBG_PRINT( if(fpUpdateIncSeg != NULL) { ftPrintIncompleteSegment(incompleteData->incomplete[inflowIdx],partOwner->node->inflow.nlinks,fpUpdateIncSeg, NULL); })
if(completing) {
	printf("");
}
				if (!partOwner->processing) {
// INVESTIGATE:
					// ONLY if we are not already updating it.  This can happen if the incomplete
					// segment follows a loop back around to itself.  This would likely indicate a
					// "bad" hydraulic solution, but this case must be handled.

					ftFreeIncompleteSegments(incompleteData, node->inflow.nlinks, &node->incompleteSegs);
////					ftCheckForInvalidSegments();
#ifdef DEBUG
				} else {
					fprintf(stdout, "");  fflush(stdout);
#endif
				}
				completed=1;
				c++;
			}
		}
		incompleteData=prevInc;
	}
	if (n>1) {
		fprintf(stdout, "Part Owner matched more than one!\n");  fflush(stdout);
	}
	return completed;
}
/*
**--------------------------------------------------------------
**   Input:   ms:The merged segment to update
**            incSeg: the incomplete segment that is now complete
**            mass: the mass in the incSeg
**   Output:  none
**   Purpose: Update the merged segment with the now known volume and mass
**--------------------------------------------------------------
*/
void ftUpdateMergedSegment(PMergedSeg ms, PIncompleteSeg incSeg, double mass, int completing) {
#ifdef ENABLE_DBG_PRINT
	if (fpAdj != NULL) {
		fprintf(fpAdj, "MERGE\t%lld\t%d\t%0.10g\t%0.10g\t%0.10g\t%0.10g\t%0.10g\t%d",
			ms->owner.source, ms->owner.sequence, incSeg->v, incSeg->v*ms->adjustment, ms->adjustment,
			ms->v-ms->knownVol-incSeg->v, ms->v-ms->knownVol-incSeg->v*ms->adjustment, 0);
		ftPrintSource(fpAdj, ms->owner.source);
		fprintf(fpAdj, "\n");
		fflush(fpAdj);
	}
#endif
    if(ms->owner.sequence==19851) {
        printf("");
    }
	DBG_PRINT(
	if (fpIncTrace!=NULL) {
		PSegOwner o=&ms->owner;
		fprintf(fpIncTrace, "UMS\t%ld\t%s\t%s\t%s\t%0.14g\t\t%0.14g\t%0.14g\n",
			o->sequence, OwnerTypes[o->type], ftGetOwnerObject(o), ftGetOwnerData(o), incSeg->v,
			ms->v, ms->knownVol);
		fflush(fpIncTrace);
	}
	)
	// check for validity - may be able to remove once there is enough confidence that this won't occur
	ftValidateMergedSegment(ms, "UMS1");
	// one of the incomplete segments is now complete.
	if (ms->knownVol+incSeg->v-FT_ATOL > ms->v) {
		fprintf(stdout, "Known Vol + seg vol > merged seg expected vol!\n");  fflush(stdout);
	}
	// update the mass and volume that is "known" (i.e. that came from normal segments or that is known like this mass)
	ms->knownMass+=mass;
	ms->knownVol+=incSeg->v;
	// remove the incomplete segment from the merged segment's incomplete segments list
	ftRemoveSegFromList((Pseg)incSeg, (Pseg *)&ms->incSegs);
	// check for validity again
	ftValidateMergedSegment(ms, "UMS2");
	// if there are no more incomplete segments, then the merged seg should be complete
	if (ms->incSegs==NULL) {
		// possibly adjust the known volume based on the tolerances
		if (ms->knownVol != ms->v) {
			if (ftWithinTolerance(ms->knownVol, ms->v, FT_ATOL, FT_RTOL)) {
				DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "UMS1\t%.10g\n", ms->knownVol-ms->v); fflush(fpDeltas); })
					ms->knownVol=ms->v;
			} else {
				fprintf(stdout, "merged segment volumes don't match\n");  fflush(stdout);
			}
		}
		// compute the concentration
		ms->c=ms->knownMass/ms->knownVol;
		switch (ms->owner.type) {
		case SO_NONE:
			break;
		case SO_THROUGH:
		{
			// for through segment types, replace the ms (will be in throughSegs list)
			// with a new normal seg with the correct volume and concentration
			PNodeTransport node=ms->owner.object.node;
			ftReplaceSeg((Pseg)ms, createseg(ms->v, ms->c), &node->throughSegs, 26);
		}
			break;
		case SO_INFLOW:
		{
			PNodeTransport node=ms->owner.object.node;
			int inflowIdx=ms->owner.data.inflowIdx;
			ftReplaceSeg((Pseg)ms, createseg(ms->v, ms->c), &node->inflowSegs.segs[inflowIdx], 260);
			/* old comment...
			// nothing to do.  In fact, there shouldn't be any at this point.
			// are we sure...
			fprintf(stdout,"In updateMergedSegment - SO_INFLOW - investigation/proper handling needed\n");  fflush(stdout);
			*/
		}
			break;
		case SO_DEMAND:
		{
			// for demand segment types, replace the ms (will be in demandSegs list)
			// with a new normal seg with the correct volume and concentration
			PNodeTransport node=ms->owner.object.node;
			ftReplaceSeg((Pseg)ms, createseg(ms->v, ms->c), &node->demandSegs, 27);
		}
			break;
		case SO_LINK:
		{
			// for link segment types, replace the ms (will be in FirstSeg/LastSeg lists)
			// with a new normal seg with the correct volume and concentration
			int linkIdx = ms->owner.object.linkIdx;
			ftReplaceLinkSeg((Pseg)ms, createseg(ms->v, ms->c), linkIdx, 22);
		}
			break;
		case SO_MERGED:
		{
			  // should never happen since merged segments have (and need) no knowledge of their parent
		}
			break;
		case SO_INCOMPLETE_PART:
		{
			// incSeg is owned by another incomplete segment
			PSegOwner owner=&ms->owner;
			PIncSegPart incpart;
			PIncompleteSeg incpartOwner;

			incpart = owner->object.incSeg;
			incpartOwner=owner->data.incSeg;
			// call updatePartSegment to handle it.
			ftUpdatePartSegment(incpart, incpartOwner, (Pseg)ms, mass,completing);
//			ftUpdateIncompleteSegs(downNode,inflowIdx,seg->c,&volMoved);
		}
			break;
		default:
		{
		   fprintf(stdout, "Unhandled merged owner type: %s (%d)\n", OwnerTypes[ms->owner.type], ms->owner.type);  fflush(stdout);
		}
			break;
		}
		FTInitNoneOwnerType(incSeg->owner, incSeg->owner.source*100+13);
	}
}
/*
**--------------------------------------------------------------
**   Input:   seg: the segment to remove from the list
**            list: the list from which to remove it.  This may be modified by this function
**   Output:  none
**   Purpose: Remove a segment from a list of segments
**--------------------------------------------------------------
*/
void ftRemoveSegFromList(Pseg seg, Pseg *list) {
	int found=0;
	Pseg p, t;
	p=*list;
	t=NULL;
	while (p!= NULL) {
		if (p==seg) { // found it
			found=1;
			if (t==NULL) { // first element
				*list=p->prev;
			} else {
				t->prev=p->prev;
			}
		}
		t=p;
		p=p->prev;
	}
	if (!found) {
		fprintf(stdout, "Segment not found!\n");  fflush(stdout);
	}
}
/*
**--------------------------------------------------------------
**   Input:   segMass: The mass to add to a segment
**            is: the incomplete segment to add it to
**   Output:  none
**   Purpose: add the specified mass to is and all of its ancestors
**--------------------------------------------------------------
*/
void ftAddMass(double segMass, PIncompleteSeg is) {
	int i;
	is->massadded += segMass;
	for (i=0; i<is->nchildren; i++) {
		ftAddMass(segMass * is->children[i].ratio, is);
	}
}
/*
**--------------------------------------------------------------
**   Input:   mergedSeg: the merged segment to be split
**            nparts: the number of parts to split it into
**            owners: the owners of the new segments (size: nparts)
**            volumes: the volumes of thge new segments (size: nparts)
**            segs: the new segments (size: nparts)
**   Output:  none
**   Purpose: Split the merged segment
**--------------------------------------------------------------
*/
void ftSplitMergedSegment(PMergedSeg mergedSeg, int nparts, PSegOwner owners, VolSplitData volSplit[], PMergedSeg *segs) {
	int i;
	PSegOwner incOwners=(PSegOwner)calloc(nparts, sizeof(SegOwner));
	PVolSplitData incVolSplit=(PVolSplitData)calloc(nparts, sizeof(VolSplitData));
	PIncompleteSeg is;
	double sumRatio=0;

	// create new segments based on the ratio of each part volume to the total volume of the segment
	for (i=0; i<nparts; i++) {
		double thisRatio=volSplit[i].volume/mergedSeg->v;
		PMergedSeg ms=(PMergedSeg)calloc(1, sizeof(SMergedSeg));

		sumRatio+=thisRatio;

		ms->v=volSplit[i].volume;
		ms->c=-2;  // identify as a merged segment
		ms->prev=NULL;
		ms->adjustment=volSplit[i].adjustment;
		ms->knownVol=mergedSeg->knownVol*thisRatio;
		ms->knownMass=mergedSeg->knownMass*thisRatio;
		memcpy(&ms->owner, &owners[i], sizeof(SegOwner));
		segs[i]=ms;
	}
	// now split the incomplete segments and add each to the corresponding incomplete segment
	for (is=mergedSeg->incSegs; is != NULL; is=is->prev) {
		sumRatio=0;
		for (i=0; i<nparts; i++) {
			double thisRatio=volSplit[i].volume/mergedSeg->v;
			sumRatio+=thisRatio;
			FTInitMergedPartOwner(incOwners[i], segs[i], 16)
			incVolSplit[i].volume=is->v*thisRatio;
			incVolSplit[i].adjustment=1.0; //??
		}
		ftSplitIncompleteSegment(is, nparts, incOwners, incVolSplit);
		for (i=0; i<nparts; i++) {
			PIncompleteSeg pseg=&is->children[i];
			pseg->prev=segs[i]->incSegs;
			segs[i]->incSegs=pseg;
		}
	}
	free(incOwners);
	free(incVolSplit);
	for (i=0; i<nparts; i++) {
		ftValidateMergedSegment(segs[i], "SMS");
	}
}

/*
**--------------------------------------------------------------
**   Input:   incSeg: The incomplete segment to split
**            nparts: the number of parts to split it into
**            owners: the owners of the new segments (size: nparts)
**            volumes: the volumes of thge new segments (size: nparts)
**   Output:  none
**   Purpose: Split the incomplete segment
**--------------------------------------------------------------
*/
void ftSplitIncompleteSegment(PIncompleteSeg incSeg, int nchildren, SegOwner owners[], VolSplitData volSplit[]) {
	int i, n;
	incSeg->nchildren = nchildren;
	// allocate memory for the children
	incSeg->children = (PIncompleteSeg)calloc(nchildren, sizeof(SIncompleteSeg));
#ifdef ENABLE_DBG_PRINT
	if (fpOwners!= NULL) {
		for (i = 0; i<nchildren; i++) {
			fprintf(fpOwners, "%s\t%d\t%d\t", OwnerTypes[owners[i].type], owners[i].source, owners[i].sequence);
			switch (owners[i].type) {
			case SO_NONE:
				fprintf(fpOwners, "\n");
				break;
			case SO_THROUGH:
			case SO_DEMAND:
			case SO_INFLOW:
			{
				SNodeTransport *node=owners[i].object.node;
				fprintf(fpOwners, "Node id\t%s\n", node->node->ID);
			}
				break;
			case SO_LINK:
			{
				int linkIdx = owners[i].object.linkIdx;
				fprintf(fpOwners, "Link id\t%s\n", Link[linkIdx].ID);
			}
				break;
			case SO_MERGED:
			{
				fprintf(fpOwners, "\n");
			}
				break;
			case SO_INCOMPLETE_PART:
			{
				PSegOwner owner=&owners[i];
				PIncSegPart incpart;
				PIncompleteSeg incpartOwner;

				incpart = owner->object.incSeg;
				incpartOwner=owner->data.incSeg;
				fprintf(fpOwners, "Incomplete part owner\t%s\n", incpartOwner->node->node->ID);
			}
				break;
			default:
				break;
			}
		}
		fflush(fpOwners);
	}
#endif
	for (i = 0; i<nchildren; i++) {
		PIncompleteSeg is=&incSeg->children[i];

		is->v=volSplit[i].volume;
		is->c=-1;  // identify as an incomplete segment
		is->prev=NULL;
		is->parent=incSeg;
		is->adjustment=1.0;//volSplit[i].adjustment;
		// calculate ratio of this child to its parent
		is->ratio=volSplit[i].volume/incSeg->v;
		is->massadded=incSeg->massadded * is->ratio;
		memcpy(&is->owner, &owners[i], sizeof(SegOwner));
		is->nparts=incSeg->nparts;
		// nparts is equal to the number of inflow links the original incomplete segment owner has
		is->parts=(PIncSegPart)calloc(is->nparts, sizeof(SIncSegPart));
		for (n=0; n<is->nparts; n++) {
			PIncSegPart isc=&is->parts[n];
			isc->ratio=incSeg->parts[n].ratio;
			isc->volume=is->v*isc->ratio;
			isc->flowed=incSeg->parts[n].flowed/incSeg->parts[n].volume*isc->volume;
			isc->conc=incSeg->parts[n].conc;
			if(isc->conc==-2) {
				// INVESTIGATE: needed??
				printf("");
			}
		}
	}
}

/*
**--------------------------------------------------------------
**   Input:   seg: the list of segments whose water volume is to be computed
**   Output:  The total water volume
**   Purpose: Compute the volume of a list of segments
**--------------------------------------------------------------
*/
double ftComputeSegmentWaterVolume(Pseg seg) {
	double v=0;
	while (seg != NULL) {
		v += seg->v;
		seg=seg->prev;
	}
	return v;
}
/*
**--------------------------------------------------------------
**   Input:   seg: the list of segments whose water volume is to be computed
**   Output:  The total water volume
**   Purpose: Compute the volume of a list of segments
**--------------------------------------------------------------
*/
double ftComputeCompleteSegmentWaterVolume(Pseg seg) {
	double v=0;
	while (seg != NULL) {
		if(seg->c >= 0) {
			v += seg->v;
		}
		seg=seg->prev;
	}
	return v;
}
/*
**--------------------------------------------------------------
**   Input:   nodeLinks: the array of link indices
**            linkIdx: the link index being searched for
**   Output:  The index into nodelinks that contains linkIdx
**   Purpose: Find the nodelink index that contains linkIdx
**--------------------------------------------------------------
*/
int ftGetLinkIndex(SNodeLinks *nodelinks, int linkIdx) {
	int i;
	for (i=0; i<nodelinks->nlinks; i++) {
		if (nodelinks->linkIdxs[i]==linkIdx) return i;
	}
	return -1;
}
/*
**--------------------------------------------------------------
**   Input:   source: the mass source
**            dt: the current timestep
**            node: the node
**   Output:  The total mass added by the source
**   Purpose: Compute the amount of mass added by the given
**            source for the specified time
**--------------------------------------------------------------
*/
// INVESTIGATE: look at comments from copied code - relevant? need to implement??
double ftMassAdded(Psource source, long dt, SNodeTransport *node) {
	//	double qout, qcutoff;
	double volout;
	double s;
	double massadded=0;
	/* Establish a flow cutoff which indicates no outflow from a node */
	//	qcutoff = 10.0*TINY;

	/* Find total flow volume leaving node */
	if (node->nodeIdx <= Njuncs) volout = node->totalVolOutflow;  /* Junctions */
	else
		volout = node->totalVolOutflow - node->demand;    /* Tanks */
	//	qout = volout / (double) dt;

	/* Evaluate source input only if node outflow > cutoff flow */
	//	if (qout > qcutoff) {

	/* Mass added depends on type of source */
	s = sourcequal(source);
	switch (source->Type) {
		/* Concen. Source: */
		/* Mass added = source concen. * -(demand) */
	case CONCEN:

		/* Only add source mass if demand is negative */
		if (NodeDemand[node->nodeIdx] < 0.0) {
			// use NodeDemand[nodeidx] here because node->demand does not contain demand for nodes with negative demand
			massadded = -s*NodeDemand[node->nodeIdx]*dt;

			/* If node is a tank then set concen. to 0. */
			/* (It will be re-set to true value in updatesourcenodes()) */
// ???				if (n > Njuncs) NodeQual[n] = 0.0;
		} else {
			massadded = 0.0;
		}
		break;
		/* Mass Inflow Booster Source: */
	case MASS:
		massadded = s*dt;
		break;

		/* Setpoint Booster Source: */
		/* Mass added is difference between source */
		/* & node concen. times outflow volume  */
		/* NEED to figure out how best to handle SETPOINT sources */
	case SETPOINT:
		//				if (s > NodeQual[n]) massadded = (s-NodeQual[n])*volout;
		//				else massadded = 0.0;
		break;

		/* Flow-Paced Booster Source: */
		/* Mass added = source concen. times outflow volume */
	case FLOWPACED:
		massadded = s*volout;
		break;
	}
	//	}
	MB_MassAdded[node->nodeIdx]+=massadded;
	return massadded;
}
/*
**--------------------------------------------------------------
**   Input:   segs: a nodesegs structure
**   Output:  1 if all inflow segs are not null, 0 otherwise
**   Purpose: Determine whether there is at least one segment
**            in each segment list within segs
**--------------------------------------------------------------
*/
int ftHaveAllInflowSegs(PNodeSegs segs) {
	int i;
	for (i=0; i<segs->nlinks; i++) {
		if (segs->segs[i].firstSeg==NULL)
			return 0;
	}
	return 1;
}
/*
**--------------------------------------------------------------
**   Input:   segs: a nodesegs structure
**   Output:  1 if any inflow segs are not null, 0 otherwise
**   Purpose: Determine whether there is at least one segment
**            in any segment list within segs
**--------------------------------------------------------------
*/
int ftHaveAnyInflowSegs(SNodeSegs *segs) {
	int i;
	for (i=0; i<segs->nlinks; i++) {
		if (segs->segs[i].firstSeg!=NULL)
			return 1;
	}
	return 0;
}
/*
**--------------------------------------------------------------
**   Input:   segvol: the known segment volume
**            segmass: the known segment mass
**            segs: The inflow segs (a mixture of incomplete
**                  and merged - normal segs are represented in
**                  segvol & segmass)
**   Output:  a new MergedSeg
**   Purpose: Create a mergedSeg from the parameters
**--------------------------------------------------------------
*/
PMergedSeg ftCreateMergedSeg(double segvol, double segmass, Pseg *segs, int n, PSegOwner owner) {
	int i;
	PMergedSeg s=(PMergedSeg)calloc(1, sizeof(SMergedSeg));

	s->c=-2;
	s->v=segvol;
	s->adjustment=1.0;
	s->knownVol=segvol;
	s->knownMass=segmass;
	s->owner=*owner;
	for (i=0; i<n; i++) {
		Pseg seg=segs[i];
		if (seg != NULL) {
			if (IsIncompleteSeg(seg)) { // incomplete segment
				PIncompleteSeg incSeg=(PIncompleteSeg)seg;
				FTInitMergedPartOwner(incSeg->owner, s, 17)
					incSeg->prev=s->incSegs;
				s->incSegs=incSeg;
				s->v+=incSeg->v;
			} else if (IsMergedSeg(seg)) { // merged segment
				PMergedSeg ms=(PMergedSeg)seg;
				PIncompleteSeg incs;
				s->knownVol+=ms->knownVol;
				s->knownMass+=ms->knownMass;
				s->v += ms->v;
				incs=ms->incSegs;
				while (incs != NULL) {
					PIncompleteSeg t=incs->prev;
					FTInitMergedPartOwner(incs->owner, s, 18)
						incs->prev=s->incSegs;
					s->incSegs=incs;
					incs=t;
				}
				// this pointer was created when splitting a merged segment in ftMergeSegments
				// after using the data to accumulate in this newly created merged segment,
				// it will never be refenced again and needs to be freed here.
				free(ms);
			}
		}
	}
	ftValidateMergedSegment(s, "CMS");
	return s;
}


/*
**--------------------------------------------------------------
**   Input:   volSegs: the segments volumes to merge
**            nodeSegs: The node's inflow segments
**            inflow: The node's inflow parameters
**   Output:  A new segment to be moved.  It will have a volume
**            of vol, and may either be a normal, incomplete or
**            merged segment.
**   Purpose: create a new segment of volume vol (from volSegs) derived from
**            the first inflow segment of each inflow link
**--------------------------------------------------------------
*/
Pseg ftMergeSegments(PNodeSegs volSegs, PNodeSegs nodeSegs, PNodeLinks inflow) {
	// segs is statically created and held here to avoid having to
	// allocate and free this memory constantly which is expensive
	static int segSize=0;
	static Pseg *segs=NULL;

	SegOwner msOwner;
	int n=inflow->nlinks;
	int i;
	double segvol=0;
	double segmass=0;
	int merged=0, incomplete=0;
	double tvol=0;

	if (n>segSize) {
		segSize=n;
		segs=(Pseg*)realloc(segs, n*sizeof(Pseg));
	}
	// for each inflow link
	for (i=0; i<n; i++) {
		// compute the volume required from this segment
		double thisv=volSegs->segs[i].firstSeg->v;
		Pseg firstSeg=nodeSegs->segs[i].firstSeg;
		tvol+=thisv;
		ftRemoveFirstSegment(&volSegs->segs[i]);

		DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f", thisv); })
		if (IsNormalSeg(firstSeg)) { // "normal" segment
			// seg segs[i] to null to signify that it is a normal segment
			segs[i]=NULL;
			// remove the volume from the segment.
			firstSeg->v-=thisv;
			// increment the new segment's volume
			segvol += thisv;
			// and the mass
			segmass += thisv*firstSeg->c;
			DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f\t%.8f\t", thisv, firstSeg->v); })
				// in some rare cases, this is < 0 (but very small ~ 1e-14).  If so, treat it as 0.
			if (firstSeg->v==0) {
				// remove the empty segment if there is no more volume
				ftRemoveFirstSegment(&nodeSegs->segs[i]);
			}
		} else {
			// incomplete or merged segment - it needs to be split
			int nparts;
			SegOwner *owners;
			VolSplitData *volSplit;
			PSegOwner owner;
			Pseg replacementSeg=NULL;

			if (thisv==firstSeg->v) {
				nparts=1;
			} else {
				nparts=2;
			}
			owners=(PSegOwner)calloc(nparts, sizeof(SegOwner));
			volSplit=(PVolSplitData)calloc(nparts, sizeof(VolSplitData));

			if (IsIncompleteSeg(firstSeg)) { // "incomplete" segment
				PIncompleteSeg incSeg=(PIncompleteSeg)firstSeg;
				incomplete++;
				owner=&incSeg->owner;
			} else if (IsMergedSeg(firstSeg)) { // "Merged" segment
				PMergedSeg mergedSeg=(PMergedSeg)firstSeg;
				merged++;
				owner=&mergedSeg->owner;
			}
			FTInitCopyOwner(owners[0], owner, 19)
				volSplit[0].volume=thisv;
			volSplit[0].adjustment=1.0;
			if (nparts==2) {
				FTInitCopyOwner(owners[1], owner, 20)
					volSplit[1].volume=firstSeg->v-thisv;
				volSplit[1].adjustment=1.0;
			}
			if (IsIncompleteSeg(firstSeg)) { // "incomplete" segment
				PIncompleteSeg incSeg=(PIncompleteSeg)firstSeg;
				ftSplitIncompleteSegment(incSeg, nparts, owners, volSplit);
				segs[i]=(Pseg)&incSeg->children[0];
				if (nparts==2) {
					replacementSeg=(Pseg)&incSeg->children[1];
				}
			} else if (IsMergedSeg(firstSeg)) { // "Merged" segment
				PMergedSeg *msegs=(PMergedSeg*)calloc(nparts, sizeof(PMergedSeg));
				ftSplitMergedSegment((PMergedSeg)firstSeg, nparts, owners, volSplit, msegs);
				segs[i]=(Pseg)msegs[0];
				if (nparts==2) {
					replacementSeg=(Pseg)msegs[1];
				}
				free(msegs);
			}
			// now that it has been split, set this type to NONE because nobody owns it anymore
			FTInitNoneOwnerType(*owner, owner->source*100+25)
				free(volSplit);
			free(owners);

			// if there was only 1 part, remove the segment from the list
			if (nparts==1) {
				ftUnlinkFirstSegment(&nodeSegs->segs[i]);
			} else {
				// otherwise replace it with the replacement seg
				ftReplaceSeg(firstSeg, replacementSeg, &nodeSegs->segs[i], 28);
			}
			DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "\t%.8f\t%.8f\t", thisv, nparts==1?0:replacementSeg->v); })
		}
	}
#ifdef ENABLE_DBG_PRINT
	if (fpi != NULL) {
		for (; i<maxIn; i++) {
			fprintf(fpi, "\t\t\t\t");
		}
	}
	if (fpi!=NULL) { fprintf(fpi, "\t%.8f\t%.8f\t%.8f\t%s\t%s\n", segvol, segmass, segmass/segvol, ftHaveAnyInflowSegs(nodeSegs)?"TRUE":"FALSE", ftHaveAllInflowSegs(nodeSegs)?"TRUE":"FALSE"); }
#endif
	// at this point, if merged and incomplete are both zero, then a normal segment can be created
	if (merged==0 && incomplete==0) {
		checkConc(segmass/segvol);
		return createseg(segvol, segmass/segvol);
	}
	if (n==1) {
		if (incomplete==1 && merged==0) {
			// the only segment is an incomplete one, simply return it
			return segs[0];
		} else if (incomplete==0 && merged==1) {
			// the only segment is a merged one, simply return it
			return segs[0];
		}
	}
	// if we are are, then we need to create a merged segment from segs array an segvol&segmass
	FTInitInflowPartOwner(msOwner, NULL, -1, 22)
		return (Pseg)ftCreateMergedSeg(segvol, segmass, segs, n, &msOwner);
}
// structs used by ftGetVolumesByTime
typedef struct SVolByTime {
	// to allow standard segment methods to operate on these...
	double v;
	double c;  // ignored
	struct SVolByTime *prev;

	Pseg inflowSeg;  // the original inflow seg
	SNodeSegs segments;   // the segments split from this segment
} SVolByTime, *PVolByTime;

typedef struct SVBTIter {
	PVolByTime vbt;
	Pseg cseg;
	Pseg pseg;
} SVBTIter, *PVBTIter;

/*
**--------------------------------------------------------------
**   Input:   node: the node whose inflow segments need to be "mixed"
**            nodeSegs: The node's inflow segments
**            inflow: The node's inflow parameters
**   Output:  volsToMerge: lists of segments from each inflow segment that can be merged
**   Purpose: create a set of segments from the node's inflow segments that can be mixed 
**--------------------------------------------------------------
*/
void ftGetVolumesByTime(PNodeTransport node, SNodeSegs *inflowSegs, SNodeLinks *inflow, SNodeSegs *volsToMerge) {
	int i, n=inflow->nlinks;
	SNodeSegs tvols;
	double mint;
//	double *totInflow=(double *)calloc(inflow->nlinks, sizeof(double));
	static int nfs=0;
	static PVolByTime *firstSegs=NULL;
	static PVBTIter iter=NULL;
	int haveAll;

	memset(&tvols, 0, sizeof(SNodeSegs));
	tvols.segs=(SegList*)calloc(n, sizeof(SegList));
	tvols.nlinks=n;
	if (volsToMerge->size < n) {
		volsToMerge->segs=(SegList*)realloc(volsToMerge->segs, n*sizeof(SegList));
		for (i=volsToMerge->size; i<n; i++) {
			volsToMerge->segs[i].firstSeg=NULL;
			volsToMerge->segs[i].lastSeg=NULL;
		}
		volsToMerge->size=n;
	}
	volsToMerge->nlinks=n;
	n=inflow->nlinks;
	// first compute total inflows for each inflow path
	for (i=0; i<n; i++) {
		Pseg s;
		for (s=inflowSegs->segs[i].firstSeg; s!=NULL; s=s->prev) {
			PVolByTime vbts = (PVolByTime)calloc(1, sizeof(SVolByTime));
			vbts->v=s->v;
			vbts->c=-3;
			vbts->inflowSeg=s;
			vbts->segments.segs=(SegList *)calloc(1, sizeof(SegList));
			vbts->segments.nlinks=1;
			ftAddNodeSeg(&tvols.segs[i], (Pseg)vbts, 99);
//			totInflow[i]+=s->v;
		}
	}
//	// now find the minimum time that will consume the total inflow
//	mint=FLT_MAX;
//	for (i=0; i<n; i++) {
//		double t=totInflow[i]/inflow->flowRate[i];
//		mint=MIN(mint, t);
//	}
	if (nfs < n) {
		firstSegs=(PVolByTime *)realloc(firstSegs, n*sizeof(PVolByTime));
		for (i=nfs; i<n; i++) {
			firstSegs[i]=NULL;
		}
		iter=(PVBTIter)realloc(iter, n*sizeof(SVBTIter));
		for (i=nfs; i<n; i++) {
			memset(&iter[i], 0, sizeof(SVBTIter));
		}
		nfs=n;
	}
	haveAll=1;
	for (i=0; i<n; i++) {
		firstSegs[i]=(PVolByTime)tvols.segs[i].firstSeg;
		haveAll&=firstSegs[i]!=NULL;
	}
	// temporary check to see if any inflow segs have more than one seg

	while (haveAll) {
		mint=FLT_MAX;
		for (i=0; i<n; i++) {
			double t=firstSegs[i]->v/inflow->flowRate[i];
#ifdef DEBUG
			if (firstSegs[i]->v < 1e-8) {
				fprintf(stdout, "");  fflush(stdout);
			}
#endif
			mint=MIN(mint, t);
		}
		for (i=0; i<n; i++) {
			double v=inflow->flowRate[i] * mint;
			PVolByTime vbt=(PVolByTime)firstSegs[i];
			Pseg s;
			vbt->v-=v;
			s=createseg(v, 0); s->c=-3;
			ftAddNodeSeg(&vbt->segments.segs[0], s, 98);
			if (vbt->v == 0) {
				firstSegs[i]=vbt->prev;
				haveAll&=firstSegs[i]!=NULL;
			} else if (ftWithinTolerance(firstSegs[i]->v, 0, 1e-8, FT_RTOL)) {
				DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "GVBT2\t%.10g\n", firstSegs[i]->v); fflush(fpDeltas); })
					//					if(vbt->v < 0) {
					vbt->segments.segs[0].lastSeg->v+=vbt->v;
				vbt->v=0;
				firstSegs[i]=vbt->prev;
				haveAll&=firstSegs[i]!=NULL;
				//					}
			}
		}
	}
	for (i=0; i<n; i++) {
		iter[i].vbt=(PVolByTime)tvols.segs[i].firstSeg;
		iter[i].cseg=iter[i].vbt->segments.segs[0].firstSeg;
		iter[i].pseg=NULL;
	}

	// find any split segments that are very small and figure out where to put them (if anywhere)
	haveAll=1;
	for (i=0; i<n; i++) {
		firstSegs[i]=(PVolByTime)tvols.segs[i].firstSeg;
		haveAll&=iter[i].vbt!=NULL;
	}
	while (haveAll) {
		int hasSmall=0;
		int hasPrev=0;
		for (i=0; i<n; i++) {
			hasSmall += iter[i].cseg->v < 1e-8;  // ? smaller? larger?
			hasPrev += iter[i].pseg!=NULL;
		}
		if (hasSmall == n && hasPrev==n) {
			for (i=0; i<n; i++) {
				iter[i].pseg->v += iter[i].cseg->v;
				iter[i].pseg->prev = iter[i].cseg->prev;
				ftReleaseSeg(iter[i].cseg);
				iter[i].cseg=iter[i].pseg;
			}
		}
		for (i=0; i<n; i++) {
			PVBTIter vbti=&iter[i];
			vbti->pseg=vbti->cseg;
			vbti->cseg=vbti->cseg->prev;
			if (vbti->cseg == NULL) {
				vbti->pseg=NULL;
				vbti->vbt=vbti->vbt->prev;
				haveAll&=vbti->vbt!=NULL;
				if (vbti->vbt != NULL) {
					vbti->cseg=vbti->vbt->segments.segs[0].firstSeg;
					vbti->pseg=NULL;
				}
			}
		}
	}
	// look through tvols to find any with v > 0 and figure out where to put them (if anywhere)
	for (i=0; i<n; i++) {
		PVolByTime vbt;
		while ((vbt=(PVolByTime)tvols.segs[i].firstSeg) != NULL) {
			if (vbt->v != 0) {
				double v=ftSumVolume(&vbt->segments.segs[0]);
				double dv=vbt->inflowSeg->v-v;
				// do something to move the volume if sufficiently small
				if (ftWithinTolerance(dv, 0, 1e-8, FT_RTOL)) {
					// add it to one of the child segments
					vbt->segments.segs[0].lastSeg->v += vbt->v;
					vbt->v=0;
				}
			}
			// now add these to volsToMerge
			while (vbt->segments.segs[0].firstSeg!=NULL) {
				Pseg s=ftUnlinkFirstSegment(&vbt->segments.segs[0]);
				ftAddNodeSeg(&volsToMerge->segs[i], s, 96);
			}
			ftUnlinkFirstSegment(&tvols.segs[i]);
			free(vbt->segments.segs);
			free(vbt);
		}
	}
//	free(totInflow);
	free(tvols.segs);
}
/*
**--------------------------------------------------------------
**   Input:   node: the node whose segs need to be moved
**   Output:  segs: the list of segments to move
**   Purpose: combine the node's inflow segments into segments and return them in segs.
**--------------------------------------------------------------
*/
void ftGetSegmentsToMove(SegList * segs, SNodeTransport *node) {
	SNodeSegs *inflowSegs=&node->inflowSegs;
	// process inflow segments first
	if (node->outflowAdjustment == 0) {
		// There is no outflow from this node but we have inflow.
		// this happens if there is a bad hydraulic solution
		// in this case, track the lost mass in MB_MassNoOutflow array
		int n=inflowSegs->nlinks;
		int i;
		for (i=0; i<n; i++) {
			Pseg firstSeg;
			while ((firstSeg=inflowSegs->segs[i].firstSeg) != NULL) {
				if (IsNormalSeg(firstSeg)) {
					if (firstSeg->c > 0) {
						double mass=firstSeg->c*firstSeg->v;
						MB_MassNoOutflow[node->nodeIdx]+=mass;
					}
					ftRemoveFirstSegment(&inflowSegs->segs[i]);
				} else if (IsIncompleteSeg(firstSeg)) {
					// INVESTIGATE:  Need to save these incomplete segments so the missing mass can be tracked
					fprintf(stdout, "No Outflow, inflow segment is an incomplete segment\n");  fflush(stdout);
					ftUnlinkFirstSegment(&inflowSegs->segs[i]);
				} else if (IsMergedSeg(firstSeg)) {
					// INVESTIGATE:  Need to save these merged segments so the missing mass can be tracked
					PMergedSeg ms = (PMergedSeg)firstSeg;
					fprintf(stdout, "No Outflow, inflow segment is a merged segment:\n");  fflush(stdout);
					fprintf(stdout, " %08x %g\n",ms,ms->v);  fflush(stdout);
					ftUnlinkFirstSegment(&inflowSegs->segs[i]);
				}
			}
		}
		return;
	}
	if (ftHaveAllInflowSegs(inflowSegs)) {
		// there are segments on all inflow links
		int i, n;
		double **vols=NULL;
		static PNodeSegs volsToMerge=NULL;

		if (volsToMerge==NULL) {
			volsToMerge=(PNodeSegs)calloc(1, sizeof(SNodeSegs));
		}
		n=inflowSegs->nlinks;
		if (n > 1) {
			int ngt1 = 0;
			for (i = 0; i < n; i++) {
				ngt1 += inflowSegs->segs[i].firstSeg->prev != NULL;
			}
			if (ngt1) {
				printf("");
			}
		}
		ftGetVolumesByTime(node, inflowSegs, &node->inflow, volsToMerge);
		while (ftHaveAllInflowSegs(volsToMerge)) {
			Pseg mergedSeg=ftMergeSegments(volsToMerge, inflowSegs, &node->inflow);
			// and add that segment to the list of segments to move.
			ftAddNodeSeg(segs, (Pseg)mergedSeg, 19);
		}
		// now look at the first segments on each of the inflow links
		for (i=0; i<n; i++) {
			Pseg firstSeg=inflowSegs->segs[i].firstSeg;
			// if it's not already null and is sufficiently close to zero
			if (firstSeg != NULL) {
				if (firstSeg->v != 0 && ftWithinTolerance(firstSeg->v, 0, .001, FT_RTOL)) {
					DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "GSTM1\t%.10g\n", firstSeg->v); fflush(fpDeltas); })
#ifdef DEBUG
				} else {
					// non-zero volume left over - verify
					fprintf(stdout, "");  fflush(stdout);
#endif
				}
				// remove the empty segment
				DBG_PRINT(if (fpi!=NULL) { fprintf(fpi, "SEGSTOMOVE-FS\t\t%.8e\n", firstSeg->v); })
				if (IsNormalSeg(inflowSegs->segs[i].firstSeg)) {
					// if it is a normal segment, it is OK to remove it
					ftRemoveFirstSegment(&inflowSegs->segs[i]);
				} else {
					// otherwise just unlink it since this segment is either an incomplete or merged segment
					// and shouldn't be placed on the FreeSeg list.
					ftUnlinkFirstSegment(&inflowSegs->segs[i]);
				}
			}
		}

	} else {
		// if we are here, it means that all remaining nodes aren't ready to process
		// likely due to a loop condition causing one or more of its inflow
		// links to not be there.
		// An incomplete segment will be created to allow this segment to be moved around
		// just like any other segment, and will ultimately be converted to a regular segment
		// when inflow to this node finally gets here
		int i, n=inflowSegs->nlinks;
		double segvol=0;
		double segmass=0;
		// find the minimum volume of all first segments
		while (ftHaveAnyInflowSegs(inflowSegs)) {
			PIncompleteSeg is;
			PIncompleteSegmentData incompleteData;
			double mint=FLT_MAX;
			double v;
			// find the time to move the largest complete leading segment
			for (i=0; i<n; i++) {
				if (inflowSegs->segs[i].firstSeg!=NULL) {
					mint=MIN(mint, inflowSegs->segs[i].firstSeg->v/node->inflow.flowRate[i]);
				}
			}
			// now compute what the total volume would be
			v=0;
			for (i=0; i<n; i++) {
				v+=node->inflow.flowRate[i] * mint;
			}
			incompleteData=(PIncompleteSegmentData)calloc(1, sizeof(SIncompleteSegmentData));
			ftAddIncompleteSegment(incompleteData, &node->incompleteSegs);
			incompleteData->incomplete=(PIncompleteSeg*)calloc(n, sizeof(SIncompleteSeg));
			// PIncompleteSegmentData instances hold an array of PIncompleteSeg of where each
			// element corresponds to an inflow.  If the value is NULL, then that inflow link
			// has a segment.  If it isn't NULL, then it will be a pointer the the incomplete
			// segment.  If there are multiple inflow links without flow, all of those array
			// elements will point to the same PIncompleteSeg.
			is=(PIncompleteSeg)calloc(1, sizeof(SIncompleteSeg));
			is->v=v;     // the volume of this segment
			is->c=-1.0;  // indicates this is a SpartialSegment
			is->parent=NULL;
			is->node=node;
			is->ratio=1.0;
			//			is->valveAdjRatio=1.0;
			is->adjustment=1.0;
			FTInitInflowPartOwner(is->owner, node, -1, 23)
			is->nchildren=0;
			is->children=NULL;
			is->nparts=n;
			is->parts=(PIncSegPart)calloc(n, sizeof(SIncSegPart));
			// there is one part for each inflow link
			for (i=0; i<n; i++) {
				Pseg firstSeg=inflowSegs->segs[i].firstSeg;
				// the volume of this part
				double thisv=v * node->inflow.flowRatio[i];

				PIncSegPart part=&is->parts[i];
				part->ratio=node->inflow.flowRatio[i];
				part->flowed=0;  // how much of this part's volume has actually been "seen" (i.e. from a normal segment)
				if (firstSeg == NULL) {
					// there is no inflow on this inflow link yet
					part->conc=-1;  // signal that this component of the incomplete segment hasn't yet been seen;
					part->volume=thisv;
					incompleteData->incomplete[i]=is;
					node->inflow.volumeFlowed[i] += thisv;
				} else {
					// there is inflow on this inflow link
					// if the segment's volume is sufficiently close to the calculated vol
					if (firstSeg->v!=thisv && ftWithinTolerance(firstSeg->v, thisv, FT_ATOL, FT_RTOL)) {
						// set the volume of the part to the segment's volume
						DBG_PRINT(if (fpDeltas!=NULL) { fprintf(fpDeltas, "GSTM2\t%.10g\n", thisv-firstSeg->v); fflush(fpDeltas); })
							thisv=firstSeg->v;
					}
					// assign part's volume & concentration
					part->volume=thisv;
					part->conc=firstSeg->c;
					if (IsNormalSeg(firstSeg)) {
						// this is a normal segment
						part->flowed=part->volume;  // all the volume has actually flowed through the link.
						firstSeg->v-=thisv;
						if (firstSeg->v==0) {
							// remove the empty segment
							ftRemoveFirstSegment(&inflowSegs->segs[i]);
						}
					} else {
						// either incomplete or merged and needs to be split into 
						// an INCOMPLETE_PART and, if there is more than thisv volume,
						// a replacement INFLOW part
						PSegOwner owners;
						PVolSplitData volSplit;
						int nparts;
						Pseg toAdd=NULL;
						Pseg replacementSeg=NULL;

						if (firstSeg->v == thisv) {
							nparts=1;
						} else {
							nparts=2;
						}
						// setup the parameters for splitting
						owners=(PSegOwner)calloc(nparts, sizeof(SegOwner));
						volSplit=(PVolSplitData)calloc(nparts, sizeof(VolSplitData));

						FTInitIncompletePartOwner(owners[0], part, is, 24)
						volSplit[0].volume=thisv;
						volSplit[0].adjustment=1.0;
						if (nparts==2) {
							volSplit[1].volume=firstSeg->v-thisv;
							volSplit[1].adjustment=1.0;
						}
						if (IsIncompleteSeg(firstSeg)) {
							PIncompleteSeg srcSeg=(PIncompleteSeg)firstSeg;
							if (nparts == 2) {
								FTInitInflowPartOwner(owners[1], srcSeg->owner.object.node, srcSeg->owner.data.inflowIdx, 26);
							}
							ftSplitIncompleteSegment(srcSeg, nparts, owners, volSplit);

							FTInitNoneOwnerType(srcSeg->owner, srcSeg->owner.source*100+27);
							toAdd=(Pseg)&srcSeg->children[0];
							if (nparts == 2) {
								replacementSeg=(Pseg)&srcSeg->children[1];
							}
						} else if (IsMergedSeg(firstSeg)) {
							PMergedSeg ms=(PMergedSeg)firstSeg;
							PMergedSeg *segs=(PMergedSeg*)calloc(nparts, sizeof(PMergedSeg));
							if (nparts == 2) {
								FTInitInflowPartOwner(owners[1], ms->owner.object.node, ms->owner.data.inflowIdx, 28);
							}
							ftSplitMergedSegment(ms, nparts, owners, volSplit, segs);
							FTInitNoneOwnerType(ms->owner, ms->owner.source*100+29);

							toAdd=(Pseg)segs[0];
							if (nparts==2) {
								replacementSeg=(Pseg)segs[1];
							}
							free(segs);
						}
						free(owners);
						free(volSplit);
						// place the incomplete part on the part's seg list
						toAdd->prev=part->segs;
						part->segs=toAdd;
						// if there is a replacement segment
						if (replacementSeg != NULL) {
							// replace the correct inflowSeg with it
							ftReplaceSeg(firstSeg, replacementSeg, &inflowSegs->segs[i], 29);
						} else {
							// otherwise unlink the correct inflowSeg
							ftUnlinkFirstSegment(&inflowSegs->segs[i]);
						}
						incompleteData->incomplete[i]=is;
					}
				}
			}
			// add the incomplete segment to the list of segments to move
			// OK to cast sp to Pseg because the first elements are the same;
			ftAddNodeSeg(segs, (Pseg)is, 20);
		}
	}
	// if there is massadded to this node for this timestep
	if (node->massadded) {
		Pseg seg;
		for (seg=segs->firstSeg; seg != NULL; seg=seg->prev) {
			// add it to each newly created segment in
			// an amount proportionate to the segment's volume
			// to the total volume inflow
			double segMass=node->massadded * seg->v / node->totalVolInflow;
			MB_MassMoved[node->nodeIdx]+=segMass;
			node->massused += segMass;
			if (IsNormalSeg(seg)) {
				// compute the new total mass and adjust the segment's concentration
				segMass += seg->v*seg->c;
				seg->c=segMass/seg->v;
				checkConc(seg->c);
			} else if (IsMergedSeg(seg)) { // merged segment
				PMergedSeg ms=(PMergedSeg)seg;
				// add the mass to the known mass for the merged segment
				ms->knownMass += segMass;
			} else if (IsIncompleteSeg(seg)) { // incomplete segment
				PIncompleteSeg is=(PIncompleteSeg)seg;
				// add the mass to the incomplete segment and all of it's children
				ftAddMass(segMass, is);
			}
		}
	}
}


/*
**--------------------------------------------------------------
**   Input:   incompleteData: The incomplete segment data to add
**   Output:  list: the list to add it to
**   Purpose: Add the incompleteSegmentData to the list
**--------------------------------------------------------------
*/
void ftAddIncompleteSegment(PIncompleteSegmentData incompleteData, PIncompleteSegmentData *list) {
	PIncompleteSegmentData p, t;
	p=*list;
	t=NULL;
	incompleteData->prev=NULL;
	if (p==NULL) {
		*list=incompleteData;
		return;
	}
	while (p->prev != NULL) { p=p->prev; }
	p->prev=incompleteData;
}

/*
**--------------------------------------------------------------
**   Input:   v1: value to compare
**            v2: value to compare
**            atol: the absolute tolerance
**            rtol: The relative tolerance
**   Output:  1 if v1 & v2 are within a tolerance, 0 otherwise
**   Purpose: Compare two values agains an absolute and relative
**            tolerance
**--------------------------------------------------------------
*/
int ftWithinTolerance(double v1, double v2, double atol, double rtol) {
	double val=v1-v2;
	if (val != 0) {
		if (ABS(val) < atol) {
			return 1;
		} else {
			double r=ABS(val)/(v1==0?v2:v1);
			if (r < rtol) {
				return 1;
			}
		}
		return 0;
	}
	// value is 0
	return 1;
}
/*
**--------------------------------------------------------------
**   Input:   linkIdx: The link index for which the first segment
**            is to be unlinked
**   Output:  The segment that was unlinked
**   Purpose: Unlink the first segment for a Link, but do not free it
**--------------------------------------------------------------
*/
Pseg ftUnlinkFirstLinkSegment(int linkIdx) {
	Pseg firstSeg=FirstSeg[linkIdx];
	FirstSeg[linkIdx]=firstSeg->prev;

	if (FirstSeg[linkIdx] == NULL) {
		LastSeg[linkIdx] = NULL;
	}
	return firstSeg;
}
/*
**--------------------------------------------------------------
**   Input:   linkIdx: The link index for which the first segment
**            is to be removed
**   Output:  none
**   Purpose: Unlink the first segment for a Link, and place it
**            back on the FreeSegs list
**--------------------------------------------------------------
*/
void ftRemoveFirstLinkSegment(int linkIdx) {
	Pseg firstSeg=ftUnlinkFirstLinkSegment(linkIdx);
	firstSeg->prev = FreeSeg;
	FreeSeg = firstSeg;
}
/*
**--------------------------------------------------------------
**   Input:   segList: The segList for which the first segment
**            is to be unlinked
**   Output:  The segment that was unlinked
**   Purpose: Unlink the first segment from the segList, but do not free it
**--------------------------------------------------------------
*/
Pseg ftUnlinkFirstSegment(PSegList segList) {
	Pseg firstSeg=segList->firstSeg;
	segList->firstSeg=firstSeg->prev;

	if (segList->firstSeg == NULL) {
		segList->lastSeg = NULL;
	}
	return firstSeg;
}
/*
**--------------------------------------------------------------
**   Input:   seg: the segment to release
**   Output:  none
**   Purpose: return the segment to the FreeSeg list for  reuse
**--------------------------------------------------------------
*/
void ftReleaseSeg(Pseg seg) {
	seg->prev=FreeSeg;
	FreeSeg=seg;
}
/*
**--------------------------------------------------------------
**   Input:   segList: The segList for which the first segment
**            is to be unlinked
**   Output:  none
**   Purpose: Unlink the first segment from the segList and place
**            it back on the FreeSegs list
**--------------------------------------------------------------
*/
void ftRemoveFirstSegment(SegList *segList) {
	Pseg firstSeg=ftUnlinkFirstSegment(segList);
	ftReleaseSeg(firstSeg);
}
/*
**--------------------------------------------------------------
**   Input:   segList: The segList for which all segments are to
**            be released
**   Output:  none
**   Purpose: Free all the segments in segList
**--------------------------------------------------------------
*/
void ftFreeAllSegs(PSegList segList) {
	Pseg seg=segList->firstSeg;
	Pseg t;

	while (seg != NULL) {
		t=seg;
		seg=seg->prev;
		ftReleaseSeg(t);
	}
	segList->firstSeg=NULL;
	segList->lastSeg=NULL;
}
/*
**--------------------------------------------------------------
**   Input:   nodesToProcess: The list of nodes that still need
**            to be processed
**   Output:  the next node to process, removed from nodesToProcess
**   Purpose: Return the next node to process
**--------------------------------------------------------------
*/
SNodeTransport *ftGetNextNodeToProcess(SNodeTransportList **nodesToProcess) {
	SNodeTransportList *p=*nodesToProcess, *t;
	t=NULL;
#ifdef ENABLE_DBG_PRINT
	if (fpTrace!=NULL && printTrace) {
		fprintf(fpTrace, "NTP: ");
		while (p!=NULL) {
			p->age++;
			fprintf(fpTrace, "%s (%d, %d)%s", p->node->node->ID, ftHaveAllInflowSegs(&p->node->inflowSegs), p->age, (p->next!=NULL?",":""));
			p=p->next;
		}
		fprintf(fpTrace, "\n");
	}
#endif
	// the next node to process is determined in the following order:
	// - the first element in the nodesToProcessList that has at least one segment in each of its inflow segments lists
	// - the first element that doesn't already have an incomplete segment.
	// - the first element in the list
	p=*nodesToProcess;
	t=NULL;
	while (p!=NULL) {
		SNodeTransport *node=p->node;
		if (ftHaveAllInflowSegs(&node->inflowSegs)) {
			// this node is good to use
			if (t==NULL) {
				*nodesToProcess=p->next;
			} else {
				t->next=p->next;
			}
			free(p);
			node->toProcess=0;
			return node;
		}
		t=p;
		p=p->next;
	}
	p=*nodesToProcess;
	t=NULL;
	// if we haven't yet found a node to process, there are no more "ready"
	// It is likely the beginning of a flow loop.
	// ftGetSegmentsToMove will have to deal with this case.
	// find the first node that doesn't already have an incomplete segment
	while (p!=NULL) {
		if (p->node->incompleteSegs == NULL) {
			SNodeTransport *node=p->node;
			// this node is good to use
			if (t==NULL) {
				*nodesToProcess=p->next;
			} else {
				t->next=p->next;
			}
			free(p);
			node->toProcess=0;
			return node;
		}
		t=p;
		p=p->next;
	}
	p=*nodesToProcess;
	t=NULL;
	// if we still haven't found one, just return the first one on the list
	if (p!=NULL) {
		SNodeTransport *node=p->node;
		*nodesToProcess=p->next;
		free(p);
		node->toProcess = 0;
		return node;
	}
	return NULL;
}
/*
**--------------------------------------------------------------
**   Input:   node: The node to be added to nodesToProcess
**            nodesToProcess: the list of nodes to process
**   Output:  none
**   Purpose: Add the node to the list of elements to process if
**            it isn't already there
**--------------------------------------------------------------
*/
void ftAddNodeToProcess(SNodeTransport *node, SNodeTransportList **nodesToProcess) {
	if (node->toProcess==0) {
		SNodeTransportList *elem=(SNodeTransportList *)calloc(1, sizeof(SNodeTransportList));
		node->toProcess=1;
		elem->node=node;

		elem->next=*nodesToProcess;
		*nodesToProcess=elem;
	}
}

/*
**--------------------------------------------------------------
**   Input:   segs: The segList to add the segment to
**            seg: the segment to add
**            loc: For debugging.  See comment from ftCombineSegments
**   Output:  none
**   Purpose: add the segment to the segs list
**--------------------------------------------------------------
*/
void ftAddNodeSeg(SegList *segs, Pseg seg, int loc) {
	if (segs->firstSeg==NULL) {
		segs->firstSeg=seg;
	}
	if (segs->lastSeg != NULL) {
		segs->lastSeg->prev=seg;
	}
	segs->lastSeg=seg;
	ftCombineSegs(&segs->firstSeg, &segs->lastSeg, loc);

}
/*
**--------------------------------------------------------------
**   Input:   firstSeg: Pointer to the pointer to the first seg
**                      in a list of segs
**            lastSeg: Pointer to the pointer to the first seg
**                      in a list of segs
**            loc: For debugging.  See comment from ftCombineSegments
**   Output:  none
**   Purpose: Go through the list of segs from firstSeg to lastSeg
**            and combine any two adjacent segments that can be combined
**
** The value of this function should be looked into to see if
** it can be optimized a bit more or even eliminated.  It was
** created to handle situations where an incomplete or merged
** segment was finally made complete to see if the number of
** segments could be reduced.  There may be smarter ways to do that.
**
**--------------------------------------------------------------
*/
void ftCombineSegs(Pseg *firstSeg, Pseg *lastSeg, int loc) {
	Pseg t=NULL, p=*firstSeg;
	while (p!=NULL) {
		if (t!=NULL) {
			if (IsNormalSeg(t) && IsNormalSeg(p) && ftCombineSegments(t->v, t->c, p->v, p->c, loc)) {
				t->c = (t->c*t->v + p->c*p->v) / (t->v + p->v);                     //(2.00.11 - LR)
				checkConc(t->c);
				t->v += p->v;
				t->prev=p->prev;
				if (*lastSeg==p) {
					*lastSeg=t;
				}
				ftReleaseSeg(p);
				p=t->prev;
			} else {
				t=p;
				p=p->prev;
			}
		} else {
			t=p;
			p=p->prev;
		}
	}
}
/*
**--------------------------------------------------------------
**   Input:   oldSeg: The segment to replace
**            newSeg: The segment to replace it with
**            firstSeg: Pointer to the pointer to the first seg
**                      in a list of segs
**            lastSeg: Pointer to the pointer to the first seg
**                      in a list of segs
**   Output:  none
**   Purpose: replace an incomplete or merged segment with a
**            normal segment
**--------------------------------------------------------------
*/
void ftReplacePseg(Pseg oldSeg, Pseg newSeg, Pseg *firstSeg, Pseg *lastSeg) {
	Pseg t=NULL, p=*firstSeg;
	while (p!=NULL) {
		if (p==oldSeg) {
			newSeg->prev=p->prev;
			if (*lastSeg==p) {
				*lastSeg=newSeg;
			}
			if (t==NULL) {
				*firstSeg=newSeg;
			} else {
				t->prev=newSeg;
			}
			oldSeg->prev=NULL;
		}
		t=p;
		p=p->prev;
	}
}
/*
**--------------------------------------------------------------
**   Input:   oldSeg: The segment to replace
**            newSeg: The segment to replace it with
**            linkIdx: the link idx containing the oldSeg
**            loc: debugging.  see comment from ftCombineSegment
**   Output:  none
**   Purpose: replace an incomplete or merged segment with a
**            normal segment
**--------------------------------------------------------------
*/
void ftReplaceLinkSeg(Pseg oldSeg, Pseg newSeg, int linkIdx, int loc) {
	ftReplacePseg(oldSeg, newSeg, &FirstSeg[linkIdx], &LastSeg[linkIdx]);
	ftCombineSegs(&FirstSeg[linkIdx], &LastSeg[linkIdx], loc);
}
/*
**--------------------------------------------------------------
**   Input:   oldSeg: The segment to replace
**            newSeg: The segment to replace it with
**            segList: the seglist containing the oldSeg
**            loc: debugging.  see comment from ftCombineSegment
**   Output:  none
**   Purpose: replace an incomplete or merged segment with a
**            normal segment
**--------------------------------------------------------------
*/
void ftReplaceSeg(Pseg oldSeg, Pseg newSeg, SegList *segs, int loc) {
	ftReplacePseg(oldSeg, newSeg, &segs->firstSeg, &segs->lastSeg);
	ftCombineSegs(&segs->firstSeg, &segs->lastSeg, loc);
}

/*
**--------------------------------------------------------------
**   Input:   c: thew concentration to check
**   Output:  none
**   Purpose: Provide a common place to check if a concentration
**            is below a certain value.  Useful for debugging...
**--------------------------------------------------------------
*/
void checkConc(double c) {
	if (c > 0 && c < 1e-5) {
		fprintf(stdout, "");  fflush(stdout);
	}
}
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  none
**   Purpose: Update Mass balance data 
**--------------------------------------------------------------
*/
void ftComputeMassBalance() {
	int i;
	Pseg s;
	double mass;
	for (i=1; i<=Nnodes; i++) {
		mass=0;
		for (s=TransportNodes[i].demandSegs.firstSeg; s!=NULL; s=s->prev) {
			mass += s->v*s->c;
		}
		MB_MassRemoved[i]+=mass;
	}
	for (i=1; i<=Nlinks; i++) {
		mass=0;
		for (s=FirstSeg[i]; s!=NULL; s=s->prev) {
			mass += s->v*s->c;
		}
		MB_MassInPipes[i]=mass;
	}
	for (i=1; i<=Ntanks; i++) {
		MB_MassInTanks[i]=Tank[i].M;
	}
}

/*
**--------------------------------------------------------------
**   Input:   s: the segment to check
**            nm: an identifier to print
**   Output:  the next segment
**   Purpose: Check if a segment has already been freed.  The
**            values it is checked against are (i believe) specific
**            to the MSVC compiler
**--------------------------------------------------------------
*/
Pseg ftCheckSegment(Pseg s, char *nm) {
	// check for free'd memory
	char *cs = (char *)s;
	if(*cs==0xee && *(cs+1)==0xef) {
		fprintf(stdout,"segment has been freed (%s)!\n",nm);
		return NULL;
	}
	return s->prev;
}
/*
**--------------------------------------------------------------
**   Input:   segs: the segments to check
**            nm: an identifier to print
**   Output:  None
**   Purpose: Check all segments for being already freed
**--------------------------------------------------------------
*/
void ftCheckSegList(Pseg segs, char *nm) {
	for(Pseg s = segs; s != NULL; ) {
		Pseg next=ftCheckSegment(s,nm);
		s=next;
	}
}
/*
**--------------------------------------------------------------
**   Input:   none
**   Output:  None
**   Purpose: Check all segments for being already freed
**--------------------------------------------------------------
*/
void ftCheckForInvalidSegments() {
	int i;

	for (i=1; i<=Nnodes; i++) {
		PNodeTransport n=&TransportNodes[i];
        if(n->inflowSegs.segs != NULL) {
    		ftCheckSegList(n->inflowSegs.segs->firstSeg,"Inflow");
        }
		ftCheckSegList(n->demandSegs.firstSeg,"Demand");
		ftCheckSegList(n->throughSegs.firstSeg,"Through");
		ftCheckSegList(n->allDemandSegs.firstSeg,"AllDemand");
		ftCheckSegList(n->allThroughSegs.firstSeg,"AllThrough");
		for(PIncompleteSegmentData isd = n->incompleteSegs; isd != NULL; isd=isd->prev) {
			for(int j=0;j<n->inflow.nlinks;j++) {
				PIncompleteSeg is = isd->incomplete[j];
				if(is != NULL) {
					int i;
					i=0;
				}
			}
		}
	}

}
