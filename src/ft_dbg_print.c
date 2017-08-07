//#include "toolkit.h"

#include <stdlib.h>
#include <string.h>
#include <direct.h>
#include <sys/stat.h>

#include "ft_dbg_print.h"
#ifdef ENABLE_DBG_PRINT
FILE *fpDeltas=NULL;
FILE *fpAdj=NULL;
FILE *fpIncTrace=NULL;
FILE *demFP;
FILE *qualFP;
FILE *massFP;
FILE *fpAllAvgConc;
FILE *fpAllInstConc;
FILE *fpMassBalance;
FILE *fpcombine;
FILE *incfp=NULL;
FILE *fpTrace = NULL;
FILE *fpi = NULL;
FILE *warnfp=NULL;
FILE *fpOwners=NULL;
FILE *fpCompleteInc=NULL;
FILE *fpUpdateIncSeg=NULL;
FILE *fpLoss=NULL;
FILE *fpTank=NULL;
FILE *fpd=NULL;
FILE *fpInOut = NULL;

FILE *fpDocNodes;
int *dbgPrintNodes;
int *dbgPrintLinks;
int printTrace;
int maxIn=0;
int maxOut=0;


int printAllLinkSegs=0;
int openMassBalanceData=0;
int *nodesToWrite;
int *linksToWrite;

double   *MB_CumMassAdded,   /* nnodes */
         *MB_CumMassRemoved; /* nnodes */

int openInfo=0;
int openBalance=0;//1;
int openFlow=0;
int openWarn=0;
int openDetail=0;
int openAllIncSegs=0;

void ftOpenDebugFiles(long stime) {
	// files for debugging
	char fn[256];
	int openTrace=0;
	int openUpdateIncSeg=0;
	int openFiles=0;
	int openCombine=0;
	int openIncomplete=0;
	int openCompleteInc=0;
	int openOwners=0;
	int openLoss=0;//1;
	int openDeltas=0;//1;
	int openAdj=0;
	int openTank=0;
	int openDemands=0;
	int openIncTrace=0;
	int openMSegs=0;
	int openDocFiles = 0;
	int openInOut = 0;

	openWarn=0;//1;
	openDetail=0;

	openFiles=0;
//	openFiles = stime == 89928;
	openInOut = stime == 89928 || stime==3600;
	openTrace = openInOut;
	openWarn = 0;
	openDetail=0;
	openAllIncSegs=0;
	printAllLinkSegs=0;
	openMassBalanceData=0;

	//	printAllLinkSegs=stime >= 1336*60 && stime <= 1360*60;
//	printAllLinkSegs|=openFiles;
	openMassBalanceData|=openFiles;
	openMassBalanceData = 1;
	openUpdateIncSeg|=openFiles;

    openUpdateIncSeg=0;

    openIncomplete|=openFiles;
	openCompleteInc|=openFiles;
	openBalance|=openFiles;
	openCombine|=openFiles;
	openTrace|=openFiles;
	openInfo|=openFiles;
	openOwners|=openFiles;
	openLoss|=openFiles;
	openDeltas|=openFiles;
	//	openAdj|=openFiles;
	openIncTrace|=openFiles;
	openWarn|=openFiles;
	openTank|=openFiles;
	openDemands|=openFiles;
	openFlow|=openFiles;
	openAllIncSegs|=openFiles;
	if (openFiles || openCombine || openWarn || openMassBalanceData ||
		openCompleteInc || openUpdateIncSeg || openTrace || openInfo ||
		openOwners || openBalance || openLoss || openDeltas || openAdj ||
		openIncTrace || openTank || openDemands || openFlow || openInOut) {
		_mkdir("trace");
	}
//	openDocFiles = stime == 3600;
	char *nodeids[] = {
		"2025.3_ND", "2025.3_NU", "28010133", "28010134", "28010151", "28010161", "28010163", "28010288", "28010333", "28010336", "28010475",
		"28010479", "28020269", "28020270", "28020273", "28020274", "28020275", "28020276", "28020277", "28020279", "28020281", "2802", "282",
		"28020283", "28020289", "28020297", "28020298", "28020299", "28020299.3_ND", "28020299.3_NU", "28020300", "28020301", "28020302",
		"28020303", "28020304", "28020304.3_ND", "28020304.3_NU", "28020319", "28020320", "28020321", "28020322", "28020323", "28020337",
		"28020341", "28020342", "28020482", "28080001", "28080002", "28080003", "28080004", "28080005", "28080006", "28080007", "28080008",
		"28080009", "2808001", "28080010", "28080011", "28080012", "28080014", "28080015", "28080016", "28080017", "280800", "7.3_ND",
		"28080017.3_NU", "28080018", "28080019", "28080021", "EQ009329_ND", "EQ009329_NU", "EQ009332_ND", "EQ009332_NU", "EQ009335_ND",
		"EQ009335_NU", "EQ009337_ND", "EQ009337_NU", "EQ009340_ND", "EQ009340_NU", "EQ009342_ND", "EQ009342_NU", "EQ009346_ND",
		"EQ009346_NU", "28080020", "28080022", NULL };
	if (nodeids != NULL) {
		if (dbgPrintNodes == NULL) {
			dbgPrintNodes = (int *)calloc(Nnodes + 1, sizeof(int));
			dbgPrintLinks = (int *)calloc(Nlinks + 1, sizeof(int));
		}
		memset(dbgPrintLinks, 0, Nnodes*sizeof(int));
		memset(dbgPrintNodes, 0, Nnodes*sizeof(int));
		for (int i = 0; nodeids[i] != NULL; i++) {
			int idx;
			PNodeTransport node;
			ENgetnodeindex(nodeids[i], &idx);
			dbgPrintNodes[idx] = 1;
			node = &TransportNodes[idx];
			for (int j = 0; j < node->inflow.nlinks; j++) {
				int idx = node->inflow.linkIdxs[j];
				if (idx > 0) {
					dbgPrintLinks[idx] = 1;
				}
			}
			for (int j = 0; j < node->outflow.nlinks; j++) {
				int idx = node->outflow.linkIdxs[j];
				if (idx > 0) {
					dbgPrintLinks[idx] = 1;
				}
			}
		}
	}
	sprintf(fn, "trace/combineseg/combineseg_%06d.txt", stime);
	if (openCombine) { _mkdir("trace/combineseg"); fpcombine=fopen(fn, "w");
	fprintf(fpcombine,"Seq\tLoc\tV1\tC1\tV2\tC2\tType\tV1+V2\tABS(C1-C2)\tNewC\tNewC-C1\tCombine\n"); }

	sprintf(fn, "trace/demands/demands_%06d.txt", stime);
	if (openDemands) { _mkdir("trace/demands"); fpd=fopen(fn, "w"); }

	sprintf(fn, "trace/loss/loss_%06d.txt", stime);
	if (openLoss) { _mkdir("trace/loss"); fpLoss=fopen(fn,"w");
	fprintf(fpLoss,"Node\tVol Added\tVol Lost\tMassLost\n"); }

	sprintf(fn, "trace/deltas/deltas_%06d.txt", stime);
	if (openDeltas) { _mkdir("trace/deltas"); fpDeltas=fopen(fn,"w");
	fprintf(fpDeltas, "Location\tdelta\n"); }

	sprintf(fn, "trace/adj/adj_%06d.txt", stime);
	if (openAdj) { _mkdir("trace/adj"); fpAdj=fopen(fn, "w");
	fprintf(fpAdj,"Type\tSource\tseq\tvol\tvol*adj\tadj\tpv-pf-vol\tpv-pf-vol*adj\tapply adj\tsources\n"); }

	sprintf(fn, "trace/tanks/tanks_%06d.txt", stime);
	if (openTank) { _mkdir("trace/tanks"); fpTank=fopen(fn, "w");
	fprintf(fpTank, "Time\tFlowDir\tTank\tType\tinit vol\tD\tdvol\tNewVol\tMassIn\toutVol\tOutConc\tOutMass\n"); }

	sprintf(fn, "trace/inc_trace/inc_trace_%06d.txt", stime);
	if (openIncTrace) { _mkdir("trace/inc_trace"); fpIncTrace=fopen(fn, "w");
	fprintf(fpIncTrace, "Loc\tSequence\tType\tobj\tData\tincseg vol\tpartIndex\tpartVol\tpartFlowed\n"); }

	// open trace file
	sprintf(fn, "trace/trace/trace_%06d.txt", stime);
	if (openTrace) { _mkdir("trace/trace"); fpTrace=fopen(fn, "w"); }
	// open incomplete segments file
	sprintf(fn, "trace/incseg/incseg_%06d.txt", stime);
	if (openIncomplete) { _mkdir("trace/incseg"); incfp=fopen(fn, "w"); }
	// open owners trace file
	sprintf(fn, "trace/owners/owners_%06d.txt", stime);
	if (openOwners) { _mkdir("trace/owners"); fpOwners=fopen(fn, "w"); }
	// open owners trace file
	sprintf(fn, "trace/comp_inc/comp_inc_%06d.txt", stime);
	if (openCompleteInc) { _mkdir("trace/comp_inc"); fpCompleteInc=fopen(fn, "w"); }
	// open incomplete segments file
	sprintf(fn, "trace/updateincseg/updateincseg_%06d.txt", stime);
	if (openUpdateIncSeg) { _mkdir("trace/updateincseg"); fpUpdateIncSeg = fopen(fn, "w"); }
	// open in/out files
	sprintf(fn, "trace/inout/inout_%06d.txt", stime);
	if (openInOut) { _mkdir("trace/inout"); fpInOut = fopen(fn, "w"); }

	if (openDocFiles) {
		sprintf(fn, "trace/doc/nodes/nodes_%06d.txt", stime);
		_mkdir("trace");
		_mkdir("trace/doc");
		_mkdir("trace/doc/nodes");
		fpDocNodes = fopen(fn, "w");

	}
}

void ftOpenInfoFile(long stime, long dt) {
	char fn[256];
	// open info file
	sprintf(fn, "trace/info/info_%06d.txt", stime);
	if (openInfo) { _mkdir("trace/info"); fpi=fopen(fn, "w"); }
	if (fpi != NULL) {
		fprintf(fpi, "dt\t%d\tATOL\t%.8e\n", dt, FT_ATOL);
		ftPrintInfoHeader(fpi, maxIn, maxOut);
		ftPrintInflowHeader(fpi, maxIn);
		ftPrintSegmentsToMoveHeader(fpi, maxIn);
		ftPrintDemandHeader(fpi);
		ftPrintFlowLinksHeader(fpi, maxOut);
		ftPrintMoveDSHeaders(fpi);
	}

}


void ftOpenMBFiles() {
	int errcode=0,i;
	MB_CumMassAdded = (double *)calloc(Nnodes+1, sizeof(double));
	MB_CumMassRemoved = (double *)calloc(Nnodes+1, sizeof(double));
	ERRCODE(MEMCHECK(MB_CumMassAdded));
	ERRCODE(MEMCHECK(MB_CumMassRemoved));
	fpMassBalance=fopen("MassBalance.txt", "w");
	fprintf(fpMassBalance, "Time\tAdded\tRemoved\ttanks\tpipes\n");

	demFP=fopen("EN_dem.txt", "w");
	qualFP=fopen("EN_qual.txt", "w");
	massFP=fopen("EN_mass.txt", "w");

	fprintf(demFP, "Time");
	fprintf(massFP, "Time");
	fprintf(qualFP, "Time");
	for (i=1; i<=Nnodes; i++) {
		fprintf(demFP, "\t%s", Node[i].ID);
		fprintf(massFP, "\t%s", Node[i].ID);
		fprintf(qualFP, "\t%s", Node[i].ID);
	}
	fprintf(massFP, "\t");

	for (i=1; i<=Nlinks; i++) {
		fprintf(massFP, "\t%s", Link[i].ID);
	}
	fprintf(demFP, "\n");
	fprintf(massFP, "\n");
	fprintf(qualFP, "\n");
}
void ftdAccumulateMassAdded(int i) {
	if (MB_CumMassAdded != NULL) {  // will be non-null if printing MB files.
		MB_CumMassAdded[i]+=MB_MassAdded[i];
		MB_CumMassRemoved[i]+=MB_MassRemoved[i];
	}
}

void appendConcFile(FILE *fp, double *c) {
	int i;
	fprintf(fp, "%d", Qtime);
	for (i=1; i<=Nnodes; i++) {
		fprintf(fp, "\t%g", NodeQual[i]*Ucf[QUALITY]);
	}
	fprintf(fp, "\n");
}
FILE *openAllConcFile(char *fn) {
	FILE *fp=fopen(fn, "w");
	int i;

	fprintf(fp, "Time");
	for (i=1; i<=Nnodes; i++) {
		fprintf(fp, "\t%s", Node[i].ID);
	}
	fprintf(fp, "\n");
	return fp;
}


// debug printing / tracing methods
void ftPrintInfoHeader(FILE *fpi, int maxIn, int maxOut) {
	int l;
	if (fpi==NULL) return;
	fprintf(fpi, "INFO\tNodeID\tdemand\ttotal vol inflow\ttotal vol outflow\tmassadded\t");
	for (l=0; l<maxIn; l++) {
		fprintf(fpi, "\tlink idx\tlink id\tflow ratio\tflow vol\tvol flowed\t");
	}
	fprintf(fpi, "\t");
	for (l=0; l<maxOut; l++) {
		fprintf(fpi, "\tlink idx\tlink id\tlink vol\tDS node\tflow ratio\tflow vol\tvol flowed\t");
	}
	fprintf(fpi, "\n");
	fflush(fpi);
}
void ftPrintNodeInfo(FILE *fpi, PNodeTransport n, int maxIn, int maxOut) {
	PNodeLinks in=&n->inflow;
	PNodeLinks out=&n->outflow;
	int l;
	if (fpi==NULL) return;
	fprintf(fpi, "\n");
	fprintf(fpi, "INFO\t%s\t%.8f\t%.8f\t%.8f\t%.8f\t", n->node->ID, n->demand, n->totalVolInflow, n->totalVolOutflow, n->massadded);
	for (l=0; l<in->nlinks; l++) {
		int linkIdx=in->linkIdxs[l];
		if (linkIdx >= 0) {
			fprintf(fpi, "\t%d\t%s\t%.8f\t%.8f\t%.8f\t", in->linkIdxs[l], Link[in->linkIdxs[l]].ID, in->flowRatio[l], in->flowVolume[l], in->volumeFlowed[l]);
		} else {
			fprintf(fpi, "\t\t\t\t\t\t");
		}
	}
	for (; l<maxIn; l++) {
		fprintf(fpi, "\t\t\t\t\t\t");
	}
	fprintf(fpi, "\t");
	for (l=0; l<out->nlinks; l++) {
		fprintf(fpi, "\t%d\t%s\t%.8f\t%s\t%.8f\t%.8f\t%.8f\t", out->linkIdxs[l], Link[out->linkIdxs[l]].ID, LINKVOL(out->linkIdxs[l]), Node[DOWN_NODE(out->linkIdxs[l])].ID, out->flowRatio[l], out->flowVolume[l], out->volumeFlowed[l]);
	}
	for (; l<maxOut; l++) {
		fprintf(fpi, "\t\t\t\t\t\t\t\t");
	}
	fprintf(fpi, "\n");
	fflush(fpi);
}
void ftPrintDemandHeader(FILE *fpi) {
	if (fpi==NULL) return;
	fprintf(fpi, "DEMAND\tRaw vol ratio\tädj vol ratio\tDem From This\tseg->v\tseg v to move\n");

}
void ftPrintFlowLinksHeader(FILE *fpi, int maxOut) {
	int l;
	if (fpi == NULL) return;
	fprintf(fpi, "FLOW2LINKS\t");
	for (l=0; l<maxOut; l++) {
		fprintf(fpi, "\tvl\tseg->c\tadj seg v\tadj seg c\tAdd or NEW\tadj seg v\tadj seg c\t");
	}
	fprintf(fpi, "\n");
	if (fpi != NULL) { fprintf(fpi, "FLOW2LINKS-(INC|MERGED)\t"); }
	for (l=0; l<maxOut; l++) {
		fprintf(fpi, "\tv moved\tc moved\t");
	}
	fprintf(fpi, "\n");
}
void ftPrintMoveDSHeaders(FILE *fpi) {
	if (fpi == NULL) return;
	fprintf(fpi, "MOVEDSINFO\t\tLink ID\tDS Node\texcess vol\tvlv type\tvlv adj ratio\ttol ratio\tadj excess vol\n");
	fprintf(fpi, "MOVEDS\t\t\texcess vol still to move\tseg v\tseg c\tvolMoved\ttol volMoved\tseg v remain\texcessVol remain\tds v\tds c\n");
}
void ftPrintInflowHeader(FILE *fpi, int maxIn) {
	int l;
	fprintf(fpi, "INFLOW\t\t\t\t");
	fprintf(fpi, "\tExt Inflow vol\tExt Inflow conc\t");
	for (l=0; l<maxIn; l++) {
		fprintf(fpi, "\tInflow vol\tInflow conc\t");
	}
	fprintf(fpi, "\n");
}
void ftPrintSegmentsToMoveHeader(FILE *fpi, int maxIn) {
	int l;
	if (fpi == NULL) return;
	fprintf(fpi, "SEGSTOMOVE\t\tvol\tconc\n");
	fprintf(fpi, "SEGSTOMOVEINT\t");
	for (l=0; l<maxIn; l++) {
		fprintf(fpi, "\tfs vol\tfs c\tfs ratio\tthisv\t");
	}
	fprintf(fpi, "\tmin vol\t");
	for (l=0; l<maxIn; l++) {
		fprintf(fpi, "\tthisv\ttol thisv\tfs vol\t");
	}
	fprintf(fpi, "\tsegvol\tsegmass\tsegc\tany\tall\n");
}
void ftPrintSegmentsToMove(FILE *fpi, PSegList segs) {
	Pseg seg=segs->firstSeg;
	if (fpi == NULL) return;
	while (seg != NULL) {
		fprintf(fpi, "SEGSTOMOVE\t\t%.8f\t%.8f\n", seg->v, seg->c);
		seg=seg->prev;
	}
	fflush(fpi);

}
void ftPrintAllLinkSegs(long stime, char *mode) {
	int i;
	static int first=1;
	FILE *fp;
	FILE *fpc;
	FILE *fpl;
	char fn[256];
	sprintf(fn, "trace/linksegs-%s", mode);
	_mkdir(fn);
	sprintf(fn, "trace/linksegs-%s/linksegs_%07d.txt", mode, stime);
	fp=fopen(fn, "w");
	sprintf(fn, "trace/linksegs-%s/clinksegs_%07d.txt", mode, stime);
	fpc=fopen(fn, "w");
	for (i=1; i<=Nlinks+Ntanks; i++) {
		Pseg s;
		char *id;
		int num;
		double vt, mass;
		vt=0;
		mass=0;
		num=0;
		id=i<=Nlinks?Link[i].ID:Node[Tank[i-Nlinks].Node].ID;
		sprintf(fn, "trace/linksegs-%s/segments_%s.txt", mode, id);
		fpl=fopen(fn, first?"w":"a");
		for (s=FirstSeg[i]; s!=NULL; s=s->prev) {
			num++;
			vt+=s->v;
			mass+=s->v*s->c;
		}
		fprintf(fp, "%d\t%s\t%d\t%f\t%f\t%f\n", i, id, num, vt, vt!=0?mass/vt:0, mass);
		fprintf(fpl, "%ld\t%d\t%s\t%d\t%g\t%g\t%g\t", stime, i, id, num, vt, vt!=0?mass/vt:0, mass);
		if (mass>0) {
			fprintf(fpc, "%d\t%s\t%d\t%f\t%f\t%f\n", i, id, num, vt, vt!=0?mass/vt:0, mass);
		}
		for (s=FirstSeg[i]; s!=NULL; s=s->prev) {
			if (mass>0) {
				fprintf(fpc, "\t\t\t\t\t\t%f\t%f\t%f", s->v, s->c, s->v*s->c);
				if (s->c>0) {
					fprintf(fpc, "\tx");
				}
				fprintf(fpc, "\n");
			}
			fprintf(fp, "\t\t\t\t\t\t%f\t%f\t%f", s->v, s->c, s->v*s->c);
			if (s->c>0) {
				fprintf(fp, "\tx");
			}
			fprintf(fp, "\n");
			fprintf(fpl, "\t\t%g\t%g\t%g", s->v, s->c, s->v*s->c);
			if (s->c>0) {
				fprintf(fpl, "\tx");
			}
		}
		fprintf(fpl, "\n");
		fclose(fpl);
		first=0;
	}
	fclose(fp);
	fclose(fpc);
}
void ftPrintInflowOutflow(long stime) {
	FILE *fpf, *fpfw;
	int i;
	char fn[128];
	if (!openBalance)
		return;
	sprintf(fn, "trace/balance/balance_%06d.txt", stime);
	_mkdir("trace/balance");
	fpf=fopen(fn, "w");
	sprintf(fn, "trace/balance/balance_w_%06d.txt", stime);
	fpfw=fopen(fn, "w");
	fprintf(fpf, "NodeID\terr\tsrc\tstatus\tinflow\tdemand\toutflow\toutflow adj\n");
	fprintf(fpfw, "NodeID\terr\tsrc\tstatus\tinflow\tdemand\toutflow\toutflow adj\n");
	for (i=1; i<=Nnodes; i++) {
		double err;
		char *emsg;
		int tol, linkIdx;
		int type, isvalve;
		double inf, outf, dem;
		int status;

		PNodeTransport n=&TransportNodes[i];
		inf=n->totalVolInflow;
		outf=n->totalVolOutflow;
		dem=n->demand;
		err=ABS(n->totalVolInflow-n->demand-n->totalVolOutflow);
		tol=ftWithinTolerance(err, 0, FT_ATOL, FT_RTOL);
		isvalve=0;
		if (n->inflow.linkIdxs != NULL) {
			linkIdx=n->inflow.linkIdxs[0];
			if (linkIdx > 0) { // linkIdx will be -1 for tanks or reservoirs with inflow, so ignore them here
				type=Link[linkIdx].Type;
				isvalve=(type==PRV || type==PSV || type==PBV || type==FCV || type==TCV || type==GPV);
			}
		}
		status=0;
		if (inf==0)  status |= 0x10;
		if (outf==0) status |= 0x08;
		if (dem==0)  status |= 0x04;
		if (err==0)  status |= 0x02;
		if (tol)     status |= 0x01;

		if (err == 0) {
			emsg="ZERO";
		} else if (inf ==0 && (outf>0 || dem > 0)) {
			emsg="NO_INF";
		} else if (outf==0 && inf-dem > 0) {
			emsg="NO_OUTF";
		} else {
			if (tol) {
				emsg="OK_TOL";
			} else {
				if (isvalve) {
					emsg="VLV";
				} else {
					emsg="OUT_TOL";
				}
			}
		}
		fprintf(fpf, "%s\t%.8f\t%s\t%d\t%.12f\t%.12f\t%.12f\t%.12f\n", n->node->ID, err, emsg, status, n->totalVolInflow, n->demand, n->totalVolOutflow, n->outflowAdjustment);
		if (fpfw != NULL && nodesToWrite && nodesToWrite[i]) {
			fprintf(fpfw, "%s\t%.12f\t%s\t%d\t%.12f\t%.12f\t%.12f\t%.12f\n", n->node->ID, err, emsg, status, n->totalVolInflow, n->demand, n->totalVolOutflow, n->outflowAdjustment);
		}
	}
	fclose(fpf);
	fclose(fpfw);
}
void ftPrintWarningsFile(long stime) {
	FILE *warnfp;
	int i, j;
	char fn[128];
	if (!openWarn)
		return;
	sprintf(fn, "trace/warn/warn_%06d.txt", stime);
	_mkdir("trace/warn");
	warnfp=fopen(fn, "w");
	for (i=1; i<=Nnodes; i++) {
		double dem, through;
		PNodeTransport n=&TransportNodes[i];

		// check inflow expected vs actual
		for (j=0; j<n->inflow.nlinks; j++) {
			double expected=n->inflow.flowVolume[j];
			double actual=n->inflow.volumeFlowed[j];
			if (expected != actual) {
				int linkIdx=n->inflow.linkIdxs[j];
				char *linkID;
				if (linkIdx < 0) {
					linkID="External";
				} else {
					linkID=Link[linkIdx].ID;
				}
				fprintf(warnfp, "INFLOWVOL\t%s\t%s\t%.10g\t%.10g\t%.10g\n", n->node->ID, linkID, expected, actual, expected-actual);
			}
		}
		// check outflow expected vs actual
		for (j=0; j<n->outflow.nlinks; j++) {
			double expected=n->outflow.flowVolume[j];
			double actual=n->outflow.volumeFlowed[j];
			if (expected != actual) {
				int linkIdx=n->outflow.linkIdxs[j];
				char *linkID;
				if (linkIdx < 0) {
					linkID="External";
				} else {
					linkID=Link[linkIdx].ID;
				}
				fprintf(warnfp, "OUTFLOWVOL\t%s\t%s\t%.10g\t%.10g\t%.10g\n", n->node->ID, linkID, expected, actual, expected-actual);
			}
		}
		if (ftHaveAnyInflowSegs(&n->inflowSegs)) {
			int l, num=0;
			double vol=0;
			for (l=0; l<n->inflow.nlinks; l++) {
				Pseg s=n->inflowSegs.segs[l].firstSeg;
				while (s!=NULL) {
					vol+=s->v;
					num++;
					s=s->prev;
				}
			}
			fprintf(warnfp, "INFLOWSEGS\t%s\t%d\t%f\n", n->node->ID, num, vol); fflush(warnfp);
		}
		//		if(n->incompleteSegs != NULL) {
		//			PIncompleteSegmentData is=n->incompleteSegs;
		//			fprintf(warnfp,"INCOMPLETE\t%s\n",n->node->ID); fflush(warnfp);
		//			while(is != NULL) {
		//				ftPrintIncompleteSegments(is->incomplete,n->inflow.nlinks,warnfp);
		//				is=is->prev;
		//			}
		//		}
		dem=ftSumVolume(&n->demandSegs);
		through=ftSumVolume(&n->throughSegs);
		if (!ftWithinTolerance(dem, n->demand, FT_ATOL, FT_RTOL)) {
			fprintf(warnfp, "DEMAND\t%s\t%g\t%g\n", n->node->ID, dem, n->demand); fflush(warnfp);
		}
		if (!ftWithinTolerance(through, n->totalVolOutflow, FT_ATOL, FT_RTOL)) {
			fprintf(warnfp, "THROUGH\t%s\t%g\t%g\n", n->node->ID, through, n->totalVolOutflow); fflush(warnfp);
		}
		if (n->visited==0 && (n->inflow.nlinks > 0 || n->outflow.nlinks > 0)) {
			fprintf(warnfp, "NOTVISITED\t%s\t%d\t%d\t%g\t%g\t%g\n", n->node->ID, n->inflow.nlinks, n->outflow.nlinks,
				n->demand, n->totalVolInflow, n->totalVolOutflow);
			fflush(warnfp);
		}
		if (!ftWithinTolerance(n->totalVolInflow, n->demand + n->totalVolOutflow, FT_ATOL, FT_RTOL)) {
			fprintf(warnfp, "FLOW_IMBALANCE\t%s\t%.5f\t%.5f\t%.5f\t%.5f\n", n->node->ID,
				n->demand, n->totalVolInflow, n->totalVolOutflow, n->totalVolInflow - n->demand - n->totalVolOutflow);
			fflush(warnfp);
		}
	}
	for (i=1; i<=Nlinks; i++) {
		double lv=LINKVOL(i);
		double lwv=ftComputeSegmentWaterVolume(FirstSeg[i]);
		if (!ftWithinTolerance(lv, lwv, FT_ATOL, FT_RTOL)) {
			fprintf(warnfp, "LINK_VOL\t%s\t%.8f\t%.8f\n", Link[i].ID, lv, lwv);
		}
	}

}
void ftCloseDebugFiles() {
	if (warnfp!=NULL) { fclose(warnfp); warnfp=NULL; }
	if (fpcombine!=NULL) { fclose(fpcombine); fpcombine=NULL; }
	if (fpTrace!=NULL) { fclose(fpTrace); fpTrace=NULL; }
	if (incfp!=NULL) { fclose(incfp); incfp=NULL; }
	if (fpOwners!=NULL) { fclose(fpOwners); fpOwners=NULL; }
	if (fpCompleteInc!=NULL) { fclose(fpCompleteInc); fpCompleteInc=NULL; }
	if (fpUpdateIncSeg!=NULL) { fclose(fpUpdateIncSeg); fpUpdateIncSeg=NULL; }
	if (fpi!=NULL) { fclose(fpi); fpi=NULL; }
	if (fpLoss!=NULL) { fclose(fpLoss); fpLoss=NULL; }
	if (fpDeltas!=NULL) { fclose(fpDeltas); fpDeltas=NULL; }
	if (fpAdj!=NULL) { fclose(fpAdj); fpAdj=NULL; }
	if (fpIncTrace!=NULL) { fclose(fpIncTrace); fpIncTrace=NULL; }
	if (fpTank!=NULL) { fclose(fpTank); fpTank=NULL; }
	if (fpDocNodes != NULL) { fclose(fpDocNodes); fpDocNodes = NULL; }
	if (fpInOut != NULL) { fclose(fpInOut); fpInOut = NULL; }
}
int getbin(double val, double bins[], int nbins) {
	int i;
	for (i=1; i<nbins; i++) {
		if (val<=bins[i]) return i-1;
	}
	return nbins-1;
}
void printSegStats(long stime, char *str) {
	FILE *fp;
	char fn[256];
	int nbins=9;
	double bins[]={ 0, 0, .001, .01, .1, 1, 10, 100, 1000 };
	int link_counts[]={ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	int demand_counts[]={ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	int through_counts[]={ 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	int i;

	sprintf(fn, "trace/segs/segs_%06d_%s.txt", stime, str);
	_mkdir("trace/segs");
	fp=fopen(fn, "w");

	for (i=1; i<=Nlinks; i++) {
		Pseg seg=FirstSeg[i];
		while (seg != NULL) {
			int bi=getbin(seg->v, bins, nbins);
			link_counts[bi]++;
			seg=seg->prev;
		}
	}
	fprintf(fp, "Link Segments\n");
	fprintf(fp, "Min\tmax\tcount\n");
	for (i=0; i<nbins; i++) {
		if (i<nbins-1) {
			fprintf(fp, "%f\t%f\t%d\n", bins[i], bins[i+1], link_counts[i]);
		} else {
			fprintf(fp, "%f\t\t%d\n", bins[i], link_counts[i]);
		}
	}
	if (TransportNodes != NULL) {
		for (i=1; i<=Nnodes; i++) {
			Pseg seg=TransportNodes[i].allDemandSegs.firstSeg;
			while (seg != NULL) {
				int bi=getbin(seg->v, bins, nbins);
				demand_counts[bi]++;
				seg=seg->prev;
			}
			seg=TransportNodes[i].allThroughSegs.firstSeg;
			while (seg != NULL) {
				int bi=getbin(seg->v, bins, nbins);
				through_counts[bi]++;
				seg=seg->prev;
			}
		}
		fprintf(fp, "Demand Segments\n");
		fprintf(fp, "Min\tmax\tcount\n");
		for (i=0; i<nbins; i++) {
			if (i<nbins-1) {
				fprintf(fp, "%f\t%f\t%d\n", bins[i], bins[i+1], demand_counts[i]);
			} else {
				fprintf(fp, "%f\t\t%d\n", bins[i], demand_counts[i]);
			}
		}
		fprintf(fp, "Through Segments\n");
		fprintf(fp, "Min\tmax\tcount\n");
		for (i=0; i<nbins; i++) {
			if (i<nbins-1) {
				fprintf(fp, "%f\t%f\t%d\n", bins[i], bins[i+1], through_counts[i]);
			} else {
				fprintf(fp, "%f\t\t%d\n", bins[i], through_counts[i]);
			}
		}
	}
	fclose(fp);
}
void ftPrintAllIncompleteSegments(PNodeTransport n) {
	PIncompleteSegmentData p=n->incompleteSegs;
	while (p != NULL) {
		ftPrintIncompleteSegments(p->incomplete, n->inflow.nlinks, fpCompleteInc);
		p=p->prev;
	}
}

void ftPrintIncompleteSegments(PIncompleteSeg *s, int n, FILE *fp) {
	int i;
	for (i=0; i<n; i++) {
		ftPrintIncompleteSegment(s[i], 0, fp, NULL);
	}
}
void ftPrintIncompleteSegment(PIncompleteSeg s, int level, FILE *fp, char *indent2) {
	int i;
    static char *indent=NULL;
    static char indentLen=0;
	if (s==NULL) return;
	if (fp==NULL) return;
    if(indentLen < level+1) {
        indent = (char *)realloc(indent,(level+1)*sizeof(char));
        indentLen=level+1;
    }
    indent[level]=0;
    if(level>0) {
        indent[level-1]='\t';
    }
//	for (i=0; i<level; i++) { strcat(indent, "\t"); }
	if (s->node != NULL) { fprintf(fp, "%sOwner:\t%s\n", indent, s->node->node->ID); }
	fprintf(fp, "%saddr\tv\tc\tratio\tvalve ratio\tmassadded\n", indent);
	fprintf(fp, "%s0x%08x\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\n", indent, (int)s,s->v, s->c, s->ratio, s->adjustment, s->massadded);
	fprintf(fp, "%s%s\t%d\t%d\t", indent, OwnerTypes[s->owner.type], s->owner.source, s->owner.sequence);
	switch (s->owner.type) {
	case SO_NONE:
		fprintf(fp, "\n");
		break;
	case SO_THROUGH:
	case SO_DEMAND:
	case SO_INFLOW:
	{
					  SNodeTransport *node=s->owner.object.node;
					  fprintf(fp, "Node id\t%s\n", node->node->ID);
	}
		break;
	case SO_LINK:
	{
					int linkIdx = s->owner.object.linkIdx;
					fprintf(fp, "Link id\t%s\n", Link[linkIdx].ID);
	}
		break;
	case SO_MERGED:
	{
					  fprintf(fp, "0x%08x\n",s->owner.object.mergedSeg);
	}
		break;
	case SO_INCOMPLETE_PART:
	{
							   PSegOwner owner=&s->owner;
							   PIncSegPart incpart;
							   PIncompleteSeg incpartOwner;

							   incpart = owner->object.incSeg;
							   incpartOwner=owner->data.incSeg;
							   fprintf(fp, "Incomplete part owner\t%s\t0x%08x\t%.10g\t%.10g\t%.10g", incpartOwner->node->node->ID, (int)incpartOwner, incpartOwner->v, incpart->volume, incpart->flowed);
							   if (incpart->segs == NULL) {
								   fprintf(fp, "\tNoSegs\n");
							   } else {
								   Pseg pp;
								   for (pp=incpart->segs; pp!=NULL; pp=pp->prev) {
									   fprintf(fp, "\t0x%08x\t%.10g", (int)pp,pp->v);
								   }
								   fprintf(fp, "\n");
							   }
	}
		break;
	default:
		break;
	}
	// print parts...
	fprintf(fp, "%sratio\tv\tflowed\tc\tSegs\tv\tc\n", indent);
	for (i=0; i<s->nparts; i++) {
		/*
		double ratio;       // ratio of this component's voluime to the total segment volume
		double volume;      // the volume of this component
		double flowed;      // how much volume has flowed
		double conc;        // the concentration of this component
		Pseg   segs;        // incomplete or merged segments that make up this part.
		*/
		Pseg seg;
		fprintf(fp, "%s%.10g\t%.10g\t%.10g\t%.10g\n", indent, s->parts[i].ratio, s->parts[i].volume, s->parts[i].flowed, s->parts[i].conc);
		fflush(fp);
		seg=s->parts[i].segs;
		if(s->parts[i].conc==-2) {
			printf("");
		}
		while (seg != NULL) {
			fprintf(fp, "%s\t\t\t\t\t%.10g\t%.10g", indent, seg->v, seg->c);
            if(IsIncompleteSeg(seg)) {
                fprintf(fp,"\t0x%08x\n",seg);
            } else if(IsMergedSeg(seg)) {
                fprintf(fp,"\t0x%08x\n",seg);
            }
			fflush(fp);
			seg=seg->prev;
		}
	}
	for (i=0; i<s->nchildren; i++) {
		ftPrintIncompleteSegment(&s->children[i], level+1, fp,indent);
	}
    if(level>0) {
        indent[level-1]=0;
//    } else {
//        free(indent);
    }
	fflush(fp);
}
void ftPrintFlow(long stime, long dt) {
	FILE *fpf=NULL;
	FILE *fpfw=NULL;
	int i;
	char fn[256];

	sprintf(fn, "trace/flow/flow_%06d.txt", stime);
	if (openFlow) { _mkdir("trace/flow"); fpf=fopen(fn, "w"); }
	if (openFlow && linksToWrite != NULL) {
		sprintf(fn, "trace/flow/flow_w_%06d.txt", stime);
		fpfw=fopen(fn, "w");
	}
	if (fpf!=NULL) { fprintf(fpf, "i\tID\tLinkVolume\tABS(Q[i]) (fv)\torigv\tvol\tus idx\tus id\tds idx\tds id\n"); }
	if (fpfw!=NULL) { fprintf(fpfw, "i\tID\tLinkVolume\tABS(Q[i]) (fv)\torigv\tvol\tus idx\tus id\tds idx\tds id\n"); }
	for (i=1; i<=Nlinks; i++) {
		double origv, v;
		int ds=DOWN_NODE(i);
		origv=ABS(Q[i])*dt;             /* Flow volume */
		if (ABS(Q[i]) < FT_QZERO) {
			v=0;
		} else {
			v=origv;
		}
		if (fpf!=NULL) { fprintf(fpf, "%d\t%s\t%f\t%.10f\t%.10f\t%.10f\t%d\t%s\t%d\t%s\n", i, Link[i].ID, LINKVOL(i), ABS(Q[i]), origv, v, UP_NODE(i), Node[UP_NODE(i)].ID, ds, Node[ds].ID); }
		if (fpfw!=NULL && linksToWrite[i]) { fprintf(fpfw, "%d\t%s\t%f\t%.10f\t%.10f\t%.10f\t%d\t%s\t%d\t%s\n", i, Link[i].ID, LINKVOL(i), ABS(Q[i]), origv, v, UP_NODE(i), Node[UP_NODE(i)].ID, ds, Node[ds].ID); }
	}
	if (fpf!=NULL) { fclose(fpf); fpf=NULL; }
	if (fpfw!=NULL) { fclose(fpfw); fpfw=NULL; }
}
void ftPrintMBArrays(long stime) {
	FILE *fpmb;
	char fn[256];
	int i;
	double tot;
	struct _stat sb;
	if (!openMassBalanceData) return;

	sprintf(fn, "trace/mb/mb_%06d.txt", stime);
	if (!(_stat("trace/mb", &sb)==0 && sb.st_mode & _S_IFDIR)) {
		_mkdir("trace/mb");
	}
	fpmb=fopen(fn, "w");
	fprintf(fpmb, "Entity ID\tAdded\tRemoved\tTanks\tPipes\tLoss-NotMoved\tLoss-NoOutflow\tLoss-NegTanvkVol\n");
	fprintf(fpmb, "Total");
	double gtot = 0;
	tot=0; for (i=1; i<=Nnodes; i++) { tot += MB_MassAdded[i]; } fprintf(fpmb, "\t%0.14g", tot); gtot+=tot;
	tot=0; for (i=1; i<=Nnodes; i++) { tot += MB_MassRemoved[i]; } fprintf(fpmb, "\t%0.14g", tot); gtot+=tot;
	tot=0; for (i=1; i<=Ntanks; i++) { tot += MB_MassInTanks[i]; } fprintf(fpmb, "\t%0.14g", tot); gtot+=tot;
	tot=0; for (i=1; i<=Nlinks; i++) { tot += MB_MassInPipes[i]; } fprintf(fpmb, "\t%0.14g", tot); gtot+=tot;
	tot=0; for (i=1; i<=Nnodes; i++) { tot += MB_MassAdded[i]-MB_MassMoved[i]; } fprintf(fpmb, "\t%0.14g", tot); gtot+=tot;
	tot=0; for (i=1; i<=Nnodes; i++) { tot += MB_MassNoOutflow[i]; } fprintf(fpmb, "\t%0.14g", tot); gtot+=tot;
	tot=0; for (i=1; i<=Ntanks; i++) { tot += MB_MassNegTankVol[i]; } fprintf(fpmb, "\t%0.14g", tot); gtot+=tot;
	fprintf(fpmb, "\t%014g\n\n",gtot);

	for (i=1; i<=Nnodes; i++) {
		if (MB_MassAdded[i] != 0) { fprintf(fpmb, "%s\t%0.14g\n", Node[i].ID, MB_MassAdded[i]); }
	}
	for (i=1; i<=Nnodes; i++) {
		if (MB_MassRemoved[i] != 0) { fprintf(fpmb, "%s\t\t%0.14g\n", Node[i].ID, MB_MassRemoved[i]); }
	}
	for (i=1; i<=Ntanks; i++) {
		if (MB_MassInTanks[i] != 0) { fprintf(fpmb, "%s\t\t\t%0.14g\n", Node[Tank[i].Node].ID, MB_MassInTanks[i]); }
	}
	for (i=1; i<=Nlinks; i++) {
		if (MB_MassInPipes[i] != 0) { fprintf(fpmb, "%s\t\t\t\t%0.14g\n", Link[i].ID, MB_MassInPipes[i]); }
	}
	for (i=1; i<=Nnodes; i++) {
		if (MB_MassAdded[i]-MB_MassMoved[i] != 0) { fprintf(fpmb, "%s\t\t\t\t\t%0.14g\n", Node[i].ID, MB_MassAdded[i]-MB_MassMoved[i]); }
	}
	for (i=1; i<=Nnodes; i++) {
		if (MB_MassNoOutflow[i] != 0) { fprintf(fpmb, "%s\t\t\t\t\t\t%0.14g\n", Node[i].ID, MB_MassNoOutflow[i]); }
	}
	for (i=1; i<=Ntanks; i++) {
		if (MB_MassNegTankVol[i] != 0) { fprintf(fpmb, "%s\t\t\t\t\t\t\t%0.14g\n", Node[Tank[i].Node].ID, MB_MassNegTankVol[i]); }
	}
	fclose(fpmb);
}
void ftPrintMassBalanceData(long stime) {

	FILE *fpmb;
	char fn[256];
	int i;
	if (!openMassBalanceData) return;

	sprintf(fn, "trace/mass/mass_%06d.txt", stime);
	_mkdir("trace/mass");
	fpmb=fopen(fn, "w");
	fprintf(fpmb, "Mass Added\nID\tmass\n");
	for (i=1; i<=Nnodes; i++) {
		if (TransportNodes[i].massadded >0) {
			fprintf(fpmb, "%s\t%f\n", TransportNodes[i].node->ID, TransportNodes[i].massadded);
		}
	}
	fprintf(fpmb, "Mass Removed\nID\tmass\n");
	for (i=1; i<=Nnodes; i++) {
		double mass=0;
		Pseg s;
		for (s=TransportNodes[i].demandSegs.firstSeg; s!=NULL; s=s->prev) {
			mass += s->v * s->c;
		}
		if (mass!=0) {
			fprintf(fpmb, "%s\t%f\n", TransportNodes[i].node->ID, mass);
		}
	}
	fprintf(fpmb, "Mass In Pipes\nID\tmass\n");
	for (i=1; i<=Nlinks; i++) {
		double mass=0;
		Pseg s;
		for (s=FirstSeg[i]; s!=NULL; s=s->prev) {
			mass += s->v * s->c;
		}
		if (mass>0) {
			fprintf(fpmb, "%s\t%f\n", Link[i].ID, mass);
		}
	}
	fprintf(fpmb, "Mass In Tanks\nID\tmass\n");
	for (i=1; i<=Ntanks; i++) {
		if (Tank[i].M>0) {
			fprintf(fpmb, "%s\t%f\n", Node[Tank[i].Node].ID, Tank[i].M);
		}
	}
	fclose(fpmb);
}

void printMassBalance(FILE *fp, long step) {
	double rem, add, tanks, pipes;
	int i;
	rem=0;
	add=0;
	for (i=1; i<=Nnodes; i++) {
		rem+=MB_CumMassRemoved[i]+MB_MassRemoved[i];
		add+=MB_CumMassAdded[i]+MB_MassAdded[i];
	}
	tanks=0; for (i=1; i<=Ntanks; i++) { tanks+=MB_MassInTanks[i]; }
	pipes=0; for (i=1; i<=Nlinks; i++) { pipes+=MB_MassInPipes[i]; }
	fprintf(fp, "%d\t%f\t%f\t%f\t%f\t\t%f\t%f\n", step, add, rem, tanks, pipes, rem+tanks+pipes, add-(rem+tanks+pipes));
	fflush(fp);
}
void printmassbalanceDEBUG(long t, long dt) {
	int i;
	double c;
	if (fpMassBalance ==NULL) return;

	fprintf(demFP, "%d", t);
	fprintf(massFP, "%d", t);
	fprintf(qualFP, "%d", t);

	for (i=1; i<=Nnodes; i++) {
		double mass=0;
		PNodeTransport n = &TransportNodes[i];
		if (D[i]>0) {

			if (n->demandSegs.lastSeg != NULL) {
				c=ftGetAvgConc(&n->demandSegs);
			} else {
				c=0;
			}
		} else {
			if (n->throughSegs.lastSeg != NULL) {
				c=ftGetAvgConc(&n->throughSegs);
			} else {
				c=0;
			}
		}
		if (i <= Njuncs) {
			mass=D[i]*dt*c;
		} else {
			mass=Tank[i-Njuncs].M;
		}
		fprintf(demFP, "\t%.12f", D[i]);
		fprintf(massFP, "\t%f", mass);
		fprintf(qualFP, "\t%f", c);
	}
	fprintf(massFP, "\t");
	for (i=1; i<=Nlinks; i++) {
		double mass=0;
		mass = ftGetMass(FirstSeg[i]);
		fprintf(massFP, "\t%f", mass);
	}
	fprintf(demFP, "\n");
	fprintf(massFP, "\n");
	fprintf(qualFP, "\n");
}

int ftGetNumSegs(Pseg s) {
	int i=0;
	while (s!= NULL) {
		i++;
		s=s->prev;
	}
	return i;
}
void ftWriteSegs(FILE *fp, Pseg *segs, char **hdrs, int num, char *nodeID) {
	int maxsegs=0;
	int i;
	fprintf(fp, "%s", nodeID);
	for (i=0; i<num; i++) {
		maxsegs=MAX(maxsegs, ftGetNumSegs(segs[i]));
		fprintf(fp, "\t%s\t", hdrs[i]);
	}
	fprintf(fp, "\n");
	for (i=0; i<maxsegs; i++) {
		int sn;
		for (sn=0; sn<num; sn++) {
			Pseg s=segs[sn];
			if (s!=NULL) {
				fprintf(fp, "\t%f\t%f", s->v, s->c);
			}
		}
		fprintf(fp,"\n");
		for (sn=0; sn<num; sn++) {
			Pseg s=segs[sn];
			if (s!=NULL) {
				if (IsIncompleteSeg(s)) {
					PIncompleteSeg is=(PIncompleteSeg)s;
					PNodeTransport tnode=NULL;
					PIncompleteSeg ts=is;
					int j;
					while (ts->parent != NULL) {
						ts=ts->parent;
					}
					tnode=ts->node;
					fprintf(fp, "\t\tInflowLink\tv\tflowed\tc\n");
					for (j=0; j<is->nparts; j++) {
						PIncSegPart part=&is->parts[j];
						fprintf(fp, "\t\t%s\t%f\t%f\t%f\n", Link[tnode->inflow.linkIdxs[j]].ID, part->volume, part->flowed, part->conc);
					}
				}
				s=s->prev;
				segs[sn]=s;
			} else {
				//				fprintf(fp,"\t\t\n");
			}
		}
	}
	fprintf(fp, "\n");
}
void ftPrintDetail(long t, char *id, PSegList segsToMove, int sequence, int loc) {
	int i;

	char fn[256];
	FILE *fp;
	int num=0;
	Pseg *segs=NULL;
	char **headers=NULL;

	if (!openDetail)
		return;

	sprintf(fn, "trace/detail/detail_%06d_%08d_%s_%02d.txt", t, sequence, id, loc);
	_mkdir("trace");
	_mkdir("trace/detail");
	fp=fopen(fn, "w");
	for (i=1; i<=Nnodes; i++) {
		int j;
		PNodeTransport n = &TransportNodes[i];
		if (num < 2+n->inflowSegs.nlinks) {
			int l;
			segs = (Pseg *)realloc(segs,(2+n->inflowSegs.nlinks)*sizeof(Pseg));
			headers = (char **)realloc(headers,(2+n->inflowSegs.nlinks)*sizeof(char*));
			for (l=num; l<2+n->inflowSegs.nlinks; l++) {
				headers[l]=(char *)calloc(80, sizeof(char));
			}
			if (num==0) {
				strcpy(headers[0], "Through");
				strcpy(headers[1], "Demand");
			}
			num=2+n->inflowSegs.nlinks;
		}
		segs[0]=n->throughSegs.firstSeg;
		segs[1]=n->demandSegs.firstSeg;
		for (j=0; j<n->inflowSegs.nlinks; j++) {
			int lidx;
			segs[2+j]=n->inflowSegs.segs[j].firstSeg;
			lidx=n->inflow.linkIdxs[j];
			if (lidx >=0) {
				sprintf(headers[2+j], "Inflow_%s", Link[lidx]);
			} else {
				sprintf(headers[2+j], "External inflow");
			}
		}
		ftWriteSegs(fp, segs, headers, 2+n->inflowSegs.nlinks, Node[n->nodeIdx].ID);
	}

	strcpy(headers[0], "LinkSegs");
	for (i=1; i<=Nlinks; i++) {
		segs[0]=FirstSeg[i];
		if (segs[0]->c == -1) {
			fprintf(stdout, "");  fflush(stdout);
		}
		ftWriteSegs(fp, segs, headers, 1, Link[i].ID);
	}

	if (segsToMove != NULL) {
		strcpy(headers[0], "SegsToMove");
		segs[0]=segsToMove->firstSeg;
		ftWriteSegs(fp, segs, headers, 1, "");
	}



	fclose(fp);
	free(segs);
	for (i=0; i<num; i++) {
		free(headers[i]);
	}
	free(headers);
}
void ftPrintAllConcentrations() {
	FILE *fpw;
	fpAllAvgConc=openAllConcFile("AllAvgConc.txt");
	fpAllInstConc=openAllConcFile("AllInstConc.txt");

	fpw=fopen("towrite.txt", "r");
	if (fpw != NULL) {
		char buf[256];
		int nodes=0;
		int links=0;
		nodesToWrite = (int *)calloc(Nnodes+1, sizeof(int));
		linksToWrite = (int *)calloc(Nlinks+1, sizeof(int));
		while (!feof(fpw)) {
			fgets(buf, 255, fpw);
			if (strlen(buf)>0) {
				buf[strlen(buf)-1]=0;
				if (strncmp(buf, "Nodes", 5)==0) {
					links=0; nodes=1;
				} else if (strncmp(buf, "Links", 5)==0) {
					links=1; nodes=0;
				} else {
					int idx;
					if (nodes) {
						ENgetnodeindex(buf, &idx);
						nodesToWrite[idx]=1;
					} else {
						ENgetlinkindex(buf, &idx);
						linksToWrite[idx]=1;
					}
				}
			}
		}
		fclose(fpw);
	}
}


#endif
void TMPwriteConcFile(long stime) {
	FILE *fp;
	char fn[128];
	int i;
	return;
	sprintf(fn, "conc_%s", transportMethod==TM_ORIGINAL?"orig":"flow");
	_mkdir(fn);
	sprintf(fn, "conc_%s/conc_%05d.txt", transportMethod==TM_ORIGINAL?"orig":"flow", stime);
	fp=fopen(fn, "w");
	fprintf(fp, "%d\n", stime);
	for (i=1; i<=Nnodes; i++) {
		fprintf(fp, "%s\t%.8f", Node[i].ID, NodeQual[i]);
		if (transportMethod==TM_FLOW) {
			fprintf(fp, "\t%.8f", AvgNodeQual[i]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

#ifdef ENABLE_DBG_PRINT
char *getSegmentType(Pseg s) {
	if (IsNormalSeg(s)) return "normal";
	if (IsIncompleteSeg(s)) return "incomplete";
	if (IsMergedSeg(s)) return "merged";
	return "UNKNOWN";
}

void ftPrintAllIncSegs(long stime) {
	FILE *fp=NULL;
	int i;
	if (!openAllIncSegs)
		return;
	if (fp==NULL) {
		char fn[256];
		_mkdir("trace");
		sprintf(fn, "trace/all_inc_segs/all_inc_segs_%06d.txt", stime);
		_mkdir("trace/all_inc_segs");
		fp=fopen(fn, "w");
	}
	for (i=1; i<=Nnodes; i++) {
		PNodeTransport n=&TransportNodes[i];
		if (n->incompleteSegs != NULL) {
			PIncompleteSegmentData p=n->incompleteSegs;
			while (p != NULL) {
				ftPrintIncompleteSegments(p->incomplete, n->inflow.nlinks, fp);
				p=p->prev;
			}
		}
	}
	fclose(fp);
}

void ftCheckForIncompleteSegments(PNodeTransport n, long stime) {
	PIncompleteSegmentData prevInc = NULL;
	PIncompleteSegmentData psd;
	DBG_PRINT(ftPrintAllIncompleteSegments(n));
	if (strcmp(n->node->ID, "309102") == 0) {
		printf("");
	}
	for (psd = n->incompleteSegs; psd != NULL; psd = psd->prev) {
		int l;
		int inc = 0, merged = 0, vol = 0, other = 0;
		PIncompleteSeg *incsegs = n->incompleteSegs->incomplete;
		for (l = 0; l<n->inflow.nlinks; l++) {
			PIncompleteSeg is = incsegs[l];
			if (is != NULL) {
				int p = l;
				PIncSegPart part = &is->parts[p];
				if (part->conc == -2) {
					merged++;
				}
				else if (part->conc == -1) {
					inc++;
				}
				else if (part->flowed < part->volume) {
					vol++;
				}
				else {
					other++;
				}
			}
		}
		if (inc + merged + vol + other) {
			fprintf(stdout, "Node:\t%s\t%d\t%d\t%d\t%d\n", n->node->ID, inc, merged, vol, other);
		}
	}
}

void ftPrintSource(FILE *fp, long long segSrc) {
	int s, i;
	int *sources = NULL;
	int maxSource = 0;

	while (segSrc > 0) {
		s = segSrc % 100;
		if (s>maxSource) {
			sources = (int *)realloc(sources, s*sizeof(int));
			for (i = maxSource; i<s; i++) {
				sources[i] = 0;
			}
			maxSource = s;
		}
		sources[s - 1] = 1;
		segSrc = segSrc / 100;
	}
	fprintf(fp, "\t%d", maxSource);
	for (i = 0; i<maxSource; i++) {
		if (sources[i] == 1) {
			fprintf(fp, "\t1");
		}
		else {
			fprintf(fp, "\t");
		}
	}
	free(sources);
}
int ftHasSource(long long segSrc, int src) {
	int s;

	while (segSrc > 0) {
		s = segSrc % 100;
		if (s == src)
			return 1;
		segSrc = segSrc / 100;
	}
	return 0;
}
char *ftGetOwnerObject(PSegOwner o) {
	switch (o->type) {
	case SO_NONE:
		return "N/A";
	case SO_THROUGH:
		return o->object.node->node->ID;
		break;
	case SO_INFLOW:
		return o->object.node->node->ID;
		break;
	case SO_DEMAND:
		return o->object.node->node->ID;
		break;
	case SO_LINK:
		return Link[o->object.linkIdx].ID;
		break;
	case SO_MERGED:
		//o->object.mergedSeg;
		return "";
		break;
	case SO_INCOMPLETE_PART:
		//o->object.incSeg;
		return "";
		break;
	}
	return "Unk";
}

char *ftGetOwnerData(PSegOwner o) {
	char tmp[40];
	switch (o->type) {
	case SO_NONE:
		return "";
	case SO_THROUGH:
		return "";
		break;
	case SO_INFLOW:
		sprintf(tmp, "%d", o->data.inflowIdx);
		return tmp;
		break;
	case SO_DEMAND:
		return "";
		break;
	case SO_LINK:
		return "";
		break;
	case SO_MERGED:
		return "";
		break;
	case SO_INCOMPLETE_PART:
		//			o->data.incSeg;
		return "";
		break;
	}
	return "Unk";
}
void ftDocPrintSegListHeader() {
	fprintf(fpDocNodes, "Volume\tConc\tSegmentID\n");
}
void ftDocPrintSegList(PSegList seglist,char *hdr,char*hdr2) {
	int first = 1;
	for (Pseg s = seglist->firstSeg; s != NULL; s = s->prev) {
		fprintf(fpDocNodes, "%s\t%.3f\t%.3f\t0x%08x",first==1?hdr:hdr2, s->v, s->c, s);
		if (s->c == -1) {
			PIncompleteSeg i = (PIncompleteSeg)s;
		} else {

		}
		if (s->c == -2) {
			PMergedSeg m = (PMergedSeg)s;
			fprintf(fpDocNodes, "\t\t%.3f\t%.3f\t%s\n", m->knownVol, m->knownMass, OwnerTypes[m->owner.type]);
			for (PIncompleteSeg is = m->incSegs; is != NULL; is = is->prev) {
				fprintf(fpDocNodes, "%s\t\t\t\t\t\t\t\t%.3f\t%.3f\t0x%08x", first == 1 ? hdr : hdr2, is->v, is->c, is);
				if (is->prev != NULL) {
					fprintf(fpDocNodes, "\n");
				}
			}
		}
		fprintf(fpDocNodes, "\n");
		first = 0;
	}

}
void ftDocPrintLinkVolumes(long stime, SNodeTransport *node, char *desc) {
	if (fpDocNodes != NULL && dbgPrintNodes[node->nodeIdx]==1) {
		fprintf(fpDocNodes, "Link Segments\n***%s\t%ld\n", desc, stime);
		fprintf(fpDocNodes, "\tLink ID\t");
		ftDocPrintSegListHeader();
		for (int l = 0; l <= Nlinks; l++) {
			if (dbgPrintLinks[l] == 1) {
				SegList sl;
				fprintf(fpDocNodes, "\t%s", Link[l].ID);
				sl.firstSeg = FirstSeg[l];
				ftDocPrintSegList(&sl, "", "\t");
			}
		}
	}
}
void ftDocPrintNodeInfo(long stime, SNodeTransport *node, char *desc) {
	if (fpDocNodes != NULL) {
		if (dbgPrintNodes[node->nodeIdx] == 1) {
			fprintf(fpDocNodes, "\n*****************************\nNode\t%s\n***%s\t%ld\n", node->node->ID,desc,stime);
			fprintf(fpDocNodes, "\tDemand\t%.3f\n", node->demand);
			fprintf(fpDocNodes, "\tTotalDemand\t%.3f\n", node->totalDemand);
			fprintf(fpDocNodes, "\tInflow\tLinkID\tLinkVolume\tFlowRate\tFlowVolume\tFlowRatio\tVolumeFlowed\n");
			for (int i = 0; i < node->inflow.nlinks; i++) {
				int linkIdx = node->inflow.linkIdxs[i];
				fprintf(fpDocNodes, "\t\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", linkIdx<0?"External":Link[linkIdx].ID,
					linkIdx<0 ? 0 : LINKVOL(linkIdx),
					node->inflow.flowRate[i], node->inflow.flowVolume[i], node->inflow.flowRatio[i],
					node->inflow.volumeFlowed[i]);
			}
			fprintf(fpDocNodes, "\tOutflow Adjustment\t%f\n", node->outflowAdjustment);
			fprintf(fpDocNodes, "\tOutflow\tLinkID\tLinkVolume\tFlowRate\tFlowVolume\tFlowRatio\tVolumeFlowed\n");
			for (int i = 0; i < node->outflow.nlinks; i++) {
				int linkIdx = node->outflow.linkIdxs[i];
				fprintf(fpDocNodes, "\t\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", Link[linkIdx].ID, LINKVOL(linkIdx),
					node->outflow.flowRate[i], node->outflow.flowVolume[i], node->outflow.flowRatio[i],
					node->outflow.volumeFlowed[i]);
			}
			fprintf(fpDocNodes, "\tInflow Segments\tLinkID\t");
			ftDocPrintSegListHeader();

			for (int i = 0; i < node->inflowSegs.nlinks; i++) {
				int linkIdx = node->inflow.linkIdxs[i];
				fprintf(fpDocNodes, "\t\t%s", linkIdx<0 ? "External" : Link[linkIdx].ID);
				ftDocPrintSegList(&(node->inflowSegs.segs[i]), "", "\t\t");
			}
			fprintf(fpDocNodes, "\tMass Added\t%.3f\n", node->massadded);
			fprintf(fpDocNodes, "\tMass Ussed\t%.3f\n", node->massused);
			fprintf(fpDocNodes, "\tTotalVolInflow\t%.3f\n", node->totalVolInflow);
			fprintf(fpDocNodes, "\tTotalVolOutflow\t%.3f\n", node->totalVolOutflow);
		}
	}
}
void ftDocPrintNode(long stime, SNodeTransport *node, char *desc) {
	if (fpDocNodes != NULL) {
		if (dbgPrintNodes[node->nodeIdx] == 1) {
			fprintf(fpDocNodes, "***%s\n", desc);
			fprintf(fpDocNodes, "\tDemandSegs\t");
			ftDocPrintSegListHeader();
			ftDocPrintSegList(&(node->demandSegs), "\t","\t");
			fprintf(fpDocNodes, "\tThroughSegs\t");
			ftDocPrintSegListHeader();
			ftDocPrintSegList(&(node->throughSegs), "\t","\t");

/*
	SegList allDemandSegs;   // the segments that left this node due to demand for the current reporting interval
	SegList allThroughSegs;  // the segments that passed through this node for the current reporing interval
	PIncompleteSegmentData incompleteSegs;
*/
		}
	}
}
void ftDocPrintSegments(long stime, SNodeTransport *node, PSegList segs, char *desc) {
	if (fpDocNodes != NULL) {
		if (dbgPrintNodes[node->nodeIdx] == 1) {
			fprintf(fpDocNodes, "***%s\n", desc);
			fprintf(fpDocNodes, "\t");
			ftDocPrintSegListHeader();
			ftDocPrintSegList(segs, "","");
		}
	}

}
void ftDocPrintLinks(long stime, SNodeTransport *node, char *desc) {
	if (fpDocNodes != NULL) {
		if (dbgPrintNodes[node->nodeIdx] == 1) {
			fprintf(fpDocNodes, "***%s\tLink Segments\n", desc);
			fprintf(fpDocNodes, "\tLink ID\t", desc);
			ftDocPrintSegListHeader();
			for (int i = 0; i < node->outflow.nlinks; i++) {
				int linkIdx = node->outflow.linkIdxs[i];
				SegList sl;
				fprintf(fpDocNodes, "\t%s", Link[linkIdx].ID);
				sl.firstSeg = FirstSeg[linkIdx];
				ftDocPrintSegList(&sl, "", "\t");
			}
		}
	}
}
void ftDocPrintConcentrations(long stime, char *desc) {
	if (fpDocNodes != NULL) {
		fprintf(fpDocNodes, "Node Concentrations\n***%s\t%ld\n", desc, stime);
		fprintf(fpDocNodes, "\tNode ID\tC\tAvgC\n");
		for (int n = 0; n <= Nnodes; n++) {
			PNodeTransport node = &TransportNodes[n];
			if (dbgPrintNodes[n] == 1) {
				double c, avgc;
				if (D[n]>0) {
					if (node->demandSegs.lastSeg != NULL) {
						double mass = ftGetMass(node->demandSegs.firstSeg);
						c = node->demandSegs.lastSeg->c;
						//  not using ftGetAvgConc here because the summed volume
						// may not be exsctly the same as the computed demand due to tolerance
						// issues, so simply sum the mass and divide it by the actaul demand.
						avgc = mass/ node->demand;
					}
					else if (n <= Njuncs) {
						c = 0;
						avgc = 0;
					}
				}
				else {
					// if the node does not have demand, set the instantaneous concentration
					// to the the most recently added through segment's concentration
					// and compute the average concetration from all the through segments
					if (node->throughSegs.lastSeg != NULL) {
						c = node->throughSegs.lastSeg->c;
						avgc = ftGetAvgConc(&node->throughSegs);
					}
					else {
						c = 0;
						avgc = 0;
					}
				}
				fprintf(fpDocNodes, "\t%s\t%.3f\t%.3f\n", Node[n].ID, c, avgc);
			}
		}
	}
}

void ftPrintInflowOutflowData(long stime) {
	if (fpInOut == NULL) return;
	fprintf(fpInOut, "Node\tadj\t\tExpInflow\tExpOutflow\texpDemand\t\tactInflow\tactOutflow\tactDemand\tactThrough\n");
	for (int i = 1; i <= Nnodes; i++) {
		PNodeTransport n = &TransportNodes[i];
		double expInflow = 0, expOutflow = 0,actInflow=0,actOutflow=0,actDemand=0,actThrough=0;
		for (int j = 0; j < n->inflow.nlinks; j++) {
			expInflow += n->inflow.flowVolume[j];
			actInflow += n->inflow.volumeFlowed[j];
		}
		for (int j = 0; j < n->outflow.nlinks; j++) {
			expOutflow += n->outflow.flowVolume[j];
			actOutflow += n->outflow.volumeFlowed[j];
		}
		for (Pseg s = n->demandSegs.firstSeg; s != NULL; s = s->prev) {
			actDemand += s->v;
		}
		for (Pseg s = n->throughSegs.firstSeg; s != NULL; s = s->prev) {
			actThrough += s->v;
		}
		fprintf(fpInOut, "%s\t%lf\t\t%lf\t%lf\t%lf\t\t%lf\t%lf\t%lf\t%lf\n", n->node->ID, n->outflowAdjustment,expInflow, expOutflow, n->demand, actInflow, actOutflow, actDemand, actThrough);
	}
}

#endif

