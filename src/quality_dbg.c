
#define QDBG_EXT 
#include "quality_dbg.h"
#include <direct.h>

extern double    *VolIn;               /* Total volume inflow to node             */
extern double    *MassIn;              /* Total mass inflow to node               */

FILE *fpMassInfo = NULL;
FILE *fpLV = NULL;
double cLostMass = 0;
double cExtraMass = 0;
double cExtraVol = 0;

void dbgPrintLinkQual(long stime, char *type) {
	int i;
	static FILE **fps = NULL;
	static Pseg *segs = NULL;
	static char *fpw = NULL;
	static int *lastNZ;
	static int maxSegs = 0;
	static long prevTime = 0;
	char fn[256];

	if (fps == NULL) {
		lastNZ = (int*)calloc(Nlinks + 1, sizeof(int));
		fps = (FILE **)calloc(Nlinks + 1, sizeof(FILE *));
		fpw = (char *)calloc(Nlinks + 1, sizeof(char));
		for (i = 1; i <= Nlinks; i++) {
			sprintf(fn, "trace/linkqual/linkqual_%s.txt", Link[i].ID);
			_unlink(fn);
		}
	}
	_mkdir("trace");
	sprintf(fn, "trace/linkqual");
	_mkdir(fn);
	for (i = 1; i <= Nlinks; i++) {
		if (LinkHasNZConcentration(FirstSeg[i])) {
			Pseg s;
			int nsegs;
			double m = 0, v = 0, avgC;
			int firstSeg, lastSeg, incr, n, stop;
			if (fps[i] == NULL) {
				char fn[256];
				sprintf(fn, "trace/linkqual/linkqual_%s.txt", Link[i].ID);
				fps[i] = fopen(fn, "w");
				if (!fpw[i]) {
					fprintf(fps[i], "time\tID\tN1\tN2\tAvgC\tC1\tC2\tv\tc\t...\n");
					fpw[i] = 1;
				}
			}
			nsegs = 0;
			for (s = FirstSeg[i]; s != NULL; s = s->prev) {
				if ((nsegs + 1) > maxSegs) {
					maxSegs = nsegs + 1;
					segs = (Pseg *)realloc(segs, maxSegs*sizeof(Pseg));
				}
				segs[nsegs] = s;
				m += s->v*s->c;
				v += s->v;
				nsegs++;
			}
			avgC = m / v;
			nsegs--;
			if (UP_NODE(i) == Link[i].N1) { // flow is N1->N2
				firstSeg = nsegs;
				lastSeg = 0;
				incr = -1;
				stop = -1;
			}
			else {
				firstSeg = 0;
				lastSeg = nsegs;
				incr = 1;
				stop = nsegs + 1;
			}
			if (lastNZ[i] == 0) {
				fprintf(fps[i], "%d\t%s\t%s\t%s\t%g\t%g\t%g\t%g\t%g\n", prevTime, Link[i].ID, Node[Link[i].N1].ID, Node[Link[i].N2].ID, 0, 0, 0, LINKVOL(i), 0);
			}
			fprintf(fps[i], "%d\t%s\t%s\t%s\t%g\t%g\t%g", stime, Link[i].ID, Node[Link[i].N1].ID, Node[Link[i].N2].ID, avgC*Ucf[QUALITY], segs[firstSeg]->c*Ucf[QUALITY], segs[lastSeg]->c*Ucf[QUALITY]);
			for (n = firstSeg; n != stop; n += incr) {
				fprintf(fps[i], "\t%g\t%g", segs[n]->v, segs[n]->c*Ucf[QUALITY]);
			}
			fprintf(fps[i], "\n");
			lastNZ[i] = 1;
		}
		else {
			if (lastNZ[i] == 1) {
				if (fps[i] == NULL) {
					char fn[256];
					sprintf(fn, "trace/linkqual/linkqual_%s.txt", Link[i].ID);
					fps[i] = fopen(fn, "w");
					if (!fpw[i]) {
						fprintf(fps[i], "time\tID\tN1\tN2\tAvgC\tC1\tC2\tv\tc\t...\n");
						fpw[i] = 1;
					}
				}
				fprintf(fps[i], "%d\t%s\t%s\t%s\t%g\t%g\t%g\t%g\t%g\n", stime, Link[i].ID, Node[Link[i].N1].ID, Node[Link[i].N2].ID, 0, 0, 0, LINKVOL(i), 0);
				lastNZ[i] = 0;
			}
		}
		if (fps[i] != NULL) {
			fclose(fps[i]);
			fps[i] = NULL;
		}
	}
	prevTime = stime;
}
void dbgPrintDetail(long stime, int loc, long dt, double *x)
{
	int i;
	char fn[256];
	FILE *fp;

	//	return;
	//	if(stime<3660 || stime >3780) return;
	if (stime>10800) return;
	sprintf(fn, "trace/detail/detail_%06d_%02d.txt", stime, loc);
	_mkdir("trace");
	_mkdir("trace/detail");
	fp = fopen(fn, "w");
	fprintf(fp, "NodeID\tMassIn(mg)\tVolIn (m^3)\tC (mg/l)\tSrcContrib (mg/l)\tD (m^3)\tMass(mg)\n");
	for (i = 1; i <= Nnodes; i++) {
		char *nodeIDs[] = { "239", "247", "249", NULL };
		char *nodeID = Node[i].ID;
		int p = 0;
		for (int k = 0; nodeIDs[k] != NULL; k++) {
			if (strcmp(nodeID, nodeIDs[k]) == 0) {
				p = 1;
			}
		}
		if (p) {
			fprintf(fp, "%s[%d]\t%g\t%g\t%g\t%g\t%g\t%g\n", Node[i].ID, i, MassIn[i], VolIn[i] * LperFT3 / 1000.0, NodeQual[i] * Ucf[QUALITY], (x != NULL ? x[i] * Ucf[QUALITY] : 0), NodeDemand[i] * LperFT3 * dt / 1000.0, NodeQual[i] * Ucf[QUALITY] * NodeDemand[i] * LperFT3 * dt);
		}
	}
	fprintf(fp, "\n");
	fprintf(fp, "Link ID\tLink Vol (m^3)\tFlow (m^3)\tVol (m^3)\tConc (mg/l)\tMass (mg)\n");
	for (i = 1; i <= Nlinks; i++) {
		char *linkIDs[] = { "273", "275", "277", "281", "283", "285", "287", "295", NULL };
		char *linkID = Link[i].ID;
		int p = 0;
		for (int k = 0; linkIDs[k] != NULL; k++) {
			if (strcmp(linkID, linkIDs[k]) == 0) {
				p = 1;
			}
		}
		if (p) {
			Pseg s = FirstSeg[i];
			fprintf(fp, "%s[%d] %s -> %s\t%g\t%g%s", Link[i].ID, i, Node[UP_NODE(i)].ID, Node[DOWN_NODE(i)].ID, LINKVOL(i) * LperFT3 / 1000, ABS(Q[i] * LperFT3 * dt / 1000.0), (s == NULL ? "\n" : ""));
			while (s != NULL) {
				if (s != FirstSeg[i]) {
					fprintf(fp, "\t\t");
				}
				fprintf(fp, "\t%g\t%g\t%g\n", s->v*LperFT3 / 1000.0, s->c*Ucf[QUALITY], s->v*LperFT3 * s->c*Ucf[QUALITY]);
				s = s->prev;
			}
		}
	}
	fclose(fp);

}
void writeLinkVolumes(FILE *fplv, long stime, char *id, long dt)
{
	FILE *fp;
	char fn[256];
	int i;
	int wroteAny = 0;

	if (fplv == NULL) {
		sprintf(fn, "linkvols/lv_%06d_%s.txt", stime, id);
		_mkdir("linkvols");
		fp = fopen(fn, "w");
		fprintf(fp, "time\ti\tID\tLinkVol\tFlowVol\tnumsegs\tsegvol\tmass\tlinkvol-segvol\n");
	}
	else {
		fp = fplv;
	}
	for (i = 1; i <= Nlinks; i++) {
		Pseg s = FirstSeg[i];
		double mass = 0;
		double segvol = 0;
		int numsegs = 0;
		for (s = FirstSeg[i]; s != NULL; s = s->prev) {
			numsegs++;
			mass += s->v*s->c;
			segvol += s->v;
		}
		if (fplv == NULL || mass>0) {
			fprintf(fp, "%d\t%d\t%s\t%f\t%f\t%d\t%f\t%f\t%f\n", stime, i, Link[i].ID, LINKVOL(i), ABS(Q[i])*dt, numsegs, segvol, mass, LINKVOL(i) - segvol);
			wroteAny = 1;
		}
	}
	if (fplv != NULL) {
		if (!wroteAny) {
			fprintf(fp, "%d\n", stime);
		}
	}
	else {
		fclose(fp);
	}

}
int LinkHasNZConcentration(Pseg seg) {
	while (seg != NULL) {
		if (seg->c > 0) return 1;
		seg = seg->prev;
	}
	return 0;
}
void dbgPrintNodeQual(long stime, char *type) {
	int i;
	static FILE **fps = NULL;
	static char *fpw = NULL;
	char fn[256];

	if (fps == NULL) {
		fps = (FILE **)calloc(Nnodes + 1, sizeof(FILE *));
		fpw = (char *)calloc(Nnodes + 1, sizeof(char));
		for (i = 1; i <= Nnodes; i++) {
			sprintf(fn, "trace/nodequal/nodequal_%s.txt", Node[i].ID);
			_unlink(fn);
		}
	}
	_mkdir("trace");
	sprintf(fn, "trace/nodequal");
	_mkdir(fn);
	for (i = 1; i <= Nnodes; i++) {
		if (NodeQual[i] > 0) {
			if (fps[i] == NULL) {
				char fn[256];
				sprintf(fn, "trace/nodequal/nodequal_%s.txt", Node[i].ID);
				fps[i] = fopen(fn, "a");
				if (!fpw[i]) {
					fprintf(fps[i], "timestep\tNode_ID\tC\tAvgC\n");
					fpw[i] = 1;
				}

			}
			fprintf(fps[i], "%d\t%s\t%g", stime, Node[i].ID, NodeQual[i] * Ucf[QUALITY]);
			if (AvgNodeQual != NULL) {
				fprintf(fps[i], "\t%g", AvgNodeQual[i] * Ucf[QUALITY]);
			}
			else {
				fprintf(fps[i], "\t%g", 0);
			}
			fprintf(fps[i], "\n");
			fclose(fps[i]);
			fps[i] = NULL;
		}
	}
}
void dbgPrintTankQual(long stime, char *type) {
	int i;
	static FILE **fps = NULL;
	char fn[256];

	if (fps == NULL) {
		fps = (FILE **)calloc(Ntanks + 1, sizeof(FILE *));
		for (i = 1; i <= Ntanks; i++) {
			sprintf(fn, "trace/tankqual/tankqual_%s.txt", Node[Tank[i].Node].ID);
			_unlink(fn);
		}
	}
	_mkdir("trace");
	sprintf(fn, "trace/tankqual");
	_mkdir(fn);
	for (i = 1; i <= Ntanks; i++) {
		if (Tank[i].C > 0) {
			if (fps[i] == NULL) {
				char fn[256];
				sprintf(fn, "trace/tankqual/tankqual_%s.txt", Node[Tank[i].Node].ID);
				fps[i] = fopen(fn, "a");
				fprintf(fps[i], "timestep\tTank_ID\tV\tC\t\n");
			}
			fprintf(fps[i], "%d\t%s\t%f\t%g\n", stime, Node[Tank[i].Node].ID, Tank[i].V, Tank[i].C*Ucf[QUALITY]);
			//			fflush(fps[i]);
		}
	}
}
