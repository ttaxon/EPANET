
#include "ft_checks.h"

int ftCheckForAlternatingSegs(PSegList segs);
/*
**--------------------------------------------------------------
**   Input:   segs: the ist of segments to check for alternating concentrations
**   Output:  the number of "flips"
**   Purpose: There have been instances where some scenarios take a significant
**     amount of time relative to others.  In some cases, this has been traced
**     back to many segments alternating between a zero concentration and a
**     very small non-zero concentration.  When this happens, the volumes of
**     those segmnents also typically are very small also.  These methods look
**     for that case and will write information to the file alt_segs.txt in the
**     current working directory
**--------------------------------------------------------------
*/

void ftAlternatingSegmentCheck(long stime) {
	static FILE *fp=NULL;
	int i;
	if(fp==NULL) {
		fp=fopen("alt_segs.txt", "w");
	}
	for (i=1; i<=Nnodes; i++) {
		int cntD, cntT;
		SNodeTransport *n=&TransportNodes[i];
		cntD=ftCheckForAlternatingSegs(&n->allDemandSegs);
		cntT=ftCheckForAlternatingSegs(&n->allThroughSegs);
		if (ftCheckForAlternatingSegs(&n->allDemandSegs) > 1 || ftCheckForAlternatingSegs(&n->allThroughSegs) > 1) {
			fprintf(fp, "%ld\t%s\t%d\t%d\n", stime, n->node->ID, cntD, cntT);  fflush(fp);
			if (cntD > 4 || cntT > 4) {
				printf("");
			}
			if (cntD) {
				Pseg s;
				fprintf(fp, "Demand Segs\n");
				for (s=n->allDemandSegs.firstSeg; s!=NULL; s=s->prev) {
					fprintf(fp, "\t%g\t%g\n", s->v, s->c);
				}
				fflush(fp);
			}
			if (cntT) {
				Pseg s;
				fprintf(fp, "Through Segs\n");
				for (s=n->allThroughSegs.firstSeg; s!=NULL; s=s->prev) {
					fprintf(fp, "\t%g\t%g\n", s->v, s->c);
				}
				fflush(fp);
			}
		}
	}
	for (i=1; i<=Nlinks; i++) {
		SegList sl;
		int cnt;
		sl.firstSeg=FirstSeg[i];
		sl.lastSeg=LastSeg[i];
		if ((cnt=ftCheckForAlternatingSegs(&sl)) > 1) {
			Pseg s;
			fprintf(fp, "%ld\t%s\t%d\n", stime, Link[i].ID, cnt);  fflush(fp);
			fprintf(fp, "Link Segs\n");
			for (s=FirstSeg[i]; s!=NULL; s=s->prev) {
				fprintf(fp, "\t%g\t%g\n", s->v, s->c);
			}
			fflush(fp);
		}

	}
}
/*
**--------------------------------------------------------------
**   Input:   segs: the ist of segments to check for alternating concentrations
**   Output:  the number of "flips"
**   Purpose: There have been instances where some scenarios take a significant
**     amount of time relative to others.  In some cases, this has been traced
**     back to many segments alternating between a zero concentration and a
**     very small non-zero concentration.  When this happens, the volumes of
**     those segmnents also typically are very small also.  These methods look
**     for that case and will write information to the file alt_segs.txt in the
**     current working directory
**--------------------------------------------------------------
*/
int ftCheckForAlternatingSegs(PSegList segs) {
	int cnt=0;
	double lastC=-1;
	Pseg seg;
	for (seg=segs->firstSeg; seg != NULL; seg=seg->prev) {
		if (lastC == -1) {
			lastC=seg->c;
		} else {
			if (seg->c==0) {
				if (lastC < 1e-6) {
					cnt++;
				}
				lastC=0;
			} else {
				if (lastC == 0) {
					if (seg->c < 1e-6 && seg->prev==NULL) {
						cnt++;
					}
					lastC=seg->c;
				}
			}
		}
	}
	return cnt;
}
