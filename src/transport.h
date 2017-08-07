#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_
#ifndef CODEVERSION // checks to see if types.h is already included
#include "types.h"
#endif
#include <stdio.h>
#include "hash.h"
#include "types.h"
#include "funcs.h"
#define  EXTERN  extern
#include "vars.h"
#include "mempool.h"

//#include "toolkit.h"
#include "ft_debug.h"

// These are used from quality.c
EXTERN Pseg FreeSeg;              /* Pointer to unused segment               */
EXTERN Pseg *FirstSeg;            /* First (downstream) segment in each pipe */
EXTERN Pseg *LastSeg;             /* Last (upstream) segment in each pipe    */
EXTERN char *FlowDir;             /* Flow direction for each pipe            */


typedef struct {
	int nlinks;           // the number of links
	int size;             // the size of the arrays.  This is so that these arrays do not have to be constantly freed/alloced. see ftAddLink for more details
	int *linkIdxs;        // the link indices (into the Link array).  Negative values (inflow only)identify where water is entering the system
	double *flowRatio;    // the inflow or outflow ratio of each links
	double *flowRate;     // the inflow or outflow rate of each link
	double *flowVolume;   // the inflow or outflow volume of each link
	double *volumeFlowed; // the actual volume flowed so far.
} SNodeLinks, *PNodeLinks;

typedef struct SegList {
	Pseg firstSeg;      // the oldest segment added - the first one to go.  Then follow prev pointers
	Pseg lastSeg;       // The last segment added
} SegList, *PSegList;

typedef struct SNodeSegs {
	int nlinks;         // number of inflow links
	int size;           // the size of the seg array.  This is so that it does not have to be constantly freed/alloced. see ftAddLink for more details
	SegList *segs;      // segments coming from upstream nodes.  array size: nlinks;
} SNodeSegs, *PNodeSegs;

typedef enum { SO_NONE,SO_THROUGH,SO_INFLOW,SO_DEMAND,SO_LINK,SO_MERGED,SO_INCOMPLETE_PART } SegOwnerType;

typedef union SegOwnerObject {
	// SO_NONE does not require an object
	struct SNodeTransport *node;     // SO_THROUGH, SO_INFLOW, SO_DEMAND
	int linkIdx;	         // SO_LINK
	struct SMergedSeg *mergedSeg;    // SO_MERGED
	struct SIncSegPart *incSeg; // SO_INOMPLETE_PART
} SegOwnerObject, *PSegOwnerObject;

typedef struct VolSplitData {
	double volume;
	double adjustment;
} VolSplitData, *PVolSplitData;

typedef union SegOwnerData {
	// SO_NONE, SO_THROUGH,SO_DEMAND,SO_LINK,SO_MERGED do not require extra data
	struct SIncompleteSeg *incSeg;    // SO_INCOMPLETE_PART
	int inflowIdx;     // SO_INFLOW
} SegOwnerData, *PSegOwnerData;

typedef struct SegOwner {
	SegOwnerType type;     // the type of object that "owns" this segment
	SegOwnerObject object; // the "object" to which this partial segment belongs.  Data type is dependent on type field
	SegOwnerData data;     // other data 
	long long source; // temp for debugging so it is known where it was created...
	long sequence;
} SegOwner, *PSegOwner;

typedef struct SIncSegPart {
	double ratio;       // ratio of this component's voluime to the total segment volume
	double volume;      // the volume of this component
	double flowed;      // how much volume has flowed
	double conc;        // the concentration of this component
	Pseg   segs;        // incomplete or merged segments that make up this part.
} SIncSegPart, *PIncSegPart;

typedef struct SIncompleteSeg { // incomplete segment
	/* The first three fields: c, v, and prev MUST be in the same order as Pseg
	    as this is used in place(temporarily) of Pseg in certain circumstances */
	double  v;             /* Segment volume      */
	double  c;             /* Water quality value - value will always be -1.0*/
	struct  SIncompleteSeg *prev; /* Record for previous segment */

	struct SIncompleteSeg *parent;// the parent SincSeg (if any)
	struct SNodeTransport *node;   // the node that owns this incomplete segment.  Only non-null at top-most level (when parent is NULL)
	double ratio;          // the ratio of the parent this segment represents
//	double valveAdjRatio;  // the volume adjustment ratio (due to valves that change the flow rate)
	double adjustment;     // the flow rate adjustment factor
	double massadded;      // mass added by a source.  Must be incorporated into c when this segment is complete
	SegOwner owner;        // the "owner" of this incomplete segment
	int nparts;            // the number of parts this incomplete segment has - the same as the number of inflow links the containing node has
	PIncSegPart parts;// the component parts
	int nchildren;         // the numnber of children this node has
	struct SIncompleteSeg *children; // the children
	int processing;        // set if this incomplete segment is currently bing processed - to detect loop conditions
} SIncompleteSeg, *PIncompleteSeg;

typedef struct SMergedSeg {
	/* The first three fields: c, v, and prev MUST be in the same order as Pseg
	    as this is used in place(temporarily) of Pseg in certain circumstances */
	double  v;                /* Segment volume      */
	double  c;                /* Water quality value - value will always be -1.0*/
	struct  SMergedSeg *prev; /* Record for previous segment */
	double adjustment;
	double knownVol;          // the volume of water with known concentration
	double knownMass;         // the mass in the volume of water with known concentration
	PIncompleteSeg incSegs;   // the list of incomplete segments
	SegOwner owner;           // the "owner" of this merged segment;
} SMergedSeg, *PMergedSeg;

typedef struct SIncompleteSegmentData {
	PIncompleteSeg *incomplete;
	struct SIncompleteSegmentData *prev;
} SIncompleteSegmentData, *PIncompleteSegmentData;

typedef struct SNodeTransport {
	Snode *node;             // the Snode data for this transport node
	double outflowAdjustment;
	double demand;           // demand at the node for the time step
	double totalDemand;      // demand at the node for the WQ time step
	SNodeLinks inflow;       // the links flowing into this node for the current hydraulic state
	SNodeLinks outflow;      // the links flowing out of this node for the current hydraulic state
	SNodeSegs inflowSegs;    // segments entering this node either from external sources or upstream nodes.
	double massadded;        // mass added from a source
	double massused;         // mass actually used
	double totalVolInflow;   // the volume that flows into the node during the current dt
	double totalVolOutflow;  // the volume that flows out of the node during the current dt
	SegList demandSegs;      // the segments that left this node due to demand for the current timestep (dt)
	SegList throughSegs;     // the segments that passed through this node for the current timestep (dt)
	SegList allDemandSegs;   // the segments that left this node due to demand for the current reporting interval
	SegList allThroughSegs;  // the segments that passed through this node for the current reporing interval
	PIncompleteSegmentData incompleteSegs;
#ifdef ENABLE_DBG_PRINT
	int visited;             // how many times this node has been visited - only used for debug print info
#endif
	int nodeIdx;             // the index into the global Node array
	int toProcess;           // whether or not this node is waiting to be processed.
} SNodeTransport, *PNodeTransport;

typedef struct SNodeTransportList {
	SNodeTransport *node;
	int age;
	struct SNodeTransportList *next;
} SNodeTransportList, *PNodeTransportList;


void ftInitQual();
void ftOpenMoveSegs();
void ftCloseMoveSegs();
void ftInitSegs();
void ftNewHydStep();

void ftMoveSegs(long dt,long stime);
void ftInitializeTransport();
void ftTransportDone();

void ftTankDirectInput(int tankIdx, double massadded);


int ftWithinTolerance(double v1, double v2,double atol, double rtol);
int ftHaveAnyInflowSegs(SNodeSegs *segs);
double ftSumVolume(PSegList segs);
double ftComputeSegmentWaterVolume(Pseg seg);
double ftComputeCompleteSegmentWaterVolume(Pseg seg);
double ftGetAvgConc(PSegList segs);
double ftGetMass(Pseg segs);

void TMPwriteConcFile(long stime);

FILE *openAllConcFile(char *fn);
void appendConcFile(FILE *fp, double *c);
// macros to help identify the type of a segment
#define IsNormalSeg(s) ((s)->c>=0)
#define IsIncompleteSeg(s) ((s)->c==-1)
#define IsMergedSeg(s) ((s)->c==-2)

#define   LINKVOL(k)   ( 0.785398*Link[(k)].Len*SQR(Link[(k)].Diam) )
extern double FT_QZERO;
// set FT_RTOL to 0 because sometimes (especially with long WQ time steps),
// segments were being ignored when they shouldn't have been.  This seems
// to be working for now.
//double FT_RTOL=1e-4;
extern double FT_RTOL;
// FT_ATOL_BASE is the basis for computing FT_ATOL based on the dt that is being
// used for any given step
extern double FT_ATOL_BASE;
extern double FT_ATOL;
// these two are used to aid in determining if two segments can be combined
// for more detail, see the comment in ftCombineSegemnts
extern double FT_CTOL;
extern double FT_VTOL;

extern SNodeTransport *TransportNodes; /* Flow transport data structures */
extern char *OwnerTypes[];

#endif
