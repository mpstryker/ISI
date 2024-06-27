/* Image analysis               */
/* Written by V.Kalatsky        */
/* 02Feb01                      */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cpgplot.h>
#include <ctype.h>
//#include <time.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
//#include <fcntl.h>
//#include <tcl.h> 
#include <errno.h>

#include "fileheader.h"

extern int errno; 

typedef struct CAR {
  int n;
  int iLock;
  struct CAR *next;
  struct CAR *prev;
  //  struct FRAMEHEADER frameheader;
  void *pFrame;
  void *pFrameHeader;
  void *pFrameImage;
  unsigned short *pimage;
} CAR; 

typedef struct FILESTRUC {
  int iVersion;
  int iExperiment;
  char *fullfilename;
  char *filename;
  FILE *fp;
  int state_of_file;
  void *pFileHeader;
  unsigned long nFileHeaderSize;
  void *pFileFooter;
  unsigned long nFileFooterSize;
  unsigned int ini_frame;
  unsigned int fin_frame;
  unsigned int offset_frame;
} FILESTRUC;

typedef struct TRAIN {

  char *filename;
  FILE *fp;

  int iGlobalVersion;
  int iGlobalExperiment;
  int n_files;
  struct FILESTRUC *files;
  unsigned int *frame_file_n;
  unsigned int max_n;
  unsigned int Xdim;
  unsigned int Ydim;
  unsigned long XYdim;
  unsigned long ulDataType;
  char *frameread;
  struct CAR **framemap;
  unsigned int n;                /* current number of items in train */
  unsigned int size_fileheader;  /* in bytes */
  unsigned int nFrameHeaderSize; /* in bytes */
  unsigned int nFrameImageSize;  /* in bytes */
  unsigned int nFrameSize;       /* sum of the two above */
  unsigned long ulNSynchChannels;
  unsigned long *pulSynchChannelMax;
  struct CAR *head;
  struct CAR *tail;
  struct CAR *replace;
  struct CAR *newest;
} TRAIN;

#define VERSION_UNKNOWN 0
#define VERSION_HEADER  1
#define VERSION_CHUNK   2

#define EXPERIMENT_CONTINUOUS_MODE_NUMBER     3
#define EXPERIMENT_CONTINUOUS_MODE_NONE	      0
#define EXPERIMENT_CONTINUOUS_MODE_CONTINUOUS 1
#define EXPERIMENT_CONTINUOUS_MODE_EPISODIC   2

#define FILE_TYPE_UNKNOWN    0
#define FILE_TYPE_GREEN      1
#define FILE_TYPE_STREAM     2
#define FILE_TYPE_COMPRESSED 3
#define FILE_TYPE_ANALYSIS   4

#define MAX_FILE_N 1

#define FILE_OPEN 1
#define FILE_CLOSED 0

#define FRAME_OUT 0
#define FRAME_IN 1

#define SCHEME_DIFF2 1
#define SCHEME_DIFFN 0

#define DIFF_RELATIVE 0
#define DIFF_ABSOLUTE 1

#define DO_SINGLE 0
#define DO_MULTIPLE 1

#define REPLACE_OLDEST 0
#define REPLACE_REQUESTED 1
#define REPLACE_RANDOM 2

#define DOUBLE 0
#define FLOAT  1

#define AVERAGE 0
#define AVERAGE_NOT 1

#define MAX_FRAMES_IN_CUE 3000
#define MAX_MEMORY_USAGE 0.5
#define NFRAMES_TO_FREE 10

#define SYNCH_MAX 64.0

#define POLYNOM_COEFF 0.1
#define POLYNOM_ADD 2

#define DISASSEMBLE 0
#define ASSEMBLE 1

#define COMPRESS_NOT 0
#define COMPRESS 1
#define DECOMPRESS 2

#define DRADIUS 0.0
#define DRADIUSBIG 0.0
#define DCRADIUS 2.5
#define TIMERADIUS 100
#define CUT_OFF_L 1
#define CUT_OFF_U 1000
#define HARMONIC 1.0
#define DEFAULT_INITFRAME 1
#define NCYCLES 1

#define DEFAULT_PERIOD 1 // frames/cycle // Val: totally meaningless, just to get it running

#define DEFAULT_SYNCH_MAX 64

#define SAVE_XY_IN_TWO_FILES 0
#define SAVE_XY_IN_ONE_FILE 1
#define SAVE_XY_DEFAULT SAVE_XY_IN_ONE_FILE

#define SAVE_TYPE_BIT_DOUBLE_FLOAT 1
#define SAVE_TYPE_BIT_ADD_DOUBLE_AT_END 2
#define SAVE_TYPE_BIT_XY_IN_FILES 4

#define DEFAULT_STATISTICS_NUM_BINS 256

#define SYNCHRONIZATION_NONE          0x0L
#define SYNCHRONIZATION_TRUE_EXTERNAL 0x1L
#define SYNCHRONIZATION_QUASI_TIME    0x2L
#define SYNCHRONIZATION_QUASI_USER    0x4L
#define SYNCHRONIZATION_EPST_INTERNAL 0x8L

#define DEFAULT_KEEP_ORIGINAL_FILES 0

#define GREEN_IMAGE_SIZE 8.0
#define GREEN_IMAGE_SIZE_SAVE 6.0

#define TRACE_INTERACTIVE_WINDOW_SIZEX 7.0
#define TRACE_PLOT_WINDOW_SIZEX 10.0

void start_graphics();
void end_graphics();

int InitializeTrain(TRAIN *tr,char *filen);
void ReleaseTrain(TRAIN *tr,int iMode);
int InitializeBuffers(TRAIN *tr);
int InitializeFitting(TRAIN *tr);
int InitializeTimeAveraging(TRAIN *tr);

CAR *AddFrame(TRAIN *tr,int fr,int a);
int RemoveFrame(TRAIN *tr,int fr);
CAR *ShrinkCue(TRAIN *tr,int fr,int a);
CAR *ReplaceFrame(TRAIN *tr,int frold,int frnew,int mode,int a);
int LoadFrame(TRAIN *tr,int fr,float *fbuffer,int a);
int LockFrame(TRAIN *tr,CAR *pcar,int iLock);
int GetRecordUS(TRAIN *tr,int fr,int iByteOffset,unsigned short *pus);
int GetRecordUL(TRAIN *tr,int fr,int iByteOffset,unsigned long *pul);
int GetFrameAverage(TRAIN *tr,int fr,unsigned long *pul,int a);

void PrintFileheader(TRAIN *tr,int fi);
void PrintFrameheader(TRAIN *tr,int fr);
int ConvertImage(TRAIN *tr,int fr,int a);
void DisplayImage(TRAIN *tr,int fr,float *imbuffer);
int ConvertImages(TRAIN *tr,int fr1,int fr2,int a);
void DisplayImages(TRAIN *tr,int fr1,int fr2,float *imbuffer);
void Pallet(int type,float contra,float bright);
void DisplayMap(TRAIN *tr,double *mx,double *my);
void DisplayVectorField2(int dx,int dy,double *mx,double *my,float *fmx,float *fmy,int mode);
void DisplaySeries(TRAIN *tr,int n,double **buf);
int CocktailAverage(TRAIN *tr);
void DisplayCocktailAverage(TRAIN *tr);
int Statistics(TRAIN *tr);
int TraceAnalysis(TRAIN *tr,int ifr,int ffr,int a);
int GreenManager(int iVersion,FILE *fp,char *pc);
int Compressor(TRAIN *tr,int mode);
int Synchronization(TRAIN *tr,int ifr,int ffr,int *iifr,int *fffr,int a);
int InterframeTime(TRAIN *tr,int ifr,int ffr,int *iifr,int *fffr,int a);
int Inverse(int n,double omega,double alpha,int h);
double PolynomialFitS(int k,int n,double x);
void PolynomialFitM(int k,int n,double x,double *y,double *z,int mode);
void options(int argc,char **argv);
int AccumulateTimeAverage(TRAIN *tr,int a,double *t_l,double **t_b);
int AdvanceTimeAverage(TRAIN *tr,int fr,int a,double *t_l,double **t_b);
int BinFrames(TRAIN *tr,int a);
int AddFrameToBins(TRAIN *tr,int fr,int a);
int Double2Float(unsigned long xy ,double *dbuf,float *fbuf);
int EPSTParser(TRAIN *tr,int ifr,int ffr,int *iifr,int *fffr);

int AskToQuit(char *pcLineOut);

char* GetString(char *str,int iMode);
int GetFloatPoint(float *pfX,float *pfY,int iMode,int iColor);
int GetIntPoint(int *piX,int *piY,int iDX,int iDY,int iMode,int iColor);


#define FRAME_RULE_N 2
#define FRAME_RULE_0 0 // single cycle per cycle
#define FRAME_RULE_1 1 // EPST frame rule 
#define FRAME_RULE_DEFAULT FRAME_RULE_0 

int FrameRule0(TRAIN *tr,int fr,CAR *pcar);
int FrameRule1(TRAIN *tr,int fr,CAR *pcar);

int (*pfunFrameRule[FRAME_RULE_N])(TRAIN *tr,int fr,CAR *pcar);

extern int imagetype;
extern char device[40];
extern int mainid;
extern float pminx,pmaxx,pminy,pmaxy;
extern float *fbuffer,*fbuffer_removed;
extern unsigned short Xdim,Ydim;
extern unsigned long XYdim;
extern float bg,fg,trans[6];
extern int ibac,ifor;
extern int n_panels,panel;
extern int scheme,difference;
extern int initframe,increment,increment1;
extern int finalframe;
extern int single_multiple;
extern unsigned short *dumbuffer;
extern double *mapx,*mapy,*map0;
extern double *timeaverage_buffer;
extern double timeaverage_buffer_frames_in;
extern float *cocktailaverage,camin,camax;
extern double phi0,dphi;
extern int radius,cradius;
extern int radiusbig;
extern int timeradius,sacrifice_boundaries;
extern double dradius,dradiusbig;
extern int *aindex_x,*aindexbig_x,*aindex_y,*aindexbig_y;
extern int aindex_n,aindexbig_n;
extern float bgg,fgg;
extern float bga,fga;
extern int image_count;
extern int average;
extern double harmonic;
extern unsigned int n_cycles;
extern int do_contours;
extern int do_verbose;
extern int do_statistics;
extern int do_trace;
extern int do_synch;
extern int do_precise,do_intervals;
extern int do_cocktailaverage;
extern int do_timeaveraging;
extern int display_removed;
extern int cut_offl,cut_offu;
extern int nfiles;
extern char *filenames[];
extern double period;
extern double *synch_phi,*synch_cos,*synch_sin;
extern double synch_phase_shift,synch_omega;
extern double cyclesupper,cycleslower,cyclesexact;
extern float *synch_trace;
extern int additionalframes;
extern double stimulus_period;
extern double omega_frame;
extern double frame_omega;
extern int find_best_cycle;
extern int do_remove_linear,do_remove_curve;
extern double *remove_sy,*remove_sxy;
extern double a_cos,a_sin,b_cos,b_sin;
extern int poly_fitn;
extern double **remove_curve;
extern double *remove_cos,*remove_sin,*remove_dum;
extern int nframes_good;
extern double timeaverage_cycles;
extern double timeaverage_blocks;
extern int timeaverage_cycles_int;
extern int replace_scheme;
extern unsigned long memory_size;
extern unsigned int max_frames_in_cue;
extern int compress;
extern int initframe_ta,finalframe_ta;
extern int n_bins;
extern double **frame_bins;
extern int *frames_in_bin;
extern int displayimageN,displayimagesN;
extern double cycles_per_frame;
extern int frame_rule;
extern int save_in_files;
extern int prompt_to_save_images;
extern int statistics_num_bins;
extern unsigned long ulSynchType;
extern int iSynchChannel;
extern const char strListOfKnownExperimentChunks[];
extern int iDoDivisionByAverage;
extern int iDivisionByAverageCount;
extern double dDivisionByAverageFactor;
extern double dBinDivisionByAverageFactor;
extern unsigned long ulBinX,ulBinY,ulBinT;
extern int iKeepOriginalFiles;
extern int iRemoveTrailingZsFromDecompressedFileName;
extern float fContrast;
extern int iDrawAxesAndLabelOnGreen;
extern float fSizeOfGreenSave;
extern float fSizeOfTraceInteractiveWindow;
extern float fSizeOfTracePlotWindow;
extern int iUseExternalFourierThreshold;
extern float fExternalFourierThreshold;
