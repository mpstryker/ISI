/* Map analysis - multiple(2)   */
/* Written by V.Kalatsky        */
/* 06Feb01                      */
/* Modified 5May01              */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cpgplot.h>
//#include <time.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
//#include <fcntl.h>
//#include <tcl.h> 
#include <errno.h>

typedef struct Pinwheel {
  double charge;
  double x,y;
} Pinwheel;

extern int errno;

#define SAFE_FREE(p) {if(p){free(p); p=NULL;}}

#define VECTPROD(x1,y1,x2,y2) ((x1)*(y2)-(x2)*(y1))
#define DOTPROD(x1,y1,x2,y2) ((x1)*(x2)+(y1)*(y2))

#define DOUBLE 0
#define FLOAT  1

#define SAVE_XY_IN_TWO_FILES 0
#define SAVE_XY_IN_ONE_FILE 1

#define SAVE_TYPE_BIT_DOUBLE_FLOAT 1
#define SAVE_TYPE_BIT_ADD_DOUBLE_AT_END 2
#define SAVE_TYPE_BIT_XY_IN_FILES 4

#define MAX_FILE_N 4

#define DEFAULT_RADIUS 0.0
#define DEFAULT_RADIUSBIG 0.0

#define MAIN_WIN_SIZE 10.0

#define MODE_INVERSE_NONE 0
#define MODE_INVERSE_Y    1
#define MODE_INVERSE_X    2
#define MODE_INVERSE_DEFAULT MODE_INVERSE_Y

#define FLIP_PHASE_SIZE_DEFAULT 50
#define FLIP_PHASE_THRESHOLD_DEFAULT 0.6

void start_graphics(int paneln);
void end_graphics();
void Pallet(int type,float contra,float bright);

int InitializeMaps();
int InitializeBuffers();

void DisplayMap(int panel,int xd,int yd,double **mp);
void DisplayOd();

// Analysis
void AverageMaps();
int AnalyzePW(int xd,int yd,double *mx,double *my);
int Rotator(int bits,int *rbits,int *position,int *angle);
int Charge(int x,int y);
double FindMapAngle(int iX,int iY,double *pdMX1,double *pdMY1,double *pdMX2,double *pdMY2,int iMode);
int FourierMax(int iNC,double *pdC,int iMode);
double FourierSum(double x,int iNC,double *pdC,int iMode);
double FourierSumDir(double x,int iNC,double *pdC,int iMode);

int PlotDominanceHistogram(int iX,int iY,double *pdMX1,double *pdMY1,double *pdMX2,double *pdMY2,int iMode);

void Interactive(int panel,int xd,int yd,double **dmp);

void options(int,char**);

extern Pinwheel *apinwheels;
extern int apinwheels_n;

extern char device[40];
extern int mainid,odid;
extern float pminx,pmaxx,pminy,pmaxy;
extern float bg,fg,trans[6];

extern int nfiles;
extern char *filenames[MAX_FILE_N];
extern char *filestr1,*filestr2;
extern char savefile[256];
extern char *modestr[];
extern int type;
extern unsigned short Xdim,Ydim;
extern unsigned long XYdim;
extern float *fbuffer,*fbufferr,*fbufferx,*fbuffery,*fbufferi,dumf;
extern int radius,radiusbig;
extern double dradius,dradiusbig;
extern int *aindex_x,*aindexbig_x,*aindex_y,*aindexbig_y;
extern int aindex_n,aindexbig_n;
extern double **maps;
extern double *dbufferx,*dbuffery,*dbufferphix,*dbufferphiy;
extern double *pdBufferRho1,*pdBufferRho2;
extern int do_od,do_pinwheels,do_contours;
extern double harmonic_d;
extern double **dbuffers;
extern int do_interactive;
extern int main_win_size;
extern int npanels;
extern int make_phi_positive;
extern int iInverseMode;
extern int iFlipPhaseSize;
extern double dFlipPhaseThreshold;
extern int iROIXLeft,iROIXRight,iROIYTop,iROIYBottom;
extern int iCropROI;
extern char *pcSaveFileSuffix;
extern int iDoDominanceHistogram;
extern double dDominanceThreshold;
extern int iTakeSQRTOfSum;
