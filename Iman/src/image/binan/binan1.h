/* Bin analysis                 */
/* Written by V.Kalatsky        */
/* 20Jul2001                    */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cpgplot.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>
#include <pthread.h> 

extern int errno;

#define SAFE_FREE(p) {if(p){free(p);p=NULL;}}

#define VECTPROD(x1,y1,x2,y2) ((x1)*(y2)-(x2)*(y1))
#define DOTPROD(x1,y1,x2,y2) ((x1)*(x2)+(y1)*(y2))
#define NORM(x1,y1,x2,y2) sqrt((double)(((x1)*(x1)+(y1)*(y1))*((x2)*(x2)+(y2)*(y2)))) 
#define HYPOT2(x1,y1) ((x1)*(x1)+(y1)*(y1)) 

#define DOUBLE 0
#define FLOAT  1

#define MAX_FILE_N 1

#define DEFAULT_RADIUS 0.0
#define DEFAULT_RADIUSBIG 0.0
#define DEFAULT_CHOP_RADIUS 0.0

#define WINDOW_SIZE 10.0

#define DISPLAY_MODE2 1
#define DISPLAY_MODE4 2
#define DISPLAY_MODE6 3
#define DEFAULT_DISPLAY_MODE DISPLAY_MODE2
#define DISPLAY_MODE_MASK 63
#define DISPLAY_MODE_SINGLE 128

#define BIG_FLOAT 10e10

#define INT_MODE_NON -1
#define INT_MODE_DEFAULT 0

#define AVERAGE_LOW 1
#define AVERAGE_HIGH 2
#define AVERAGE_CHOP 3

#define MAX_SPATIAL_BINNING 16

#define DEFAULT_ANIMATION_DELAY 10000L

#define DEFAULT_BIN_SIZE_SEC 1.0

void start_graphics(int mode,char *str,int *id);
void end_graphics(int id);
void Pallet(int type,float contra,float bright);

int InitializeMaps();
int InitializeBuffers();
void AverageMaps(int xd,int yd,double **dm,float **fm,int nn,double dr,double drb,int mode);
void Double2FloatShort(unsigned long nn,double *pd,float *pf);
void Double2FloatLong(unsigned long nn,double *pd,float *pf,float *min,float *max);
void FindExtrema(float *fb,float *fmin,float *fmax,unsigned long *imin,unsigned long *imax);
void FindExtremaEx(unsigned long n,unsigned long nn,float *fb,float *fmin,float *fmax,unsigned long *imin,unsigned long *imax);
void FindMax(unsigned long n,float *fb,float *fmax,unsigned long *imax);
void FindMin(unsigned long n,float *fb,float *fmin,unsigned long *imin);
int FindFTComponents(int iN,float *pfT,int iC,float *pfRe,float *pfIm);
int FindFTComponents_DOUBLE(int iN,double *pdT,int iC,double *pdRe,double *pdIm);
double FindExtremaFComponents(int iN,float *pfFC,int iNC,int iCMin,int iSign);
double FourierSum(int iN,float *pfFC,int iNC,int iCMin,float fT);
double FourierSumD(int iN,float *pfFC,int iNC,int iCMin,float fT);
double FindExtremaFComponents_DOUBLE(int iN,double *pdFC,int iNC,int iCMin,int iSign);
double FourierSum_DOUBLE(int iN,double *pdFC,int iNC,int iCMin,double dT);
double FourierSumD_DOUBLE(int iN,double *pdFC,int iNC,int iCMin,double dT);
void MakeMaps(int iNB,double **ppdMaps,int iXdim,int iYdim,double *pdMapX,double *pdMapY);

float FindAverage(unsigned long n,float *fb,float *faverage);
void RemoveBias(int xd,int yd,double *db,float *fb,int mode);
void DisplayBins(int mode,int xd,int yd,int n,float **bins);
void AnimateBins(int mode,int xd,int yd,int n,float **bins,unsigned long delay);
void DisplayMap(int xd,int yd,double **mp);

void Interactive(int xd,int yd,float **mp,int *mode);
void Stencil(int xd,int yd,int n,double *dx,double *dy,unsigned long *sn,unsigned long **s);
void options(int,char**);

void AnimateBinsThread(void *arg);

typedef struct AnimateBinsThreadArgs {
  int mode;
  int xd;
  int yd;
  int n;
  float **bins;
  unsigned long delay;
} AnimateBinsThreadArgs;

extern pthread_t animate_thread;
extern pthread_attr_t animate_thread_attr;
extern AnimateBinsThreadArgs animatebinsthreadargs;
extern int thread_return_value;
extern int stop_animation;
extern unsigned long animation_delay;
extern int animate_forward;

extern char device[40],device_ask[40];
extern int mainid,saveid,fieldid;
extern float pminx,pmaxx,pminy,pmaxy;
extern float bg,fg,trans[6];
extern float transs[6];
extern float wedge_width,wedge_pos;

extern int nfiles;
extern char **filenames;
extern char *filestr;
extern int type;
extern unsigned short Xdim,Ydim;
extern unsigned long XYdim;
extern float *fbufferi;
extern float *fbuffer;
extern int radius,radiusbig,cradius;
extern double dradius,dradiusbig,dcradius;
extern int *aindex_x,*aindex_y;
extern int aindex_n;
extern double **maps,**maps_orig;
extern float **fmaps,**fmaps_orig;
extern int do_marking;
extern int do_verbose;
extern int do_wedge;
extern int do_animation,do_fixed_background;
extern int removebias;
extern int harmonic;
extern double harmonic_d;
extern int int_mode;
extern int display_mode;
extern unsigned long *stencil,stenciln;
extern int smaxx,sminx,smaxy,sminy;
extern float window_size;
extern int spatial_binning;
extern int save_on;
extern int nbins;
extern int panelnx,panelny;
extern int ibac,ifor;
extern float *bin_max,*bin_min;
extern int shift_max_to_zero;
extern int shift_average_to_zero;
extern int iNormalizeBinning;
extern float fBinSizeSec;
extern int iPlotFComponents;
extern int iNFComponents;
extern int iMinFComponent;
extern double *pdMapX,*pdMapY;
extern int iGenerateMaps;
extern int iMakeMaps;
