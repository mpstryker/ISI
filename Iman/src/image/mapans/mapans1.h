/* Map analysis                 */
/* Written by V.Kalatsky        */
/* 06Feb2001                    */
/* Modified 11May2001           */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cpgplot.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>

typedef struct Pinwheel {
  double charge;
  double x,y;
} Pinwheel;

typedef struct PAIR {
  unsigned long n;
  float x1,y1;
  float x2,y2;
} PAIR;

typedef struct POINT {
  float x,y;
} POINT;

typedef struct POINTD {
  double x,y;
} POINTD;

typedef struct RECT {
  float x1,y1;
  float x2,y2;
} RECT;

#define RECORD_HIDE_NOT      0
#define RECORD_HIDE_INT      1
#define RECORD_HIDE_EXT      2
#define RECORD_HIDE_EXT_FAKE 4

typedef struct RECORD{
  int iIndex;
  double dX,dY,dZ;
  double dZR;
  double dPhi,dRho;
  double dPhiR;
  int iHide;
} RECORD;

typedef struct DOMAIN2{
  int iIndex;
  int iNPoints;
  struct POINT *pPoints;
} DOMAIN2;

#define FILTER_KIND_NONE   0
#define FILTER_KIND_GAUSS  1
#define FILTER_KIND_COS    2
#define FILTER_KIND_RANDOM 4

typedef struct FILTER{
  int iKind;
  int iNPoints;
  double dRadius;
  double dDepth;
  int *piX;
  int *piY;
  double *pdValues;
} FILTER;

typedef struct MARKRECORD{
  int iIndex;
  int iSymbol;
  int iSelect;
  float x,y;
}MARKRECORD;

extern int errno;

#define SAFE_FREE(p) {if(p){free(p);p=NULL;}}

#define VECTPROD(x1,y1,x2,y2) ((double)(x1)*(double)(y2)-(double)(x2)*(double)(y1))
#define DOTPROD(x1,y1,x2,y2) ((double)(x1)*(double)(x2)+(double)(y1)*(double)(y2))
#define NORM(x1,y1,x2,y2) sqrt((double)(((x1)*(x1)+(y1)*(y1))*((x2)*(x2)+(y2)*(y2)))) 
#define HYPOT2(x1,y1) ((x1)*(x1)+(y1)*(y1)) 
#define SQUARED(x) ((x)*(x)) 
#define INSIDE(x,l,r) (x>=l && x<r)

#define DOUBLE 0
#define FLOAT  1

#define SAVE_XY_IN_TWO_FILES 0
#define SAVE_XY_IN_ONE_FILE 1

#define SAVE_TYPE_BIT_DOUBLE_FLOAT 1
#define SAVE_TYPE_BIT_ADD_DOUBLE_AT_END 2
#define SAVE_TYPE_BIT_XY_IN_FILES 4
#define SAVE_TYPE_BIT_GREEN 8

#define MAX_FILE_N 2

#define DEFAULT_RADIUS 0.0
#define DEFAULT_RADIUSBIG 0.0
#define DEFAULT_CHOP_RADIUS 0.0
#define DEFAULT_SMOOTH_RADIUS 0.0
#define DEFAULT_SECTION_RADIUS 1.0

#define WINDOW_SIZE 4.0
#define INTERACTIVE_WINDOW_SIZEX 5.0

#define DISPLAY_MODE2 1
#define DISPLAY_MODE4 2
#define DISPLAY_MODE6 3
#define DEFAULT_DISPLAY_MODE DISPLAY_MODE2
#define DISPLAY_MODE_MASK 63
#define DISPLAY_MODE_SINGLE 128
#define DISPLAY_MODE_SINGLE2 256
#define DISPLAY_MODE_NOT 512

#define DEFAULT_FIELD_MODE DISPLAY_MODE4

#define DEFAULT_SAVE_MODE DISPLAY_MODE2

#define BIG_FLOAT 10e10

#define INT_MODE_NON 0
#define INT_MODE_RHO 1
#define INT_MODE_PHI 2
#define INT_MODE_DEFAULT INT_MODE_PHI //will be switched anyway

#define CONTOUR_NON 0
#define CONTOUR_PHI 1
#define CONTOUR_RHO 2

#define CONTOUR_NON_NON 0
#define CONTOUR_PHI_NON 1
#define CONTOUR_RHO_NON 2
#define CONTOUR_NON_PHI 4
#define CONTOUR_PHI_PHI 5
#define CONTOUR_RHO_PHI 6
#define CONTOUR_NON_RHO 8
#define CONTOUR_PHI_RHO 9
#define CONTOUR_RHO_RHO 10
#define CONTOUR_DEFAULT CONTOUR_PHI_RHO

#define CONTOUR_TYPE_DUMB 0
#define CONTOUR_TYPE_SMART 1
#define CONTOUR_TYPE_DEFAULT CONTOUR_TYPE_SMART

#define COLOR_N 2
#define COLOR_RGB 0 
#define COLOR_HLS 1 
#define COLOR_DEFAULT COLOR_RGB

#define AVERAGE_LOW 1
#define AVERAGE_HIGH 2
#define AVERAGE_CHOP 3

#define DEFAULT_SMART_CONTOURS_THRESHOLD 0.1
#define DEFAULT_CONTOURS_INI_ANGLE_RAD (0.0)

#define MAX_SPATIAL_BINNING 16

#define PALLET_FUNNY -2
#define PALLET_HEAT -1
#define PALLET_RGB1 0
#define PALLET_RGB2 1
#define PALLET_HSL6 2
#define PALLET_HSL12 3
#define PALLET_HSL18 4
#define PALLET_HSL24 5
#define DEFAULT_PALLET PALLET_HSL6
//#define DEFAULT_PALLET PALLET_FUNNY
#define PALLET_SHIFT_SHIFT 1000
#define PALLET_SPECIAL 1001

#define PALLET_REDUCED 2000
#define PALLET_NULL_OUTSIDERS 2001

#define INVERSION_NONE 0
#define INVERSION_X 1
#define INVERSION_Y 2
#define DEFAULT_INVERSION INVERSION_Y

#define DEFAULT_PIXEL_SIZE_MICRON 8.889

#define LINE_SEGMENTS_STYLE_SIMPLE     0
#define LINE_SEGMENTS_STYLE_ARROW      1
#define LINE_SEGMENTS_STYLE_MARK_ENDS  2
#define LINE_SEGMENTS_STYLE_ENUMERATE  4

#define DEFAULT_BLEND_BORDER 30

#define MODE_LOAD_BITMAP_REPLACE 0
#define MODE_LOAD_BITMAP_AND     1
#define MODE_LOAD_BITMAP_OR      2
#define MODE_LOAD_BITMAP_DEFAULT MODE_LOAD_BITMAP_REPLACE

#define HIDDEN_COLOR_N 90
#define HIDDEN_COLOR_R 0.19
#define HIDDEN_COLOR_G 0.19
#define HIDDEN_COLOR_B 0.19

#define MARK_SELECTION_MODE_N       3
#define MARK_SELECTION_MODE_INCLUDE 0
#define MARK_SELECTION_MODE_EXCLUDE 1
#define MARK_SELECTION_MODE_SPECIAL 2

#define DISTRIBUTION_PHI_BINS_N 100
#define DISTRIBUTION_PHI_OD_BINS_N 70
#define DISTRIBUTION_RHO_BINS_N 100

void start_graphics(int mode,char *str,int *id);
void end_graphics(int id);
void Pallet(int type,float contra,float bright);
void MakePallet(int type,float fS,float fPr,float fPo);
void MakeImagePallet(int type,float fS,float fPr,float fPo);

int InitializeMaps();
int InitializeBuffers();
int SetFilter(FILTER *pFilter,int xd,int yd);

void AverageMaps(int xd,int yd,double **dm,float **fm,int nn,double dr,double drb,int mode);
int AnalyzePW(int xd,int yd,double *mx,double *my);
int SmartContoursPolar(int xd,int yd,float *rho,float *phi,PAIR **sc,int *scn,float value,float threshold);
int SmartContoursCartesian(int xd,int yd,float *x,float *y,PAIR **sc,int *scn,float val,float threshold);
void Cartesian2Polar(unsigned long nn,float *x,float *y,float *rho,float *phi);
void Polar2Cartesian(unsigned long nn,float *rho,float *phi,float *x,float *y);
void Double2Float(unsigned long nn,double *pd,float *pf);
void FindExtrema(float *fb,float *fmin,float *fmax,unsigned long *imin,unsigned long *imax);
void RemoveBias(int xd,int yd,double *mx,double *my);
void DisplayMap(int mode,int xd,int yd,double *mx,double *my,double *gx,double *gy);
void DisplayFields(int mode,int xd,int yd,double *mx,double *my,double *gx,double *gy);
void Interactive(int xd,int yd,double **mp,int *mode);
void DrawGrayImage(float *im,int xd,int yd,float max,float min,float *tran);
void DrawColorImage(float *im,int xd,int yd,float max,float min,float *tran);
void DrawNiceColorImagePolar(float *phi,float *rho,int xd,int yd,float max,float min,float *tran);
void DrawNiceColorImageCartesian(float *mx,float *my,int xd,int yd,float *tran);
void DrawContours(int mode,int scheme,float **sb,int xd,int yd,float fmax,float fmin,float threshold,float *tran);
void DrawColorWedge(float fWedgePos,float fWedgeWidth,float fMin,float fMax,char *strUnits);
void PlotSection(int mode,int xd,int yd,float *pfMapX,float *pfMapY);
void MarkPoints(int iDX,int iColor,int iMode);
void PlotLineSegments(int iColor,int iWidth,int Style);


int Rotator(int bits,int *rbits,int *position,int *angle);
int Charge(int x,int y);
void Stencil(int xd,int yd,int n,POINTD *pDom,unsigned long *sn,unsigned long **s);
void StencilMap(int xd,int yd,float *pfX,float *pfY,int n,double *dx,double *dy,unsigned long *sn,unsigned long **s);
int LoadStencilBitmap(int xd,int yd,unsigned long *sn,unsigned long **s,int iMode);
int SaveStencilBitmap(int xd,int yd,unsigned long sn,unsigned long *s);
int AddPixelToStencil(int x,int y,int xd,int yd,unsigned long *sn,unsigned long **s,int iMode);


void options(int,char**);

int RectSegmentIntersect(POINT *pIniSegment,POINT *pFinSegment,RECT rect);
int SegmentSegmentContinuationIntersect(POINT *pSegment1,POINT *pSegment2,POINT *pPoint);
int SegmentSegmentIntersect(POINT *pSegment1,POINT *pSegment2,POINT *pPoint);
int IsUnitSquareInsideCircle(int iX,int iY,double dX,double dY,double dRadius2);
int IsPartOfUnitSquareInsideCircle(int iX,int iY,double dX,double dY,double dRadius2);
int IsPointInsideDomain(double dX,int dY,int n,POINTD *pDom,RECT *pRect);

double PolynomialFitS(int k,int n,double x);
void PolynomialFitM(int k,int n,double x,double *y,double *z,int mode);

double ran3(long *idum);
int CorrelationToExternalSource(int iN,RECORD *pRec,int xd,int yd,double **mp,int iMode);

char* GetString(char *str);
int GetStringFromFile(FILE *pF,char *str,int iNChars,int iMode);


extern Pinwheel *apinwheels;
extern int apinwheels_n;

extern char device[40],devices[40];
extern int mainid,saveid,fieldid;
extern float pminx,pmaxx,pminy,pmaxy;
extern float bg,fg,trans[6];
extern float transs[6];
extern float wedge_width,wedge_pos;

extern int nfiles;
extern char *filenames[2];
extern char *filestr;
extern int type;
extern unsigned short Xdim,Ydim;
extern unsigned long XYdim;
extern float *fbufferx,*fbuffery;
extern double *pdBufferX,*pdBufferY;
extern float *fbufferphi,*fbufferrho;
extern float *fbufferi;
extern float *fbufferdiv,*fbufferrot,*fbufferdivrot_phi,*fbufferdivrot_rho;
extern int radius,radiusbig,cradius,sradius;
extern double dradius,dradiusbig,dcradius,dsradius;
extern int *aindex_x,*aindex_y;
extern int aindex_n;
extern double **maps,**maps_orig;
extern float **fmaps_orig;
extern double **dbuffers,**dbuffersa;
extern float **smoothfbuffers;
extern int do_pinwheels,do_fields,do_profiles,do_dump,do_marking;
extern int do_zerolines,do_rescaling,do_verbose;
extern int do_smartcontours,depth_lines,depth_lines_rho,depth_lines_phi;
extern int do_wedge,do_polarcontours;
extern int do_sectionplot;
extern int do_pixel_size_transform;
extern int do_discontinuous_phase;
extern int removebias;
extern int harmonic;
extern double harmonic_d;
extern int int_mode;
extern int save_mode;
extern int display_mode;
extern int field_mode;
extern int iNDomainPoints;
extern POINTD *pDomain;

extern int iShiftDomainForIntegerPoints;
extern unsigned long *stencil,stenciln;
extern int smaxx,sminx,smaxy,sminy;
extern float *contours;
extern float contours_iniangle,emphasize_contours_iniangle;
extern float fContourSpacingAngle;
extern double alpha,inversion;
extern int contour_scheme;
extern int contour_scheme_phi,contour_scheme_rho;
extern double xshift,yshift;
extern float window_size;
extern float fInteractiveWindowWidth;
extern int do_phase_transform;
extern double phase_transform;
extern int do_distribution;
extern int pallet;
extern int color_scheme;
extern PAIR *smart_contours;
extern int smart_contoursn;
extern PAIR **ppPAIRSmartContours;
extern int *piSmartContoursN;
extern PAIR **ppPAIRSmartContoursExternal;
extern int *piSmartContoursNExternal;
extern int iNDepthLinesExternal;
extern int iContoursColorExternal;
extern float amplitude_chop_l,amplitude_chop_u;
extern float *chopfbuffer;
extern float smartcontours_threshold;
extern int spatial_binning;
extern int inversion_scheme;
extern int save_on;
extern POINT section[2];
extern RECT mainRect;
extern double dPixelSizeU;
extern double dSectionRadius;
extern int iDoWedgeShift;
extern float fWedgeShift;
extern float fWedgeMin,fWedgeMax;
extern char *strWedgeUnits;
extern int iDoWedgeAnnotation;
extern int iRotateMaps;
extern int iFlipMaps;
extern int iDoLargeWedge;
extern int iReplotSavedContours;
extern int iDrawColorWedge;
extern int iDrawBWWedge;
extern int iWhiteBackgroundOnPolarPlot;
extern float fWedgeVeryLargeCH;
extern float fWedgeLargeCH;
extern float fWedgeSmallCH;
extern int iWedgeFont;
extern int iPlotExternalContours;
extern int iEmptySource;
extern int iOverlayContoursColor;
extern int iDoNotPlotPhase;
extern int iColorContours;
extern int iContourLinesWidth;
extern int iContoursColorExternalSave;
extern int iContoursColorExternalInitial;
extern int iContoursColorExternalInitialSave;
extern int iBlankSource;
extern double dRhoNoise;
extern int iGetNoiseLevel;
extern int iSectionSinglePixel;
extern int iMarkMinMax;
extern int iCutBlanksOnWedge;
extern float fStim;
extern float fPostBlank;
extern float fPreBlank;
extern int iReducedColorWedge;
extern int iNullOutsiders;
extern int iMarkPoints;
extern int iMarkPointsSymbol;
extern float fMarkPointsRadius; 
extern float fMarkPointsCharSize; 
extern int iMarkPointsColor;
extern int iMarkPointsN;
extern MARKRECORD *pMarkPoints;
extern int iMarkPointsHide;
extern int iMarkPointsEnumerate;
extern int iMarkPointsLoadN;
extern int iMarkPointsSelectN;
extern int iMarkPointsSelectMode;
extern int iFullSizeViewPort;
extern int iUseRhoMinOnSectionAmplitudePlot;
extern long lSeed;
extern int iDoPhaseScale;
extern double dPhaseScale;
extern double dStimPeriod;
extern double dPreStimTime;
extern double dPostStimTime;
extern int iPrePostStimTimeRescale;
extern double dCorrelationThreshold;
extern int iNCorrRecords;
extern RECORD *pCorrRecords;
extern int iCorrelationModeShift;
extern int iNLineSegments;
extern PAIR *pLineSegments;
extern int iPlotLineSegments;
extern int iPlotLineSegmentsStyle;
extern int iSaveCropBorderAsLineSegments;
extern int iBlendBorder;
extern int iBlendBorderPixel;
extern int iSaverPinwheels;
extern int iLoadBitmapMode;
extern int iFindMarkedInsideSelected;
extern int iMarkCorrelationPoints;
extern int iAdditiveCropChanges;
extern FILTER stFilter;
extern double dLensRadius;
extern int iLensCorrection;
extern int iFixedGrayScaleLo,iFixedGrayScaleHi;
extern float fFixedGrayScaleLo,fFixedGrayScaleHi;
extern int iDoAmplitudeNorm;
extern double dAmplitudeNormCoeff;
extern int iDoODAnalysis;
extern int iODAnalysisScheme;
extern int iODAnalysisNBins;
extern char *pcSaveFileSuffix;
