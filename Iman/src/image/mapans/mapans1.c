/* Map analysis                 */
/* Written by V.Kalatsky        */
/* 06Feb2001                    */
/* Modified 10May2001           */

#include "mapans1.h"

Pinwheel *apinwheels=(Pinwheel *)0;
int apinwheels_n=0;

char device[40]="/xw\0",devices[40]="?\0";
int mainid,saveid,fieldid;
float pminx,pmaxx,pminy,pmaxy;
float bg,fg,trans[6]={-1.0,1.0,0.0,-1.0,0.0,1.0};
float transs[6]={-1.0,1.0,0.0,-1.0,0.0,1.0};
//float wedge_width=0.4,wedge_pos=0.0;
// Large werge
//Overlays image
float wedge_width=0.37,wedge_pos=-0.11;
//Does not overlay image
//float wedge_width=0.35,wedge_pos=0.01;
//Does not overlay image, use to plot fake images with wedge without annotation
//float wedge_width=0.28,wedge_pos=0.01;
// Val: use to plot wedge only
//float wedge_width=9.0,wedge_pos=-3.5;

int nfiles=0;
char *filenames[2];
char *filestr;
int type;
unsigned short Xdim,Ydim;
unsigned long XYdim;
float *fbufferx=NULL,*fbuffery=NULL;
double *pdBufferX=NULL,*pdBufferY=NULL;
float *fbufferphi=NULL,*fbufferrho=NULL;
float *fbufferi;
float *fbufferdiv,*fbufferrot,*fbufferdivrot_phi,*fbufferdivrot_rho;
int radius,radiusbig,cradius,sradius;
double dradius=DEFAULT_RADIUS,dradiusbig=DEFAULT_RADIUSBIG,dcradius=DEFAULT_CHOP_RADIUS,dsradius=DEFAULT_SMOOTH_RADIUS;
double **maps,**maps_orig;
float **fmaps_orig;
double **dbuffers;
float **smoothfbuffers;
int do_pinwheels=0,do_fields=0,do_profiles=0,do_dump=0;
int do_zerolines=0,do_rescaling=0,do_verbose=0;
int do_phase_transform=0,do_distribution=0;
int do_wedge=1,do_polarcontours=0;
int do_sectionplot=0;
int do_pixel_size_transform=0;
int do_discontinuous_phase=0;
int depth_lines=0,depth_lines_rho=0,depth_lines_phi=0,do_smartcontours=CONTOUR_TYPE_DEFAULT;
int iSaverPinwheels=0;
int iDoWedgeShift=0;
int iDoWedgeAnnotation=0;
float fWedgeShift;
float fWedgeMin,fWedgeMax;
char *strWedgeUnits=NULL;
int removebias=0;
int harmonic=1;
double harmonic_d=1.0;
int int_mode=INT_MODE_NON;
int save_mode=DEFAULT_SAVE_MODE;
int display_mode=DEFAULT_DISPLAY_MODE;
int field_mode=DEFAULT_FIELD_MODE;
int iNDomains;
DOMAIN2 *pDomains=NULL;
int iNDomainPoints=0;
POINTD *pDomain=NULL;
int iShiftDomainForIntegerPoints=1;
unsigned long *stencil=(unsigned long *)0,stenciln=0;
int smaxx,sminx,smaxy,sminy;
float *contours=(float *)0;
float contours_iniangle=DEFAULT_CONTOURS_INI_ANGLE_RAD,emphasize_contours_iniangle=1;
float fContourSpacingAngle=0.0;
double alpha=0.0,inversion=1.0;
int contour_scheme=CONTOUR_DEFAULT;
int contour_scheme_phi=(CONTOUR_DEFAULT & 3),contour_scheme_rho=(CONTOUR_DEFAULT >> 2);
double xshift=0.0,yshift=0.0;
float window_size=WINDOW_SIZE;
float fInteractiveWindowWidth=INTERACTIVE_WINDOW_SIZEX;
double phase_transform;
int pallet=DEFAULT_PALLET;
int color_scheme=COLOR_DEFAULT;
PAIR *smart_contours=(PAIR *)0;
int smart_contoursn=0;
PAIR **ppPAIRSmartContours=(PAIR **)0;
int *piSmartContoursN=NULL;
PAIR **ppPAIRSmartContoursExternal=(PAIR **)0;
int *piSmartContoursNExternal=NULL;
int iNDepthLinesExternal=0;
float amplitude_chop_l=0.0,amplitude_chop_u=1.0;
float *chopfbuffer=(float*)0;
float smartcontours_threshold=DEFAULT_SMART_CONTOURS_THRESHOLD;
int spatial_binning=1;
int inversion_scheme=DEFAULT_INVERSION;
int save_on=0;
POINT section[2];
RECT mainRect;
double dPixelSizeU=DEFAULT_PIXEL_SIZE_MICRON;
double dSectionRadius=DEFAULT_SECTION_RADIUS;
int iRotateMaps=0;
int iFlipMaps=0;
int iDoLargeWedge=0;
int iReplotSavedContours=0;
int iDrawColorWedge=0;
int iDrawBWWedge=1;
int iWhiteBackgroundOnPolarPlot=0;
float fWedgeVeryLargeCH=15.5;
float fWedgeLargeCH=10.5;//10.5
float fWedgeSmallCH=5.5;
int iWedgeFont=1;
int iPlotExternalContours=1;
int iContoursColorExternal=0;
int iContoursColorExternalSave=0;
int iContoursColorExternalInitial=2;//1;
int iContoursColorExternalInitialSave=2;//1;
int iOverlayContoursColor=0;
int iDoNotPlotPhase=0;
int iColorContours=0;
int iContourLinesWidth=2;
int iBlankSource=0;
double dRhoNoise=-1.0;
int iGetNoiseLevel=0;
int iSectionSinglePixel=0;
int iMarkMinMax=0;
int iCutBlanksOnWedge=0;
float fStim=-1.0;
float fPostBlank=-1.0;
float fPreBlank=-1.0;
int iReducedColorWedge=0;
int iNullOutsiders=0;
int iMarkPoints=0;
int iMarkPointsSymbol=-3; //1-square, 3-*, 4-circle, 9-circle+dot
int iMarkPointsColor=1;
//Large
//float fMarkPointsRadius=4.0; // if iMarkPointsSymbol=-3
//float fMarkPointsCharSize=1.4;
//Small 
float fMarkPointsRadius=1.3; // if iMarkPointsSymbol=-3
float fMarkPointsCharSize=1.0;//0.6; 
int iMarkPointsN=0;
MARKRECORD *pMarkPoints=NULL;
int iMarkPointsHide=0;
int iMarkPointsEnumerate=0;
int iMarkPointsLoadN=0;
int iMarkPointsSelectN=0;
int iMarkPointsSelectMode=MARK_SELECTION_MODE_SPECIAL;
int iFullSizeViewPort=0;
int iUseRhoMinOnSectionAmplitudePlot=0;
long lSeed=666;
int iDoPhaseScale=0;
double dPhaseScale=1.0;
double dStimPeriod=0.0;
double dPreStimTime=0.0;
double dPostStimTime=0.0;
int iPrePostStimTimeRescale=0;
double dCorrelationThreshold=0.2;
int iNCorrRecords=0;
RECORD *pCorrRecords=NULL;
int iCorrelationModeShift=1;
int iNLineSegments=0;
PAIR *pLineSegments=NULL;
int iPlotLineSegments=3; // 0-not, 1-on /xw, 2-on save /ps
int iPlotLineSegmentsStyle=LINE_SEGMENTS_STYLE_SIMPLE;
int iSaveCropBorderAsLineSegments=0;
int iBlendBorder=0;
int iBlendBorderPixel=DEFAULT_BLEND_BORDER;
int iLoadBitmapMode=MODE_LOAD_BITMAP_DEFAULT;
int iFindMarkedInsideSelected=0;
int iMarkCorrelationPoints=1;
int iAdditiveCropChanges=1;
FILTER stFilter;
double dLensRadius=2048.0;
int iLensCorrection=0;
int iFixedGrayScaleLo=0,iFixedGrayScaleHi=0;
float fFixedGrayScaleLo=0.0,fFixedGrayScaleHi=1.0;
int iDoAmplitudeNorm=0;
double dAmplitudeNormCoeff=1.0;
int iDoODAnalysis=0;
int iODAnalysisScheme=0;
int iODAnalysisNBins=DISTRIBUTION_PHI_OD_BINS_N;
char *pcSaveFileSuffix=NULL;

int main(int argc,char *argv[]){
  int i,j;
  char str[256],strDum[8];
  FILE *fp;

  options(argc,argv);
  if(nfiles>MAX_FILE_N){
    printf("Accepts two(X and Y) or one (XY) files\n");
    exit(1);
  }
  if(nfiles<=0){
    printf("Missing file name?\n");
    exit(1);
  }

  if((i=InitializeMaps())){
    printf("Initialization failed %i\n",i);
    exit(1);
  }
  if((i=InitializeBuffers())){
    printf("Buffer initialization failed %i\n",i);
    exit(2);
  }

  //  i=XYdim-1;
  //  printf("%i %i %f %f %f\n",i%Xdim,i/Xdim,maps[0][i],maps[1][i],hypot(maps[0][i],maps[1][i]));

  if(removebias && !radiusbig) RemoveBias(Xdim,Ydim,maps[0],maps[1]);

  AverageMaps(Xdim,Ydim,maps,(float**)0,2,dradius,dradiusbig,(int)DOUBLE);

  start_graphics(display_mode,device,&mainid);

  if(do_pinwheels) AnalyzePW(Xdim,Ydim,maps[0],maps[1]);

  DisplayMap(display_mode,Xdim,Ydim,maps[0],maps[1],dbuffers[0],dbuffers[1]);

  if(do_fields){ 
    start_graphics(field_mode,device,&fieldid);
    DisplayFields(field_mode,Xdim,Ydim,maps[0],maps[1],dbuffers[0],dbuffers[1]);
  }

  if(int_mode) Interactive(Xdim,Ydim,maps,&int_mode);

  if(do_fields) end_graphics(fieldid);
  end_graphics(mainid);

  printf("Save image?(yY/n): ");
  fgets(strDum,4,stdin);
  if(*strDum == 'y' || *strDum == 'Y'){
    save_on=1;
    if(*strDum == 'Y'){
      save_mode|=DISPLAY_MODE_SINGLE;
      if(strDum[1] == '2'){
	save_mode|=DISPLAY_MODE_SINGLE2;
	sprintf(str,"%s%s",filestr,".mapRho.ps/VCPS");
      }
      else{
	if(color_scheme==COLOR_HLS) sprintf(str,"%s%s",filestr,".mapp.ps/VCPS");
	else sprintf(str,"%s%s",filestr,".mapPhi.ps/VCPS");
      }
    }
    else sprintf(str,"%s%s",filestr,".map.ps/VCPS");
    printf("%s\n",str);
    start_graphics(save_mode,devices,&saveid);
    DisplayMap(save_mode,Xdim,Ydim,maps[0],maps[1],dbuffers[0],dbuffers[1]);
    end_graphics(saveid);
    if(strDum[1] == '3'){
      save_mode|=DISPLAY_MODE_SINGLE2;
      sprintf(str,"%s%s",filestr,".mapRho.ps/VCPS");
      printf("%s\n",str);
      start_graphics(save_mode,devices,&saveid);
      DisplayMap(save_mode,Xdim,Ydim,maps[0],maps[1],dbuffers[0],dbuffers[1]);
      end_graphics(saveid);
    }
    save_on=0;
  }

  if(do_dump){
    sprintf(str,"dump.dat");
    if(!(fp=fopen(str,"w"))){
      printf("Cannot open file %s\n",str);
    }
    for(j=Ydim-1;j>=0;j--){
      for(i=0;i<Xdim;i++){
	//	fprintf(fp,"%3i %3i %f\n",i,j,fbufferrho[i+Xdim*j]);
	fprintf(fp,"%f ",fbufferrho[i+Xdim*j]);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }

  return(0);
}
