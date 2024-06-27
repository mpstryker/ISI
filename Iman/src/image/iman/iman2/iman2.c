/* Image analysis               */
/* Written by V.Kalatsky        */
/* 02Feb01                      */

#include "iman1.h"

int imagetype;
char device[40]="/xw\0";
int mainid=0;
float pminx,pmaxx,pminy,pmaxy;
float *fbuffer=(float*)0,*fbuffer_removed=(float*)0;
unsigned short Xdim,Ydim;
unsigned long XYdim;
float bg,fg,trans[6]={-1.0,1.0,0.0,-1.0,0.0,1.0};
int ibac,ifor;
int n_panels=1,panel=0;
int scheme=SCHEME_DIFF2,difference=DIFF_RELATIVE;
int initframe=-1;
int finalframe=-1;
int increment=1,increment1=1;
int single_multiple=DO_SINGLE;
unsigned short *dumbuffer;
double *mapx=NULL,*mapy=NULL,*map0=NULL;
double *timeaverage_buffer=NULL;
double timeaverage_cycles=0.0;
double timeaverage_blocks=0.0;
double timeaverage_buffer_frames_in=-1.0;
int timeaverage_cycles_int;
float *cocktailaverage,camin,camax;
double phi0,dphi;
int radius,radiusbig;
double dradius=DRADIUS,dradiusbig=DRADIUSBIG;
int timeradius=TIMERADIUS,sacrifice_boundaries=1;
int *aindex_x,*aindexbig_x,*aindex_y,*aindexbig_y;
int aindex_n,aindexbig_n;
float bgg=0.0,fgg=20.0;
float bga=0.0,fga=0.0;
int image_count=0;
int average=AVERAGE_NOT;
double harmonic=HARMONIC;
unsigned int n_cycles=NCYCLES;
int do_contours=0,do_verbose=0,do_statistics=0,do_trace=0,do_synch=0,do_precise=0,do_intervals=0;
int find_best_cycle=0,do_remove_linear=0,do_remove_curve=0;
int display_removed=1,do_cocktailaverage=0,do_timeaveraging=0;
int print_fileheader=1;
double *synch_phi=(double *)0,*synch_cos=(double *)0,*synch_sin=(double *)0;
double synch_phase_shift=0.0,synch_omega=-1.0;
double frame_omega=-1.0;
int cut_offl=CUT_OFF_L,cut_offu= CUT_OFF_U;
int nfiles=0;
char *filenames[2];
double period=-1.0;
double cyclesupper,cycleslower,cyclesexact=NCYCLES;
float *synch_trace=(float *)0;
int additionalframes=0;
double stimulus_period=-1.0;
double *remove_sy=(double *)0,*remove_sxy=(double *)0;
double **remove_curve=(double **)0;
double *remove_cos=(double *)0,*remove_sin=(double *)0,*remove_dum=(double *)0;
int poly_fitn=0;
double a_cos,a_sin,b_cos,b_sin;
int nframes_good=0;
int replace_scheme=/*REPLACE_OLDEST*/REPLACE_REQUESTED/*REPLACE_RANDOM*/;
unsigned long memory_size=0;
unsigned int max_frames_in_cue=MAX_FRAMES_IN_CUE;
int compress=COMPRESS_NOT;
int initframe_ta=0,finalframe_ta=0;
int n_bins=0;
double **frame_bins=(double**)0;
int *frames_in_bin=(int*)0;
int displayimageN=100,displayimagesN=1;
double cycles_per_frame=0.0;
int save_in_files=SAVE_XY_DEFAULT;
int prompt_to_save_images=0;                                                    
int statistics_num_bins=DEFAULT_STATISTICS_NUM_BINS;
unsigned long ulSynchType=SYNCHRONIZATION_NONE;
int iSynchChannel=-1;
int iPlayNonStop=1;
int iDoDivisionByAverage=0;
int iDivisionByAverageCount=0;
double dDivisionByAverageFactor=10000.0;
double dBinDivisionByAverageFactor=1.0;
unsigned long ulBinX,ulBinY,ulBinT;
int iKeepOriginalFiles=DEFAULT_KEEP_ORIGINAL_FILES;
int iRemoveTrailingZsFromDecompressedFileName=1;
float fContrast=1.0;
int iDrawAxesAndLabelOnGreen=1;
float fSizeOfGreenSave=GREEN_IMAGE_SIZE_SAVE;
float fSizeOfTraceInteractiveWindow=TRACE_INTERACTIVE_WINDOW_SIZEX;
float fSizeOfTracePlotWindow=TRACE_PLOT_WINDOW_SIZEX;
int iUseExternalFourierThreshold=0;
float fExternalFourierThreshold=0.0;

int (*pfunFrameRule[FRAME_RULE_N])(TRAIN *tr,int fr,CAR *pcar)={FrameRule0,FrameRule1};
int frame_rule=FRAME_RULE_DEFAULT;

int main(int argc,char *argv[]){
  TRAIN train;
  long duml;
  int i,j,nfr,stop_frame=1;
  double dumd;
  double s,sx,sxx,delta;
  double sy_cos,sy_sin,sxy_cos,sxy_sin;
  CAR *pcar;
  char str[5];
  int iExitCode=0;

  memset(&train,0,sizeof(TRAIN));

  options(argc,argv);

  if(nfiles != MAX_FILE_N){
    printf("MAN Usage: iman2 <filename>(%i)\n",nfiles);
    exit(1);
  }

  if(InitializeTrain(&train,filenames[0])){
    printf("MAN Train initialization failed\n");
    exit(3);
  }
  //  printf("*** I=%i F=%i\n",iframe,fframe);
  if(InitializeBuffers(&train)){
    printf("MAN Buffer initialization failed\n");
    iExitCode=4;
    goto bailout2;
  }

  if(print_fileheader){ 
    for(i=0;i<train.n_files;i++) PrintFileheader(&train,i);
  }

  if(compress!=COMPRESS_NOT){
    iExitCode=Compressor(&train,compress);
    goto bailout2;
  }

  printf("MAN Xdim Ydim: %i %i\n",Xdim,Ydim);
  printf("MAN Frames %i Ini frame %i Cycles %i\n",train.max_n,initframe,n_cycles);

  if(!ulSynchType){
    if(train.iGlobalExperiment==EXPERIMENT_CONTINUOUS_MODE_EPISODIC){
      if(train.iGlobalVersion==VERSION_CHUNK){
	if(EPSTParser(&train,initframe,finalframe,&initframe,&finalframe)){
	  printf("MAN EPSTParser failed\n");
	  iExitCode=11;
	  goto bailout2;
	}
      }
      else{
	n_bins=((FILEHEADER*)(train.files[0].pFileHeader))->sumNconds;
	if(do_timeaveraging){
	  sacrifice_boundaries=0;
	  initframe_ta=initframe;
	  finalframe_ta=finalframe;
	  timeradius=(int)rint(timeaverage_blocks*(((FILEHEADER*)(train.files[0].pFileHeader))->sumNframeiti+((FILEHEADER*)(train.files[0].pFileHeader))->sumNframestim)/((FILEHEADER*)(train.files[0].pFileHeader))->sumTempBinning);
	  if(InitializeTimeAveraging(&train)){ 
	    printf("BIN Time averaging initialization failed\n");
	    iExitCode=13;
	    goto bailout2;
	  }
	}
      }
    }
    else{
      printf("MAN No synchronization chosen\n");
      printf("MAN Running with default meaningless period = %i Fr/cycle\n",DEFAULT_PERIOD);
      period=DEFAULT_PERIOD;
    }
  }

  pminx=-0.5;
  pmaxy=-0.5;
  pmaxx=(float)Xdim-0.5;
  pminy=(float)Ydim-0.5;

  if(period>0.0){ 
    dphi=harmonic*2.0*M_PI/period;
    finalframe=(int)ceil(period*(double)n_cycles)+initframe-1;
    if(!(pcar=AddFrame(&train,initframe,AVERAGE_NOT))){
      printf("MAN AddFrame return 0 fr=%i\n",initframe);
    }
    phi0=harmonic*(2.0*M_PI*(double)(((FRAMEHEADER*)pcar->pFrameHeader)->synch_in>>2)/DEFAULT_SYNCH_MAX);
    printf("MAN Frames range %i-%i phi0=%f dphi=%f\n",initframe,finalframe,phi0,dphi);
  }

  if(do_cocktailaverage){  
    CocktailAverage(&train);
    DisplayCocktailAverage(&train); 
    getchar();
  }

  if(finalframe >= train.max_n || finalframe<=0){ 
    printf("MAN Fixing finalframe from %i to %i\n",finalframe,train.max_n-1);
    finalframe=train.max_n-1;
  }

  start_graphics();

  if(do_synch){
    if(Synchronization(&train,initframe,finalframe,&initframe,&finalframe,average)){
      printf("MAN Synchronization failed\n");
      if(do_intervals%2){
	if((i=InterframeTime(&train,initframe,finalframe,&initframe,&finalframe,average))){
	  printf("MAN Frame time averaging failed (error=%i) Enter 'y' to continue: ",i);
	  fgets(str,4,stdin);
	  if(*str != 'y'){
	    iExitCode=2;
	    goto bailout1;
	  }	  
	}
	else{
	  //	  omega=frame_omega;
	  if(Inverse(finalframe-initframe+1,frame_omega,0.0,(int)(n_cycles*harmonic))){
	    printf("MAN Inversion failed (frame)\n");
	    iExitCode=1;
	    goto bailout1;
	  }
	  else{
	    printf("MAN Inverted frame basis\n");
	  }
	}
	do_intervals=0;
      }
    }
    else{
      //      omega=synch_omega;
      if(Inverse(finalframe-initframe+1,synch_omega,synch_phase_shift,(int)(n_cycles*harmonic))){
	printf("MAN Inversion failed (synch)\n");
	iExitCode=2;
	goto bailout1;
      }
      else{
	printf("MAN Inverted synch basis\n");
      }
    }
  }
  if(do_intervals){
    if((i=InterframeTime(&train,initframe,finalframe,&initframe,&finalframe,average))){
      printf("MAN Frame time averaging failed (error=%i) Enter 'y' to continue: ",i);
      fgets(str,4,stdin);
      if(*str != 'y'){
	iExitCode=2;
	goto bailout1;
      }	  
    }
    else{
      if(!synch_cos){
	if(Inverse(finalframe-initframe+1,frame_omega,0.0,(int)(n_cycles*harmonic))){
	  printf("MAN Inversion failed (frame)\n");
	  //	  getchar();
	  iExitCode=3;
	  goto bailout1;
	}
	else{
	  printf("MAN Inverted frame basis\n");
	}
      }
    }
  }

  printf("MAN Harmonic %.2f (%i) Cycles of %i\n",harmonic,(int)harmonic,(int)(n_cycles*harmonic));

  if(do_statistics) Statistics(&train);

  if(do_trace) TraceAnalysis(&train,initframe,finalframe,average);

  //  printf("do_timeaveraging=%i\n",do_timeaveraging);
  //  printf("do_remove_curve=%i\n",do_remove_curve);
  //  printf("do_remove_linear=%i\n",do_remove_linear);

  if(single_multiple == DO_SINGLE){
    if(do_timeaveraging){
      // Accumulate time averaged
      if(AccumulateTimeAverage(&train,average,&timeaverage_buffer_frames_in,&timeaverage_buffer)){
	printf("MAN Cannot accumulate time average\n");
	iExitCode=5;
	goto bailout1;
      }
    }
    if(n_bins){
      if(BinFrames(&train,average)){
	printf("MAN Frame binning failed\n");
	iExitCode=1;
      }
      iExitCode=0;
      goto bailout1;
    }
    
    printf("MAN ");
    for(nfr=0,i=initframe;i<=finalframe;i+=increment1,nfr++){
      ConvertImage(&train,i,average);
      if(!(nfr%displayimageN)){
	LoadFrame(&train,i,fbuffer,average);
	DisplayImage(&train,i,fbuffer);
      }
      if(!(nfr%100)){ 
	if(!(nfr%1000)){ 
	  printf("|"); fflush(stdout);
	}
	else printf(".");
	fflush(stdout);
      }
    }
    printf("\n");
  }
  if(single_multiple == DO_MULTIPLE){
    for(nfr=0,i=initframe,j=0;i<=finalframe;i+=increment1,j=(j+1)%stop_frame,nfr++){
      ConvertImages(&train,(i+increment)%train.max_n,i%train.max_n,average);
      if(!(nfr%displayimagesN)){
	DisplayImages(&train,(i+increment)%train.max_n,i%train.max_n,fbuffer);
      }
      if(j==stop_frame-1 && !iPlayNonStop) getchar();
    }
  }
  
  if(do_remove_curve){
    for(i=0;i<XYdim;i++){
      for(sy_cos=sy_sin=0.0,j=0;j<poly_fitn;j++){
	sy_cos+=remove_curve[j][i]*remove_cos[j];
	sy_sin+=remove_curve[j][i]*remove_sin[j];
      }
      mapx[i]-=sy_cos;
      mapy[i]-=sy_sin;
    
      if(display_removed){ 
	fbuffer[i]=sy_cos;
	fbuffer_removed[i]=sy_sin;
      }      
    }
    if(display_removed) DisplayVectorField2(Xdim,Ydim,(double*)0,(double*)0,fbuffer,fbuffer_removed,FLOAT);
  }
  else{
    if(do_remove_linear){
      s=(double)nfr;
      sx=0.5*s*(double)(nfr-1);
      sxx=sx*(double)(2*nfr-1)/3.0;
      delta=sx*sx-s*sxx;
      
      sy_cos=(-sxx*a_cos+sx*b_cos)/delta;
      sxy_cos=(sx*a_cos-s*b_cos)/delta;
      sy_sin=(-sxx*a_sin+sx*b_sin)/delta;
      sxy_sin=(sx*a_sin-s*b_sin)/delta;
      
      for(i=0;i<XYdim;i++){
	mapx[i]-=sy_cos*remove_sy[i]+sxy_cos*remove_sxy[i];
	mapy[i]-=sy_sin*remove_sy[i]+sxy_sin*remove_sxy[i];
      }
      if(display_removed){ 
	for(i=0;i<XYdim;i++){
	  dumd=sy_cos*remove_sy[i]+sxy_cos*remove_sxy[i];
	  remove_sy[i]=sy_sin*remove_sy[i]+sxy_sin*remove_sxy[i];
	  remove_sxy[i]=dumd;
	}      
	DisplayVectorField2(Xdim,Ydim,remove_sxy,remove_sy,(float *)0,(float *)0,DOUBLE);
      }
    }
  }
  
  if(!synch_sin){
    dumd=2.0/(double)(nfr);
    for(i=0;i<XYdim;i++){
      mapx[i]*=dumd;
      mapy[i]*=dumd;
    }
  }
  
  if(iDoDivisionByAverage){
    if(iDivisionByAverageCount!=nfr) 
      printf("MAN Division Frame Count missmatch %i %i\n",iDivisionByAverageCount,nfr);
    dumd=1.0/(dDivisionByAverageFactor*(double)iDivisionByAverageCount);
    for(i=0;i<XYdim;i++){
      s=(map0[i]*=dumd);
      mapx[i]/=s;
      mapy[i]/=s;
    }
  }
  
  DisplayMap(&train,mapx,mapy);

 bailout1:
  cpgslct(mainid);
  cpgask(0);
  end_graphics();

 bailout2:

  ReleaseTrain(&train,iExitCode);
  return(iExitCode);
}
