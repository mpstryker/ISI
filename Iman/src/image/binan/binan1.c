/* Bin analysis                 */
/* Written by V.Kalatsky        */
/* 20Jul2001                    */

#include "binan1.h"

char device[40]="/xw\0",device_ask[40]="?\0";
int mainid;
float pminx,pmaxx,pminy,pmaxy;
float bg,fg,trans[6]={-1.0,1.0,0.0,-1.0,0.0,1.0};
int ibac=0,ifor=0;
float wedge_width=0.4,wedge_pos=0.0;

int nfiles=0;
char **filenames;
char *filestr;
int type;
unsigned short Xdim,Ydim;
unsigned long XYdim;
float *fbufferi;
float *fbuffer;
int radius,radiusbig,cradius;
double dradius=DEFAULT_RADIUS,dradiusbig=DEFAULT_RADIUSBIG,dcradius=DEFAULT_CHOP_RADIUS;
double **maps,**maps_orig;
double *pdMapX=NULL,*pdMapY=NULL;
float **fmaps,**fmaps_orig;
float *bin_max,*bin_min;
int do_marking=0;
int do_verbose=0;
int do_wedge=1;
int do_animation=0,do_fixed_background=0;
int removebias=0;
int harmonic=1;
double harmonic_d=1.0;
int int_mode=INT_MODE_NON;
int display_mode=DEFAULT_DISPLAY_MODE;
unsigned long *stencil=(unsigned long *)0,stenciln=0;
int smaxx,sminx,smaxy,sminy;
float window_size=WINDOW_SIZE;
int spatial_binning=1;
int nbins=0;
int panelnx=0,panelny=0;
int shift_max_to_zero=0;
int shift_average_to_zero=0;

pthread_t animate_thread;
pthread_attr_t animate_thread_attr;
AnimateBinsThreadArgs animatebinsthreadargs;
int thread_return_value=0;
int stop_animation=0;
unsigned long animation_delay=DEFAULT_ANIMATION_DELAY;
int animate_forward=1;
int iNormalizeBinning=1;
float fBinSizeSec=DEFAULT_BIN_SIZE_SEC;
int iPlotFComponents=0;
int iNFComponents=2;
int iMinFComponent=1;
int iGenerateMaps=1;
int iMakeMaps=0;

int main(int argc,char *argv[]){
  int i;
  double *ppd[2];

  options(argc,argv);
  /*
  if(nfiles!=MAX_FILE_N){
    printf("Accepts %i file(s)\n",MAX_FILE_N);
    exit(1);
  }
  */

  if(nfiles<1){
    printf("requires at least 1 file\n");
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

  if(do_animation) AnimateBins(display_mode,Xdim,Ydim,nbins,fmaps,animation_delay);

  if(iMakeMaps){
    MakeMaps(nbins,maps_orig,Xdim,Ydim,pdMapX,pdMapY);
    ppd[0]=pdMapX;
    ppd[1]=pdMapY;
    DisplayMap(Xdim,Ydim,ppd);
  }

  start_graphics(display_mode,device,&mainid);

  DisplayBins(display_mode,Xdim,Ydim,nbins,fmaps);

  if(int_mode!=INT_MODE_NON) Interactive(Xdim,Ydim,fmaps,&int_mode);

  end_graphics(mainid);

  return(0);
}
