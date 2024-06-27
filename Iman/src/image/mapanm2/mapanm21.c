/* Map analysis - multiple(2)   */
/* Written by V.Kalatsky        */
/* 06Feb01                      */
/* Modified 5May01              */

#include "mapanm21.h"

Pinwheel *apinwheels=(Pinwheel *)0;
int apinwheels_n=0;

char device[40]="/xw\0";
int mainid,odid;
float pminx,pmaxx,pminy,pmaxy;
float bg,fg,trans[6]={0.0,1.0,0.0,0.0,0.0,1.0};

int nfiles;
char *filenames[MAX_FILE_N];
char *filestr1,*filestr2;
char savefile[256];
char *modestr[]={"0","1","dif","sum","average","true_sum","true_dif","dif_sqrt","skew_dif","fourier_sum","0s","1s"};
int type;
unsigned short Xdim,Ydim;
unsigned long XYdim;
float *fbuffer,*fbufferr,*fbufferx,*fbuffery,*fbufferi,dumf;
int radius=0,radiusbig=0;
double dradius=DEFAULT_RADIUS,dradiusbig=DEFAULT_RADIUS;
int *aindex_x,*aindexbig_x,*aindex_y,*aindexbig_y;
int aindex_n,aindexbig_n;
double **maps,**mapsp;
double *dbufferx,*dbuffery,*dbufferphix,*dbufferphiy;
double *pdBufferRho1=NULL,*pdBufferRho2=NULL;
double *pdDum1,*pdDum2,*pdDum3,*pdDum0;
double harmonic_d;
int do_od=0,do_pinwheels=0,do_contours=0;
int make_phi_positive=1;
double **dbuffers;
int do_interactive=0;
int main_win_size=MAIN_WIN_SIZE;
int npanels=5;
int iInverseMode=MODE_INVERSE_DEFAULT;
int iFlipPhaseSize=FLIP_PHASE_SIZE_DEFAULT;
double dFlipPhaseThreshold=FLIP_PHASE_THRESHOLD_DEFAULT;
int iROIXLeft,iROIXRight,iROIYTop,iROIYBottom;
int iCropROI=0;
char *pcSaveFileSuffix=NULL;
int iDoDominanceHistogram=0;
double dDominanceThreshold=0.0;
int iTakeSQRTOfSum=0;

int main(int argc,char *argv[]){
  double dumd,dumd0,dumd1,dumd2,dumd3,dump,dumr,dumdx,dumdy;
  double phi,phi2,rho;
  double *dmp[2];
  double co,si;
  float dumf;
  int i,panel;
  double adComp[6];

  options(argc,argv);
  if(nfiles>MAX_FILE_N || nfiles<1){
    printf("MAP Accepts up to %i files\n",MAX_FILE_N);
    exit(1);
  }

  if(i=InitializeMaps()){
    printf("MAP Initialization failed %i\n",i);
    exit(1);
  }
  if(i=InitializeBuffers()){
    printf("MAP Buffer initialization failed %i\n",i);
    exit(2);
  }
  AverageMaps();

  FindMapAngle(Xdim,Ydim,maps[0],maps[1],maps[2],maps[3],iInverseMode);

  start_graphics(npanels);

  if(do_pinwheels) AnalyzePW(Xdim,Ydim,maps[0],maps[1]);

  if(iDoDominanceHistogram) PlotDominanceHistogram(Xdim,Ydim,maps[0],maps[1],maps[2],maps[3],0);


  adComp[0]=1.0;
  adComp[1]=-0.11;
  adComp[2]=1.0;
  adComp[3]=-0.22;
  FourierMax(2,adComp,1);
  printf("Maximum %f %f\n",adComp[4],adComp[5]);

  panel=0;

  // 0 First map
  if(panel<npanels){
    DisplayMap(++panel,Xdim,Ydim,maps);
    if(do_interactive & 1<<(panel-1)) Interactive(panel,Xdim,Ydim,maps);
  }

  // 1 Second map
  if(panel<npanels){
    DisplayMap(++panel,Xdim,Ydim,maps+2);
    if(do_interactive & 1<<(panel-1)) Interactive(panel,Xdim,Ydim,maps+2);
  }

  // 2 Difference
  dumdx=dumdy=0.0;
  dumd2=dumd3=0.0;
  for(i=0;i<XYdim;i++){
    dumdx+=(dumd0=dbufferx[i]= DOTPROD(maps[0][i],maps[1][i],maps[2][i],maps[3][i]));
    dumdy+=(dumd1=dbuffery[i]=-VECTPROD(maps[0][i],maps[1][i],maps[2][i],maps[3][i]));
    if((dumd=hypot(dumd0,dumd1))!=0.0){
      dumd2+=dumd0/dumd;
      dumd3+=dumd1/dumd;
    }
    // Division
    /*
    dumd*=dumd;
    dbufferx[i]/=dumd;
    dbuffery[i]/=dumd;
    */
  }
  phi=atan2(dumdy,dumdx)*180.0/M_PI;
  printf("MAP DIF Average P=%f R=%f NAverage P=%f\n",phi,hypot(dumdy,dumdx)/(double)(XYdim),atan2(dumd3,dumd2)*180.0/M_PI);
  if(panel<npanels){
    dmp[0]=dbufferx;
    dmp[1]=dbuffery;
    DisplayMap(++panel,Xdim,Ydim,dmp);
    if(do_interactive & 1<<(panel-1)) Interactive(panel,Xdim,Ydim,dmp);
  }

  // 3 Sum
  dumdx=dumdy=0.0;
  dumd2=dumd3=0.0;
  for(i=0;i<XYdim;i++){
    dumdx+=(dumd0=dbufferx[i]= DOTPROD(maps[0][i],-maps[1][i],maps[2][i],maps[3][i]));
    dumdy+=(dumd1=dbuffery[i]=VECTPROD(maps[0][i],-maps[1][i],maps[2][i],maps[3][i]));
    if((dumd=hypot(dumd0,dumd1))!=0.0){
      dumd2+=dumd0/dumd;
      dumd3+=dumd1/dumd;
    }
  }
  dumdx/=(double)(XYdim);
  dumdy/=(double)(XYdim);
  rho=hypot(dumdy,dumdx);
  phi=atan2(dumdy,dumdx)*180.0/M_PI;
  dumd0=-dumdy/rho;
  dumd1=dumdx/rho;
  phi2=0.0;
  for(i=0;i<XYdim;i++){
    dumd=DOTPROD(dumd0,dumd1,dbufferx[i],dbuffery[i]);
    phi2+=dumd*dumd;
  }
  phi2=atan2(sqrt(phi2/(double)(XYdim)),rho)*180.0/M_PI;
  printf("MAP SUM Average P=%f(Psigma=%f) R=%f NAverage P=%f\n",phi,phi2,rho,atan2(dumd3,dumd2)*180.0/M_PI);

  //Val take sqrt of amplitude and half the phase
  if(iTakeSQRTOfSum){
    for(i=0;i<XYdim;i++){
      phi=atan2(dbuffery[i],dbufferx[i]);
      rho=sqrt(hypot(dbuffery[i],dbufferx[i]));
      //dumd0=atan2(maps[1][i],maps[0][i]);
      //if(dumd0<0.0) dumd0+=2*M_PI;
      //dumd1=atan2(maps[3][i],maps[2][i]);
      //if(dumd1<0.0) dumd1+=2*M_PI;
      //phi=0.5*(dumd0+dumd1);
      //rho=sqrt(hypot(maps[1][i],maps[0][i])*hypot(maps[3][i],maps[2][i]));
      dbufferx[i]=rho*cos(phi);
      dbuffery[i]=rho*sin(phi);
    }
  }
  if(panel<npanels){
    dmp[0]=dbufferx;
    dmp[1]=dbuffery;
    DisplayMap(++panel,Xdim,Ydim,dmp);
    if(do_interactive & 1<<(panel-1)) Interactive(panel,Xdim,Ydim,dmp);
  }

  // 4 Average, first corrected, second corrected
  for(i=0;i<XYdim;i++){
    phi=0.5*atan2(dbuffery[i],dbufferx[i]);
    if(make_phi_positive && phi<0.0) phi+=M_PI;
    co=cos(phi);
    si=sin(phi);
    dbuffers[0][i]= -DOTPROD(maps[0][i],maps[1][i],co,si);
    dbuffers[1][i]=VECTPROD(maps[0][i],maps[1][i],co,si);
    dbuffers[2][i]= -DOTPROD(maps[2][i],maps[3][i],co,si);
    dbuffers[3][i]=-VECTPROD(maps[2][i],maps[3][i],co,si);
    dbufferphix[i]=0.5*(dbuffers[0][i]+dbuffers[2][i]);
    dbufferphiy[i]=0.5*(dbuffers[1][i]+dbuffers[3][i]);
  }
  if(panel<npanels){
    dmp[0]=dbufferphix;
    dmp[1]=dbufferphiy;
    DisplayMap(++panel,Xdim,Ydim,dmp);
    if(do_interactive & 1<<(panel-1)) Interactive(panel,Xdim,Ydim,dmp);
  }

  // 5  True sum
  dumdx=dumdy=0.0;
  dumd2=dumd3=0.0;
  for(i=0;i<XYdim;i++){
    dumdx+=(dumd0=dbufferx[i]= 0.5*(maps[0][i]+maps[2][i]));
    dumdy+=(dumd1=dbuffery[i]= 0.5*(maps[1][i]+maps[3][i]));
    if((dumd=hypot(dumd0,dumd1))!=0.0){
      dumd2+=dumd0/dumd;
      dumd3+=dumd1/dumd;
    }
  }
  phi=atan2(dumdy,dumdx)*180.0/M_PI;
  printf("MAP TRUE SUM Average P=%f R=%f NAverage P=%f\n",phi,hypot(dumdy,dumdx)/(double)(XYdim),atan2(dumd3,dumd2)*180.0/M_PI);
  if(panel<npanels){
    dmp[0]=dbufferx;
    dmp[1]=dbuffery;
    DisplayMap(++panel,Xdim,Ydim,dmp);
    if(do_interactive & 1<<(panel-1)) Interactive(panel,Xdim,Ydim,dmp);
  }

  // 6 True dif
  dumdx=dumdy=0.0;
  dumd2=dumd3=0.0;
  for(i=0;i<XYdim;i++){
    dumdx+=(dumd0=dbufferx[i]= (maps[0][i]-maps[2][i]));
    dumdy+=(dumd1=dbuffery[i]= (maps[1][i]-maps[3][i]));
    if((dumd=hypot(dumd0,dumd1))!=0.0){
      dumd2+=dumd0/dumd;
      dumd3+=dumd1/dumd;
    }
  }
  phi=atan2(dumdy,dumdx)*180.0/M_PI;
  printf("MAP TRUE DIF Average P=%f R=%f NAverage P=%f\n",phi,hypot(dumdy,dumdx)/(double)(XYdim),atan2(dumd3,dumd2)*180.0/M_PI);
  if(panel<npanels){
    dmp[0]=dbufferx;
    dmp[1]=dbuffery;
    DisplayMap(++panel,Xdim,Ydim,dmp);
    if(do_interactive & 1<<(panel-1)) Interactive(panel,Xdim,Ydim,dmp);
  }

  // 7 dif_sqrt
  dumd2=dumd3=0.0;
  for(i=0;i<XYdim;i++){
    dumd0=dbufferx[i]= DOTPROD(maps[0][i],maps[1][i],maps[2][i],maps[3][i]);
    dumd1=dbuffery[i]=-VECTPROD(maps[0][i],maps[1][i],maps[2][i],maps[3][i]);
    if((dumd=sqrt(hypot(dumd0,dumd1)))!=0.0){
      dbufferx[i]/=dumd;
      dbuffery[i]/=dumd;
    }
  }
  printf("MAP DIF SQRT\n");
  if(panel<npanels){
    dmp[0]=dbufferx;
    dmp[1]=dbuffery;
    DisplayMap(++panel,Xdim,Ydim,dmp);
    if(do_interactive & 1<<(panel-1)) Interactive(panel,Xdim,Ydim,dmp);
  }

  
  // 8 Skew dif 
  pdDum0=maps[0];
  pdDum1=maps[1];
  pdDum2=maps[2];
  pdDum3=maps[3];
  dumd=1.0/sqrt(2.0);
  for(i=0;i<XYdim;i++){
    dumd2=hypot(pdDum0[i],pdDum1[i]);
    dumd3=hypot(pdDum2[i],pdDum3[i]);
    // 1st power
    // Rotated by 45deg
    //pdBufferRho1[i]= dumd*(dumd3+dumd2);
    //pdBufferRho2[i]= dumd*(dumd3-dumd2);
    pdBufferRho1[i]= dumd2;
    pdBufferRho2[i]= dumd3;
    // Higher powers
    //    dumr=hypot(dumd3,dumd2);
    //    dump=atan2(dumd3,dumd2);    
    // 2nd power
    //    pdBufferRho1[i]= dumr*sin(2.0*dump);
    //    pdBufferRho2[i]= -dumr*cos(2.0*dump);
    // 4th power
    //    pdBufferRho1[i]= -dumr*cos(4.0*dump);
    //    pdBufferRho2[i]= -dumr*sin(4.0*dump);
  }
  printf("MAP SKEW DIF \n");
  if(panel<npanels){
    dmp[0]=pdBufferRho1;
    dmp[1]=pdBufferRho2;
    DisplayMap(++panel,Xdim,Ydim,dmp);
    if(do_interactive & 1<<(panel-1)) Interactive(panel,Xdim,Ydim,dmp);
  }

  // 9 Fourier sum 
  pdDum0=maps[0];
  pdDum1=maps[1];
  pdDum2=maps[2];
  pdDum3=maps[3];
  for(i=0;i<XYdim;i++){
    adComp[0]=hypot(pdDum0[i],pdDum1[i]);
    adComp[1]=atan2(pdDum1[i],pdDum0[i]);
    adComp[2]=hypot(pdDum2[i],pdDum3[i]);
    adComp[3]=atan2(pdDum3[i],pdDum2[i]);
    FourierMax(2,adComp,1);
    pdBufferRho1[i]=adComp[4]*cos(adComp[5]);
    pdBufferRho2[i]=adComp[4]*sin(adComp[5]);
  }
  printf("MAP FOURIER SUM \n");
  if(panel<npanels){
    dmp[0]=pdBufferRho1;
    dmp[1]=pdBufferRho2;
    DisplayMap(++panel,Xdim,Ydim,dmp);
    if(do_interactive & 1<<(panel-1)) Interactive(panel,Xdim,Ydim,dmp);
  }

  // 10 First + shift
  if(panel<npanels){
    DisplayMap(++panel,Xdim,Ydim,dbuffers);
    if(do_interactive & 1<<(panel-1)) Interactive(panel,Xdim,Ydim,dbuffers);
  }

  // 11 Second + shift
  if(panel<npanels){
    DisplayMap(++panel,Xdim,Ydim,dbuffers+2);
    if(do_interactive & 1<<(panel-1)) Interactive(panel,Xdim,Ydim,dbuffers+2);
  }

  if(do_od) DisplayOd();

  end_graphics();

  return(0);
}
