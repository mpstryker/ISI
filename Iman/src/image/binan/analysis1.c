/* Bin analysis                 */
/* Written by V.Kalatsky        */
/* 20Jul2001                    */

#include "binan1.h"

void Stencil(int xd,int yd,int n,double *dx,double *dy,unsigned long *sn,unsigned long **s){
  static unsigned long *pul=(unsigned long *)0;
  unsigned long k;
  int j;
  unsigned long i,xy;
  double x,y;
  double angle;

  xy=(unsigned long)xd*(unsigned long)yd;
  pul=(unsigned long *)realloc(*s,xy*sizeof(unsigned long));

  if(n==0){
    *sn=xy;
    for(i=0;i<stenciln;i++) pul[i]=i;
    *s=pul;
    return;
  }

  smaxx=smaxy=0;
  sminx=xd-1;
  sminy=yd-1;
  for(k=i=0;i<xy;i++){
    x=(double)(i%xd);
    y=(double)(i/xd);
    for(j=0,angle=0.0;j<n;j++){
      angle+=atan2(VECTPROD(dx[j]-x,dy[j]-y,dx[j+1]-x,dy[j+1]-y),DOTPROD(dx[j]-x,dy[j]-y,dx[j+1]-x,dy[j+1]-y));
    }
    if(fabs(angle)>M_PI){
      pul[k]=i;
      k++;
      j=i%xd;
      if(j>smaxx) smaxx=j;
      if(j<sminx) sminx=j;
      j=i/xd;
      if(j>smaxy) smaxy=j;
      if(j<sminy) sminy=j;
    }
  }
  *sn=k;
  pul=(unsigned long *)realloc(pul,k*sizeof(unsigned long));
  *s=pul;
}

void Double2FloatShort(unsigned long nn,double *pd,float *pf){
  unsigned long i;
  for(i=0;i<nn;i++) pf[i]=(float)pd[i];
}

void Double2FloatLong(unsigned long nn,double *pd,float *pf,float *min,float *max){
  unsigned long i;
  float dumf;

  *max=0.0;
  *min=(float)(1<<16);
  for(i=0;i<nn;i++){
    dumf=pf[i]=(float)pd[i];
    if(i>=Xdim){ // Do not take into account first screwed row of each frame
      if(*max<dumf){
        *max=dumf;
        ifor=i;
      }
      if(*min>dumf){
        *min=dumf;
        ibac=i;
      }
    }
  }  return;  
}

void FindExtrema(float *fb,float *fmin,float *fmax,unsigned long *imin,unsigned long *imax){
  unsigned long i,j;
  *fmax=-(*fmin=BIG_FLOAT);
  for(j=0;j<stenciln;j++){
    i=stencil[j];
    if(fb[i]>*fmax){ 
      *fmax=fb[i];
      *imax=i;
    }
    if(fb[i]<*fmin){ 
      *fmin=fb[i];
      *imin=i;
    }
  }
}

void FindExtremaEx(unsigned long n,unsigned long nn,float *fb,float *fmin,float *fmax,unsigned long *imin,unsigned long *imax){
  unsigned long i;
  float dumf=0;

  *fmax=-(*fmin=BIG_FLOAT);
  for(i=0;i<n;i++){
    if(i>=nn){ // Do not take into account first screwed row of each frame
      dumf=fb[i];
      if(dumf>*fmax){ 
	*fmax=dumf;
	*imax=i;
      }
      if(dumf<*fmin){
	*fmin=dumf;
	*imin=i;
      }
    }
  }
}

void FindMax(unsigned long n,float *fb,float *fmax,unsigned long *imax){
  unsigned long i;
  for(*fmax=-BIG_FLOAT,i=0;i<n;i++){
    if(fb[i]>*fmax){ 
      *fmax=fb[i];
      *imax=i;
    }
  }
}

void FindMin(unsigned long n,float *fb,float *fmin,unsigned long *imin){
  unsigned long i;
  for(*fmin=BIG_FLOAT,i=0;i<n;i++){
    if(fb[i]<*fmin){ 
      *fmin=fb[i];
      *imin=i;
    }
  }
}

int FindFTComponents(int iN,float *pfT,int iC,float *pfRe,float *pfIm){
  double re,im;
  double co,si,coN,siN;
  double dumd;
  int i;
  if(iC>iN/2){
    printf("ANA Bad Fourier component requested\n");
    return(1);
  }
  if(iC==0){
    re=im=0.0;
    for(i=0;i<iN/2;i++){
      re+=(double)pfT[2*i]+(double)pfT[2*i+1];
      im+=(double)pfT[2*i]-(double)pfT[2*i+1];
    }
    *pfRe=re/iN;
    *pfIm=2.0*im/iN;
  }
  else{
    co=cos(2.0*M_PI*(double)iC/(double)iN);
    si=sin(2.0*M_PI*(double)iC/(double)iN);
    coN=1.0; 
    siN=0.0; 
    re=(double)pfT[0];
    im=0.0;
    for(i=1;i<iN;i++){
      dumd=coN*co-siN*si;
      siN=siN*co+coN*si;
      coN=dumd;
      dumd=(double)pfT[i];
      re+=dumd*coN;
      im+=dumd*siN;
    }
    *pfRe=2.0*re/iN;
    *pfIm=2.0*im/iN;
  }
  return(0);
}

int FindFTComponents_DOUBLE(int iN,double *pdT,int iC,double *pdRe,double *pdIm){
  double re,im;
  double co,si,coN,siN;
  double dumd;
  int i;
  if(iC>iN/2){
    printf("ANA Bad Fourier component requested\n");
    return(1);
  }
  if(iC==0){
    re=im=0.0;
    for(i=0;i<iN/2;i++){
      re+=pdT[2*i]+pdT[2*i+1];
      im+=pdT[2*i]-pdT[2*i+1];
    }
    *pdRe=re/iN;
    *pdIm=2.0*im/iN;
  }
  else{
    co=cos(2.0*M_PI*(double)iC/(double)iN);
    si=sin(2.0*M_PI*(double)iC/(double)iN);
    coN=1.0; 
    siN=0.0; 
    re=pdT[0];
    im=0.0;
    for(i=1;i<iN;i++){
      dumd=coN*co-siN*si;
      siN=siN*co+coN*si;
      coN=dumd;
      dumd=pdT[i];
      re+=dumd*coN;
      im+=dumd*siN;
    }
    *pdRe=2.0*re/iN;
    *pdIm=2.0*im/iN;
  }
  return(0);
}

#define RESAMPLING 4
#define N_ITERATIONS 6
#define BIG_DOUBLE 1.0e10

// iSign=0 - min, else - max
double FindExtremaFComponents(int iN,float *pfFC,int iNC,int iCMin,int iSign){
  double dDum=0,dEx=0,dIncr;
  double dX,dXL,dXR;
  int i,iNS,iDum=0;
  iNS=iN*RESAMPLING;
  dIncr=1.0/RESAMPLING;
  if(iSign){
    dEx=-BIG_DOUBLE;
    for(i=0;i<iNS;i++){
      if(dEx<(dDum=FourierSum(iN,pfFC,iNC,iCMin,dIncr*i))){
	dEx=dDum;
	iDum=i;
      }
    }
    dX=iDum*dIncr;
    dXL=dX-dIncr;
    dXR=dX+dIncr;
    for(i=0;i<N_ITERATIONS;i++){
      dX=0.5*(dXL+dXR);
      if(FourierSumD(iN,pfFC,iNC,iCMin,dX)>0.0) dXL=dX;
      else dXR=dX;
    }
  }
  else{
    dEx=BIG_DOUBLE;
    for(i=0;i<iNS;i++){
      if(dEx>(dDum=FourierSum(iN,pfFC,iNC,iCMin,dIncr*i))){
	dEx=dDum;
	iDum=i;
      }
    }
    dX=iDum*dIncr;
    dXL=dX-dIncr;
    dXR=dX+dIncr;
    for(i=0;i<N_ITERATIONS;i++){
      dX=0.5*(dXL+dXR);
      if(FourierSumD(iN,pfFC,iNC,iCMin,dX)<0.0) dXL=dX;
      else dXR=dX;
    }
  }
  return(0.5*(dXR+dXL));
}

double FourierSum(int iN,float *pfFC,int iNC,int iCMin,float fT){
  double dOmega;
  int i,iMin;
  double dDum=0;
  if(iCMin==0) dDum=pfFC[0];
  iMin=iCMin>0 ? iCMin : 1;
  for(i=iMin;i<iMin+iNC;i++){
    dOmega=fT*2.0*M_PI*(double)i/(double)iN;
    dDum+=pfFC[2*i]*cos(dOmega)+pfFC[2*i+1]*sin(dOmega);
  }
  return(dDum);
}

double FourierSumD(int iN,float *pfFC,int iNC,int iCMin,float fT){
  double dOmega;
  int i,iMin;
  double dDum=0;
  iMin=iCMin>0 ? iCMin : 1;
  for(i=iMin;i<iMin+iNC;i++){
    dOmega=2.0*M_PI*(double)i/(double)iN;
    dDum+=dOmega*(-pfFC[2*i]*sin(fT*dOmega)+pfFC[2*i+1]*cos(fT*dOmega));
  }
  return(dDum);
}

double FindExtremaFComponents_DOUBLE(int iN,double *pdFC,int iNC,int iCMin,int iSign){
  double dDum=0,dEx=0,dIncr;
  double dX,dXL,dXR;
  int i,iNS,iDum=0;
  iNS=iN*RESAMPLING;
  dIncr=1.0/RESAMPLING;
  if(iSign){
    dEx=-BIG_DOUBLE;
    for(i=0;i<iNS;i++){
      if(dEx<(dDum=FourierSum_DOUBLE(iN,pdFC,iNC,iCMin,dIncr*i))){
	dEx=dDum;
	iDum=i;
      }
    }
    dX=iDum*dIncr;
    dXL=dX-dIncr;
    dXR=dX+dIncr;
    for(i=0;i<N_ITERATIONS;i++){
      dX=0.5*(dXL+dXR);
      if(FourierSumD_DOUBLE(iN,pdFC,iNC,iCMin,dX)>0.0) dXL=dX;
      else dXR=dX;
    }
  }
  else{
    dEx=BIG_DOUBLE;
    for(i=0;i<iNS;i++){
      if(dEx>(dDum=FourierSum_DOUBLE(iN,pdFC,iNC,iCMin,dIncr*i))){
	dEx=dDum;
	iDum=i;
      }
    }
    dX=iDum*dIncr;
    dXL=dX-dIncr;
    dXR=dX+dIncr;
    for(i=0;i<N_ITERATIONS;i++){
      dX=0.5*(dXL+dXR);
      if(FourierSumD_DOUBLE(iN,pdFC,iNC,iCMin,dX)<0.0) dXL=dX;
      else dXR=dX;
    }
  }
  return(0.5*(dXR+dXL));
}

double FourierSum_DOUBLE(int iN,double *pdFC,int iNC,int iCMin,double dT){
  double dOmega;
  int i,iMin;
  double dDum=0;
  if(iCMin==0) dDum=pdFC[0];
  iMin=iCMin>0 ? iCMin : 1;
  for(i=iMin;i<iMin+iNC;i++){
    dOmega=dT*2.0*M_PI*(double)i/(double)iN;
    dDum+=pdFC[2*i]*cos(dOmega)+pdFC[2*i+1]*sin(dOmega);
  }
  return(dDum);
}

double FourierSumD_DOUBLE(int iN,double *pdFC,int iNC,int iCMin,double dT){
  double dOmega;
  int i,iMin;
  double dDum=0;
  iMin=iCMin>0 ? iCMin : 1;
  for(i=iMin;i<iMin+iNC;i++){
    dOmega=2.0*M_PI*(double)i/(double)iN;
    dDum+=dOmega*(-pdFC[2*i]*sin(dT*dOmega)+pdFC[2*i+1]*cos(dT*dOmega));
  }
  return(dDum);
}

void MakeMaps(int iNB,double **ppdMaps,int iXdim,int iYdim,double *pdMapX,double *pdMapY){
  int i;
  unsigned long ul,ulXY;
  double dPhi,dRho;
  double *pdB=NULL,*pdFB=NULL;
  if(!(pdB=(double*)malloc(iNB*sizeof(double)))){
    printf("ANA Cannot allocate for pdB\n");
    return;
  }
  if(!(pdFB=(double*)malloc(iNB*sizeof(double)))){
    printf("ANA Cannot allocate for pdFB\n");
    return;
  }
  ulXY=(unsigned long)iXdim*(unsigned long)iYdim;
  for(ul=0;ul<ulXY;ul++){
    for(i=0;i<iNB;i++) pdB[i]=ppdMaps[i][ul];
    for(i=iMinFComponent;i<iMinFComponent+iNFComponents;i++){
      FindFTComponents_DOUBLE(iNB,pdB,i,pdFB+(2*i),pdFB+(2*i+1));
    }
    dPhi=FindExtremaFComponents_DOUBLE(iNB,pdFB,iNFComponents,iMinFComponent,1);
    dRho=FourierSum_DOUBLE(iNB,pdFB,iNFComponents,iMinFComponent,dPhi);
    dPhi*=2*M_PI/iNB;
    pdMapX[ul]=dRho*cos(dPhi);
    pdMapY[ul]=dRho*sin(dPhi);
  }
  if(pdB) free(pdB);
  if(pdFB) free(pdFB);
  printf("ANA Maps Created\n");
}
