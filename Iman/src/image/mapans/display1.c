/* Map analysis                 */
/* Written by V.Kalatsky        */
/* 06Feb2001                    */
/* Modified 10May2001           */

#include "mapans1.h"

#define DISTRIBUTION_WIN_SIZE 8.0
#define HISTOGRAM_PHI_SIZE 7.0
#define HISTOGRAM_RHO_SIZE 7.0

#define ELLIPSE_N 100

#define SECTION_SAVE_NONE 0
#define SECTION_SAVE_PHI 1
#define SECTION_SAVE_PHI_FIT 2
#define SECTION_SAVE_RHO 4
#define SECTION_SAVE_RHO_FIT 8
#define SECTION_SAVE_PHI_DER 16
#define SECTION_SAVE_PHI_DER_FIT 32

#define OD_SCORE_N 7

//#define ENUMERATE_SQUARES

int iHistogramFontSize=0;

int BinPhi(int iNBins,int *piBins,float *pfX,float *pfY,double *pdX,double *pdY,int *piMaxY,int *piMinY,int *piMaxX,int *piMinX);
int BinPhiOD(int iNBins,int *piBins,float *pfX,float *pfY,double *pdX,double *pdY,int *piMaxY,int *piMinY,int *piMaxX,int *piMinX,int iODScoreN,float *pfODScore,float *pfCBI,float *pfMI);
void PlotPhi(int iID,int iNBins,int *piBins,int iMaxY,float fMaxX,int iMode);
void PlotPhiOD(int iID,int iNBins,int *piBins,int iMaxY,float fMaxX,int iODScoreN,float *pfODScore,float fCBI,float fMI,int iMode);
void PlotRho(int iID,int iNBins,int *piBins,int iMaxY,float fMaxX,double dRhoMax,double dTLo,double dTHi,int iMode);

void DisplayMap(int mode,int xd,int yd,double *mx,double *my,double *gx,double *gy){
  float fmaxx,fminx,fmaxy,fminy,fmin,fmax,fmaxg;
  float fminc,fmaxc;
  float fmins,fmaxs;
  unsigned long imax=0,imin=0;
  double arm,arm2,srm,arg,spmv,axm,aym,axxm,ayym,axym,axmn,aymn,apmn,armn,apmv,armv;
  double normal2averagex,normal2averagey,normal2averager,normal2averagep;
  double sigmaxx,sigmayy;
  double dumd,dumdx,dumdy;
  double dRhoMin,dRhoMax;
  unsigned long i,j,ul,xy;
  double co,si,r,phi,dphi;
  int distributionid=0,callerid,iCallerID1;
  float dumf,dumfx,dumfy;
  char ch;
  double inversionx,inversiony;
  int save_single_bit=0,save_single_bit2=0;
  int iNPhiBins=0,iNRhoBins=0,iNPhiODBins=0;
  int *piPhiBins=NULL,*piRhoBins=NULL,*piPhiODBins=NULL;
  int iMaxPhiBin=0,iMaxRhoBin=0,iLocalMaxPhiBin=0;
  int iMaxPhiBinX=0;
  int iMaxPhiODBinX=0,iMaxPhiODBinY=0;
  float fAngleBinCoeff;
  double dThreshold;
  double dMaxResponseHeightSmooth,dMaxResponseHeightRaw;
  double dAreaAboveThreshold;
  double dVolumeAboveThreshold;
  double *pdDomainX=(double*)0,*pdDomainY =(double*)0;
  int iDomainN=0;
  float x,y;
  char str[256];
  unsigned long *pulSaveStencil=NULL,ulSaveStencilNP=0;
  int iRhoID,iPhiID,iPhiODID,iSaveID,iCallerID2;
  FILE *pF=NULL;
  char strFileName[256],strDum[128];
  double dRhoThresholdLo=-1.0,dRhoThresholdHi=-1.0;
  int iRhoThreshold=0;
  int iNormalizePhiHist=1;
  unsigned long ulPhiHistTotal=0;
  float fGrayScaleLo,fGrayScaleHi;
  int iRhoFlat;
  float afODScore[OD_SCORE_N],fCBI,fMI;
  int aiODScoreCount[OD_SCORE_N];
  int iODScoreN=OD_SCORE_N;
  int iSaveNBins,*piSaveBins;

  cpgqid(&callerid);

  xy=(unsigned long)xd*(unsigned long)yd;

  save_single_bit=mode&DISPLAY_MODE_SINGLE;
  save_single_bit2=mode&DISPLAY_MODE_SINGLE2;

  mode&=DISPLAY_MODE_MASK;

  fminx=fminy=BIG_FLOAT;
  fmaxx=fmaxy=-fminx;
  fmax=0.0;
  fmin=BIG_FLOAT;
  dRhoMin=BIG_FLOAT;
  dRhoMax=-BIG_FLOAT;
  arm=arm2=0.0;
  axm=aym=0.0;
  axxm=ayym=axym=0.0;
  axmn=aymn=0.0;
  memset((void *)fbufferx,0,xy*sizeof(float));
  memset((void *)fbuffery,0,xy*sizeof(float));
  memset((void *)pdBufferX,0,xy*sizeof(double));
  memset((void *)pdBufferY,0,xy*sizeof(double));
  memset((void *)fbufferphi,0,xy*sizeof(float));
  memset((void *)fbufferrho,0,xy*sizeof(float));
  co=cos(alpha);
  si=sin(alpha);
  if(inversion < -0.5){
    if(inversion_scheme==INVERSION_X){
      inversionx=-1.0;
      inversiony=1.0;
    }
    else{
      inversionx=1.0;
      inversiony=-1.0;
    }
  }
  else{
    inversionx=1.0;
    inversiony=1.0;    
  }
  for(j=0;j<stenciln;j++){
    i=stencil[j];
    //    axm+=(dumdx=mx[i]+xshift);
    //    aym+=(dumdy=my[i]+yshift);
    axm+=(dumdx=inversionx*(mx[i]+xshift)*co-inversiony*(my[i]+yshift)*si);
    aym+=(dumdy=inversiony*(my[i]+yshift)*co+inversionx*(mx[i]+xshift)*si);
    axxm+=dumdx*dumdx;
    ayym+=dumdy*dumdy;
    axym+=dumdx*dumdy;
    fbufferphi[i]=(float)atan2(dumdy,dumdx);
    dumf=fbufferx[i]=(float)(dumdx);
    //pdBufferX[i]=dumdx;
    if(fmaxx<dumf) fmaxx=dumf;
    if(fminx>dumf) fminx=dumf;
    dumf=fbuffery[i]=(float)(dumdy);
    //pdBufferY[i]=dumdy;
    if(fmaxy<dumf) fmaxy=dumf;
    if(fminy>dumf) fminy=dumf;
    fbufferrho[i]=(float)(dumd=hypot(dumdx,dumdy));
    arm+=dumd;
    arm2+=dumd*dumd;
    if(dumd!=0){
      axmn+=dumdx/dumd;
      aymn+=dumdy/dumd;
    }
    /*
    if(fmax<dumf){ 
      fmax=dumf;
      imax=i;
    }
    if(fmin>dumf){ 
      fmin=dumf;
      imin=i;
    }
    */
    if(dRhoMax<dumd){ 
      dRhoMax=dumd;
      imax=i;
    }
    if(dRhoMin>dumd){ 
      dRhoMin=dumd;
      imin=i;
    }
    fmax=dRhoThresholdHi=dRhoMax;
    fmin=dRhoThresholdLo=dRhoMin;
    dRhoThresholdHi=dRhoMax*1.001;
    dRhoThresholdLo=dRhoMin*0.999;
    iRhoFlat=0;
    if(dRhoThresholdHi==dRhoThresholdLo){
      iRhoFlat=1;
      dRhoThresholdHi=1.01*dRhoMax;
      dRhoThresholdLo=0.99*dRhoMin;
    }
  }

  if(!mode) return;

  if(fmaxx==fminx){
    if(fmaxy==fminy){
      fmaxx+=0.5;
      fminx-=0.5;
      fmaxy+=0.5;
      fminy-=0.5;
    }
    else{
      fmaxx+=(fmaxy-fminy)*0.5;
      fminx-=(fmaxy-fminy)*0.5;
    }
  }
  if(fmaxy==fminy){
    fmaxy+=(fmaxx-fminx)*0.5;
    fminy-=(fmaxx-fminx)*0.5;
  }

  axm/=(double)stenciln;
  aym/=(double)stenciln;
  axxm/=(double)stenciln;
  ayym/=(double)stenciln;
  axym/=(double)stenciln;
  axxm-=axm*axm;
  ayym-=aym*aym;
  axym-=axm*aym;

  if(iGetNoiseLevel){ 
    dRhoNoise=sqrt(axxm+ayym);
    iGetNoiseLevel=0;
    printf("DIS Noise level %f\n",dRhoNoise);
  }

  if(do_verbose){
    printf("DIS Mean X=%f Y=%f\n",axm,aym);
    printf("DIS Dev  XX=%f YY=%f XY=%f sqrt(XX+YY)=%f\n",axxm,ayym,axym,sqrt(axxm+ayym));
  }

  dumd=sqrt((axxm-ayym)*(axxm-ayym)+4.0*axym*axym);
  sigmaxx=0.5*(axxm+ayym+dumd);
  sigmayy=0.5*(axxm+ayym-dumd);

  axmn/=(double)stenciln;
  aymn/=(double)stenciln;
  if((armn=hypot(aymn,axmn))!=0.0) apmn=atan2(aymn,axmn);
  else  apmn=0.0;
  if(do_verbose) printf("DIS armn=%e apmn=%f(x=%e y=%e)\n",armn,apmn,axmn,aymn);

  arm/=(double)stenciln;
  arm2/=(double)stenciln;
  srm=sqrt(arm2-arm*arm);

  apmv=atan2(aym,axm);
  if((armv=hypot(aym,axm))!=0.0){
    normal2averagex=-aym/armv;
    normal2averagey=axm/armv;
    normal2averager=sqrt(aym*aym*axxm-2.0*axm*aym*axym+axm*axm*ayym)/armv;
    normal2averagep=atan2(normal2averager,armv);
  }
  else{
    normal2averagex=1.0;
    normal2averagey=0.0;
    normal2averager=sqrt((axxm-2.0*axym+ayym)/2.0);
    normal2averagep=atan2(normal2averager,armv);
  }
  spmv=sqrt(axm*axm*axxm+aym*aym*ayym)/(armv*armv);

  // X(min max) and Y  (min max)
  //  printf("DIS X(%f %f) Y(%f %f)\n",fminx,fmaxx,fminy,fmaxy);
  printf("DIS Rho: Min=%.3f(%lu %lu) Max=%.3f(%lu %lu)\n",fmin,imin%xd,imin/xd,fmax,imax%xd,imax/xd);
  printf("DIS Average Rho: R=%.3f sqrt(R2)=%.3f sqrt(R2-R*R)=%.3f\n",arm,sqrt(arm2),srm);
  printf("DIS Average X,Y: N=%li P=%.3f(%.3f) R=%.3f(%.3f) X/Y=%.3f\n",stenciln,apmv*180.0/M_PI,normal2averagep*180.0/M_PI,armv,normal2averager,(axm-aym)/(axm+aym));
  printf("DIS Normalized Average PN=%.3f RN=%.3f\n",apmn*180.0/M_PI,armn);

  if(do_phase_transform){
    if(apmv<0.0) apmv+=2.0*M_PI;
    printf("DIS Time=%f(%f) %f\n",apmv*phase_transform,spmv*phase_transform,apmn*phase_transform);
  }


  if(do_distribution){
    iNPhiBins=DISTRIBUTION_PHI_BINS_N;
    //iNPhiBins=16;
    if(!(piPhiBins=(int *)calloc(iNPhiBins+1,sizeof(int)))){
      printf("INI Cannot allocate for piPhiBins\n");
    }
    iNRhoBins=DISTRIBUTION_RHO_BINS_N;
    if(!(piRhoBins=(int *)calloc(iNRhoBins,sizeof(int)))){
      printf("INI Cannot allocate for piRhoBins\n");
    }
    if(iDoODAnalysis){
      if(iODAnalysisNBins>0) iNPhiODBins=iODAnalysisNBins;
      else iNPhiODBins=DISTRIBUTION_PHI_BINS_N;
      if(!(piPhiODBins=(int *)calloc(iNPhiODBins+1,sizeof(int)))){
	iDoODAnalysis=0;
	printf("INI Cannot allocate for iNPhiODBins\n");
      }
    }
    if(!(pulSaveStencil=(unsigned long *)malloc(stenciln*sizeof(unsigned long)))){
      printf("INI Cannot allocate for pulSaveStencil\n");
    }
    else{
      ulSaveStencilNP=stenciln;
      memcpy(pulSaveStencil,stencil,ulSaveStencilNP*sizeof(unsigned long));
    }

    distributionid=cpgopen("/xw");
    cpgask(0);
    cpgpap(DISTRIBUTION_WIN_SIZE,(fmaxy-fminy)/(fmaxx-fminx));
    cpgsci(1);
    cpgenv(fminx,fmaxx,fminy,fmaxy,1,1);

  label_start_distribution:

    cpgsci(3);
    cpgbbuf();

    for(j=0;j<stenciln;j++){
      i=stencil[j];
      cpgpt1(fbufferx[i],fbuffery[i],1);
    }
    BinPhi(iNPhiBins,piPhiBins,fbufferx,fbuffery,NULL,NULL,&iMaxPhiBin,NULL,&iMaxPhiBinX,NULL);

    fAngleBinCoeff=(fmaxy-fminy)<(fmaxx-fminx) ? (fmaxy-fminy)*0.5 : (fmaxx-fminx)*0.5;
    fAngleBinCoeff/=1.5*iMaxPhiBin;

    if(fmax>0.0){
      for(j=0;j<stenciln;j++){
	i=stencil[j];
	r=hypot((double)(fbuffery[i]),(double)(fbufferx[i]))/fmax;
	piRhoBins[(int)floor(r*iNRhoBins)]++;
      }
      for(iMaxRhoBin=piRhoBins[0],j=1;j<iNRhoBins;j++){
	if(iMaxRhoBin<piRhoBins[j]) iMaxRhoBin=piRhoBins[j];
      }
    }

    cpgsch(1.2);
    cpgsci(4);
    cpgarro(0.0,0.0,(float)(arm*axmn/armn),(float)(arm*aymn/armn));
    cpgsci(2);
    cpgarro(0.0,0.0,(float)axm,(float)aym);
    if(axym!=0.0) phi=atan2((sigmaxx-axxm)/axym,1.0);
    else phi=0.0;
    co=cos(phi);
    si=sin(phi);
    dphi=2.0*M_PI/(double)(ELLIPSE_N);
    dumdx=sqrt(sigmaxx);
    dumdy=0.0;
    cpgmove((float)(axm+co*dumdx-si*dumdy),(float)(aym+co*dumdy+si*dumdx));
    for(j=1;j<=ELLIPSE_N;j++){
      phi=dphi*(double)j;
      dumdx=cos(phi);
      dumdy=sin(phi);
      r=sqrt(sigmaxx*sigmayy/(sigmayy*dumdx*dumdx+sigmaxx*dumdy*dumdy));
      dumdx*=r;
      dumdy*=r;
      cpgdraw((float)(axm+co*dumdx-si*dumdy),(float)(aym+co*dumdy+si*dumdx));
    }
    cpgsci(7);
    cpgarro((float)axm,(float)aym,(float)(axm+normal2averagex*normal2averager),(float)(aym+normal2averagey*normal2averager));

    cpgsci(1);
    //cpgmove((float)(fAngleBinCoeff*piPhiBins[0]),0.0);
    cpgmove((float)(fAngleBinCoeff*piPhiBins[0]*cos((2*M_PI*(0.5))/iNPhiBins)),
	    (float)(fAngleBinCoeff*piPhiBins[0]*sin((2*M_PI*(0.5))/iNPhiBins)));
    for(j=1;j<iNPhiBins+1;j++){
      i=j%iNPhiBins;
      cpgdraw((float)(fAngleBinCoeff*piPhiBins[i]*cos((2*M_PI*(i+0.5))/iNPhiBins)),
	      (float)(fAngleBinCoeff*piPhiBins[i]*sin((2*M_PI*(i+0.5))/iNPhiBins)));
    }
    cpgebuf();
    
  label_loop:

    cpgband(0,0,dumfx,dumfy,&dumfx,&dumfy,&ch);
    if(ch=='c'){
      
      cpgsci(1);
  
    label_start:

      cpgband(0,0,x,y,&x,&y,&ch);
      if(ch=='q'){
	goto label_loop;
      }
      if(ch=='m'){
	if(GetString(str)) x=atof(str);
	else goto label_start;
	if(GetString(str)) y=atof(str);
	else goto label_start;
      }
      iDomainN=1;
      pdDomainX=(double *)malloc(iDomainN*sizeof(double));
      pdDomainY=(double *)malloc(iDomainN*sizeof(double));
      pdDomainX[iDomainN-1]=(double)x;
      pdDomainY[iDomainN-1]=(double)y;
      printf("%i %f %f\n",iDomainN,x,y);
      cpgmove(x,y);
      cpgsci(2);
      cpgpt1(x,y,2);
      cpgsci(1);
      while(1){
	cpgband(1,0,pdDomainX[iDomainN-1],pdDomainY[iDomainN-1],&x,&y,&ch);
	if(ch=='q'){
	  iDomainN=0;
	  if(pdDomainX) free(pdDomainX);
	  if(pdDomainY) free(pdDomainY);
	  goto label_loop;
	}
	if(ch=='A' || ch=='m'){
	  if(ch=='m'){
	    if(GetString(str)) x=atof(str);
	    else continue;
	    if(GetString(str)) y=atof(str);
	    else continue;
	  }
	  iDomainN++;
	  pdDomainX=(double *)realloc((void *)pdDomainX,iDomainN*sizeof(double));
	  pdDomainY=(double *)realloc((void *)pdDomainY,iDomainN*sizeof(double));
	  pdDomainX[iDomainN-1]=(double)x;
	  pdDomainY[iDomainN-1]=(double)y;
	  cpgdraw(x,y);
	  cpgsci(2);
	  cpgpt1(x,y,2);
	  cpgsci(1);
	  printf("%i %f %f\n",iDomainN,x,y);
	}
	if(ch=='D'){
	  iDomainN=0;
	  free(pdDomainX);
	  free(pdDomainY);
	  goto label_start_distribution;
	}
	if(ch=='X'){
	  if(iDomainN>2){ 
	    pdDomainX=(double *)realloc((void *)pdDomainX,(iDomainN+1)*sizeof(double));
	    pdDomainY=(double *)realloc((void *)pdDomainY,(iDomainN+1)*sizeof(double));
	    pdDomainX[iDomainN]=pdDomainX[0];
	    pdDomainY[iDomainN]=pdDomainY[0];
	    cpgdraw((float)pdDomainX[iDomainN],(float)pdDomainY[iDomainN]);	    
	    StencilMap(xd,yd,fbufferx,fbuffery,iDomainN,pdDomainX,pdDomainY,&stenciln,&stencil);	
	    iDomainN=0;
	    if(pdDomainX) free(pdDomainX);
	    if(pdDomainY) free(pdDomainY);

	    do_distribution=0;
	    DisplayMap(mode,xd,yd,mx,my,gx,gy);
	    do_distribution=1;

	    goto label_loop;
	  }
	  else printf("INT Needs three points to complete selection\n");
	}
      }      
    } 

    if(ch=='H'){
      cpgqid(&iCallerID1);
      cpgsave();
      iPhiID=cpgopen("/xw");
      cpgask(0);
      cpgpap(HISTOGRAM_PHI_SIZE,0.6);
      PlotPhi(iPhiID,iNPhiBins,piPhiBins,iMaxPhiBin,(float)iMaxPhiBinX,1);

      iRhoID=cpgopen("/xw");
      cpgask(0);
      cpgpap(HISTOGRAM_RHO_SIZE,0.6);
      cpgsci(1);
      PlotRho(iRhoID,iNRhoBins,piRhoBins,iMaxRhoBin,0,dRhoMax,dRhoThresholdLo,dRhoThresholdHi,1);

      if(iDoODAnalysis){
	BinPhiOD(iNPhiODBins,piPhiODBins,fbufferx,fbuffery,NULL,NULL,&iMaxPhiODBinY,NULL,&iMaxPhiODBinX,NULL,iODScoreN,afODScore,&fCBI,&fMI);
	iPhiODID=cpgopen("/xw");
	cpgask(0);
	cpgpap(HISTOGRAM_PHI_SIZE,0.6);
	PlotPhiOD(iPhiODID,iNPhiODBins,piPhiODBins,(float)iMaxPhiODBinY,(float)iMaxPhiODBinX,iODScoreN,afODScore,fCBI,fMI,1);
      }
      cpgslct(iRhoID);
      while(1){
	cpgband(6,0,x,y,&x,&y,&ch);
	if(ch=='q' || ch=='Q') break;

	if(ch=='h'){
	  printf("\n");
	  printf("Left mouse button or X - select Low Rho threshold\n");
	  printf("Right mouse button or A - select High Rho threshold\n");
	  printf("F - change font size switch, cycles among 3 values: small, normal, large\n");
	  printf("S - save Phi, OD if in OD mode, histogram data, file name is auto generated\n");
	  printf("h - print this help\n");
	  printf("q or Q - quit interactive dialog\n");
	  printf("s - save Phi, OD if in OD mode, histogram as a picture (PS) file, user prompted for file name\n");
	  continue;
	}
	if(ch=='s' && iDoODAnalysis){
	  printf("%sODHist.ps\n",filenames[0]);
	  fflush(stdout);
	  if(GetStringFromSTDIN(str,255,"DIS Enter file name for OD historgram: ")) continue;
	  snprintf(strFileName,255,"%s/CPS",str);
	  cpgsave();
	  cpgqid(&iCallerID2);
	  iSaveID=cpgopen(strFileName);
	  cpgask(0);
	  cpgpap(HISTOGRAM_PHI_SIZE,0.6);
	  PlotPhiOD(iSaveID,iNPhiODBins,piPhiODBins,(float)iMaxPhiODBinY,(float)iMaxPhiODBinX,iODScoreN,afODScore,fCBI,fMI,3);
	  cpgslct(iCallerID2);
	  cpgunsa();
	  continue;
	}
	if(ch=='S'){
	  if(iDoODAnalysis){
	    snprintf(strFileName,255,"%s_HistODLo%0.0fHi%0.0f",filenames[0],iNRhoBins*dRhoThresholdLo/fmax,iNRhoBins*dRhoThresholdHi/fmax);
	    printf("INT Saving OD Histogram in %s\n",strFileName);
	    iSaveNBins=iNPhiODBins;
	    piSaveBins=piPhiODBins;
	  }
	  else{
	    snprintf(strFileName,255,"%s_HistPhiLo%0.0fHi%0.0f",filenames[0],iNRhoBins*dRhoThresholdLo/fmax,iNRhoBins*dRhoThresholdHi/fmax);
	    printf("INT Saving Phi Histogram in %s\n",strFileName);
	    iSaveNBins=iNPhiBins;
	    piSaveBins=piPhiBins;
	  }
	  if((pF=fopen(strFileName,"w"))==NULL){
	    printf("INT Cannot open file %s\n",strFileName);
	    goto bailout3;
	  }
	  sprintf(strDum,"#BEGIN SAMPLE SECTION\n");
	  if(fprintf(pF,strDum) != strlen(strDum)){
	  printf("INT Cannot write %s to file %s\n",strDum,strFileName);
	  goto bailout3;
	  }
	  sprintf(strDum,"# Bins %i\n",iSaveNBins);
	  if(fprintf(pF,strDum) != strlen(strDum)){
	    printf("INT Cannot write iNPhiBins to file %s\n",strFileName);
	    goto bailout3;
	  }
	  for(i=0,ulPhiHistTotal=0;i<iSaveNBins;i++) ulPhiHistTotal+=piSaveBins[i];
	  sprintf(strDum,"# Total Samples %li\n",ulPhiHistTotal);
	  if(fprintf(pF,strDum) != strlen(strDum)){
	    printf("INT Cannot write ulPhiHistTotal to file %s\n",strFileName);
	    goto bailout3;
	  }

	  if(iNormalizePhiHist){
	    for(i=0;i<iSaveNBins;i++){
	      sprintf(strDum,"%li %0.5f\n",i,piSaveBins[i]/(float)ulPhiHistTotal);
	      if(fprintf(pF,strDum) != strlen(strDum)){
		printf("INT Cannot write sample %li to file %s\n",i,strFileName);
		goto bailout3;
	      }
	    }
	  }
	  else{
	    for(i=0;i<iSaveNBins;i++){
	      sprintf(strDum,"%li %i\n",i,piSaveBins[i]);
	      if(fprintf(pF,strDum) != strlen(strDum)){
		printf("INT Cannot write sample %li to file %s\n",i,strFileName);
		goto bailout3;
	      }
	    }
	  }
	  sprintf(strDum,"#END SAMPLE SECTION\n");
	  if(fprintf(pF,strDum) != strlen(strDum)){
	    printf("INT Cannot write %s to file %s\n",strDum,strFileName);
	    goto bailout3;
	  }
	  if(iDoODAnalysis){
	    sprintf(strDum,"#BEGIN CBI SECTION\n");
	    if(fprintf(pF,strDum) != strlen(strDum)){
	      printf("INT Cannot write %s to file %s\n",strDum,strFileName);
	      goto bailout3;
	    }
	    sprintf(strDum,"# CBI Bins %i\n",iODScoreN);
	    if(fprintf(pF,strDum) != strlen(strDum)){
	      printf("INT Cannot write iODScoreN to file %s\n",strFileName);
	      goto bailout3;
	    }
	    for(i=0;i<iODScoreN;i++){
	      sprintf(strDum,"# %li %f\n",i,afODScore[i]);
	      if(fprintf(pF,strDum) != strlen(strDum)){
		printf("INT Cannot write sample %li to file %s\n",i,strFileName);
		goto bailout3;
	      }
	    }
	    sprintf(strDum,"# CBI %f\n",fCBI);
	    if(fprintf(pF,strDum) != strlen(strDum)){
	      printf("INT Cannot write CBI to file %s\n",strFileName);
	      goto bailout3;
	    }
	    sprintf(strDum,"# MI %f\n",fMI);
	    if(fprintf(pF,strDum) != strlen(strDum)){
	      printf("INT Cannot write MI to file %s\n",strFileName);
	      goto bailout3;
	    }
	    sprintf(strDum,"#END CBI SECTION\n");
	    if(fprintf(pF,strDum) != strlen(strDum)){
	      printf("INT Cannot write %s to file %s\n",strDum,strFileName);
	      goto bailout3;
	    }
	  }
	bailout3:
	  fclose(pF);
	  continue;
	}
	
	if(ch=='A'){ 
	  if(dRhoThresholdHi<fmax*x/(iNRhoBins-1)){
	    printf("DIS Rho Threshold Hi<Lo");
	  }
	  else{
	    iRhoThreshold=1;
	    dRhoThresholdLo=fmax*x/(iNRhoBins-1);
	    if(dRhoThresholdLo<0.0){
	      if(iRhoFlat) dRhoThresholdLo=dRhoMin*0.99;
	      else dRhoThresholdLo=dRhoMin*0.999;
	    }
	  }
	}
	if(ch=='X'){
	  if(dRhoThresholdLo>fmax*x/(iNRhoBins-1)){
	    printf("DIS Rho Threshold Lo>Hi");
	  }
	  else{
	    iRhoThreshold=2;
	    dRhoThresholdHi=fmax*x/(iNRhoBins-1);
	    if(dRhoThresholdHi>dRhoMax){
	      if(iRhoFlat) dRhoThresholdHi=dRhoMax*1.01;
	      else dRhoThresholdHi=dRhoMax*1.001;
	    }
	  }
	}

	if(ch=='F'){ 
	  iHistogramFontSize=(iHistogramFontSize+1)%3;
	  iRhoThreshold=4;
	}

	if(iRhoThreshold){
	  printf("Threshold Lo %f Hi %f\n",dRhoThresholdLo,dRhoThresholdHi);
	  
	  if(pulSaveStencil){
	    stenciln=ulSaveStencilNP;
	    if(!(stencil=(unsigned long *)realloc(stencil,stenciln*sizeof(unsigned long)))){
	      printf("INI Cannot realloc for stencil\n");
	    }
	    memcpy(stencil,pulSaveStencil,ulSaveStencilNP*sizeof(unsigned long));
	  }
	  
	  DisplayMap(DISPLAY_MODE_NOT,xd,yd,mx,my,gx,gy);
	  
	  StencilMap(xd,yd,fbufferx,fbuffery,1,&dRhoThresholdLo,&dRhoThresholdHi,&stenciln,&stencil);	
	  
	  BinPhi(iNPhiBins,piPhiBins,fbufferx,fbuffery,NULL,NULL,&iLocalMaxPhiBin,NULL,&iMaxPhiBinX,NULL);
	  PlotPhi(iPhiID,iNPhiBins,piPhiBins,iLocalMaxPhiBin,(float)iMaxPhiBinX,1);
	  if(iDoODAnalysis){
	    BinPhiOD(iNPhiODBins,piPhiODBins,fbufferx,fbuffery,NULL,NULL,&iMaxPhiODBinY,NULL,&iMaxPhiODBinX,NULL,iODScoreN,afODScore,&fCBI,&fMI);
	    PlotPhiOD(iPhiODID,iNPhiODBins,piPhiODBins,(float)iMaxPhiODBinY,(float)iMaxPhiODBinX,iODScoreN,afODScore,fCBI,fMI,1);
	  }
	  PlotRho(iRhoID,iNRhoBins,piRhoBins,iMaxRhoBin,0,dRhoMax,dRhoThresholdLo,dRhoThresholdHi,1);
	  
	  do_distribution=0;
	  DisplayMap(mode,xd,yd,mx,my,gx,gy);
	  do_distribution=1;
	}
      }

      cpgslct(iPhiODID);
      cpgclos();
      cpgslct(iRhoID);
      cpgclos();
      cpgslct(iPhiID);
      cpgclos();
      cpgslct(iCallerID1);
      cpgunsa();
    }

    cpgclos();
    do_distribution=0;
    cpgslct(callerid);

    if(pulSaveStencil) free(pulSaveStencil);
    if(piRhoBins) free(piRhoBins);
    if(piPhiBins) free(piPhiBins);
    if(piPhiODBins) free(piPhiODBins);
    return;
  }

  if(do_sectionplot){
    PlotSection(0/*mode*/,xd,yd,fbufferx,fbuffery);
    return;
  }
 
  if(cradius!=0){
    memcpy(chopfbuffer,fbufferrho,xy*sizeof(float));
    AverageMaps(xd,yd,(double**)0,&chopfbuffer,1,dcradius,0,(int)FLOAT);
    FindExtrema(chopfbuffer,&fminc,&fmaxc,&imin,&imax);
  }
  else{
    fmaxc=fmax;
    fminc=fmin;
  }
  
  if(fmaxc==fminc){
    fmaxc+=1.0;
    fminc-=1.0;
  }
  else{
    dMaxResponseHeightRaw=0.0;
    dMaxResponseHeightSmooth=fmaxc;
    dThreshold=dMaxResponseHeightSmooth*amplitude_chop_l;
    dAreaAboveThreshold=0.0;
    dVolumeAboveThreshold=0.0;
    for(j=0;j<stenciln;j++){
      i=stencil[j];
      if((dumd=hypot(mx[i],my[i]))>dThreshold){
	dAreaAboveThreshold+=1.0;
	dVolumeAboveThreshold+=(dumd-dThreshold);
      }
      if(dMaxResponseHeightRaw<dumd) dMaxResponseHeightRaw=dumd;
    }
    dumd=dPixelSizeU*dPixelSizeU/1000000.0;
    // Val: uncomment if need all sorts of averages
    /*
    printf("DIS -> Threshold %0.1f%% -> %.4f\n",amplitude_chop_l*100.0,dThreshold);
    printf("DIS -> Height Max Smooth %.4f Raw %.4f\n",dMaxResponseHeightSmooth,dMaxResponseHeightRaw);
    printf("DIS -> H %.3f A %.3f V %.3f mm2 HA %.3f\n",dMaxResponseHeightRaw,dAreaAboveThreshold*dumd,dVolumeAboveThreshold*dumd,dVolumeAboveThreshold/dAreaAboveThreshold);
    */
  }

  if(do_verbose) printf("DIS MinChop=%f MaxChop=%f\n",fminc,fmaxc);

  cpgslct(mainid);
  cpgsch(5.5);
  cpgpanl(1,1);

  if(save_single_bit2) goto label_rho_map;

  if(iDoPhaseScale){
    for(j=0;j<stenciln;j++){
      i=stencil[j];
      fbufferphi[i]*=dPhaseScale;
    }
  }
  /*
  for(j=0,i=0;i<xy;i++){
    if(j<stenciln){
      if(i!=stencil[j]) fbufferphi[i]=100;//M_PI*1.01;
      else j++;
    }
    else fbufferphi[i]=M_PI*1.1;
  }
  */
  //  cpgsch(10.5);
  cpgeras();

  switch(color_scheme){
  case COLOR_RGB:  
    //Val: to draw just the wedge or for blank image comment out next line
    if(!iDoNotPlotPhase) DrawColorImage(fbufferphi,xd,yd,M_PI,-M_PI,transs);
    if(iDrawColorWedge) 
      DrawColorWedge(wedge_pos, wedge_width, -180.0/harmonic_d, 180.0/harmonic_d,"");
    break;
  case COLOR_HLS:
    //    DrawNiceColorImagePolar(fbufferphi,fbufferrho,xd,yd,fmaxc,fminc,transs);
    DrawNiceColorImageCartesian(fbufferx,fbuffery,xd,yd,transs);
  //Moved from DrawNiceColorImageCartesian
    if(iDrawColorWedge) 
      DrawColorWedge(wedge_pos, wedge_width, -180.0/harmonic_d, 180.0/harmonic_d,"");
    break;
  default:
    printf("DIS Bad color scheme %i\n",color_scheme);
    return;
  }

  if(depth_lines && contours){
    for(i=0;i<xy;i++){
      smoothfbuffers[0][i]=(float)(inversionx*(maps_orig[0][i]+xshift)*co-inversiony*(maps_orig[1][i]+yshift)*si);
      smoothfbuffers[1][i]=(float)(inversiony*(maps_orig[1][i]+yshift)*co+inversionx*(maps_orig[0][i]+xshift)*si);
    }
    if(sradius) AverageMaps(xd,yd,(double**)0,smoothfbuffers,2,dsradius,0,(int)FLOAT);
    Cartesian2Polar(xy,smoothfbuffers[0],smoothfbuffers[1],smoothfbuffers[2],smoothfbuffers[3]);
    FindExtrema(smoothfbuffers[2],&fmins,&fmaxs,&imin,&imax);
    DrawContours(CONTOUR_PHI,contour_scheme_phi,smoothfbuffers,xd,yd,fmaxs,fmins,smartcontours_threshold,transs);
  }
  else{
    if(iPlotExternalContours && iNDepthLinesExternal) DrawContours(CONTOUR_PHI,contour_scheme_phi,smoothfbuffers,xd,yd,fmaxs,fmins,smartcontours_threshold,transs);
  }

  if(iMarkPoints && !iMarkPointsHide) MarkPoints(xd,iMarkPointsColor-2*save_on,INT_MODE_PHI);
  if(iPlotLineSegments){
    if(save_on) PlotLineSegments(0,4,iPlotLineSegmentsStyle);// 0 - white, 1 - black on save
    else PlotLineSegments(0,4,iPlotLineSegmentsStyle);
  }

  if(save_single_bit) return;

  cpgpanl(1,2);

 label_rho_map:

  cpgeras();

  if(iFixedGrayScaleLo) fGrayScaleLo=fFixedGrayScaleLo;
  else fGrayScaleLo=fminc;
  if(iFixedGrayScaleHi) fGrayScaleHi=fFixedGrayScaleHi;
  else fGrayScaleHi=fmaxc;

  DrawGrayImage(fbufferrho,xd,yd,
		fGrayScaleLo+amplitude_chop_u*(fGrayScaleHi-fGrayScaleLo),
		fGrayScaleLo+amplitude_chop_l*(fGrayScaleHi-fGrayScaleLo),transs);

  if(depth_lines && contours){
    DrawContours(CONTOUR_RHO,contour_scheme_rho,smoothfbuffers,xd,yd,fmaxs,fmins,smartcontours_threshold,transs);
  }

  // imax and imin are obtained at smoothing not chopping
  if(iMarkMinMax){
    cpgslw(2);
    cpgsci(2);
    cpgpt1((float)(imax%xd)*transs[1]+transs[0],(float)(imax/xd)*transs[5]+transs[3],2);
    cpgsci(4);
    cpgpt1((float)(imin%xd)*transs[1]+transs[0],(float)(imin/xd)*transs[5]+transs[3],2);
    cpgslw(1);
  }

  if(iMarkPoints && !iMarkPointsHide) MarkPoints(xd,iMarkPointsColor+1,INT_MODE_RHO);
  if(iPlotLineSegments){
   if(save_on) PlotLineSegments(3,4,iPlotLineSegmentsStyle);
   else PlotLineSegments(3,4,iPlotLineSegmentsStyle);
  }

  if(save_single_bit2) return;

  if(mode>DISPLAY_MODE2){
    cpgpanl(2,1);
    cpgeras();
    DrawGrayImage(fbufferx,xd,yd,fmaxx,fminx,transs);
    cpgpanl(2,2);
    cpgeras();
    DrawGrayImage(fbuffery,xd,yd,fmaxy,fminy,transs);
  }

  if(mode>DISPLAY_MODE4){
    fminx=fminy=10000000000.0;
    fmaxx=fmaxy=-fminx;
    fmaxg=0.0;
    arg=0.0;
    for(ul=0;ul<stenciln;ul++){
      i=stencil[ul];
      dumf=fbufferx[i]=(float)gx[i];
      if(fmaxx<dumf) fmaxx=dumf;
      if(fminx>dumf) fminx=dumf;
      dumf=fbuffery[i]=(float)gy[i];
      if(fmaxy<dumf) fmaxy=dumf;
      if(fminy>dumf) fminy=dumf;
      arg+=(dumf=(float)hypot(gx[i],gy[i]));
      if(fmaxg<dumf) fmaxg=dumf;
    }
    arg/=(double)stenciln;
    printf("X(%f %f) Y(%f %f) Max=%f Ave=%f\n",fminx,fmaxx,fminy,fmaxy,fmaxg,arg);
    printf("MaxR=%f AveR=%f\n",fmaxg/fmax,arg/arm);

    cpgpanl(3,1);
    cpgeras();
    DrawGrayImage(fbufferx,xd,yd,fmaxx,fminx,transs);
    cpgpanl(3,2);
    cpgeras();
    DrawGrayImage(fbuffery,xd,yd,fmaxy,fminy,transs);
  }
  cpgslct(callerid);
}

int BinPhi(int iNBins,int *piBins,float *pfX,float *pfY,double *pdX,double *pdY,int *piMaxY,int *piMinY,int *piMaxX,int *piMinX){
  int i,j,iDum,iTotal=0;
  double phi,dTwoPiOver;
  int iMaxX,iMaxY,iOutBoundCount=0;

  dTwoPiOver=1.0/(2.0*M_PI);
  memset(piBins,0,sizeof(int)*iNBins);
  if(pdX && pdY){
    for(j=0;j<stenciln;j++){
      i=stencil[j];
      phi=atan2(pdY[i],pdX[i])*dTwoPiOver;
      if(phi<0.0) phi+=1.0;
      iDum=(int)floor(phi*iNBins);
      if(iDum>=iNBins){
	printf("DIS BinPhi Out of bound\n");
	iOutBoundCount++;
	iDum=iNBins-1;
      }
      piBins[iDum]++;
      iTotal++;
    }
  }
  else{
    for(j=0;j<stenciln;j++){
      i=stencil[j];
      phi=atan2((double)pfY[i],(double)pfX[i])*dTwoPiOver;
      if(phi<0.0) phi+=1.0;
      iDum=(int)floor(phi*iNBins);
      if(iDum>=iNBins){
	printf("DIS BinPhi Out of bound %i(%i)\n",iDum,iNBins);
	continue;
	iOutBoundCount++;
	iDum=iNBins-1;
      }
      piBins[iDum]++;
      iTotal++;
    }
  }
  piBins[iNBins]=piBins[0];

  for(iMaxY=piBins[0],iMaxX=0,j=1;j<iNBins;j++){
    if(iMaxY<piBins[j]) iMaxY=piBins[(iMaxX=j)];
  }
  if(piMaxY) *piMaxY=iMaxY;
  if(piMaxX) *piMaxX=iMaxX;
  printf("DIS BinPhi Total %i Max X %i Y %i\n",iTotal,iMaxX,iMaxY);
  if(iOutBoundCount) printf("DIS BinPhi Out of bound count %i\n",iOutBoundCount);

  return(0);
}

int BinPhiOD(int iNBins,int *piBins,float *pfX,float *pfY,double *pdX,double *pdY,int *piMaxY,int *piMinY,int *piMaxX,int *piMinX,int iODScoreN,float *pfODScore,float *pfCBI,float *pfMI){
  int i,j,iDum,iTotal=0;
  double phi,dTwoPiOver;
  int iMaxX,iMaxY,iOutBoundCount=0;
  float fDum,fCBI,fMI;

  dTwoPiOver=1.0/(2.0*M_PI);
  memset(piBins,0,sizeof(int)*iNBins);
  memset(pfODScore,0,sizeof(float)*iODScoreN);
  if(pdX && pdY){
    for(j=0;j<stenciln;j++){
      i=stencil[j];
      //phi=atan2(pdY[i],pdX[i])*dTwoPiOver;
      phi=pdY[i]/(pdY[i]+pdX[i]);
      if(phi<0.0) phi=0.0;
      //iDum=(int)floor(4.0*phi*iNBins);
      iDum=(int)floor(phi*iNBins);
      if(iDum==iNBins) iDum=iNBins-1;
      if(iDum>iNBins){
	printf("DIS BinPhiOD Out of bound\n");
	iOutBoundCount++;
	continue;
      }
      piBins[iDum]++;
      iTotal++;
    }
  }
  else{
    for(j=0;j<stenciln;j++){
      i=stencil[j];
      //phi=atan2((double)pfY[i],(double)pfX[i])*dTwoPiOver;
      phi=(double)pfY[i]/((double)pfY[i]+(double)pfX[i]);
      iDum=(int)floor(phi*iNBins);
      if(iDum<0) iDum=0;
      if(iDum==iNBins) iDum=iNBins-1;
      if(iDum>iNBins){
	printf("DIS BinPhiOD Out of bound %i(%i)\n",iDum,iNBins);
	iOutBoundCount++;
	continue;
      }
      piBins[iDum]++;
      iTotal++;
    }
  }
  if(iTotal==0) iTotal=1;
 
  for(iMaxY=piBins[0],iMaxX=0,j=1;j<iNBins;j++){
    if(iMaxY<piBins[j]) iMaxY=piBins[(iMaxX=j)];
  }
  if(piMaxY) *piMaxY=iMaxY;
  if(piMaxX) *piMaxX=iMaxX;

  fDum=iODScoreN/(float)iNBins;
  for(j=0;j<iNBins;j++){
    pfODScore[(int)floor(j*fDum)]+=piBins[j];
  }
  fMI=fCBI=0.0;
  for(j=0;j<iODScoreN;j++){
    fCBI+=(pfODScore[j]/=iTotal)*(iODScoreN-1-j);
    fMI+=abs(j-iODScoreN/2)*pfODScore[j];
  }
  if(iODScoreN>1){
    fCBI/=iODScoreN-1;
    fMI/=(iODScoreN/2);
  }
  if(pfCBI) *pfCBI=fCBI;
  if(pfMI) *pfMI=fMI;

  printf("DIS BinPhiOD Total %i Max X %i Y %i CBI %0.3f MI %0.3f\n",iTotal,iMaxX,iMaxY,fCBI,fMI);
  if(iOutBoundCount) printf("DIS BinPhiOD Out of bound count %i\n",iOutBoundCount);

  return(0);
}


void PlotPhi(int iID,int iNBins,int *piBins,int iMaxY,float fMaxX,int iMode){
   int iCallerID;
  int j;
  char str[128];

  cpgqid(&iCallerID);
  cpgsave();
  cpgslct(iID);
  cpgsci(1);
  cpgsls(1);
  if(iMode&1){
    cpgsch(0.5*iHistogramFontSize+1);
    cpgenv(0,(float)iNBins,0.0,(float)iMaxY,0,1);
  }
  sprintf(str,"Phi");
  cpglab("","",str);
  /*
  cpgsls(4);
  cpgmove(iNBins*0.5,0.0);
  cpgdraw(iNBins*0.5,(float)iMaxY);
  cpgsls(1);
  cpgsci(3);
  cpgmove(fMaxX,0.0);
  cpgdraw(fMaxX,(float)iMaxY);
  cpgsci(1);
  */
  cpgmove(-0.5,piBins[iNBins-1]);
  for(j=0;j<=iNBins;j++){
    cpgdraw((float)(j+0.5),piBins[j]);
  }
  cpgslct(iCallerID);
  cpgunsa();
}


void PlotPhiOD(int iID,int iNBins,int *piBins,int iMaxY,float fMaxX,int iODScoreN,float *pfODScore,float fCBI,float fMI,int iMode){
  int iCallerID;
  int j;
  char str[128];
  float fDum,fWidth,fPlotScaleY;
  int iTotal;

  fDum=(float)iODScoreN/(float)iNBins;
  for(iTotal=j=0;j<iNBins;j++){
    iTotal+=piBins[j];
  }
  if(iTotal==0) iTotal=1;
  fPlotScaleY=((float)iNBins)/((float)iTotal*(float)iODScoreN);
  cpgqid(&iCallerID);
  cpgsave();
  cpgslct(iID);
  cpgsci(1);
  cpgsls(1);
  if(iMode&2) cpgslw(2);
  else cpgslw(1);
  if(iMode&1){
    cpgsch(0.5*iHistogramFontSize+1);
    cpgenv(0,(float)iNBins,0.0,1.05*iMaxY*fPlotScaleY,0,1);
  }
  sprintf(str,"OD Max @ %0.3f CBI %0.3f MI %0.3f",fMaxX/iNBins,fCBI,fMI);
  cpglab("","",str);
  cpgmove(0.5,piBins[0]*fPlotScaleY);
  for(j=1;j<iNBins;j++){
    cpgdraw((float)(j+0.5),piBins[j]*fPlotScaleY);
  }
  fDum=(float)iNBins/(float)iODScoreN;
  fWidth=0.8*fDum;
  cpgsci(2);
  for(j=0;j<iODScoreN;j++){
    cpgmove(fDum*(j+0.5)-fWidth*0.5,0.0);
    cpgdraw(fDum*(j+0.5)-fWidth*0.5,pfODScore[j]);
    cpgdraw(fDum*(j+0.5)+fWidth*0.5,pfODScore[j]);
    cpgdraw(fDum*(j+0.5)+fWidth*0.5,0.0);
  }
  cpgslct(iCallerID);
  cpgunsa();
}

void PlotRho(int iID,int iNBins,int *piBins,int iMaxY,float fMaxX,double dRhoMax,double dTLo,double dTHi,int iMode){
  int iCallerID;
  int j;
  char str[128];

  cpgqid(&iCallerID);
  cpgsave();
  cpgslct(iID);
  cpgsci(1);
  cpgsls(1);
  if(iMode&1){
    cpgsch(0.5*iHistogramFontSize+1);
    cpgenv(0,(float)iNBins,0.0,(float)iMaxY,0,1);
  }
  sprintf(str,"Rho Treshold Lo %0.2f Hi %0.2f",dTLo,dTHi);
  cpglab("","",str);
  cpgmove(0.0,piBins[0]);
  for(j=1;j<iNBins;j++){
    cpgdraw((float)j,piBins[j]);
  }
  if(dRhoNoise){
    cpgsci(4);
    cpgmove((float)(dRhoNoise*(iNBins-1)/dRhoMax),0.0);
    cpgdraw((float)(dRhoNoise*(iNBins-1)/dRhoMax),(float)iMaxY);
    cpgsci(1);
  }

  cpgsci(3);
  cpgmove((float)(dTLo*(iNBins-1)/dRhoMax),0.0);
  cpgdraw((float)(dTLo*(iNBins-1)/dRhoMax),(float)iMaxY);
  cpgsci(2);
  cpgmove((float)(dTHi*(iNBins-1)/dRhoMax),0.0);
  cpgdraw((float)(dTHi*(iNBins-1)/dRhoMax),(float)iMaxY);

  cpgslct(iCallerID);
  cpgunsa();
}


void DisplayFields(int mode,int xd,int yd,double *mx,double *my,double *gx,double *gy){
  int x,y;
  int i,j,k,l,m;
  int callerid;
  float dumf,arm;
  float fmaxx,fminx,fmaxy,fminy,fmax,fminc,fmaxc;
  double rx,ry,phix,phiy;
  unsigned long imax=0,imin=0;

  cpgqid(&callerid);
  cpgslct(fieldid);

  for(j=0;j<yd-1;j++){
    for(i=0;i<xd-1;i++){
      k=i+j*xd;
      l=i+j*(xd-1);
      fbufferdiv[l]=(float)(0.5*(((mx[k+1]-mx[k])+(mx[k+xd+1]-mx[k+xd]))+ 
				 ((my[k+xd]-my[k])+(my[k+xd+1]-my[k+1]))));
      fbufferrot[l]=(float)(0.5*(((my[k+1]-my[k])+(my[k+xd+1]-my[k+xd]))- 
				 ((mx[k+xd]-mx[k])+(mx[k+xd+1]-mx[k+1]))));
    }
  }

  fminx=fminy=10000000000.0;
  fmaxx=fmaxy=-fminx;
  fmax=0.0;
  m=0;
  arm=0.0;
  memset((void *)fbufferdivrot_phi,0,(size_t)(xd-1)*(size_t)(yd-1)*sizeof(float));
  memset((void *)fbufferx,0,(size_t)xd*(size_t)yd*sizeof(float));
  memset((void *)fbuffery,0,(size_t)xd*(size_t)yd*sizeof(float));
  memset((void *)fbufferdivrot_rho,0,(size_t)(xd-1)*(size_t)(yd-1)*sizeof(float));
  for(j=0;j<stenciln;j++){
    i=stencil[j];
    x=i%xd;
    y=i/xd;
    if(x<xd-1 && y<yd-1){
      m++;
      l=x+y*(xd-1);
      fbufferdivrot_phi[l]=(float)atan2(fbufferdiv[l],fbufferrot[l]);
      dumf=fbufferx[l]=(float)fbufferdiv[l];
      if(fmaxx<dumf) fmaxx=dumf;
      if(fminx>dumf) fminx=dumf;
      dumf=fbuffery[l]=(float)fbufferrot[l];
      if(fmaxy<dumf) fmaxy=dumf;
      if(fminy>dumf) fminy=dumf;
      arm+=(fbufferdivrot_rho[l]=dumf=(float)hypot(fbufferdiv[l],fbufferrot[l]));
      if(fmax<dumf){ 
	fmax=dumf;
	imax=i;
      }
    }
  }
  arm/=(float)m;
  printf("DIV(%f %f) ROT(%f %f) %f\n",fminx,fmaxx,fminy,fmaxy,fmax);
  
  cpgsch(5.5);
  cpgpanl(1,1);
  cpgeras();
  if(color_scheme==COLOR_RGB){  
    DrawColorImage(fbufferdivrot_phi,xd-1,yd-1,M_PI,-M_PI,transs);
    DrawColorWedge(wedge_pos, wedge_width, -180.0/harmonic_d, 180.0/harmonic_d,"");
  }
  else{
    DrawNiceColorImagePolar(fbufferdivrot_phi,fbufferdivrot_rho,xd,yd,fmax,0.0,transs);
  }
  cpgpanl(1,2);
  cpgeras();
  DrawGrayImage(fbufferdivrot_rho,xd-1,yd-1,fmax,0.0,transs);

  if(mode>DISPLAY_MODE2){
    cpgpanl(2,1);
    cpgeras();
    DrawGrayImage(fbufferx,xd-1,yd-1,fmaxx,fminx,transs);
    cpgpanl(2,2);
    cpgeras();
    DrawGrayImage(fbuffery,xd-1,yd-1,fmaxy,fminy,transs);
  }

  if(radiusbig!=0){
    fminx=fminy=10000000000.0;
    fmaxx=fmaxy=-fminx;
    fmax=0.0;
    for(j=0;j<yd-1;j++){
      for(i=0;i<xd-1;i++){
	k=i+j*xd;
	l=i+j*(xd-1);
	dumf=fbufferdiv[l]=(float)(0.5*(((gx[k+1]-gx[k])+(gx[k+xd+1]-gx[k+xd]))+ 
					((gy[k+xd]-gy[k])+(gy[k+xd+1]-gy[k+1]))));
	if(fmaxx<dumf) fmaxx=dumf;
	if(fminx>dumf) fminx=dumf;
	dumf=fbufferrot[l]=(float)(0.5*(((gy[k+1]-gy[k])+(gy[k+xd+1]-gy[k+xd]))- 
					((gx[k+xd]-gx[k])+(gx[k+xd+1]-gx[k+1]))));
	if(fmaxy<dumf) fmaxy=dumf;
	if(fminy>dumf) fminy=dumf;
      }
    }
    printf("DIV(%f %f) ROT(%f %f)\n",fminx,fmaxx,fminy,fmaxy);
  
    cpgpanl(1,2);
    DrawGrayImage(fbufferdiv,xd-1,yd-1,fmaxx,fminx,transs);
    cpgpanl(2,2);
    DrawGrayImage(fbufferrot,xd-1,yd-1,fmaxy,fminy,transs);
  }

  if(mode>DISPLAY_MODE4){
    // Calculate max d\phi
    fmaxy=-(fminy=10000000000.0);
    for(j=0;j<yd-1;j++){
      for(i=0;i<xd-1;i++){
	k=i+j*xd;
	l=i+j*(xd-1);

	rx=0.25*(mx[k]+mx[k+1]+mx[k+xd]+mx[k+xd+1]);
	ry=0.25*(my[k]+my[k+1]+my[k+xd]+my[k+xd+1]);
	phix=0.5*(rx*((my[k+1]-my[k])+(my[k+xd+1]-my[k+xd]))-
		  ry*((mx[k+1]-mx[k])+(mx[k+xd+1]-mx[k+xd])));
	phiy=0.5*(rx*((my[k+xd]-my[k])+(my[k+xd+1]-my[k+1]))-
		  ry*((mx[k+xd]-mx[k])+(mx[k+xd+1]-mx[k+1])));

	dumf=fbufferdiv[l]=(float)(sqrt(phix*phix+phiy*phiy)/(rx*rx+ry*ry));
	if(fmaxy<dumf) fmaxy=dumf;
	if(fminy>dumf) fminy=dumf;
      }
    }

    if(cradius!=0){
      memset((void *)chopfbuffer,0,(size_t)xd*(size_t)yd*sizeof(float));
      for(j=0;j<yd-1;j++){
	for(i=0;i<xd-1;i++){
	  chopfbuffer[i+j*xd]=fbufferdiv[i+j*(xd-1)];
	}
      }
      AverageMaps(xd,yd,(double**)0,&chopfbuffer,1,dcradius,0,(int)FLOAT);
      FindExtrema(chopfbuffer,&fminc,&fmaxc,&imin,&imax);
    }
    else{
      fmaxc=fmaxy;
      fminc=fminy;
    }

    printf("D Phi Min %f @ %lu (%f) Max %f @ %lu (%f)\n",fminc,imin,fminy,fmaxc,imax,fmaxy);

    if(xd*yd!=stenciln){
      memset((void *)chopfbuffer,0,(size_t)xd*(size_t)yd*sizeof(float));
      for(j=0;j<yd-1;j++){
	for(i=0;i<xd-1;i++){
	  chopfbuffer[i+j*xd]=fbufferdiv[i+j*(xd-1)];
	}
      }
      memset((void *)fbufferdiv,0,(size_t)(xd-1)*(size_t)(yd-1)*sizeof(float));
           
      //      fmaxy=-(fminy=10000000000.0);
      for(j=0;j<stenciln;j++){
	i=stencil[j];
	x=i%xd;
	y=i/xd;
	if(x<xd-1 && y<yd-1){
	  //	  dumf=
	  fbufferdiv[x+y*(xd-1)]=chopfbuffer[x+y*xd];
	  //	  if(fmaxy<dumf) fmaxy=dumf;
	  //	  if(fminy>dumf) fminy=dumf;
	}
      }
      //      printf("D Phi %f %f\n",fminy,fmaxy);      
    }

    cpgpanl(3,1);
    cpgeras();
    i=pallet;
    pallet=-1;

    // DrawGrayImage(fbufferdiv,xd-1,yd-1,fmaxy*0.01,fminy,transs);
    DrawColorImage(fbufferdiv,xd-1,yd-1,fmaxc,fminc,transs);
    DrawColorWedge(wedge_pos, wedge_width, fminc, fmaxc,"");
    pallet=i;
  }

  cpgslct(callerid);
}

void DrawGrayImage(float *im,int xd,int yd,float max,float min,float *tran){
  int i,x,y;
  float contour;
  float scale,xb,yb;

  if(iDoLargeWedge) cpgsch(fWedgeLargeCH);
  else cpgsch(fWedgeSmallCH);
  cpggray(im,xd,yd,1,xd,1,yd,max,min,tran);
  cpgsci(1);
  cpgscf(iWedgeFont);
  if(iDrawBWWedge) cpgwedg("B", wedge_pos, wedge_width, max, min,"");
  if(do_zerolines){
    cpgsci(0);
    contour=0.0;
    cpgcont(im,xd,yd,1,xd,1,yd,&contour,-1,tran);
  }
  cpgsch(2.5);
  for(i=0;i<apinwheels_n;i++){
    cpgsci(3-(int)apinwheels[i].charge);
    cpgpt1((float)apinwheels[i].x,(float)apinwheels[i].y,2);
  }
  if(do_profiles){
    x=xd/2;
    y=yd/2;
    scale=0.4*fabs((double)(pmaxx-pminx))/(max-min);
    xb=pminx;
    yb=pminy+(pmaxy-pminy)*0.5;
    cpgsci(4);
    cpgmove(xb,yb);
    cpgdraw(pmaxx,yb);
    cpgsci(2);
    cpgmove(xb,yb-scale*im[xd*y]);
    for(i=1;i<xd;i++){
      cpgdraw(xb+i,yb-scale*im[xd*y+i]);    
    }
    scale=0.4*fabs((double)(pmaxy-pminy))/(max-min);
    xb=pminx+(pmaxx-pminx)*0.5;
    yb=pmaxy;
    cpgsci(4);
    cpgmove(xb,yb);
    cpgdraw(xb,pminy);
    cpgsci(3);
    cpgmove(xb+scale*im[x],yb);
    for(i=1;i<yd;i++){
      cpgdraw(xb+scale*im[xd*i+x],yb+i);    
    }
  }
}

void DrawColorImage(float *im,int xd,int yd,float max,float min,float *tran){
  int i;
  float contour;
  float a,b;
  //  printf("DrawColorImage\n");

  if(iNullOutsiders && iCutBlanksOnWedge) 
    MakeImagePallet(PALLET_NULL_OUTSIDERS,fStim,fPreBlank,fPostBlank);
  else Pallet(pallet,1.0,0.5);
  //normal color pallet
  cpgimag(im,xd,yd,1,xd,1,yd,min,max,tran);
  // Streched color pallet
  a=(max-min)/4.0; //fix 4 appropriately
  b=min-a; //fix here too
  //  cpgimag(im,xd,yd,1,xd,1,yd,a*0+b,a*5.58496+b,tran);
  //  cpgimag(im,xd,yd,1,xd,1,yd,a*0.3+b,a*5.3+b,tran);
  cpgsci(1);
  cpgsch(5.5);

  //  DrawColorWedge(wedge_pos, wedge_width, 180.0*min/M_PI/harmonic_d, 180.0*max/M_PI/harmonic_d,"");

  if(do_zerolines){
    cpgsci(0);
    contour=0.0;
    cpgcont(fbufferx,xd,yd,1,xd,1,yd,&contour,-1,tran);
    cpgsci(1);
    contour=0.0;
    cpgcont(fbuffery,xd,yd,1,xd,1,yd,&contour,-1,tran);
  }
  
  cpgsch(2.5);
  for(i=0;i<apinwheels_n;i++){
    cpgsci(3-(int)apinwheels[i].charge);
    cpgpt1((float)apinwheels[i].x,(float)apinwheels[i].y,2);
  }
}

void DrawNiceColorImagePolar(float *phi,float *rho,int xd,int yd,float max,float min,float *tran){
  int i,xi,yi;
  unsigned long xy,ul;
  float contour;
  float dif;
  
  //  printf("DrawNiceColorImagePolar\n");
  //  Pallet(pallet,1.0,0.5);
  //  cpgimag(im,xd,yd,1,xd,1,yd,min,max,tran);

  xy=(unsigned long)xd*(unsigned long)yd;
  dif=2.0*(max-min);
  //  cpgbbuf();
  cpgsch(1);
  for(ul=0;ul<xy;ul++){
    xi=ul%xd;
    yi=ul/xd;
    cpgshls(64,180.0*phi[ul]/M_PI,(rho[ul]-min)/dif,1.0);
    cpgsci(64);
    cpgrect((float)xi-0.5,(float)xi+0.5,(float)yi-0.5,(float)yi+0.5);
  }
  //  cpgebuf();
  Pallet(pallet,1.0,0.5);
  cpgsci(1);
  cpgsch(5.5);

  DrawColorWedge(wedge_pos, wedge_width, -180.0/harmonic_d, 180.0/harmonic_d,"");
  //  cpgwedg("BI", wedge_pos, wedge_width, -180.0/harmonic_d, 180.0/harmonic_d,"");

  if(do_zerolines){
    cpgsci(0);
    contour=0.0;
    cpgcont(fbufferx,xd,yd,1,xd,1,yd,&contour,-1,tran);
    cpgsci(1);
    contour=0.0;
    cpgcont(fbuffery,xd,yd,1,xd,1,yd,&contour,-1,tran);
  }
  
  cpgsch(2.5);
  for(i=0;i<apinwheels_n;i++){
    cpgsci(3-(int)apinwheels[i].charge);
    cpgpt1((float)apinwheels[i].x,(float)apinwheels[i].y,2);
  }
}

void DrawNiceColorImageCartesian(float *mx,float *my,int xd,int yd,float *tran){
  int i,xi,yi;
  unsigned long xy,ul;
  float contour;
  float hexh;
  float x,y;
  float dumf,sqrt3,sqrt2,hls;
  float r,g,b,h,l;
  float *pmx=NULL,*pmy=NULL;
  double dChopMin,dChopMax;
  double dR,dP,dRMax;
#ifdef ENUMERATE_SQUARES
  char str[5];
#endif

  //  printf("DrawNiceColorImageCartesian\n");

  xy=(unsigned long)xd*(unsigned long)yd;
  sqrt3=(float)sqrt(3.0);
  sqrt2=(float)sqrt(2.0);



  if(cradius){
    pmx=(float*)malloc(xy*sizeof(float));
    pmy=(float*)malloc(xy*sizeof(float));
    memcpy(pmx,mx,xy*sizeof(float));
    memcpy(pmy,my,xy*sizeof(float));
    AverageMaps(xd,yd,(double**)0,&pmx,1,dcradius,0,(int)FLOAT);
    AverageMaps(xd,yd,(double**)0,&pmy,1,dcradius,0,(int)FLOAT);
  }
  else{
    pmx=mx;
    pmy=my;
  }

  hexh=0.0;
  if(amplitude_chop_l!=0.0 || amplitude_chop_u!=1.0){
    dRMax=0.0;
    for(ul=0;ul<xy;ul++){
      if(dRMax<(dR=(double)(pmx[ul])*(double)(pmx[ul])+(double)(pmy[ul])*(double)(pmy[ul])))
	dRMax=dR;
    }
    dRMax=sqrt(dRMax);
    dChopMin=amplitude_chop_l*dRMax;
    dChopMax=(amplitude_chop_u-amplitude_chop_l)*dRMax;

    for(ul=0;ul<xy;ul++){
      if((dR=hypot((double)(mx[ul]),(double)(my[ul]))-dChopMin)<=0.0){
	mx[ul]=0.0;
	my[ul]=0.0;
      }
      else{
	if(dR>dChopMax) dR=dChopMax;
	dP=atan2((double)(my[ul]),(double)(mx[ul]));
	mx[ul]=dR*cos(dP);
	my[ul]=dR*sin(dP);
      }
    }
    for(ul=0;ul<xy;ul++){
      x=sqrt3*(float)fabs((double)(mx[ul]));
      y=(float)fabs((double)(my[ul]));
      if(y<x) dumf=0.5*(x+y);
      else dumf=y;
      if(hexh<dumf) hexh=dumf;
    }
  }
  else{
    for(ul=0;ul<xy;ul++){
      x=sqrt3*(float)fabs((double)(pmx[ul]));
      y=(float)fabs((double)(pmy[ul]));
      if(y<x) dumf=0.5*(x+y);
      else dumf=y;
      if(hexh<dumf) hexh=dumf;
    }
    dChopMax=0.0;
    for(ul=0;ul<xy;ul++){
      if(dChopMax<(dR=(double)(pmx[ul])*(double)(pmx[ul])+(double)(pmy[ul])*(double)(pmy[ul])))
	dChopMax=dR;
    }
    dChopMax=sqrt(dChopMax);
    if(cradius){
      free(pmx);
      free(pmy);
    }
  }
  //  printf("DIS HEXH=%f\n",hexh);

  hexh=1.0/(sqrt2*hexh);
  hls=180.0/M_PI;
  dChopMax*=2.0;

  //  cpgbbuf();
  cpgsch(1.5);
  for(ul=0;ul<xy;ul++){
    xi=ul%xd;
    yi=ul/xd;
  
    // Val: hexagon
    /*
    x=mx[ul]*hexh;
    y=my[ul]*hexh;
    if(y>0.0){
      if((dumf=sqrt3*x+y)>0.0){
	// No green
	r=sqrt2*y;
	g=0.0;
	b=dumf/sqrt2;
      }
      else{
	// No blue
	r=(y-x*sqrt3)/sqrt2;
	g=-dumf/sqrt2;
	b=0.0;
      }
    }
    else{
      if((dumf=sqrt3*x-y)>0.0){
	// No red
	r=0.0;
	g=-sqrt2*y;
	b=dumf/sqrt2;
      }
      else{
	// No blue
	r=-dumf/sqrt2;
	g=-(y+x*sqrt3)/sqrt2;
	b=0.0;
      }
    }
    //    printf("R%f G%f B%f\n",r,g,b);
    if(r>1.0) r=1.0;
    if(g>1.0) g=1.0;
    if(b>1.0) b=1.0;

    //Val: white background
    if(iWhiteBackgroundOnPolarPlot){
      if(mx[ul]==0.0 && my[ul]==0.0){
	r=1.0;
	g=1.0;
	b=1.0;
      }
    }
    cpgscr(64,r,g,b);
    */

  // HLS
    
    h=atan2(my[ul],mx[ul])*hls;
    if(h<0.0) h+=360.0;
    l=hypot(my[ul],mx[ul])/dChopMax;
    cpgshls(64,h,l,1.0);
    

    cpgsci(64);
    cpgrect((float)xi-0.5,(float)xi+0.5,(float)yi-0.5,(float)yi+0.5);
#ifdef ENUMERATE_SQUARES
    cpgsci(1);
    sprintf(str,"%lu",ul);
    cpgtext((float)xi,(float)yi,str);
#endif
  }
  //  cpgebuf();
  
  Pallet(PALLET_HSL6,1.0,0.5); //This routine does not use pallet, color scheme=PALLET_HSL6
  cpgsci(1);
  cpgsch(5.5);

  //Moved to DisplayMap
  //DrawColorWedge(wedge_pos, wedge_width, -180.0/harmonic_d, 180.0/harmonic_d,"");

  //  cpgwedg("BI", wedge_pos, wedge_width, -180.0/harmonic_d, 180.0/harmonic_d,"");
  
  if(do_zerolines){
    cpgsci(0);
    contour=0.0;
    cpgcont(fbufferx,xd,yd,1,xd,1,yd,&contour,-1,tran);
    cpgsci(1);
    contour=0.0;
    cpgcont(fbuffery,xd,yd,1,xd,1,yd,&contour,-1,tran);
  }
  
  cpgsch(2.5);
  for(i=0;i<apinwheels_n;i++){
    cpgsci(3-(int)apinwheels[i].charge);
    cpgpt1((float)apinwheels[i].x,(float)apinwheels[i].y,2);
  }
}

void DrawContours(int mode,int scheme,float **sb,int xd,int yd,float fmax,float fmin,float threshold_pcnt,float *tran){
  unsigned long i,j,xy;
  float fminc;
  static int iReplotCount=1;//Val: change to 0 later
  PAIR *pPAIR;
  int iSaveContourLinesWidth;
  float fContourSpacingAngleRad;
  int iOldLineStyle;

  xy=(unsigned long)xd*(unsigned long)yd;
  fminc=fmin+threshold_pcnt*(fmax-fmin);
  cpgqlw(&iSaveContourLinesWidth);
  cpgslw(iContourLinesWidth);
  switch(scheme){
  case CONTOUR_NON:
    break;
  case CONTOUR_PHI:      
    if(fContourSpacingAngle==0.0) fContourSpacingAngleRad=2.0*M_PI/depth_lines;
    else fContourSpacingAngleRad=fContourSpacingAngle*M_PI/180.0;

    for(i=0;i<depth_lines;i++){ 
      //      contours[i]=contours_iniangle+(float)i*fContourSpacingAngleRad;
      contours[i]=contours_iniangle+(2*((int)i&1)-1)*(((int)i+1)/2)*fContourSpacingAngleRad;
      contours[i]-=2.0*M_PI*(float)floor(contours[i]/(2.0*M_PI));
      if(contours[i]>M_PI) contours[i]-=2.0*M_PI;
    }

    if(do_smartcontours){
      /*
      // Plots one contour at a time (memory efficient)
      for(i=0;i<depth_lines;i++){
	SmartContoursCartesian(xd,yd,sb[0],sb[1],&smart_contours,&smart_contoursn,contours[i],fminc);
	if(do_verbose) printf("DIS Contour %li %f %i\n",i,contours[i]*180.0/M_PI,smart_contoursn);
	if(mode==CONTOUR_RHO){
	  if(emphasize_contours_iniangle && !i) cpgsci(4);
	  else cpgsci(2);
	}
	else{
	  if(emphasize_contours_iniangle && !i) cpgsci(0);
	  else cpgsci(1);
	  if(save_on) cpgsci(0);
	}
	for(j=0;j<smart_contoursn;j++){
	  if(do_verbose) printf("DIS %li %lu (%f %f) (%f %f)\n",j,smart_contours[j].n,smart_contours[j].x1,smart_contours[j].y1,smart_contours[j].x2,smart_contours[j].y2);
	  cpgmove(smart_contours[j].x1,smart_contours[j].y1);
	  cpgdraw(smart_contours[j].x2,smart_contours[j].y2);
	}
      }
      */

      // Plots all contours at once, saves all contours too (flexible)
      if(iReplotSavedContours!=1 || !iReplotCount){
	piSmartContoursN=(int*)realloc(piSmartContoursN,depth_lines*sizeof(int));
	memset(piSmartContoursN,0,depth_lines*sizeof(int));
	ppPAIRSmartContours=(PAIR**)realloc(ppPAIRSmartContours,depth_lines*sizeof(PAIR*));
	memset(ppPAIRSmartContours,0,depth_lines*sizeof(PAIR*));
	
	for(i=0;i<depth_lines;i++){
	  SmartContoursCartesian(xd,yd,sb[0],sb[1],ppPAIRSmartContours+i,piSmartContoursN+i,contours[i],fminc);
	  if(do_verbose) printf("DIS Contour %li %f %i\n",i,contours[i]*180.0/M_PI,piSmartContoursN[i]);
	}
      }
      
      for(i=0;i<depth_lines;i++){
	cpgqls(&iOldLineStyle);
	if(!(pPAIR=ppPAIRSmartContours[i])) continue;
	if(mode==CONTOUR_RHO){
	  if(emphasize_contours_iniangle && !i) cpgsci(4);
	  else cpgsci(2);
	}
	else{
	  if(emphasize_contours_iniangle){
	    if(!i){
	      cpgsci(0);//0
	      cpgsls(1);//1
	    }
	    else{
	      cpgsci(1);//1
	      cpgsls(1);//2
	    }
	  }
	  else{
	    cpgsci(0);
	    cpgsls(1);
	  }
	  //	  if(save_on) cpgsci(1);
	  if(iOverlayContoursColor) cpgsci(3);
	}
	if(iColorContours) cpgsci(i+1);
	for(j=0;j<piSmartContoursN[i];j++){
	  if(do_verbose) 
	    printf("DIS %li %lu (%f %f) (%f %f)\n",j,pPAIR[j].n,
		   pPAIR[j].x1,pPAIR[j].y1,
		   pPAIR[j].x2,pPAIR[j].y2);
	  cpgmove(pPAIR[j].x1,pPAIR[j].y1);
	  cpgdraw(pPAIR[j].x2,pPAIR[j].y2);
	}
	cpgsls(iOldLineStyle);
      }

      if(iPlotExternalContours){
	for(i=0;i<iNDepthLinesExternal;i++){
	  cpgqls(&iOldLineStyle);
	  if(!(pPAIR=ppPAIRSmartContoursExternal[i])) continue;
	  if(mode==CONTOUR_RHO){
	    if(emphasize_contours_iniangle && !i) cpgsci(4);
	    else cpgsci(2);
	  }
	  else{
	    if(emphasize_contours_iniangle){
	      if(!i){
		cpgsci(iContoursColorExternalInitial);
		cpgsls(1);
	      }
	      else{ 
		cpgsci(iContoursColorExternal);
		cpgsls(1);
	      }
	    }
	    else{
	      cpgsci(iContoursColorExternal);
	      cpgsls(1);
	    }
	    // if(save_on) cpgsci(iContoursColorExternalSave);
	    if(iOverlayContoursColor) cpgsci(2);
	  }
	  for(j=0;j<piSmartContoursNExternal[i];j++){
	    if(do_verbose) 
	      printf("DIS %li %lu (%f %f) (%f %f)\n",j,pPAIR[j].n,
		     pPAIR[j].x1,pPAIR[j].y1,
		     pPAIR[j].x2,pPAIR[j].y2);
	    cpgmove(pPAIR[j].x1,pPAIR[j].y1);
	    cpgdraw(pPAIR[j].x2,pPAIR[j].y2);
	  }
	  cpgsls(iOldLineStyle);
	}
      }
    }
    else{
      cpgsci(1);
      cpgcont(sb[3],xd,yd,1,xd,1,yd,contours,depth_lines,tran);
      if(emphasize_contours_iniangle){ 
	cpgsci(0);
	cpgcont(sb[3],xd,yd,1,xd,1,yd,contours,-1,tran);
      }
    }
    
    break;
  case CONTOUR_RHO:
    // Val: uncomment later
    
    for(i=0;i<depth_lines;i++){ 
      contours[i]=fminc+(float)(i)*(fmax-fminc)/(float)(depth_lines);
    }
    cpgsci(1);
    cpgcont(smoothfbuffers[2],xd,yd,1,xd,1,yd,contours,depth_lines,tran);
    
    //Val: draws only one contour at 50%
    //    contours[0]=fminc+(fmax-fminc)*0.5;
    /* 
   contours[0]=fmax*0.34;
    cpgsci(1);
    cpgslw(1);
    cpgcont(smoothfbuffers[2],xd,yd,1,xd,1,yd,contours,1,tran);
    cpgslw(1);
    */
    break;
  default:
    printf("DIS Bad contour scheme %i for mode %i\n",scheme,mode);
  }    
  iReplotCount++; 
  cpgslw(iSaveContourLinesWidth);
}

#define SECTION_WIN_SIZE_X 10.0
#define SECTION_WIN_SIZE_Y 8.0
//#define N_SAMPLES_PER_PIXEL 0.707
#define N_SAMPLES_PER_PIXEL 2.1
#define DEFAULT_N_FIT_POLYNOMIALS 4
#define N_MC_TRIALS 5000

void PlotSection(int iMode,int iXD,int iYD,float *pfMapX,float *pfMapY){
  double ddx,ddy,ddl,dGradX,dGradY;
  unsigned long nNSamples;
  float *pfPhiSamples=NULL,*pfRhoSamples=NULL;
  float *pfDerivative=NULL;
  int iNPolynomials=DEFAULT_N_FIT_POLYNOMIALS;
  double *pdPolyFitCoeff=(double*)0;
  double *pdDerPolyFitCoeff=(double*)0;
  unsigned long i,j;
  double dDerCoeff=1.0,dPhiCoeff=1.0,dPhiOffset=0.0;
  float fX,fY;
  float fPhiMin,fPhiMax;
  float fRhoMin,fRhoMax;
  float fDum;
  double dDum,dDum1;
  char ch;
  int iCallerID,iSectionID=0;
  unsigned long ulLeftExtent,ulTopExtent,ulRightExtent,ulBottomExtent;
  unsigned long ulNX,ulNY;
  unsigned long ulArea;
  long lXi,lYi;
  long lX,lY;
  int iX,iY;
  double dX,dY;
  double dXc,dYc;
  double dWeight,dFieldX,dFieldY,dFieldR;
  double dRadius2,dArea;
  int iDum;
  double dRX,dRY;
  double dRhoA2;
  double dPhiA;
  double dXA,dXXA,dYA,dYYA,dXYA;
  double dLA1,dLA2,dLAA;
  double dPolyA1=0.0,dPolyA2=0.0;
  FILE *pF=NULL;
  char strFileName[256],strDum[128];
  int iSaveOptions=SECTION_SAVE_RHO|SECTION_SAVE_PHI;
  static float *pfASamples=NULL;
  static int iNASamples=0,iNA=0;
  float *pf;
  int iSegmentSaved=0;
  int iUseFitMinMaxForPhiDirPlot=0;

  cpgqid(&iCallerID);

  printf("DIS Section B %f %f E %f %f\n",section[0].x,section[0].y,section[1].x,section[1].y);
    
  if(do_phase_transform){
    dPhiCoeff=phase_transform;
    dPhiOffset=0.0;
    dDerCoeff=phase_transform*M_PI/180.0;
  }
  else{ 
    if(iDoWedgeAnnotation){
      dPhiCoeff=(fWedgeMax-fWedgeMin)/(2*M_PI);
      dPhiOffset=fWedgeMin;
      dDerCoeff=dPhiCoeff;
    }
    else{
      dPhiCoeff=1.0;
      dPhiOffset=0.0;
      dDerCoeff=1.0;
    }
  }
  

  dGradX=section[1].x-section[0].x;
  dGradY=section[1].y-section[0].y;
  nNSamples=(unsigned long)floor(hypot(dGradX,dGradY)*N_SAMPLES_PER_PIXEL);
  ddx=dGradX/nNSamples;
  ddy=dGradY/nNSamples;
  ddl=hypot(ddx,ddy);
  if(do_pixel_size_transform) ddl*=dPixelSizeU*0.001;
  nNSamples++;
    
  if(!(pfPhiSamples=(float*)malloc(nNSamples*sizeof(float)))){
    printf("DIS Cannot allocate for pfPhiSamples");
    return;
  }
  if(!(pfRhoSamples=(float*)malloc(nNSamples*sizeof(float)))){
    printf("DIS Cannot allocate for pfRhoSamples");
    return;
  }
  if(!(pfDerivative=(float*)malloc((nNSamples-1)*sizeof(float)))){
    printf("DIS Cannot allocate for pfDerivative");
    return;
  }
  
  dRadius2=dSectionRadius*dSectionRadius;

  dXA=dXXA=dYA=dYYA=dXYA=0.0;
  fPhiMax=-(fPhiMin=BIG_FLOAT);
  fRhoMax=-(fRhoMin=BIG_FLOAT);
  for(i=0,dGradX=section[0].x,dGradY=section[0].y;i<nNSamples;i++,dGradX+=ddx,dGradY+=ddy){
    dX=dGradX-floor(dGradX);
    dY=dGradY-floor(dGradY);
    lXi=(long)floor(dGradX);
    lYi=(long)floor(dGradY);

    ulLeftExtent=1+(unsigned long)floor(dSectionRadius-dX);
    ulTopExtent=1+(unsigned long)floor(dSectionRadius-dY);
    ulRightExtent=1+(unsigned long)floor(dSectionRadius-1.0+dX);
    ulBottomExtent=1+(unsigned long)floor(dSectionRadius-1.0+dY);

    //  printf("DIS %lu-1-%lu  %lu-1-%lu\n",ulLeftExtent,ulRightExtent,ulTopExtent,ulBottomExtent);

    ulNX=ulLeftExtent+ulRightExtent+1;
    ulNY=ulTopExtent+ulBottomExtent+1;

    if(ulNX*ulNY==1 || iSectionSinglePixel){
      j=lXi+iXD*lYi;
      dFieldX=pfMapX[j];
      dFieldY=pfMapY[j]; 
      dFieldR=hypot(pfMapX[j],pfMapY[j]);
      goto skiplabel;
    }

    lXi-=(long)ulLeftExtent;
    lYi-=(long)ulTopExtent;

    dXc=dX+(double)ulLeftExtent;
    dYc=dY+(double)ulTopExtent;

    dArea=dFieldX=dFieldY=dFieldR=0.0;
    for(iY=0;iY<ulNY;iY++){
      lY=lYi+iY;
      if(lY<0 || lY>=iYD) continue;
      lY*=iXD;
      dY=(double)iY-dYc;
      for(iX=0;iX<ulNX;iX++){
	lX=lXi+iX;
	if(lX<0 || lX>=iXD) continue;
	
	if((iDum=IsUnitSquareInsideCircle(iX,iY,dXc,dYc,dRadius2))==4){
	  dArea+=1.0;
	  dFieldX+=pfMapX[lY+lX];
	  dFieldY+=pfMapY[lY+lX];
	  dFieldR+=hypot(pfMapX[lY+lX],pfMapY[lY+lX]);
	}
	else{
	  if(iDum || IsPartOfUnitSquareInsideCircle(iX,iY,dXc,dYc,dRadius2)){
	    dX=(double)iX-dXc;

	    for(j=0,ulArea=0;j<N_MC_TRIALS;j++){
	      dRX=dX+ran3(&lSeed);
	      dRY=dY+ran3(&lSeed);
	      if(dRX*dRX+dRY*dRY<=dRadius2) ulArea++;
	    }
	    dArea+=(dWeight=(double)ulArea/N_MC_TRIALS);
	    dFieldX+=dWeight*pfMapX[lY+lX];
	    dFieldY+=dWeight*pfMapY[lY+lX]; 
	    dFieldR+=dWeight*hypot(pfMapX[lY+lX],pfMapY[lY+lX]);
	  }
	}
      }
    }
    dFieldX/=dArea;
    dFieldY/=dArea;   
    dFieldR/=dArea;   

  skiplabel:
    dXA+=dFieldX;
    dXXA+=dFieldX*dFieldX;
    dYA+=dFieldY;
    dYYA+=dFieldY*dFieldY;
    dXYA+=dFieldX*dFieldY;
    pfPhiSamples[i]=fDum=(float)atan2(dFieldY,dFieldX);
    if(fDum<0.0) pfPhiSamples[i]=fDum=2*M_PI+fDum;
    if(fPhiMax<fDum) fPhiMax=fDum;
    if(fPhiMin>fDum) fPhiMin=fDum;
    //    pfRhoSamples[i]=fDum=(float)(hypot(dFieldY,dFieldX));
    pfRhoSamples[i]=fDum=(float)(dFieldR);
    if(fRhoMax<fDum) fRhoMax=fDum;
    if(fRhoMin>fDum) fRhoMin=fDum;
  }
  dXA/=nNSamples;
  dXXA/=nNSamples;
  dYA/=nNSamples;
  dYYA/=nNSamples;
  dXYA/=nNSamples;
  dXXA-=dXA*dXA;
  dYYA-=dYA*dYA;
  dXYA-=dXA*dYA;
  dLA1=sqrt((dXXA-dYYA)*(dXXA-dYYA)+4.0*dXYA*dXYA);
  dLA2=0.5*(dXXA+dYYA-dLA1);
  dLA1=0.5*(dXXA+dYYA+dLA1);
  dLAA=sqrt(sqrt(dLA1*dLA2));
  dRhoA2=sqrt(dXXA+dYYA);

  printf("DIS X Mean %f Dev %f Y Mean %f Dev %f XY %f\n",dXA,dXXA,dYA,dYYA,dXYA);
  printf("DIS Lambda %f(%f) %f(%f) Aver %f\n",dLA1,sqrt(dLA1),dLA2,sqrt(dLA2),dRhoA2);

  if(dRhoNoise>0.0){
    dPhiA=0.0;
    for(i=0;i<nNSamples;i++){
      dPhiA+=asin(dRhoNoise/pfRhoSamples[i]);
    }
    dPhiA/=nNSamples;
    printf("DIS Shot noise Phi %f(%f) RhoShot %f\n",dPhiA,dPhiA*dPhiCoeff,dRhoNoise);
  }

  if(!(pdPolyFitCoeff=(double *)calloc(iNPolynomials,sizeof(double)))){
    printf("DIS Cannot allocate for pdPolyFitCoeff\n");
  }

  if(!(pdDerPolyFitCoeff=(double *)calloc(iNPolynomials,sizeof(double)))){
    printf("DIS Cannot allocate for pdDerPolyFitCoeff\n");
  }
    
  dDum=(double)(fPhiMax=fPhiMin=pfPhiSamples[0]);
  PolynomialFitM(nNSamples,iNPolynomials,0.0,&dDum,pdPolyFitCoeff,0);
  for(i=1;i<nNSamples;i++){
    if(!do_discontinuous_phase && abs(pfPhiSamples[i]-pfPhiSamples[i-1])>M_PI){
      if(pfPhiSamples[i]>pfPhiSamples[i-1]) pfPhiSamples[i]-=2*M_PI;
      else  pfPhiSamples[i]+=2*M_PI;
    }
    dDum=(double)(fDum=pfPhiSamples[i]);
    if(fPhiMax<fDum) fPhiMax=fDum;
    if(fPhiMin>fDum) fPhiMin=fDum;
    PolynomialFitM(nNSamples,iNPolynomials,(double)(i),&dDum,pdPolyFitCoeff,0);
  }
    /*    
    for(i=0;i<nNSamples-1;i++){
      dDum=dDerCoeff*((double)(pfPhiSamples[i+1])-(double)(pfPhiSamples[i]))/ddl;
      PolynomialFitM(nNSamples-1,iNPolynomials,(double)(i),&dDum,pdDerPolyFitCoeff,0);
    }
    */

  // Deviation from fit
  if(pdPolyFitCoeff){
    for(i=0;i<nNSamples;i++){
      PolynomialFitM(nNSamples,iNPolynomials,(double)(i),&dDum,pdPolyFitCoeff,1);
      dDum-=pfPhiSamples[i];
      dPolyA1+=dDum;
      dPolyA2+=dDum*dDum;
    }
    dPolyA1/=nNSamples;
    dPolyA2/=nNSamples;
    dPolyA2-=dPolyA1*dPolyA1;
    dPolyA2=sqrt(dPolyA2);
    printf("DIS Phi Dev %f(%f)\n",dPolyA2,dPolyA2*dPhiCoeff);
    if(fabs(dPolyA1)>0.01*dPolyA2){
      printf("DIS Large mean %g\n",dPolyA1);
    }
  }
  

  iSectionID=cpgopen("/xw");
  cpgask(0);
  cpgsubp(1,3);
  cpgpap(SECTION_WIN_SIZE_X,SECTION_WIN_SIZE_Y/SECTION_WIN_SIZE_X);
  cpgsch(1.0);


  // Plot phase
  if(fPhiMin==fPhiMax){
    fPhiMin-=0.5;
    fPhiMax+=0.5;
  }
  cpgsci(1);
  cpgenv(0,(nNSamples-1)*ddl,fPhiMin*dPhiCoeff+dPhiOffset,fPhiMax*dPhiCoeff+dPhiOffset,0,1);
  cpglab("","","Phi");
  cpgbbuf();
  
  if(pdPolyFitCoeff){
    cpgsci(1);
    PolynomialFitM(nNSamples,iNPolynomials,0.0,&dDum,pdPolyFitCoeff,1);
    cpgmove(0.0,(float)(dDum*dPhiCoeff+dPhiOffset));
    for(i=1;i<nNSamples;i++){
      PolynomialFitM(nNSamples,iNPolynomials,(double)(i),&dDum,pdPolyFitCoeff,1);
      cpgdraw((float)(i*ddl),(float)(dDum*dPhiCoeff+dPhiOffset));
    }
  }
  
  cpgsci(4);
  cpgmove(0.0,pfPhiSamples[0]*dPhiCoeff+dPhiOffset);
  for(i=1;i<nNSamples;i++){
    cpgdraw((float)(i*ddl),pfPhiSamples[i]*dPhiCoeff+dPhiOffset);
  }
  //  printf("DIS Coe=%f Off=%f\n",dPhiCoeff,dPhiOffset);
  cpgsci(2);
  for(i=0;i<nNSamples;i++){
    cpgpt1((float)(i*ddl),pfPhiSamples[i]*dPhiCoeff+dPhiOffset,-1);
  }
  //Check
  /*
  for(i=0,dGradX=section[0].x,dGradY=section[0].y;i<nNSamples;i++,dGradX+=ddx,dGradY+=ddy){
    j=(unsigned long)floor(dGradX)+iXD*(unsigned long)floor(dGradY);
    pfPhiSamples[i]=fbufferphi[j]*dCoeff;
  }
  for(i=1;i<nNSamples;i++){
    if(abs(pfPhiSamples[i]-pfPhiSamples[i-1])>180.0){
      if(pfPhiSamples[i]>pfPhiSamples[i-1]) pfPhiSamples[i]-=360.0;
      else  pfPhiSamples[i]+=360.0;
    }
  }
  
  cpgsci(6);
  cpgmove(0.0,pfPhiSamples[0]);
  for(i=1;i<nNSamples;i++){
    cpgdraw((float)(i*ddl),pfPhiSamples[i]);
  }
  */
  cpgebuf();

  // Calculate derivative of the polynomial fit
    
  fPhiMax=-(fPhiMin=BIG_FLOAT);

  PolynomialFitM(nNSamples,iNPolynomials,0.0,&dDum1,pdPolyFitCoeff,1);

  for(i=1;i<nNSamples;i++){
    PolynomialFitM(nNSamples,iNPolynomials,(double)(i),&dDum,pdPolyFitCoeff,1);
    fDum=pfDerivative[i-1]=dDerCoeff*(dDum-dDum1)/ddl;
    if(fPhiMax<fDum) fPhiMax=fDum;
    if(fPhiMin>fDum) fPhiMin=fDum;
    dDum1=dDum;
  }
  if(!iUseFitMinMaxForPhiDirPlot){
    fPhiMax=-(fPhiMin=BIG_FLOAT);
    for(i=1;i<nNSamples;i++){
      fDum=(float)(dDerCoeff*(pfPhiSamples[i]-pfPhiSamples[i-1])/ddl);
      if(fPhiMax<fDum) fPhiMax=fDum;
      if(fPhiMin>fDum) fPhiMin=fDum;
    }
  }


  // Plot phase derivative
  if(fPhiMin==fPhiMax){
    fPhiMin-=0.5;
    fPhiMax+=0.5;
  }
  cpgsci(1);
  cpgenv(0.0,(float)(nNSamples*ddl),fPhiMin,fPhiMax,0,1);
  cpglab("Cortex mm","Magnification factor deg/mm","dPhi/dl");
  cpgbbuf();
  /*
  if(pdDerPolyFitCoeff){
    cpgsci(8);
    PolynomialFitM(nNSamples-1,iNPolynomials,0.0,&dDum,pdDerPolyFitCoeff,1);
    cpgmove((float)(0.5*ddl),(float)(dDum));
    for(i=1;i<nNSamples-1;i++){
      PolynomialFitM(nNSamples-1,iNPolynomials,(double)(i),&dDum,pdDerPolyFitCoeff,1);
      cpgdraw((float)(((double)(i)+0.5)*ddl),(float)(dDum));
    }
  }
  */

  // Plot derivative of the polynomial fit
  cpgsci(5);
  cpgmove((float)(0.5*ddl),pfDerivative[0]);
  for(i=1;i<nNSamples-1;i++){
    cpgdraw( (float)(((double)(i)+0.5)*ddl),pfDerivative[i]);
  }
  
  // Plot derivative 
  cpgsci(7);
  cpgmove((float)(0.5*ddl),(float)(dDerCoeff*(pfPhiSamples[1]-pfPhiSamples[0])/ddl));
  for(i=1;i<nNSamples-1;i++){
    cpgdraw((float)(((double)(i)+0.5)*ddl),(float)(dDerCoeff*(pfPhiSamples[i+1]-pfPhiSamples[i])/ddl));
  }
  cpgsci(2);
  for(i=0;i<nNSamples-1;i++){
    cpgpt1((float)(((double)(i)+0.5)*ddl),(float)(dDerCoeff*(pfPhiSamples[i+1]-pfPhiSamples[i])/ddl),-1);
  }
  
  cpgebuf(); 
    
  // Plot absolute value   
  if(fRhoMin==fRhoMax){
    fRhoMin-=0.5;
    fRhoMax+=0.5;
  }
  cpgsci(1);
  if(iUseRhoMinOnSectionAmplitudePlot) cpgenv(0,(nNSamples-1)*ddl,fRhoMin,fRhoMax,0,1);
  else cpgenv(0,(nNSamples-1)*ddl,0.0/*fRhoMin*/,fRhoMax,0,1);
  cpglab("","","Rho");
  cpgbbuf();
  // Draw deviation
  cpgsci(7);
  cpgmove(0.0,(float)dRhoA2);
  cpgdraw((float)((nNSamples-1)*ddl),(float)dRhoA2);
  
  cpgsci(3);
  cpgmove(0.0,pfRhoSamples[0]);
  for(i=1;i<nNSamples;i++){
    cpgdraw((float)(i*ddl),pfRhoSamples[i]);
  }
  cpgsci(2);
  for(i=0;i<nNSamples;i++){
    cpgpt1((float)(i*ddl),pfRhoSamples[i],-1);
  }
  //Check
  /*
  for(i=0,dGradX=section[0].x,dGradY=section[0].y;i<nNSamples;i++,dGradX+=ddx,dGradY+=ddy){
    j=(unsigned long)floor(dGradX)+iXD*(unsigned long)floor(dGradY);
    pfRhoSamples[i]=fbufferrho[j];
  }
  cpgsci(4);
  cpgmove(0.0,pfRhoSamples[0]);
  for(i=1;i<nNSamples;i++){
    cpgdraw((float)(i*ddl),pfRhoSamples[i]);
  }
  */
  cpgebuf();

  while(1){  
    cpgband(0,0,fX,fY,&fX,&fY,&ch); 
    if(ch=='Q' || ch=='q' || ch=='D') break;
    if(ch=='N'){
      dRhoNoise=dRhoA2;
    }

    if(ch=='s'){
      if(iSegmentSaved) continue;
      if(!(pLineSegments=(PAIR*)realloc(pLineSegments,(iNLineSegments+1)*sizeof(PAIR)))){
	printf("DIS Cannot allocate for pLineSegments\n");
	continue;
      }
      pLineSegments[iNLineSegments].n=iNLineSegments;
      pLineSegments[iNLineSegments].x1=section[0].x;
      pLineSegments[iNLineSegments].y1=section[0].y;
      pLineSegments[iNLineSegments].x2=section[1].x;
      pLineSegments[iNLineSegments].y2=section[1].y;
      iNLineSegments++;
      iSegmentSaved=1;
    }

    // Save section
    if(ch=='S'){
      if(iSaveOptions&SECTION_SAVE_PHI){
	sprintf(strFileName,"%s_sectionPhi%i-%i_%i-%i",filenames[0],
		(int)section[0].x,(int)section[0].y,(int)section[1].x,(int)section[1].y);
	if(!iSectionSinglePixel){
	  sprintf(strDum,"_%.1f",dSectionRadius);
	  strcat(strFileName,strDum);
	}
	printf("INT Saving samples in %s\n",strFileName);
	if((pF=fopen(strFileName,"w"))==NULL){
	  printf("INT Cannot open file %s\n",strFileName);
	  goto bailout1;
	}
	sprintf(strDum,"#BEGIN SAMPLE SECTION\n");
	if(fprintf(pF,strDum) != strlen(strDum)){
	  printf("INT Cannot write %s to file %s\n",strDum,strFileName);
	  goto bailout1;
	}
	sprintf(strDum,"#%li\n",nNSamples);
	if(fprintf(pF,strDum) != strlen(strDum)){
	  printf("INT Cannot write nNSamples to file %s\n",strFileName);
	  goto bailout1;
	}
	for(i=0;i<nNSamples;i++){
	  sprintf(strDum,"%f %f\n",ddl*i,pfPhiSamples[i]*dPhiCoeff+dPhiOffset);
	  if(fprintf(pF,strDum) != strlen(strDum)){
	    printf("INT Cannot write sample %li to file %s\n",i,strFileName);
	    goto bailout1;
	  }
	}
	sprintf(strDum,"#END SAMPLE SECTION\n");
	if(fprintf(pF,strDum) != strlen(strDum)){
	  printf("INT Cannot write %s to file %s\n",strDum,strFileName);
	  goto bailout1;
	}
      bailout1:
	fclose(pF);
      }
      
      if(iSaveOptions&SECTION_SAVE_PHI_FIT){
	sprintf(strFileName,"%s_sectionPhiFit%i-%i_%i-%i",filenames[0],
		(int)section[0].x,(int)section[0].y,(int)section[1].x,(int)section[1].y);
	if(!iSectionSinglePixel){
	  sprintf(strDum,"_%.1f",dSectionRadius);
	  strcat(strFileName,strDum);
	}
	printf("INT Saving samples in %s\n",strFileName);
	if((pF=fopen(strFileName,"w"))==NULL){
	  printf("INT Cannot open file %s\n",strFileName);
	  goto bailout2;
	}
	sprintf(strDum,"#BEGIN SAMPLE SECTION\n");
	if(fprintf(pF,strDum) != strlen(strDum)){
	  printf("INT Cannot write %s to file %s\n",strDum,strFileName);
	  goto bailout2;
	}
	sprintf(strDum,"#%li\n",nNSamples);
	if(fprintf(pF,strDum) != strlen(strDum)){
	  printf("INT Cannot write nNSamples to file %s\n",strFileName);
	  goto bailout2;
	}
	for(i=0;i<nNSamples;i++){
	  PolynomialFitM(nNSamples,iNPolynomials,(double)(i),&dDum,pdPolyFitCoeff,1);
	  sprintf(strDum,"%f %f\n",ddl*i,dDum*dPhiCoeff+dPhiOffset);
	  if(fprintf(pF,strDum) != strlen(strDum)){
	    printf("INT Cannot write sample %li to file %s\n",i,strFileName);
	    goto bailout2;
	  }
	}
	sprintf(strDum,"#END SAMPLE SECTION\n");
	if(fprintf(pF,strDum) != strlen(strDum)){
	  printf("INT Cannot write %s to file %s\n",strDum,strFileName);
	  goto bailout2;
	}      
      bailout2:
	if(pF) fclose(pF);
      }

      if(iSaveOptions&SECTION_SAVE_RHO){
	sprintf(strFileName,"%s_sectionRho%i-%i_%i-%i",filenames[0],
		(int)section[0].x,(int)section[0].y,(int)section[1].x,(int)section[1].y);
	if(!iSectionSinglePixel){
	  sprintf(strDum,"_%.1f",dSectionRadius);
	  strcat(strFileName,strDum);
	}
	printf("INT Saving samples in %s\n",strFileName);
	if((pF=fopen(strFileName,"w"))==NULL){
	  printf("INT Cannot open file %s\n",strFileName);
	  goto bailout3;
	}
	sprintf(strDum,"#BEGIN SAMPLE SECTION\n");
	if(fprintf(pF,strDum) != strlen(strDum)){
	  printf("INT Cannot write %s to file %s\n",strDum,strFileName);
	  goto bailout3;
	}
	sprintf(strDum,"#%li\n",nNSamples);
	if(fprintf(pF,strDum) != strlen(strDum)){
	  printf("INT Cannot write nNSamples to file %s\n",strFileName);
	  goto bailout3;
	}
	for(i=0;i<nNSamples;i++){
	  sprintf(strDum,"%f %f\n",ddl*i,pfRhoSamples[i]);
	  if(fprintf(pF,strDum) != strlen(strDum)){
	    printf("INT Cannot write sample %li to file %s\n",i,strFileName);
	    goto bailout3;
	  }
	}
	sprintf(strDum,"#END SAMPLE SECTION\n");
	if(fprintf(pF,strDum) != strlen(strDum)){
	  printf("INT Cannot write %s to file %s\n",strDum,strFileName);
	  goto bailout3;
	}
      bailout3:
	fclose(pF);
      }
    }

    if(ch=='a'){
      // Val: change pointer here
      pf=pfRhoSamples;
      if(!iNA){
	iNASamples=nNSamples;
	if(!(pfASamples=(float*)malloc(nNSamples*sizeof(float)))){
	  printf("DIS Cannot allocate for pfASamples");
	  return;
	}
	for(i=0;i<iNASamples;i++) pfASamples[i]=pf[i];
      }
      else{
	iNASamples = iNASamples<nNSamples ? iNASamples:nNSamples;
	for(i=0;i<iNASamples;i++) pfASamples[i]+=pf[i];
      }
      iNA++;
    }
    
    if(ch=='d' && iNA){
      sprintf(strFileName,"%s_sectionARho%i-%i_%i-%i_a%i",filenames[0],
	      (int)section[0].x,(int)section[0].y,(int)section[1].x,(int)section[1].y,iNA);
      if(!iSectionSinglePixel){
	sprintf(strDum,"_%.1f",dSectionRadius);
	strcat(strFileName,strDum);
      }
      printf("INT Saving accumulated samples in %s\n",strFileName);
      if((pF=fopen(strFileName,"w"))==NULL){
	printf("INT Cannot open file %s\n",strFileName);
	goto bailout4;
      }
      sprintf(strDum,"#BEGIN SAMPLE SECTION\n");
      if(fprintf(pF,strDum) != strlen(strDum)){
	printf("INT Cannot write %s to file %s\n",strDum,strFileName);
	goto bailout4;
      }
      sprintf(strDum,"#%i\n",iNASamples);
      if(fprintf(pF,strDum) != strlen(strDum)){
	printf("INT Cannot write nNSamples to file %s\n",strFileName);
	goto bailout4;
      }
      for(i=0;i<iNASamples;i++){
	sprintf(strDum,"%f %f\n",ddl*i,pfASamples[i]/iNA);
	if(fprintf(pF,strDum) != strlen(strDum)){
	  printf("INT Cannot write sample %li to file %s\n",i,strFileName);
	  goto bailout4;
	}
      }
      sprintf(strDum,"#END SAMPLE SECTION\n");
      if(fprintf(pF,strDum) != strlen(strDum)){
	printf("INT Cannot write %s to file %s\n",strDum,strFileName);
	goto bailout4;
      }
    bailout4:
      fclose(pF);
      if(pfASamples) free(pfASamples);
      iNA=iNASamples=0;
    }
  }
  
  cpgclos();
  do_sectionplot=0;
  cpgslct(iCallerID);
  
  if(pfDerivative) free(pfDerivative);
  if(pdDerPolyFitCoeff) free(pdDerPolyFitCoeff);
  if(pdPolyFitCoeff) free(pdPolyFitCoeff);
  free(pfPhiSamples);
  free(pfRhoSamples);
  return;  
}

void DrawColorWedge(float fWedgePos,float fWedgeWidth,float fMin,float fMax,char *strUnits){
  float fPre=0.0,fPost=0.0;
  cpgsci(1);
  if(iDoLargeWedge) cpgsch(fWedgeLargeCH);
  else cpgsch(fWedgeSmallCH);
  if(iDoWedgeShift){
    Pallet(pallet+PALLET_SHIFT_SHIFT,1.0,0.5);
  }
  else Pallet(pallet,1.0,0.5);
  if(iDoWedgeAnnotation){
    if(iCutBlanksOnWedge){
      if(iReducedColorWedge){
	MakePallet(PALLET_REDUCED,fStim,fPreBlank,fPostBlank);
      }
      else{
	fPre=(fWedgeMax-fWedgeMin)*fPreBlank/fStim;
	fPost=(fWedgeMax-fWedgeMin)*fPostBlank/fStim;
      }
    }
    cpgwedg("BIL", fWedgePos, fWedgeWidth, fWedgeMin-fPre, fWedgeMax+fPost,strWedgeUnits);
  }
  else{
    cpgwedg("BI", fWedgePos, fWedgeWidth, fMin, fMax,strUnits);
  }
  if(iDoWedgeShift) Pallet(pallet,1.0,0.5);
}

int SaveSamples(int iN,int iNS,void *pvS,char *pcFileName,char *pcBegin,char *pcEnd){

  return(0);
}

#define TEXT_SHIFT_X_PCT 0.005
#define TEXT_SHIFT_Y_PCT 0.005

#define MARK_POINTS_SELECT_MAX 4
int aiMarkPointsSelectPhiColor[MARK_POINTS_SELECT_MAX]={1,1,0,0};
int aiMarkPointsSelectPhiColorSave[MARK_POINTS_SELECT_MAX]={0,0,1,1};
int aiMarkPointsSelectRhoColor[MARK_POINTS_SELECT_MAX]={2,3,4,4};
int aiMarkPointsSelectPhiSymbol[MARK_POINTS_SELECT_MAX]={-3,-4,-3,-4};
int aiMarkPointsSelectRhoSymbol[MARK_POINTS_SELECT_MAX]={-3,-4,-3,-4};
//int aiMarkPointsSelectPhiSymbol[MARK_POINTS_SELECT_MAX]={4,2,4,2};
//int aiMarkPointsSelectRhoSymbol[MARK_POINTS_SELECT_MAX]={4,2,4,2};

// -4 - cross 

void MarkPoints(int iDX,int iColor,int iMode){
  float fX1,fY1,fX2,fY2;
  float fXTextShift,fYTextShift;
  char str[128];
  int i;
  RECORD *pR;
  int iColorR,iColorE,iColorI;
  int iSymbolR,iSymbolE,iSymbolI,iSymbol;
  int iSelect;
  int *piColor,*piSymbol;
  float fMarkPointsCrossSize=3.0;
  cpgsave();
  cpgslw(4);
  cpgsch(1);
  //cpgsci(iColor);
  fMarkPointsRadius=3.5;// used for correlation
  //fMarkPointsRadius=3.0;

  cpgqwin(&fX1,&fX2,&fY1,&fY2);
  fXTextShift=(fX2-fX1)*TEXT_SHIFT_X_PCT;
  fYTextShift=(fY2-fY1)*TEXT_SHIFT_Y_PCT;
 
  if(iMarkPointsHide) goto label_correlation;

  if(iMode==INT_MODE_PHI){
    if(save_on){
      printf("DIS Mark: Save on\n");
      //weird, next line should Save
      piColor=aiMarkPointsSelectPhiColor;
    }
    else{
      //weird, next line shuld be no Save
      piColor=aiMarkPointsSelectPhiColorSave;
    }
    piSymbol=aiMarkPointsSelectPhiSymbol;
  }
  else{
    piColor=aiMarkPointsSelectRhoColor;
    piSymbol=aiMarkPointsSelectRhoSymbol;
  }
  cpgsfs(1);
  for(i=0;i<iMarkPointsN;i++){
    iSelect=pMarkPoints[i].iSelect<MARK_POINTS_SELECT_MAX ? pMarkPoints[i].iSelect : MARK_POINTS_SELECT_MAX-1;
    if(iSelect<0) continue;
    cpgsci(piColor[iSelect]);
    iSymbol=piSymbol[iSelect];
    switch(iSymbol){
    case -4:
      cpgmove(pMarkPoints[i].x+fMarkPointsCrossSize,pMarkPoints[i].y+fMarkPointsCrossSize);
      cpgdraw(pMarkPoints[i].x-fMarkPointsCrossSize,pMarkPoints[i].y-fMarkPointsCrossSize);
      cpgmove(pMarkPoints[i].x-fMarkPointsCrossSize,pMarkPoints[i].y+fMarkPointsCrossSize);
      cpgdraw(pMarkPoints[i].x+fMarkPointsCrossSize,pMarkPoints[i].y-fMarkPointsCrossSize);
      break;
    case -3:
      cpgcirc(pMarkPoints[i].x,pMarkPoints[i].y,fMarkPointsRadius);
      break;
    default:
      cpgpt1(pMarkPoints[i].x,pMarkPoints[i].y,iSymbol);
    }
  }
  cpgsfs(1);

  /*
  if(iMarkPointsSymbol!=-3){
    for(i=0;i<iMarkPointsN;i++)
      cpgpt1(pMarkPoints[i].x,pMarkPoints[i].y,iMarkPointsSymbol);
  }
  else{
    cpgsfs(1);
    for(i=0;i<iMarkPointsN;i++)
      cpgcirc(pMarkPoints[i].x,pMarkPoints[i].y,fMarkPointsRadius);
    cpgsfs(1);
  }
  */

  if(iMarkPointsEnumerate){
    cpgsch(fMarkPointsCharSize);
    for(i=0;i<iMarkPointsN;i++){
      sprintf(str,"%i",i+1);
      cpgtext(fXTextShift+pMarkPoints[i].x,-fYTextShift+pMarkPoints[i].y,str);
    }   
  }

 label_correlation:

  // Should figure out a decent color scheme
  // Phi: 1(white) - regular, 0(black) - HIDE_EXT, HIDDEN_COLOR_N(custom gray) - HIDE_INT
  // Phi: 0(black) - regular, 1(white) - HIDE_EXT, 0(black) - HIDE_INT
  // Rho: 2(red)   - regular, 4(blue)  - HIDE_EXT, 3(green)        - HIDE_INT
  if(iMarkCorrelationPoints && iNCorrRecords){
    if(iMode==INT_MODE_PHI){
      iColorR=1;
      iColorE=0;
      iColorI=HIDDEN_COLOR_N;
      cpgscr(HIDDEN_COLOR_N,HIDDEN_COLOR_R,HIDDEN_COLOR_G,HIDDEN_COLOR_B);
      iColorR=0;
      iColorE=1;
      iColorI=0;
      iSymbolR=-3;
      iSymbolE=2;
      iSymbolI=2;
     }
    else{
      iColorR=2;
      iColorE=4;
      iColorI=3;
      iSymbolR=2;
      iSymbolE=4;
      iSymbolI=2;
    }
    cpgsfs(2);
    cpgslw(3);
    cpgsch(fMarkPointsCharSize);
    for(i=0,pR=pCorrRecords;i<iNCorrRecords;i++,pR++){
      if(pR->iHide&RECORD_HIDE_EXT){
	if(pR->iHide&RECORD_HIDE_EXT_FAKE){
	  cpgsci(iColorI);
	  iSymbol=iSymbolI;
	}
	else{
	  cpgsci(iColorE);
	  iSymbol=iSymbolE;
	}
      }
      else{
	if(pR->iHide&RECORD_HIDE_INT){
	  cpgsci(iColorI);
	  iSymbol=iSymbolI;
	}
	else{
	  cpgsci(iColorR);
	  iSymbol=iSymbolR;
	}
      }
      if(iSymbol==-3) cpgcirc((float)(pR->dX),(float)(pR->dY),fMarkPointsRadius);
      else cpgpt1((float)(pR->dX),(float)(pR->dY),iSymbol);
      
      if(iMarkPointsEnumerate){
	sprintf(str,"%i",pR->iIndex);
	cpgtext(fXTextShift+(float)(pR->dX),-fYTextShift+(float)(pR->dY),str);
      }
    }
    cpgsfs(1);
    
  }
  
  cpgunsa();
}

void PlotLineSegments(int iColor,int iWidth,int Style){
  float fX1,fY1,fX2,fY2;
  float fXTextShift,fYTextShift;
  char str[128];
  int i;
  int iColorB=3,iColorE=2;

  cpgsave();
  cpgslw(iWidth);
  cpgsci(iColor);
 
  if(Style&LINE_SEGMENTS_STYLE_ARROW){
    cpgsah(1,30.0,0.3);
    cpgsch(0.7);
    for(i=0;i<iNLineSegments;i++){
      //cpgarro(pLineSegments[i].x1*transs[1]+transs[0],pLineSegments[i].y1*transs[5]+transs[3],pLineSegments[i].x2*transs[1]+transs[0],pLineSegments[i].y2*transs[5]+transs[3]);
      cpgarro(pLineSegments[i].x1,pLineSegments[i].y1,pLineSegments[i].x2,pLineSegments[i].y2);
   }
  }
  else{
    for(i=0;i<iNLineSegments;i++){
      //      cpgmove(pLineSegments[i].x1*transs[1]+transs[0],pLineSegments[i].y1*transs[5]+transs[3]);
      //      cpgdraw(pLineSegments[i].x2*transs[1]+transs[0],pLineSegments[i].y2*transs[5]+transs[3]);
      cpgmove(pLineSegments[i].x1,pLineSegments[i].y1);
      cpgdraw(pLineSegments[i].x2,pLineSegments[i].y2);
    }
  }

  if(Style&LINE_SEGMENTS_STYLE_MARK_ENDS){
    for(i=0;i<iNLineSegments;i++){
      cpgsci(iColorB);
      //cpgpt1(pLineSegments[i].x1*transs[1]+transs[0],pLineSegments[i].y1*transs[5]+transs[3],2);
      cpgpt1(pLineSegments[i].x1,pLineSegments[i].y1,2);
      cpgsci(iColorE);
      //cpgpt1(pLineSegments[i].x2*transs[1]+transs[0],pLineSegments[i].y2*transs[5]+transs[3],2);
      cpgpt1(pLineSegments[i].x2,pLineSegments[i].y2,2);
    }
  }


  if(Style&LINE_SEGMENTS_STYLE_ENUMERATE){
    cpgqwin(&fX1,&fX2,&fY1,&fY2);
    fXTextShift=(fX2-fX1)*TEXT_SHIFT_X_PCT;
    fYTextShift=(fY2-fY1)*TEXT_SHIFT_Y_PCT;
    cpgsch(fMarkPointsCharSize);
    for(i=0;i<iNLineSegments;i++){
      sprintf(str,"%i",i+1);
      //cpgtext(fXTextShift+pLineSegments[i].x1*transs[1]+transs[0],-fYTextShift+pLineSegments[i].y1*transs[5]+transs[3],str);
      cpgtext(fXTextShift+pLineSegments[i].x1,-fYTextShift+pLineSegments[i].y1,str);
    }   
  }

  cpgunsa();
}
