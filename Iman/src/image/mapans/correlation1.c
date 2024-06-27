/* Map analysis                 */
/* Written by V.Kalatsky        */
/* 06Feb2001                    */
/* Modified 10May2001           */
/* Made separate file 24Dec2003 */

#include "mapans1.h"

int GetValuesFromMap(int iN,RECORD *pRec,int xd,int yd,double **mp,double *pdTMax,double *pdTMin,int iMode);
double Chi2ShiftRec(int iN,RECORD *pRec,double dS,double *pdAver1,double *pdAver2,double *pdAver4);
double Chi2Rec(int iN,RECORD *pRec,double *pdAver1,double *pdAver2,double *pdAver4);
int ShiftWrapRec(int iN,RECORD *pRec,double dS,int iMode);
double MinimizeChi2NoWrap(int iN,RECORD *pRec,double *pdMinChi2,double *pdMinS,int iMode);
int PlotRecords(int iId,int iN,RECORD *pRec,float fS,float fMinY,float fMaxY,int iMode);
int ThresholdRecords(int iNR,RECORD *pRec,double dT,int *piNR,int iMode);
int SampleChi2(int iN,RECORD *pRec,int iNSamples,float *pfS,float *pfMin,float *pfMinX,float *pfMax,double *pdAver1,double *pdAver2,double *pdAver4,int iMode);
int PlotChi2(int iId,int iN,float *pfS,float fMin,float fMax,float fMinX,float fT);
void SwapChunks(void *p1,void *p2, int iSize,void *p);

#define N_SAMPLES 200

#define CORRELATION_MODE_ENUMERATE   1
#define CORRELATION_MODE_PLOT_HIDDEN 2
#define CORRELATION_MODE_NOSHIFT     4

#define THRESHOLD_MODE_VERBOSE 1

#define THRESHOLD_NSAMPLES 20
#define CORRELATION_TEST_NTRIALS 10000
#define CORRELATION_DISTRIBUTION_NBINS 500
#define CORRELATION_NSWAPS 50


//#define MAKE_PLOTS_HACK
#define PLOT_FILES_N 6
/*
// Run 7
char *strsPlotFiles[]={
  "r045CorrelationOptical-7L3_CFkHz5ExtraZeros.txt",
  "r045CorrelationOptical-7L3_BWAtIHiLoBFSQRT.txt",
  "r045CorrelationOptical-7L3_BFCMkHz.txt",
  "r045CorrelationOptical-7L3_BFCMRkHz.txt",
  "r045CorrelationOptical-7L3_FourierPhi.txt"
};
*/
/*
char *strsPlotFiles[]={
  "r045CorrelationOptical-7L3_CFkHz110Below40.txt",
  "r045CorrelationOptical-7L3_BWAtIHiLoBFSQRT110Below40.txt",
  "r045CorrelationOptical-7L3_BFCMkHz110Below40.txt",
  "r045CorrelationOptical-7L3_BFCMRkHz110Below40.txt",
  "r045CorrelationOptical-7L3_FourierPhi110Below40.txt",
  "r045CorrelationOptical-7L3_CFkHz120.txt"
};
*/
/*
char *strsPlotFiles[]={
  "r045CorrelationOptical-8iL3_CFkHz5ExtraZeros.txt",
  "r045CorrelationOptical-8iL3_BWAtIHiLoBFSQRT.txt",
  "r045CorrelationOptical-8iL3_BFCMkHz.txt",
  "r045CorrelationOptical-8iL3_BFCMRkHz.txt",
  "r045CorrelationOptical-8iL3_FourierPhi.txt"
};
*/
/*
char *strsPlotFiles[]={
  "r045CorrelationOptical-8iL3_CFkHz110Below40.txt",
  "r045CorrelationOptical-8iL3_BWAtIHiLoBFSQRT110Below40.txt",
  "r045CorrelationOptical-8iL3_BFCMkHz110Below40.txt",
  "r045CorrelationOptical-8iL3_BFCMRkHz110Below40.txt",
  "r045CorrelationOptical-8iL3_FourierPhi110Below40.txt",
  "r045CorrelationOptical-8iL3_CFkHz120.txt"
};
*/
  /*
char *strsPlotFiles[]={
  "r045CorrelationOptical-78averageL3_CFkHz5ExtraZeros.txt",
  "r045CorrelationOptical-78averageL3_BWAtIHiLoBFSQRT.txt",
  "r045CorrelationOptical-78averageL3_BFCMkHz.txt",
  "r045CorrelationOptical-78averageL3_BFCMRkHz.txt",
  "r045CorrelationOptical-78averageL3_FourierPhi.txt"
};
  */

char *strsPlotFiles[]={
  "r045CorrelationOptical-78averageL3_CFkHz110Below40.txt",
  "r045CorrelationOptical-78averageL3_BWAtIHiLoBFSQRT110Below40.txt",
  "r045CorrelationOptical-78averageL3_BFCMkHz110Below40.txt",
  "r045CorrelationOptical-78averageL3_BFCMRkHz110Below40.txt",
  "r045CorrelationOptical-78averageL3_FourierPhi110Below40.txt",
  "r045CorrelationOptical-78averageL3_CFkHz120.txt"
};

#define PLOT_FILES_SPECIAL_N 1

char *strsPlotFilesSpecial[]={
  "r045CorrelationOptical-7iL3_FourierPhi110Below40.txt"
};

/*
char *strsPlotFilesSpecial[]={
  "r045CorrelationOptical-8L3_FourierPhi110Below40.txt"
};
*/

char *strsPlotText[]={"CF","BF","BFCM","BFCMR","Fourier","CF All"};

int aiColor[]={2,3,4,5,1,2};
int aiLineStyle[]={1,1,1,1,1,2};

int CorrelationToExternalSource(int iN,RECORD *pRec,int xd,int yd,double **mp,int iMode){
  int i,k,l;
  unsigned long j;
  RECORD *pR;
  int iHE=0;
  int iNR=0;
  int iPlotID,iChiID;
  float afS[N_SAMPLES];
  float fMin,fMax,fMinX,fDum,fShift;
  float fMinY,fMaxY;
  float fX,fY;
  char ch,str[256];
  int iSaveID;
  int iPlotMode=CORRELATION_MODE_ENUMERATE;;
  double dThresholdMin,dThresholdMax;
  double dThresholdMinDum,dThresholdMaxDum;
  double dDum,dT,dStep;
  double adChi21[THRESHOLD_NSAMPLES],adChi22[THRESHOLD_NSAMPLES];
  float afShift1[THRESHOLD_NSAMPLES];
  double adChi2Aver1[THRESHOLD_NSAMPLES],adChi2Aver2[THRESHOLD_NSAMPLES];
  double adShiftAver1[THRESHOLD_NSAMPLES],adShiftAver2[THRESHOLD_NSAMPLES];
  int iSaveID1;
  double dAver1,dAver2,dAver4;
  float fDistributionCoeff=4.0*CORRELATION_DISTRIBUTION_NBINS;
  float fDistributionPlotCoef=0.005*CORRELATION_DISTRIBUTION_NBINS/CORRELATION_TEST_NTRIALS;
  unsigned int aauiDBins[THRESHOLD_NSAMPLES][CORRELATION_DISTRIBUTION_NBINS];
  int iNSamplesSignif,iSignif;
  int iDoShuffling=0;
  int iSwap1,iSwap2;
  //char strDevice[]={"CorrPlot8-110Random.ps/vcps"};
  char strDevice[]={"/xw"};
  float fOpticalShift=0.23;
  float fNoiseLevel=0.19;//run 7
  //float fNoiseLevel=0.12;//run 8; average 7 and 8
  //float fNoiseLevel=0.15;//run 8
  float fCorreleationPlotMaxY=0.4;
  double dShift,dExtraShift,dChi2Min;
  int iPlotFileSpecial=1;
  FILE *pF=NULL;
  char *pc,strDum[256];
  int iSlice=12;
  int iNSamplesOut=0;

  if(!iN || !pRec) return(1); 
  if(!iCorrelationModeShift) iMode|=CORRELATION_MODE_NOSHIFT;

  for(iHE=i=0,pR=pRec;i<iN;i++,pR++){
    pR->iHide=0;
    if(pR->dZ<=0.0){
      pR->iHide|=RECORD_HIDE_EXT;
      iHE++;
    }
    else{
      if(iDoWedgeAnnotation) pR->dZR=(log(pR->dZ)/M_LN2-fWedgeMin)/(fWedgeMax-fWedgeMin);
      else pR->dZR=pR->dZ;
    }
  }
  if(iHE){
    printf("INT Removed %i External Records: ",iHE);
    for(i=0,pR=pRec;i<iN;i++,pR++) if(pR->iHide&RECORD_HIDE_EXT) printf("%i ",pR->iIndex);
    printf("\n");
  }

  GetValuesFromMap(iN,pRec,xd,yd,mp,&dThresholdMax,&dThresholdMin,0);

  if(iDoWedgeAnnotation) fDum=fWedgeMax-fWedgeMin;
  else fDum=1.0;
  printf("#NRecs %i(%i) Rho Min %f Max %f Noise %f Muliplier %f Oct\n",iN,iN-iHE,dThresholdMin,dThresholdMax,dRhoNoise,fDum);
  printf("#Rho\tsqr(C2)\tChi2\tError\tShift\tExtVal\tPcnt\t#Index\n");
  for(i=0,pR=pRec;i<iN;i++,pR++){
    if(!(pR->iHide&RECORD_HIDE_EXT)){
      ThresholdRecords(iN,pRec,pR->dRho,NULL,0);
      SampleChi2(iN,pRec,N_SAMPLES,afS,&fMin,&fMinX,&fMax,&dAver1,&dAver2,&dAver4,0);
      printf("%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.3f\t%0.1f%%\t#%i\n",
	     pR->dRho,sqrt(dAver2),dAver2,sqrt(dAver4-dAver2*dAver2),
	     fMinX,pR->dZ,100.0*(pR->dRho-dThresholdMin)/(dThresholdMax-dThresholdMin),pR->iIndex);
    }
    else printf("#%0.4f\tNA\tNA\tNA\tNA\t%0.3f\t%0.1f%%\t#%i\n",
		pR->dRho,pR->dZ,100.0*(pR->dRho-dThresholdMin)/(dThresholdMax-dThresholdMin),pR->iIndex);
  }


  ThresholdRecords(iN,pRec,dCorrelationThreshold,&iNR,THRESHOLD_MODE_VERBOSE);
  if(!iNR){
    printf("INT No useful records\n");
    return(2);
  }

  fMaxY=-(fMinY=1.0e20);
  for(i=0,j=0,pR=pRec;i<iN;i++,pR++){
    if(!(pR->iHide&RECORD_HIDE_EXT)){
      fDum=pR->dZR;
      if(fMinY>fDum) fMinY=fDum;
      if(fMaxY<fDum) fMaxY=fDum;
    }
  }

  SampleChi2(iN,pRec,N_SAMPLES,afS,&fMin,&fMinX,&fMax,&dAver1,&dAver2,&dAver4,THRESHOLD_MODE_VERBOSE);

  cpgqid(&iSaveID);
  cpgsave();

  if((iPlotID=cpgopen("/xw"))<1){
    printf("ANA Cannot open iPlotID\n");
    goto bailout;
  }
  cpgask(0);
  cpgpap(7,1);
  cpgscr(HIDDEN_COLOR_N,HIDDEN_COLOR_R,HIDDEN_COLOR_G,HIDDEN_COLOR_B);

  ShiftWrapRec(iN,pRec,(double)(fShift=fMinX),0);
  dChi2Min=MinimizeChi2NoWrap(iN,pRec,&dAver2,&dExtraShift,0);
  printf("COR ExtrShift %f Chi %f\n",dExtraShift,sqrt(dChi2Min));
  PlotRecords(iPlotID,iN,pRec,fShift,fMinY,fMaxY,iPlotMode);
  //getchar();
  if((iChiID=cpgopen("/xw"))<1){
    printf("ANA Cannot open iChiID\n");
    goto bailout;
  }
  cpgask(0);
  cpgpap(7,1);
  PlotChi2(iChiID,N_SAMPLES,afS,fMin,fMax,fMinX,(float)dCorrelationThreshold);

  ThresholdRecords(iN,pRec,dCorrelationThreshold,&iNR,THRESHOLD_MODE_VERBOSE);

  while(1){
    cpgband(0,0,fX,fY,&fX,&fY,&ch);
    if(ch=='Q') break;
    if(ch=='0'){
      ShiftWrapRec(iN,pRec,(double)fMinX,0);
      PlotRecords(iPlotID,iN,pRec,fMinX,fMinY,fMaxY,iPlotMode);
      continue;
    }
    if(ch=='e'){
      if(iPlotMode&CORRELATION_MODE_ENUMERATE) iPlotMode&=~CORRELATION_MODE_ENUMERATE;
      else iPlotMode|=CORRELATION_MODE_ENUMERATE;
      ShiftWrapRec(iN,pRec,(double)fShift,0);
      PlotRecords(iPlotID,iN,pRec,fShift,fMinY,fMaxY,iPlotMode);
      continue;
    }
    if(ch=='t'){
      cpgqid(&iSaveID1);
      cpgsave();
      if((cpgopen("/xw"))<1){
	printf("ANA Cannot open Interactive Threshold\n");
      }
      cpgask(0);
      cpgpap(9,1);
      cpgscr(HIDDEN_COLOR_N,HIDDEN_COLOR_R,HIDDEN_COLOR_G,HIDDEN_COLOR_B);
      cpgenv(0.0,(float)dThresholdMax,0,(float)iN,0,1);
      sprintf(str,"Threshold %0.3f Shift %0.3f Chi %0.3f",dCorrelationThreshold,fMinX,fMin);
      cpglab("Rho","# Penetrations",str);
      for(i=0,pR=pRec;i<iN;i++,pR++){
	if(pR->iHide&RECORD_HIDE_EXT) cpgsci(4);
	else{
	  if(pR->iHide&RECORD_HIDE_INT) cpgsci(3);
	  else cpgsci(2);
	}
	cpgpt1(pR->dRho,iN-i,2);
      }
      cpgsci(1);
      while(1){
	cpgband(6,0,fX,fY,&fX,&fY,&ch);
	if(ch=='q' || ch=='Q') break;

	if(ch=='d'){
	  dDum=fabs(pRec->dRho-fX);
	  k=0;
	  pR=pRec+1;
	  for(i=1;i<iN;i++,pR++){
	    if(dDum>fabs(pR->dRho-fX)){
	      dDum=fabs(pR->dRho-fX);
	      k=i;
	    }
	  }
	  printf("NEAR record %i\n",pRec[k].iIndex);
	  if(pRec[k].iHide&RECORD_HIDE_EXT){
	    if(pRec[k].iHide&RECORD_HIDE_EXT_FAKE) pRec[k].iHide&=~(RECORD_HIDE_INT|RECORD_HIDE_EXT|RECORD_HIDE_EXT_FAKE);
	  }
	  else{
	    pRec[k].iHide|=(RECORD_HIDE_INT|RECORD_HIDE_EXT|RECORD_HIDE_EXT_FAKE);
	  }
	  goto label;
	}
	if(ch=='0'){
	  for(i=0,pR=pRec;i<iN;i++,pR++){
	    if(pR->iHide&RECORD_HIDE_EXT){
	      if(pR->iHide&RECORD_HIDE_EXT_FAKE) pR->iHide&=~(RECORD_HIDE_INT|RECORD_HIDE_EXT|RECORD_HIDE_EXT_FAKE);
	    } 
	  } 
	  goto label;
	}

	dCorrelationThreshold=fX;

      label:

	ThresholdRecords(iN,pRec,dCorrelationThreshold,NULL,THRESHOLD_MODE_VERBOSE);
	SampleChi2(iN,pRec,N_SAMPLES,afS,&fMin,&fMinX,&fMax,&dAver1,&dAver2,&dAver4,THRESHOLD_MODE_VERBOSE);
	PlotChi2(iChiID,N_SAMPLES,afS,fMin,fMax,fMinX,(float)dCorrelationThreshold);

	fShift=fMinX;
	ShiftWrapRec(iN,pRec,(double)fShift,0);
	PlotRecords(iPlotID,iN,pRec,fShift,fMinY,fMaxY,iPlotMode);

	DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
	
	cpgenv(0.0,(float)dThresholdMax,0,(float)iN,0,1);
	sprintf(str,"Threshold %0.3f Shift %0.3f Chi %0.3f",dCorrelationThreshold,fMinX,fMin);
	cpglab("Rho","# Penetrations",str);
	for(i=0,pR=pRec;i<iN;i++,pR++){
	  if(pR->iHide&RECORD_HIDE_EXT){
	    if(pR->iHide&RECORD_HIDE_EXT_FAKE) cpgsci(1);
	    else cpgsci(4);
	  }
	  else{
	    if(pR->iHide&RECORD_HIDE_INT) cpgsci(3);
	    else cpgsci(2);
	  }
	  cpgpt1(pR->dRho,iN-i,2);
	}
	cpgsci(1);
	dStep=(dThresholdMax-dThresholdMin)/THRESHOLD_NSAMPLES;
	for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	  ThresholdRecords(iN,pRec,dT,&iNR,0);
	  SampleChi2(iN,pRec,N_SAMPLES,afS,&fMin,&fMinX,&fMax,&dAver1,adChi21+i,adChi22+i,iMode);
	  afShift1[i]=fMinX;
	  adChi22[i]=sqrt(adChi22[i]-adChi21[i]*adChi21[i]);
	  if(i) cpgdraw((float)dT,400*sqrt(fMin));
	  else  cpgmove((float)dT,400*sqrt(fMin));
	}

	cpgsci(6);
	for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	  cpgpt1((float)dT,400*sqrt(adChi21[i]),2);
	  cpgmove((float)dT,400*sqrt(adChi21[i]-adChi22[i]>0.0?adChi21[i]-adChi22[i]:0.0));
	  cpgdraw((float)dT,400*sqrt(adChi21[i]+adChi22[i]));
	}

	cpgsci(7);
	cpgmove((float)dThresholdMin,400*(afShift1[0]>0.0?afShift1[0]:-afShift1[0]));
	for(i=1,dT=dThresholdMin+dStep;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	  cpgdraw((float)dT,400*(afShift1[i]>0.0?afShift1[i]:-afShift1[i]));
	}
	cpgsci(1);
     }
      cpgclos();
      cpgslct(iSaveID1);
      cpgunsa();
      continue;
    }

    if(ch=='T'){
      dStep=(dThresholdMax-dThresholdMin)/THRESHOLD_NSAMPLES;
      if(iDoWedgeAnnotation) fDum=fWedgeMax-fWedgeMin;
      else fDum=1.0;
      printf("#NRecs %i(%i) Rho Min %f Max %f Noise %f Muliplier %fOct\n",iN,iN-iHE,dThresholdMin,dThresholdMax,dRhoNoise,fDum);
      printf("#Treshold\tsqrt(Chi2)\tShift\t\tNSampls\tThreshold%%\n");
      for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	ThresholdRecords(iN,pRec,dT,&iNR,0);
	SampleChi2(iN,pRec,N_SAMPLES,afS,&fMin,&fMinX,&fMax,&dAver1,&dAver2,&dAver4,iMode);
	adChi21[i]=fMin;
	afShift1[i]=fMinX;
	printf("%f\t%f\t%f\t%i\t%0.0f%%\n",dT,sqrt((double)fMin),fMinX, iNR,i*100.0/THRESHOLD_NSAMPLES);
      }
      ThresholdRecords(iN,pRec,dCorrelationThreshold,NULL,0);
      SampleChi2(iN,pRec,N_SAMPLES,afS,&fMin,&fMinX,&fMax,&dAver1,&dAver2,&dAver4,0);
      //PlotChi2(iChiID,N_SAMPLES,afS,fMin,fMax,fMinX);
      //fShift=fMinX;
      //ShiftAndPlotRec(iPlotID,iN,pRec,fShift,fMinY,fMaxY,iPlotMode);
      //ShiftWrapRec(iN,pRec,(double)fShift,0);
      //PlotRecords(iPlotID,iN,pRec,fShift,fMinY,fMaxY,iPlotMode);
      continue;
    }


    if(ch=='Z'){
      dStep=(dThresholdMax-dThresholdMin)/THRESHOLD_NSAMPLES;

      cpgqid(&iSaveID1);
      cpgsave();
      if((cpgopen(strDevice))<1){
	printf("ANA Cannot open iCorrelationLumpSum\n");
      }
      cpgask(0);
      cpgpap(7,0.7);
      cpgscr(HIDDEN_COLOR_N,HIDDEN_COLOR_R,HIDDEN_COLOR_G,HIDDEN_COLOR_B);
      cpgslw(2);
      // 4 correspond to # octaves
      cpgscf(2);
#ifdef MAKE_PLOTS_HACK
      fCorreleationPlotMaxY=0.2;
#endif
      cpgenv(0.0,(float)dThresholdMax,0,4*fCorreleationPlotMaxY,0,1);
      cpgswin(0.0,(float)dThresholdMax,0,fCorreleationPlotMaxY);

      //Draw optical shift and noise level
      cpgsci(3);
      cpgsls(4);
      cpgmove(fNoiseLevel,0.0);
      cpgdraw(fNoiseLevel,fCorreleationPlotMaxY);	  
      cpgsls(1);


#ifdef MAKE_PLOTS_HACK
      cpgslw(4);
      fDum=0.04;
      for(i=0;i<PLOT_FILES_N;i++){
	if(!(pF=fopen(strsPlotFiles[i],"r"))){
	  printf("ANA Cannot open file %s\n",strsPlotFiles[i]);
	  continue;
	}
	cpgsci(aiColor[i]);
	cpgsls(aiLineStyle[i]);
	k=0;
	dT=dThresholdMin;
	while(!GetStringFromFile(pF,strDum,255,1)){
	  pc=strtok(strDum," \t\n");
	  if(pc) fX=atof(pc);
	  else{
	    printf("INT Incomplete X record\n");
	    break;
	  }
	  pc=strtok(NULL," \t\n");
	  if(pc) fY=atof(pc);
	  else{
	    printf("INT Incomplete Y record\n");
	    break;
	  }
	  if(k) cpgdraw(fX,fY);
	  else cpgmove(fX,fY);
	  k++;
	  dT+=dStep;
	}
	//plot legend
	cpgmove((float)dThresholdMax*0.73,fCorreleationPlotMaxY*(1-fDum*(i+2)));
	cpgdraw((float)dThresholdMax*0.84,fCorreleationPlotMaxY*(1-fDum*(i+2)));
	cpgsci(1);
	cpgsls(1);
	cpgslw(2);
	cpgtext((float)dThresholdMax*0.85,fCorreleationPlotMaxY*(1-fDum*(i+2.25)),strsPlotText[i]);
	cpgslw(4);
	fclose(pF);
      }

      cpgband(0,0,fX,fY,&fX,&fY,&ch);
      cpgclos();
      cpgslct(iSaveID1);
      cpgunsa();

      continue;
#endif
      cpgslw(2);
      cpgsci(3);
      cpgsls(4);
      cpgmove(0,fOpticalShift);
      cpgdraw((float)dThresholdMax,fOpticalShift);
      cpgsls(1);

      memset(adChi21,0,THRESHOLD_NSAMPLES*sizeof(double));

      for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	ThresholdRecords(iN,pRec,dT,&iNR,0);
	SampleChi2(iN,pRec,N_SAMPLES,afS,&fMin,&fMinX,&fMax,&dAver1,&dAver2,&dAver4,iMode);
	adChi21[i]=fMin;
	afShift1[i]=fMinX;
      }

      // Random optical penetrations
      
      memset(adChi2Aver1,0,THRESHOLD_NSAMPLES*sizeof(double));
      memset(adChi2Aver2,0,THRESHOLD_NSAMPLES*sizeof(double));
      memset(adShiftAver1,0,THRESHOLD_NSAMPLES*sizeof(double));
      memset(adShiftAver2,0,THRESHOLD_NSAMPLES*sizeof(double));
      for(i=0;i<THRESHOLD_NSAMPLES;i++) 
	memset(aauiDBins[i],0,sizeof(unsigned int)*CORRELATION_DISTRIBUTION_NBINS);

      iNSamplesOut=0;
      for(k=0;k<CORRELATION_TEST_NTRIALS;k++){
	GetValuesFromMap(iN,pRec,xd,yd,mp,&dThresholdMaxDum,&dThresholdMinDum,1);
	//cpgbbuf();
	for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	  //Only one thresh level
	  if(i!=iSlice) continue;
	  ThresholdRecords(iN,pRec,dT,&iNR,0);
	  SampleChi2(iN,pRec,N_SAMPLES,afS,&fMin,&fMinX,&fMax,&dAver1,&dAver2,&dAver4,0);
	  adChi2Aver1[i]+=fMin;
	  adChi2Aver2[i]+=(double)fMin*(double)fMin;
	  adShiftAver1[i]+=fMinX;
	  adShiftAver2[i]+=(double)fMinX*(double)fMinX;
	  if((l=(int)floor(fDistributionCoeff*fMin))<CORRELATION_DISTRIBUTION_NBINS){
	    aauiDBins[i][l]++;
	  }
	  else{
	    printf("Sample %i out of bound %i\n",i,l);
	    iNSamplesOut++;
	  }
	 
	  //if(k%100){
	  //  if(i) cpgdraw((float)dT,sqrt(fMin));
	  //  else cpgmove((float)dT,sqrt(fMin));
	  //}
	}
	//cpgebuf();
      }
      if(iNSamplesOut) printf("WARNING %i samples out\n",iNSamplesOut);

      printf("#Random Penetration\n");
      //      printf("#Tresh\tSR(C2)\tSR(C2-)\tSR(C2+)\tShift\tShift-\tShift+\tPcnt\n");
      printf("#Tresh\tSR(C2)\tSR(95%%)\tShift\tShift-\tShift+\tPcnt\n");
      if(iDoWedgeAnnotation) fDum=fWedgeMax-fWedgeMin;
      else fDum=1.0;

      for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	adChi2Aver1[i]/=CORRELATION_TEST_NTRIALS;
	adChi2Aver2[i]=sqrt(adChi2Aver2[i]/CORRELATION_TEST_NTRIALS-adChi2Aver1[i]*adChi2Aver1[i]);
	adShiftAver1[i]/=CORRELATION_TEST_NTRIALS;
	adShiftAver2[i]=sqrt(adShiftAver2[i]/CORRELATION_TEST_NTRIALS-adShiftAver1[i]*adShiftAver1[i]);
	/*
	printf("%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.0f%%\n",dT,
	       fDum*sqrt(adChi2Aver1[i]),
	       adChi2Aver1[i]-adChi2Aver2[i]>0.0?fDum*sqrt(adChi2Aver1[i]-adChi2Aver2[i]):0.0,
	       fDum*sqrt(adChi2Aver1[i]+adChi2Aver2[i]),
	       fDum*adShiftAver1[i],
	       fDum*(adShiftAver1[i]-adShiftAver2[i]),
	       fDum*(adShiftAver1[i]+adShiftAver2[i]),
	       i*100.0/THRESHOLD_NSAMPLES);
	*/
      }

      iNSamplesSignif=CORRELATION_TEST_NTRIALS;
      GetValuesFromMap(iN,pRec,xd,yd,mp,&dThresholdMaxDum,&dThresholdMinDum,0);
      for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	ThresholdRecords(iN,pRec,dT,&iNR,0);
	dDum=adChi2Aver1[i]*fDistributionCoeff;
	for(iNSamplesSignif=k=0;k<CORRELATION_DISTRIBUTION_NBINS;k++){
	  if(k<dDum) iNSamplesSignif+=aauiDBins[i][k];
	  else break;
	}
	//printf("%i\n",iNSamplesSignif);
	iNSamplesSignif*=0.05;
	// Plot distributions
	if((l=aauiDBins[i][0])>=iNSamplesSignif) iSignif=0;
	else iSignif=-1;
	cpgsci(5);
	cpgmove((float)(dT+aauiDBins[i][0]*fDistributionPlotCoef),0);
	for(k=1;k<CORRELATION_DISTRIBUTION_NBINS;k++){
	  cpgdraw((float)(dT+aauiDBins[i][k]*fDistributionPlotCoef),sqrt(k/fDistributionCoeff));	  
	  if(iSignif<0) if((l+=aauiDBins[i][k])>=iNSamplesSignif) iSignif=k;
	}

	//Print out distribution
	if(i==iSlice){
	  printf("# Random penetation Trials %i BinsN %i Bin size %f\n",CORRELATION_TEST_NTRIALS,CORRELATION_DISTRIBUTION_NBINS,1.0/fDistributionCoeff);
	  printf("# line %i Thresh %f\n",i,dT);
	  printf("# Average %f Absolute 95%% level %f\n",4*sqrt(adChi2Aver1[i]),4*sqrt(iSignif/fDistributionCoeff));
	  for(k=0;k<CORRELATION_DISTRIBUTION_NBINS;k++){
	    printf("%f %u\n",4*sqrt(k/fDistributionCoeff),aauiDBins[i][k]);	  
	  }
	}

	// Plot number of penetrations
	cpgsci(1);
	cpgscf(2);
	sprintf(str,"%i",iNR);
	cpgtext((float)dT,fCorreleationPlotMaxY*0.97,str);
	// Plot error bars
	cpgsci(2);
	cpgmove((float)(dT-0.01),sqrt(iSignif/fDistributionCoeff));
	cpgdraw((float)(dT+0.01),sqrt(iSignif/fDistributionCoeff));	  
	cpgmove((float)(dT),sqrt(iSignif/fDistributionCoeff));
	cpgdraw((float)(dT),sqrt(adChi2Aver1[i]));
 	printf("%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.0f%%\n",dT,
	       fDum*sqrt(adChi2Aver1[i]),
	       fDum*sqrt(iSignif/fDistributionCoeff),
	       fDum*adShiftAver1[i],
	       fDum*(adShiftAver1[i]-adShiftAver2[i]),
	       fDum*(adShiftAver1[i]+adShiftAver2[i]),
	       i*100.0/THRESHOLD_NSAMPLES);
      }
      /*
      // PLot ShiftAverage
      cpgsci(3);
      for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	dDum=adShiftAver1[i]>0 ? adShiftAver1[i]:-adShiftAver1[i];
	if(i) cpgdraw((float)dT,dDum);
	else cpgmove((float)dT,dDum);
	cpgpt1((float)dT,dDum,2);
      }
      // PLot ShiftAverage +/- ShiftSD
       for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	 dDum=adShiftAver1[i]>0 ? adShiftAver1[i]:-adShiftAver1[i];
	 cpgmove((float)dT,dDum+adShiftAver2[i]);
	 cpgdraw((float)dT,dDum-adShiftAver2[i]);
	 cpgmove((float)dT+0.01,dDum+adShiftAver2[i]);
	 cpgdraw((float)dT-0.01,dDum+adShiftAver2[i]);
	 cpgmove((float)dT+0.01,dDum-adShiftAver2[i]);
	 cpgdraw((float)dT-0.01,dDum-adShiftAver2[i]);
      }
      */
     // PLot sqrt(Chi2Average)
      cpgsci(2);
      cpgslw(4);
      for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	dDum=sqrt(adChi2Aver1[i]);
	if(i) cpgdraw((float)dT,dDum);
	else cpgmove((float)dT,dDum);
	cpgpt1((float)dT,dDum,2);
      }
      cpgslw(2);
     /*
      // PLot sqrt(Chi2Average-SD)
      cpgsci(3);
      for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	fDum=adChi2Aver1[i]-adChi2Aver2[i]>0.0?sqrt(adChi2Aver1[i]-adChi2Aver2[i]):0.0;
	if(i) cpgdraw((float)dT,fDum);
	else cpgmove((float)dT,fDum);
      }      
      */

      if(iDoShuffling){

	// Random swaps of External frequency
	memset(adChi2Aver1,0,sizeof(double)*THRESHOLD_NSAMPLES);
	memset(adChi2Aver2,0,sizeof(double)*THRESHOLD_NSAMPLES);
	memset(adShiftAver1,0,sizeof(double)*THRESHOLD_NSAMPLES);
	memset(adShiftAver2,0,sizeof(double)*THRESHOLD_NSAMPLES);
	for(i=0;i<THRESHOLD_NSAMPLES;i++) 
	  memset(aauiDBins[i],0,sizeof(unsigned int)*CORRELATION_DISTRIBUTION_NBINS);
	
	for(k=0;k<CORRELATION_TEST_NTRIALS;k++){
	  for(l=0;l<CORRELATION_NSWAPS;l++){
	    for(;pRec[(iSwap1=(int)floor((iN*ran3(&lSeed))))].iHide&RECORD_HIDE_EXT;);
	    for(;pRec[(iSwap2=(int)floor((iN*ran3(&lSeed))))].iHide&RECORD_HIDE_EXT;);
	    dDum=pRec[iSwap1].dZR;
	    pRec[iSwap1].dZR=pRec[iSwap2].dZR;
	    pRec[iSwap2].dZR=dDum;
	  }
	  //cpgbbuf();
	  for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	    ThresholdRecords(iN,pRec,dT,&iNR,0);
	    SampleChi2(iN,pRec,N_SAMPLES,afS,&fMin,&fMinX,&fMax,&dAver1,&dAver2,&dAver4,0);
	    adChi2Aver1[i]+=fMin;
	    adChi2Aver2[i]+=(double)fMin*(double)fMin;
	    adShiftAver1[i]+=fMinX;
	    adShiftAver2[i]+=(double)fMinX*(double)fMinX;
	    if((l=(int)floor(fDistributionCoeff*fMin))<CORRELATION_DISTRIBUTION_NBINS){
	      aauiDBins[i][l]++;
	    }
	    else{
	      printf("Sample %i out of bound %i\n",i,l);
	    }
	    
	  //if(k%100){
	  //if(i) cpgdraw((float)dT,sqrt(fMin));
	  // else cpgmove((float)dT,sqrt(fMin));
	  //}
	  
	  }
	  //cpgebuf();
	}
	
	printf("#Shuffled\n");
	//      printf("#Tresh\tSR(C2)\tSR(C2-)\tSR(C2+)\tShift\tShift-\tShift+\tPcnt\n");
	printf("#Tresh\tSR(C2)\tSR(95%%)\tShift\tShift-\tShift+\tPcnt\n");
	if(iDoWedgeAnnotation) fDum=fWedgeMax-fWedgeMin;
	else fDum=1.0;
	
	for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	  adChi2Aver1[i]/=CORRELATION_TEST_NTRIALS;
	  adChi2Aver2[i]=sqrt(adChi2Aver2[i]/CORRELATION_TEST_NTRIALS-adChi2Aver1[i]*adChi2Aver1[i]);
	  adShiftAver1[i]/=CORRELATION_TEST_NTRIALS;
	  adShiftAver2[i]=sqrt(adShiftAver2[i]/CORRELATION_TEST_NTRIALS-adShiftAver1[i]*adShiftAver1[i]);
	  /*
	    printf("%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.0f%%\n",dT,
	    fDum*sqrt(adChi2Aver1[i]),
	    adChi2Aver1[i]-adChi2Aver2[i]>0.0?fDum*sqrt(adChi2Aver1[i]-adChi2Aver2[i]):0.0,
	    fDum*sqrt(adChi2Aver1[i]+adChi2Aver2[i]),
	    fDum*adShiftAver1[i],
	    fDum*(adShiftAver1[i]-adShiftAver2[i]),
	    fDum*(adShiftAver1[i]+adShiftAver2[i]),
	    i*100.0/THRESHOLD_NSAMPLES);
	  */
	}
	
	iNSamplesSignif=0.05*CORRELATION_TEST_NTRIALS;
	for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	  if((l=aauiDBins[i][0])>=iNSamplesSignif)iSignif=0;
	  else iSignif=-1;
	  cpgsci(8);
	  cpgmove((float)(dT+0.005+aauiDBins[i][0]*fDistributionPlotCoef),0);
	for(k=1;k<CORRELATION_DISTRIBUTION_NBINS;k++){
	  cpgdraw((float)(dT+0.005+aauiDBins[i][k]*fDistributionPlotCoef),sqrt(k/fDistributionCoeff));	  
	  if(iSignif<0) if((l+=aauiDBins[i][k])>=iNSamplesSignif) iSignif=k;
	}
	cpgsci(2);
	cpgmove((float)(dT-0.015),sqrt(iSignif/fDistributionCoeff));
	cpgdraw((float)(dT+0.015),sqrt(iSignif/fDistributionCoeff));	  
 	printf("%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.0f%%\n",dT,
	       fDum*sqrt(adChi2Aver1[i]),
	       fDum*sqrt(iSignif/fDistributionCoeff),
	       fDum*adShiftAver1[i],
	       fDum*(adShiftAver1[i]-adShiftAver2[i]),
	       fDum*(adShiftAver1[i]+adShiftAver2[i]),
	       i*100.0/THRESHOLD_NSAMPLES);
	}
	cpgsci(7);
	for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	  if(i) cpgdraw((float)dT,sqrt(adChi2Aver1[i]));
	  else cpgmove((float)dT,sqrt(adChi2Aver1[i]));
	}
	cpgsci(4);
	for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	  fDum=adChi2Aver1[i]-adChi2Aver2[i]>0.0?sqrt(adChi2Aver1[i]-adChi2Aver2[i]):0.0;
	  if(i) cpgdraw((float)dT,fDum);
	  else cpgmove((float)dT,fDum);
	}      
      }
      

      // PLot Shift
      cpgsci(4);
      cpgslw(4);
      for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	dDum=afShift1[i]>0 ? afShift1[i]:-afShift1[i];
	if(i) cpgdraw((float)dT,dDum);
	else cpgmove((float)dT,dDum);
	if(i==iSlice) printf("#Shift %f ",dDum);
      }  
      // PLot sqrt(Chi2)
      cpgsci(1);
      for(i=0,dT=dThresholdMin;i<THRESHOLD_NSAMPLES;i++,dT+=dStep){
	dDum=sqrt(adChi21[i]);
	if(i) cpgdraw((float)dT,(float)dDum);
	else cpgmove((float)dT,(float)dDum);
	if(i==iSlice) printf(" sqrt(Chi2) %f\n",dDum);
      }

      if(iPlotFileSpecial){
	cpgsls(2);

	for(i=0;i<PLOT_FILES_SPECIAL_N;i++){
	if(!(pF=fopen(strsPlotFilesSpecial[i],"r"))){
	  printf("ANA Cannot open file %s\n",strsPlotFilesSpecial[i]);
	  continue;
	}
	k=0;
	dT=dThresholdMin;
	while(!GetStringFromFile(pF,strDum,255,1)){
	  pc=strtok(strDum," \t\n");
	  if(pc) fX=atof(pc);
	  else{
	    printf("INT Incomplete X record\n");
	    break;
	  }
	  pc=strtok(NULL," \t\n");
	  if(pc) fY=atof(pc);
	  else{
	    printf("INT Incomplete Y record\n");
	    break;
	  }
	  if(k) cpgdraw(fX,fY);
	  else cpgmove(fX,fY);
	  k++;
	  dT+=dStep;
	}
	//plot legend
	//cpgmove((float)dThresholdMax*0.73,fCorreleationPlotMaxY*(1-fDum*(i+2)));
	//cpgdraw((float)dThresholdMax*0.84,fCorreleationPlotMaxY*(1-fDum*(i+2)));
	//cpgsci(1);
	//cpgsls(1);
	//cpgslw(2);
	//cpgtext((float)dThresholdMax*0.85,fCorreleationPlotMaxY*(1-fDum*(i+2.25)),strsPlotText[i]);
	//cpgslw(4);
	fclose(pF);
      }

 	cpgsls(1);
      }

      cpgslw(2);
    
      cpgband(0,0,fX,fY,&fX,&fY,&ch);
      cpgclos();
      cpgslct(iSaveID1);
      cpgunsa();

      GetValuesFromMap(iN,pRec,xd,yd,mp,&dThresholdMax,&dThresholdMin,0);
      ThresholdRecords(iN,pRec,dCorrelationThreshold,NULL,0);
      SampleChi2(iN,pRec,N_SAMPLES,afS,&fMin,&fMinX,&fMax,&dAver1,&dAver2,&dAver4,0);
      continue;
    }

    fShift=fX;
    if(fShift>1.0) fShift=1.0;
    if(fShift<0.0) fShift=0.0;
    ShiftWrapRec(iN,pRec,(double)fShift,0);
    PlotRecords(iPlotID,iN,pRec,fShift,fMinY,fMaxY,iPlotMode);
  }
  cpgslct(iChiID);
  cpgclos();
  cpgslct(iPlotID);
  cpgclos();
  cpgslct(iSaveID);

 bailout:

  cpgunsa();
  return(0);
}

#define VALUES_FROM_MAP_SPECIFIED 0
#define VALUES_FROM_MAP_RANDOM    1

int GetValuesFromMap(int iN,RECORD *pRec,int xd,int yd,double **mp,double *pdTMax,double *pdTMin,int iMode){
  double co,si;
  int i;
  unsigned long j;
  RECORD *pR,*pR1,stRec;
  double dX,dY;
  RECT stRect;

  if(iMode!=0){
    if(!stenciln || !stencil){
      printf("ANA EMPTY sstencil\n");
      return(1);
    }
  }

  stRect.x1=0;stRect.x2=xd-1;
  stRect.y1=0;stRect.y2=yd-1;

  co=cos(alpha);
  si=sin(alpha);
  *pdTMax=-(*pdTMin=1.0e20);
  for(i=0,pR=pRec;i<iN;i++,pR++){
    if(iMode==0){
      j=(unsigned long)pR->dX+xd*(unsigned long)pR->dY;
    }
    else{
      //for(;!IsPointInsideDomain((dX=(xd-1)*ran3(&lSeed)),(dY=(yd-1)*ran3(&lSeed)),iNDomainPoints,pDomain,&stRect););
      //j=(unsigned long)dX+xd*(unsigned long)dY;
      j=stencil[(int)((stenciln-1)*ran3(&lSeed))];
      //pR->dX=j%xd;
      //pR->dY=j/xd;      
    }

    dX=(mp[0][j]+xshift)*co-inversion*(mp[1][j]+yshift)*si;
    dY=inversion*(mp[1][j]+yshift)*co+(mp[0][j]+xshift)*si;
    pR->dRho=hypot(dY,dX);
    pR->dPhi=atan2(dY,dX)/(2*M_PI);
    if(pR->dPhi<0.0) pR->dPhi+=1.0;
    if(*pdTMax<pR->dRho) *pdTMax=pR->dRho;
    if(*pdTMin>pR->dRho) *pdTMin=pR->dRho;
  }
  //DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);

  for(i=1,pR1=pRec;i<iN;i++,pR1++){
    pR=pR1;
    while(pR->dRho>(pR+1)->dRho){
      SwapChunks(pR,pR+1,sizeof(RECORD),&stRec);
      if(pR==pRec) break;
      else pR--;
    }
  }
  return(0);
}

void SwapChunks(void *p1,void *p2, int iSize,void *p){
  memcpy(p,p2,iSize);
  memcpy(p2,p1,iSize);
  memcpy(p1,p,iSize);
}

int ThresholdRecords(int iN,RECORD *pRec,double dT,int *piNR,int iMode){
  int i,iNR,iNRPass;
  RECORD *pR;
  if(!iN || !pRec) return(-1);
  for(i=0,iNRPass=iNR=0,pR=pRec;i<iN;i++,pR++){
    if(!(pR->iHide&RECORD_HIDE_EXT)){
      if(pR->dRho<dT){
	pR->iHide|=RECORD_HIDE_INT;
	iNR++;
      }
      else{
	pR->iHide&=~RECORD_HIDE_INT;
	iNRPass++;
      }
    }
  }
  if(iMode&THRESHOLD_MODE_VERBOSE && iNR){
    printf("ANA Records below threshold %i: ",iNR);
    for(i=0,pR=pRec;i<iN;i++,pR++) if(pR->iHide&RECORD_HIDE_INT) printf("%i ",pR->iIndex);
    printf("\n");
  }
  if(piNR) *piNR=iNRPass;
  return(0);
}

#define REFINE_MIN_N 10

int SampleChi2(int iN,RECORD *pRec,int iNSamples,float *pfS,float *pfMin,float *pfMinX,float *pfMax,double *pdAver1,double *pdAver2,double *pdAver4,int iMode){
  double dX,dY;
  double dAver1,dAver2,dAver4;
  double dAver1Min,dAver2Min,dAver4Min;
  float fDum;
  int i;

  *pfMax=*pfMin=pfS[0]=Chi2ShiftRec(iN,pRec,*pfMinX=0.0,&dAver1,&dAver2,&dAver4);
  dAver1Min=dAver1;
  dAver2Min=dAver2;
  dAver4Min=dAver4;
  if(!(iMode&CORRELATION_MODE_NOSHIFT)){
    dY=1.0/iNSamples;
    for(i=1,dX=dY;i<iNSamples;i++,dX+=dY){
      if((fDum=pfS[i]=Chi2ShiftRec(iN,pRec,dX,&dAver1,&dAver2,&dAver4))<*pfMin){
	*pfMin=fDum;
	*pfMinX=dX;
	dAver1Min=dAver1;
	dAver2Min=dAver2;
	dAver4Min=dAver4;
      }
      if(*pfMax<fDum) *pfMax=fDum;
    }
    for(i=1,dX=(*pfMinX)-dY,dY/=REFINE_MIN_N;i<2*REFINE_MIN_N;i++,dX+=dY){
      if((fDum=Chi2ShiftRec(iN,pRec,dX,&dAver1,&dAver2,&dAver4))<*pfMin){
	*pfMin=fDum;
	*pfMinX=dX;
	dAver1Min=dAver1;
	dAver2Min=dAver2;
	dAver4Min=dAver4;
      }
    }
  }

  if(*pfMinX>0.5) *pfMinX-=1.0;

  if(pdAver1) *pdAver1=dAver1Min;
  if(pdAver2) *pdAver2=dAver2Min;
  if(pdAver4) *pdAver4=dAver4Min;

  if(iMode&THRESHOLD_MODE_VERBOSE) printf("ANA Min %f @ %f\n",*pfMin,*pfMinX);


  return(0);
}

double Chi2ShiftRec(int iN,RECORD *pRec,double dS,double *pdAver1,double *pdAver2,double *pdAver4){
  int i,iNR;
  double dDum;
  double dSum1=0,dSum2=0,dSum4=0;
  RECORD *pR;
  if(!iN || !pRec) return(-1);
  for(i=0,iNR=0,pR=pRec;i<iN;i++,pR++){
    if(!pR->iHide){
      dDum=pR->dPhi+dS;
      dSum1+=(dDum=dDum-floor(dDum)-pR->dZR);
      dSum2+=(dDum=dDum*dDum);
      dSum4+=dDum*dDum;
      iNR++;
    }
  }
  if(iNR){
    dSum1/=iNR;
    dSum2/=iNR;
    dSum4/=iNR;
  }
  if(pdAver1) *pdAver1=dSum1;
  if(pdAver2) *pdAver2=dSum2;
  if(pdAver4) *pdAver4=dSum4;
  return(dSum2);
}

double Chi2Rec(int iN,RECORD *pRec,double *pdAver1,double *pdAver2,double *pdAver4){
  int i,iNR;
  double dDum;
  double dSum1=0,dSum2=0,dSum4=0;
  RECORD *pR;
  if(!iN || !pRec) return(-1);
  for(i=0,iNR=0,pR=pRec;i<iN;i++,pR++){
    if(!pR->iHide){
      dSum1+=(dDum=pR->dPhiR);
      dSum2+=(dDum=dDum*dDum);
      dSum4+=dDum*dDum;
      iNR++;
    }
  }
  if(iNR){
    dSum1/=iNR;
    dSum2/=iNR;
    dSum4/=iNR;
  }
  if(pdAver1) *pdAver1=dSum1;
  if(pdAver2) *pdAver2=dSum2;
  if(pdAver4) *pdAver4=dSum4;
  return(dSum2);
}

int ShiftWrapRec(int iN,RECORD *pRec,double dS,int iMode){
  RECORD *pR;
  int i;
  double dDum;
  for(i=0,pR=pRec;i<iN;i++,pR++){
    if(!(pR->iHide&RECORD_HIDE_EXT)){
      dDum=pR->dPhi+dS;
      pR->dPhiR=dDum-floor(dDum);
    }
  }
  return(0);
}

double MinimizeChi2NoWrap(int iN,RECORD *pRec,double *pdMinChi2,double *pdMinS,int iMode){
  int i,iNR;
  double dDum;
  double dSum1=0,dSum2=0,dMin=0;
  RECORD *pR;
  if(!iN || !pRec) return(-1);
  for(i=0,iNR=0,pR=pRec;i<iN;i++,pR++){
    if(!pR->iHide){
      dSum1+=(dDum=pR->dPhiR-pR->dZR);
      dSum2+=dDum*dDum;
      iNR++;
    }
  }
  if(iNR){
    dSum1/=iNR;
    dSum2/=iNR;
    dMin=dSum2-dSum1*dSum1;
    for(i=0,iNR=0,pR=pRec;i<iN;i++,pR++){
      if(!pR->iHide){
	pR->dPhiR-=dSum1;
      }
    }
  }
  if(pdMinChi2) *pdMinChi2=dSum1;
  if(pdMinS) *pdMinS=dMin;
  return(dMin);
}


#define TEXT_SHIFT_X_PCT 0.005
#define TEXT_SHIFT_Y_PCT 0.005

int PlotRecords(int iId,int iN,RECORD *pRec,float fS,float fMinY,float fMaxY,int iMode){
  int i,iSaveID;
  char str[256];
  double dDum;
  RECORD *pR;
  float fX1,fX2,fY1,fY2;
  float fXTextShift,fYTextShift;
  int iNR,iNRH;
  float fMinX,fMaxX,fDum,fPlotMinX,fPlotMaxX;
  int iUseExtremaX=1;

  fMaxX=-(fMinX=BIG_FLOAT);
  for(iNR=iNRH=i=0,pR=pRec;i<iN;i++,pR++){
    if(!(pR->iHide&RECORD_HIDE_EXT)){
      fDum=pR->dPhiR;
      if(fMinX>fDum) fMinX=fDum;
      if(fMaxX<fDum) fMaxX=fDum;
      if(pR->iHide&RECORD_HIDE_INT) iNRH++;
      else iNR++;
    }
  }
  if(iUseExtremaX){
    fPlotMinX=fMinX;
    fPlotMaxX=fMaxX;
  }
  else{
    fPlotMinX=0.0;
    fPlotMaxX=1.0;
  }
  cpgqid(&iSaveID);
  cpgsave();
  cpgslct(iId);
  cpgsci(1);
  cpgenv(fPlotMinX,fPlotMaxX,fMinY<0.0 ? fMinY:0.0,fMaxY>1.0 ? fMaxY:1.0,0,1);

  dDum=Chi2ShiftRec(iN,pRec,fS,NULL,NULL,NULL);
  sprintf(str,"Chi2 %.4f(%.4f) Shift %.3f Items %i H %i T %i",dDum,sqrt(dDum),fS,iNR,iNRH,iN);
  cpglab("Local","External",str);
  cpgsci(4);
  cpgmove(0.0,0.0);
  cpgdraw(1.0,1.0);
  cpgsci(3);
  for(i=0,pR=pRec;i<iN;i++,pR++){
    if(!(pR->iHide&RECORD_HIDE_EXT)){
      if(pR->iHide&RECORD_HIDE_INT) cpgsci(HIDDEN_COLOR_N);
      else cpgsci(3);
      cpgpt1((float)pR->dPhiR,(float)pR->dZR,2); 
    }
    else{
      if(pR->iHide&RECORD_HIDE_EXT_FAKE){
	cpgsci(1);
	cpgpt1((float)pR->dPhiR,(float)pR->dZR,2); 
      }
    }
  }

  if(iMode&CORRELATION_MODE_ENUMERATE){
    cpgqwin(&fX1,&fX2,&fY1,&fY2);
    fXTextShift=(fX2-fX1)*TEXT_SHIFT_X_PCT;
    fYTextShift=(fY2-fY1)*TEXT_SHIFT_Y_PCT;
    cpgsch(0.65);
    for(i=0,pR=pRec;i<iN;i++,pR++){
      if(!(pR->iHide&RECORD_HIDE_EXT)){
	if(pR->iHide&RECORD_HIDE_INT) cpgsci(HIDDEN_COLOR_N);
	else cpgsci(5);
	sprintf(str,"%i",pR->iIndex);
	cpgtext((float)pR->dPhiR+fXTextShift,(float)pR->dZR+fYTextShift,str);
      }
      else{
	if(pR->iHide&RECORD_HIDE_EXT_FAKE){
	  cpgsci(1);
	  sprintf(str,"%i",pR->iIndex);
	  cpgtext((float)pR->dPhiR+fXTextShift,(float)pR->dZR+fYTextShift,str);
	}
      }
    }   
  }
  
  cpgunsa();
  cpgslct(iSaveID);
  return(0);
}

int PlotChi2(int iId,int iN,float *pfS,float fMin,float fMax,float fMinX,float fT){
  char str[256];
  int i,iSaveId;
  if(!iN || !pfS) return(1);
  cpgqid(&iSaveId);
  cpgsave();
  cpgslct(iId);
  cpgsci(1);
  cpgenv(0.0,1.0,fMin,fMax,0,1);
  sprintf(str,"Chi2 Min=%.5f(%.5f) @ %.3f T %0.3f N %0.3f",fMin,sqrt(fMin),fMinX,fT,dRhoNoise);
  cpglab("","",str);
  cpgsci(3);
  for(i=0;i<iN;i++){
    cpgpt1((float)i/iN,pfS[i],2);    
  }
  cpgunsa();
  cpgslct(iSaveId);
  return(0);
}
