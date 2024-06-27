#include "iman1.h"
#include <functions.h>

#define GAUSSF(x,sigma) ((float)exp(-(double)(0.5*(x)*(x)/((sigma)*(sigma)))))
#define GAUSS(x,sigma) ((float) (exp(-(double)(0.5*(x)*(x)/((sigma)*(sigma))))/(sqrt(2.0*M_PI)*(sigma))) )

//#define INCLUDE_HEADER

int Statistics(TRAIN *tr){
  unsigned long *bins=(unsigned long *)0;
  unsigned long integral,duml;
  double dumd,dumd2,integral2,integral4,integral6,integral8;
  float maxy,miny,maxd,sigma,dumf,sigma2,sigma4,sigma6,sigma8;
  float inv2,inv4,inv6,inv8;
  float ppminx,ppmaxx,coeff1=0.9,int1,int2;
  int binmax,binmin;
  int n_bins;
  int iframe,fframe,framen;
  int *extrabins=(int *)0,n_extrabins=0;
  int i,j,dumi;
  int frame_print_count1=10,frame_print_count2=100,frame_print_count3=500;
  CAR *pcar1,*pcar2;
  double lambdaa;
  int lambdamin,lambdamax;
  double *probability=(double *)0,*trueprobability=(double *)0,multiplier=122.0;
  int *lambda_stat=(int *)0;
  unsigned short *pus1,*pus2;
  void *pv1,*pv2;
  int average_frames=AVERAGE;
  int iRecordSize=sizeof(unsigned short),iNRecords;
  int iReplaceSchemeOrig;
  unsigned long ulBinTotal;
  double *pdM1=NULL,*pdM2=NULL,*pd1,*pd2;
  double dWellDepth,dWellDepthNorm,dWellCoeff;

  iReplaceSchemeOrig=replace_scheme;
  replace_scheme=REPLACE_RANDOM;
  ulBinTotal=ulBinX*ulBinY*ulBinT;

  n_bins= statistics_num_bins;
  binmin=n_bins;
  binmax=-n_bins;
  iframe=1;
  fframe=tr->max_n;
  //  fframe=2500;
  if(!(bins=(unsigned long *)calloc(1+2*n_bins,sizeof(unsigned long)))){
    printf("Cannot allocate for bins\n");
    return(1);
  }
  if(!(lambda_stat=(int *)calloc(1<<16,sizeof(int)))){
    printf("Cannot allocate for lambda_stat\n");
    return(1);
  }

  iNRecords=tr->nFrameImageSize/iRecordSize;
#ifdef INCLUDE_HEADER
  iNRecords+=tr->nFrameHeaderSize/iRecordSize;
#endif

  if(!(pdM1=(double *)calloc(iNRecords,sizeof(double)))){
    printf("Cannot allocate for pdM1\n");
    return(1);
  }
  if(!(pdM2=(double *)calloc(iNRecords,sizeof(double)))){
    printf("Cannot allocate for pdM2\n");
    return(1);
  }

  if(!(pcar1=AddFrame(tr,iframe,average_frames))){
    printf("STA AddFrame return 0 for iframe %i\n",iframe);     
    return(2);
  }
  pv1=pcar1->pFrame;

  lambdamin=1L<<16;
  lambdamax=0;
  for(j=0,lambdaa=0.0,pus1=(unsigned short*)(pv1+tr->nFrameHeaderSize);j<XYdim;j++){
    lambdaa+=(double)(dumi=(int)*(pus1++));
    lambda_stat[dumi]++;
    if(lambdamin>dumi) lambdamin=dumi;
    if(lambdamax<dumi) lambdamax=dumi;
  }
  lambdaa/=(double)XYdim;

  
  for(j=0,pd1=pdM1,pd2=pdM2,pus1=(USHORT*)(pv1+tr->nFrameSize);j<iNRecords;j++){
    *(pd1++)=(dumd=(double)*(--pus1));
    *(pd2++)=dumd*dumd;
  }
  

  /*
  maxy=0;
  miny=1<<16;
  for(j=lambdamin;j<=lambdamax;j++){
    dumf=(float)lambda_stat[j];
    if(miny>dumf) miny=dumf;
    if(maxy<dumf) maxy=dumf;
  }    
  printf("Lstat min=%f max=%f\n",miny,maxy);
  cpgopen("/xw");
  cpgpap(8.0,0.6);
  cpgenv((float)lambdamin,(float)lambdamax,miny,maxy,0,1);
  cpgsci(2);
  cpgmove((float)lambdamin,(float)lambda_stat[lambdamin]);
  for(i=lambdamin+1;i<=lambdamax;i++){
    cpgdraw((float)i,(float)lambda_stat[i]);
  }
  cpgclos();
  if(mainid) cpgslct(mainid);
  */
  maxy=5000.0;
  miny=0.0;
  printf("STA ");
  for(i=1+iframe,framen=0;i<fframe;i++,framen++){
    if(!(i%frame_print_count1)){
      printf(".");
      fflush(stdout);
      if(!(i%frame_print_count2)){
	printf("|");
	fflush(stdout);
	if(!(i%frame_print_count3)){
	  printf("%i\n",i);
	  fflush(stdout);
	}
      }
    }
    
    if(LockFrame(tr,pcar1,TRUE)){
      printf("STA Cannot lock frame %i\n",i-1);     
    }

    if(!(pcar2=AddFrame(tr,i,average_frames))){
      printf("STA AddFrame return 0 for frame %i\n",i);     
      return(2);
    }
    pv2=pcar2->pFrame;

    pd1=pdM1;
    pd2=pdM2;

    for(j=0,pus1=(USHORT*)(pv1+tr->nFrameSize),pus2=(USHORT*)(pv2+tr->nFrameSize);j<iNRecords;j++){

      dumi=(int)*(--pus2)-(int)*(--pus1);
      *(pd1++)+=(dumd=(double)(*pus2));
      *(pd2++)+=dumd*dumd;

      if(binmin>dumi) binmin=dumi;
      if(binmax<dumi) binmax=dumi;
      if(dumi>=-n_bins && dumi<=n_bins){
	bins[dumi+n_bins]++;
      }
      else{
	if(!(extrabins=(int *)realloc(extrabins,(++n_extrabins)*sizeof(int)))){
	  printf("Cannot reallocate for extrabins\n");
	  return(1);
	}
	extrabins[n_extrabins-1]=dumi;
	if(do_verbose) printf("Extra: %i %i\n",n_extrabins,dumi);
      }
    }
    pv1=pv2;
    if(LockFrame(tr,pcar1,FALSE)){
      printf("STA Cannot unlock frame %i\n",i-1);     
    }
    pcar1=pcar2;
  }
  printf("\n");

  dumd=1.0/(framen+1);
  dWellCoeff=0.0;
  dWellDepthNorm=0.0;
  for(j=0,pd1=pdM1,pd2=pdM2;j<iNRecords;j++,pd1++,pd2++){
    *pd1*=dumd;
    *pd2*=dumd;
    *pd2-=(*pd1)*(*pd1);
    dWellCoeff+=(*pd1)*(*pd1)/(*pd2);
    dWellDepthNorm+=(*pd1)/(*pd2);
       if(!(j%100)) printf("%f\n",sqrt(*pd2));
  }
  dWellCoeff/=(double)iNRecords*(double)ulBinTotal;
  dWellDepth=dWellCoeff*4096.0;
  dWellDepthNorm*=4096.0/(double)iNRecords;
  printf("STA Well depth %f coeff %f (%i)\n",dWellDepth,dWellCoeff,framen+1);
  printf("STA Normalized Well depth %f coeff %f\n",dWellDepthNorm,dWellDepthNorm/4096.0);

  //   DisplayMap(tr,pdM1,pdM1);
  //  DisplayMap(tr,pdM2,pdM2);


  if(binmin<-n_bins) binmin=-n_bins;
  if(binmax>n_bins) binmax=n_bins;
  integral=0;
  integral2=integral4=integral6=integral8=0.0;
  maxd=0.0;
  for(i=binmin;i<=binmax;i++){
    integral+=(duml=bins[i+n_bins]);
    dumd=(double)i*(double)i;
    integral2+=(dumd2=(double)duml*dumd);
    dumd2*=dumd;
    integral4+=dumd2;
    dumd2*=dumd;
    integral6+=dumd2;
    dumd2*=dumd;
    integral8+=dumd2;
    dumf=(float)duml;
    if(maxd<dumf) maxd=dumf;
  }
  sigma2=(float)sqrt(integral2/(double)integral);
  sigma4=(float)sqrt(sqrt(integral4/(double)integral/3.0));
  sigma6=(float)cbrt(sqrt(integral6/(double)integral/15.0));
  sigma8=(float)sqrt(sqrt(sqrt(integral8/(double)integral/105.0)));

  inv2=integral2/(double)integral;
  inv4=(float)sqrt((integral4-integral2)/(double)integral/3.0);
  inv6=(float)cbrt((integral6-5.0*integral4+4.0*integral2)/(double)integral/15.0);
  inv8=(float)sqrt(sqrt((integral8-14.0*integral6+49.0*integral4-36.0*integral2)/(double)integral/105.0));

  sigma=integral/maxd/sqrt(2.0*M_PI);
  maxy=maxd;
  printf("STA Integral=%li(%i %f)\n",integral,framen,(double)XYdim*(double)framen-n_extrabins);
  printf("STA Bin: xmin=%i xmax=%i maxd=%f sigma=%f(%f %f %f %f)\n",binmin,binmax,maxd,sigma,sigma2,sigma4,sigma6,sigma8);
  printf("STA Inv: %f %f %f %f\n",sqrt(inv2),sqrt(inv4),sqrt(inv6),sqrt(inv8));
  if(!(probability=(double *)calloc(2*n_bins+1,sizeof(double)))){
    printf("STA Cannot allocate for probability\n");
    return(1);
  }
  if(!(trueprobability=(double *)calloc(2*n_bins+1,sizeof(double)))){
    printf("STA Cannot allocate for trueprobability\n");
    return(1);
  }
  printf("STA Lmin=%i Lmax=%i La=%f\n",lambdamin,lambdamax,lambdaa);
 
  for(i=binmin;i<=binmax;i++){
    probability[i+n_bins]=0.6*multiplier*besselIne((int)((double)i*multiplier),2.0*lambdaa*multiplier);
    /*    
    for(dumd=0.0,j=lambdamin;j<=lambdamax;j++){
      dumd+=(double)lambda_stat[j]*besselIne((int)(i*multiplier),2.0*(double)j*multiplier);
    }
    trueprobability[i+n_bins]=dumd/(double)XYdim;
    */
    //    printf("P %4i (%i %f) %e %e\n",i,(int)((double)i*multiplier),2.0*lambdaa*multiplier,trueprobability[i+n_bins],probability[i+n_bins]);
    probability[i+n_bins]*=(double)integral;
    trueprobability[i+n_bins]*=(double)integral;
  }

  ppminx=binmin*coeff1;
  ppmaxx=binmax*coeff1;
  cpgopen("/xw");
  cpgpap(10.0,0.6);
  cpgenv(ppminx,ppmaxx,miny,maxy,0,1);
  cpgsci(2);
  cpgmove((float)binmin,(float)bins[binmin]);
  for(i=binmin+1;i<=binmax;i++){
    cpgdraw((float)i,(float)bins[i+n_bins]);
  }
  cpgsci(1);
  cpgmove((float)binmin,integral*GAUSS((float)binmin,sigma2));
  for(i=binmin+1;i<=binmax;i++){
    cpgdraw((float)i,integral*GAUSS((float)i,sigma2));
  }

  cpgsci(3);
  cpgmove((float)binmin,(float)probability[binmin]);
  for(i=binmin+1;i<=binmax;i++){
    cpgdraw((float)i,(float)probability[i+n_bins]);
  }
  cpgsci(4);
  cpgmove((float)binmin,(float)trueprobability[binmin]);
  for(i=binmin+1;i<=binmax;i++){
    cpgdraw((float)i,(float)trueprobability[i+n_bins]);
  }

  miny=maxy;
  maxy=-maxy;
  for(i=binmin;i<=binmax;i++){
    dumf=bins[i+n_bins]-integral*GAUSS((float)i,sigma2);
    if(miny>dumf) miny=dumf;
    if(maxy<dumf) maxy=dumf;
  }

  cpgsci(1);
  cpgenv(ppminx,ppmaxx,miny,maxy,0,1);
  cpgsci(2);
  int1=(float)fabs(bins[binmin]-integral*GAUSS((float)binmin,sigma2));
  cpgmove((float)binmin,bins[binmin]-integral*GAUSS((float)binmin,sigma2));
  for(i=binmin+1;i<=binmax;i++){
    int1+=(float)fabs(dumf=bins[i+n_bins]-integral*GAUSS((float)i,sigma2));
    cpgdraw((float)i,dumf);
  }
  int1/=(float)(binmax-binmin+1);

  cpgsci(3);
  int2=(float)fabs((double)(bins[binmin]-probability[binmin]));
  cpgmove((float)binmin,(float)(bins[binmin]-probability[binmin]));
  for(i=binmin+1;i<=binmax;i++){
    int2+=(float)fabs(dumf=(float)(bins[i+n_bins]-probability[i+n_bins]));
    cpgdraw((float)i,dumf);
  }
  int2/=(float)(binmax-binmin+1);
  
  printf("I1=%f I2=%f\n",int1,int2);
  cpgsci(4);
  cpgmove((float)binmin,probability[binmin]-integral*GAUSS((float)binmin,sigma2));
  for(i=binmin+1;i<=binmax;i++){
    cpgdraw((float)i,probability[i+n_bins]-integral*GAUSS((float)i,sigma2));
  }

  cpgclos();
  if(mainid) cpgslct(mainid);

  if(pdM1) free(pdM1);
  if(pdM2) free(pdM2);
  free(lambda_stat);
  free(trueprobability);
  free(probability);
  free(extrabins);
  free(bins);

  replace_scheme=iReplaceSchemeOrig;
  return(0);
}

