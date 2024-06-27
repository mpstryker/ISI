#include "mapanm21.h"

void DisplayMap(int panel,int xd,int yd,double **mp){
  int i;
  float contour,dumf;
  float fmaxx,fminx,fmaxy,fminy,fmax,fmin;
  float wedge_offset=0.05,wedge_width=1.4;

  fmax=fmaxy=fmaxx=-(fmin=fminy=fminx=1.0e10);
  for(i=0;i<xd*yd;i++){
    fbuffer[i]=(float)(atan2(mp[1][i],mp[0][i]));
    dumf=fbufferx[i]=(float)(mp[0][i]);
    if(fmaxx<dumf) fmaxx=dumf;
    if(fminx>dumf) fminx=dumf;
    dumf=fbuffery[i]=(float)(mp[1][i]);
    if(fmaxy<dumf) fmaxy=dumf;
    if(fminy>dumf) fminy=dumf;
    dumf=fbufferr[i]=hypot(mp[0][i],mp[1][i]);
    if(fmax<dumf) fmax=dumf;
    if(fmin>dumf) fmin=dumf;
  }
  printf("DIS X(%f %f) Y(%f %f) R(%f %f)\n",fminx,fmaxx,fminy,fmaxy,fmin,fmax);
  cpgpanl(panel,1);
  cpgsch(1.0);
  cpglab("","",modestr[panel-1]);
  printf("%s\n",modestr[panel-1]);
  Pallet(2,1.0,0.5);
  cpgimag(fbuffer,xd,yd,1,xd,1,yd,-M_PI,M_PI,trans);
  cpgsch(2.6);
  cpgwedg("BI",wedge_offset,wedge_width, -180.0, 180.0,"");
  if(do_contours){
    cpgsci(0);
    contour=0.0;
    cpgcont(fbufferx,xd,yd,1,xd,1,yd,&contour,-1,trans);
    cpgsci(1);
    contour=0.0;
    cpgcont(fbuffery,xd,yd,1,xd,1,yd,&contour,-1,trans);
  }
  
  for(i=0;i<apinwheels_n;i++){
    cpgsci(3-(int)apinwheels[i].charge);
    cpgpt1((float)apinwheels[i].x,(float)apinwheels[i].y,2);
  }
    
  cpgpanl(panel,2);
  if(fmax==fmin){ fmin=0.0;fmax=2.0;}
  cpggray(fbufferr,xd,yd,1,xd,1,yd,fmax,fmin,trans);
  cpgsci(1);
  cpgsch(2.6);
  cpgwedg("B",wedge_offset,wedge_width, fmax, fmin,"");
  if(do_contours){
    cpgsci(0);
    contour=0.0;
    cpgcont(fbufferx,xd,yd,1,xd,1,yd,&contour,-1,trans);
    cpgsci(1);
    contour=0.0;
    cpgcont(fbuffery,xd,yd,1,xd,1,yd,&contour,-1,trans);
  }

  /*
  cpgpanl(panel,2);
  cpggray(fbufferx,xd,yd,1,xd,1,yd,fmaxx,fminx,trans);
  cpgsci(1);
  cpgsch(2.5);
  cpgwedg("B", -0.3, 1.3, fmaxx, fminx,"");
  if(do_contours){
    cpgsci(0);
    contour=0.0;
    cpgcont(fbufferx,xd,yd,1,xd,1,yd,&contour,-1,trans);
  }

  cpgpanl(panel,3);
  cpggray(fbuffery,xd,yd,1,xd,1,yd,fmaxy,fminy,trans);
  cpgsci(1);
  cpgsch(2.5);
  cpgwedg("B", -0.3, 1.3, fmaxy, fminy,"");
  if(do_contours){
    cpgsci(1);
    contour=0.0;
    cpgcont(fbuffery,xd,yd,1,xd,1,yd,&contour,-1,trans);
  }
  */
}

#define OD_WIN_SIZE 10.0

void DisplayOd(){
  int i,saveodid;
  float max,min,dumf;
  char str[128];
  int iODIndex1,iODIndex2,iTotal;
  float fODIntegral1,fODIntegral2;

  odid=cpgopen(device);
  cpgsubp(1,1);
  if(fabs((double)((pmaxy-pminy)/(pmaxx-pminx)))<1.0) 
    cpgpap(OD_WIN_SIZE,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
  else 
    cpgpap(OD_WIN_SIZE/fabs((double)((pmaxy-pminy)/(pmaxx-pminx))),(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
  //  cpgsch(1.0);
  //  cpgenv(pminx,pmaxx,pminy,pmaxy,1,1);
  cpgenv(pminx+1,pmaxx-1,pminy+1,pmaxy-1,1,-2);
  max=-(min=1.0e10);
  fODIntegral1=fODIntegral2=0.0;
  for(iODIndex1=iODIndex2=iTotal=0,i=0;i<XYdim;i++){
    dumf=fbuffer[i]=log(hypot(maps[0][i],maps[1][i])/hypot(maps[2][i],maps[3][i]));
    //    dumf=fbuffer[i]=hypot(maps[0][i],maps[1][i])-hypot(maps[2][i],maps[3][i]);
    if(dumf>0.0){
      iODIndex1++;
      fODIntegral1+=dumf;
    }
    if(dumf<0.0){
      iODIndex2++;
      fODIntegral2-=dumf;
    }
    iTotal++;
    if(min>dumf) min=dumf;
    if(max<dumf) max=dumf;
  }
  printf("DIS OD Index %f Integral %f\n",(float)(iODIndex1-iODIndex2)/(float)(iODIndex1+iODIndex2),(fODIntegral1-fODIntegral2)/(fODIntegral1+fODIntegral2));

  min/=2.0;
  max/=2.0;
  cpggray(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,max,min,trans);  
  cpgsci(2);
  dumf=0.0;
  cpgcont(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,&dumf,-1,trans);
  cpgsch(2.0);
  cpgwedg("B", -0.3, 1.3, max, min,"");

  printf("Save image?(y/n): ");
  fgets(str,4,stdin);
  if(*str == 'y'){
    sprintf(str,"%s%s",filestr1,".od.ps/CPS");
    printf("%s\n",str);
    sprintf(str,"%s%s",filestr2,".od.ps/CPS");
    printf("%s\n",str);
    saveodid=cpgopen("?");
    cpgsubp(1,1);
    cpgpap(OD_WIN_SIZE,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
    cpgenv(pminx+1,pmaxx-1,pminy+1,pmaxy-1,1,-2);
    cpggray(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,max,min,trans);  
    cpgsci(2);
    dumf=0.0;
    cpgcont(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,&dumf,-1,trans);
    cpgsch(2.0);
    cpgwedg("B", -0.3, 1.3, max, min,"");
    cpgclos();
    cpgslct(odid);
  }
}

void DisplayCombo(){
  int i,saveodid,localid,callerid;
  float max,min,dumf;
  char str[128];

  cpgqid(&callerid);
  localid=cpgopen(device);
  cpgsubp(1,1);
  cpgpap(OD_WIN_SIZE,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
  //  cpgsch(1.0);
  //  cpgenv(pminx,pmaxx,pminy,pmaxy,1,1);
  cpgenv(pminx+1,pmaxx-1,pminy+1,pmaxy-1,1,-2);
  for(i=0;i<XYdim;i++){
    dumf=fbuffer[i]=log(hypot(maps[0][i],maps[1][i])/hypot(maps[2][i],maps[3][i]));
  }
  cpggray(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,max,min,trans);  
  cpgsci(2);
  dumf=0.0;
  cpgcont(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,&dumf,-1,trans);
  cpgsch(2.0);
  cpgwedg("B", -0.3, 1.3, max, min,"");

  printf("Save image?(y/n): ");
  fgets(str,4,stdin);
  if(*str == 'y'){
    sprintf(str,"%s%s",filestr1,".delay.ps/CPS");
    printf("%s\n",str);
    sprintf(str,"%s%s",filestr2,".delay.ps/CPS");
    printf("%s\n",str);
    saveodid=cpgopen("?");
    cpgsubp(1,1);
    cpgpap(OD_WIN_SIZE,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
    cpgenv(pminx+1,pmaxx-1,pminy+1,pmaxy-1,1,-2);
    cpggray(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,max,min,trans);  
    cpgsci(2);
    dumf=0.0;
    cpgcont(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,&dumf,-1,trans);
    cpgsch(2.0);
    cpgwedg("B", -0.3, 1.3, max, min,"");
    cpgclos();
    cpgslct(odid);
  }
  cpgslct(callerid);
}

#define DOMINANCE_HISTOGRAM_WIN_SIZE 7
#define DOMINANCE_HISTOGRAM_WIN_RATIO 0.6
#define BIN_MIN_N 50
#define BIN_MAX_N 400
#define BIN_STEP 0.02

int PlotDominanceHistogram(int iX,int iY,double *pdMX1,double *pdMY1,double *pdMX2,double *pdMY2,int iMode){
  int iSaveID,iID;
  float fX,fY;
  char ch;
  double *pdMR1=NULL,*pdMR2=NULL;
  double *pd1,*pd2;
  unsigned long ulXY;
  double dMax,dMin,dDum,dStep;
  float *pfBins=NULL;
  int iNBins;
  int i,j;
  float fBinMax,fBinMin;

  ulXY=(unsigned long)iX*(unsigned long)iY;
  if(!(pdMR1=(double*)malloc(ulXY*sizeof(double)))){
    printf("DIS Cannot allocate for pdMR1\n");
    goto bailout;
  }
  if(!(pdMR2=(double*)malloc(ulXY*sizeof(double)))){
    printf("DIS Cannot allocate for pdMR2\n");
    goto bailout;
  }
  dMax=dMin=hypot(pdMX1[0],pdMY1[0])-hypot(pdMX2[0],pdMY2[0]);
  for(i=0,pd1=pdMR1,pd2=pdMR2;i<ulXY;i++,pd1++,pd2++){
    *pd1=hypot(pdMX1[i],pdMY1[i]);
    *pd2=hypot(pdMX2[i],pdMY2[i]);
    if(dMin>(dDum=*pd1-*pd2)) dMin=dDum;
    if(dMax<dDum) dMax=dDum;
  }
  if(dMax==dMin){
    printf("DIS Flat images\n");
    goto bailout;
  }

  iNBins=(dMax-dMin)/BIN_STEP;
  if(iNBins>BIN_MAX_N) iNBins=BIN_MAX_N;
  if(iNBins<BIN_MIN_N) iNBins=BIN_MIN_N;
  dStep=(dMax-dMin)/iNBins;

  if(!(pfBins=(float*)calloc(iNBins,sizeof(float)))){
    printf("DIS Cannot callocate for pfBins\n");
    goto bailout;
  }
  for(i=0,pd1=pdMR1,pd2=pdMR2;i<ulXY;i++,pd1++,pd2++){
    if(*pd1<dDominanceThreshold && *pd2<dDominanceThreshold) continue;
    j=(int)(((*pd1-*pd2)-dMin)/dStep);
    if(j>=0 && j<iNBins) pfBins[j]++;
  }
  fBinMin=0.0;
  fBinMax=0.0;
  for(i=0;i<iNBins;i++){
    if(fBinMax<pfBins[i]) fBinMax=pfBins[i];
  }

  cpgsave();
  cpgqid(&iSaveID);

  iID=cpgopen("/xw");
  cpgsubp(1,1);
  cpgpap(DOMINANCE_HISTOGRAM_WIN_SIZE,DOMINANCE_HISTOGRAM_WIN_RATIO);
  cpgsci(1);
  cpgenv((float)dMin,(float)dMax,fBinMin,fBinMax,0,1);

  cpgmove((float)dMin,pfBins[0]);
  for(i=1;i<iNBins;i++){
    cpgdraw((float)(dMin+i*dStep),pfBins[i]);
  }  

  cpgband(0,0,fX,fY,&fX,&fY,&ch);
  cpgask(0);
  cpgclos();
  cpgslct(iSaveID);
  cpgunsa();

 bailout:
  SAFE_FREE(pfBins);
  SAFE_FREE(pdMR1);
  SAFE_FREE(pdMR2);
  return(0);
}
