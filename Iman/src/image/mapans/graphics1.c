/* Map analysis                 */
/* Written by V.Kalatsky        */
/* 06Feb2001                    */
/* Modified 10May2001           */

#include "mapans1.h"

void start_graphics(int mode,char *str,int *id){
  int i,lmode;

  lmode=DISPLAY_MODE_MASK&mode;
  *id=cpgopen(str);
  if((mode&DISPLAY_MODE_SINGLE)){
    cpgpap(window_size,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
    if(iFullSizeViewPort){
      cpgsvp(0.0,1.0,0.0,1.0);
      cpgswin(pminx,pmaxx,pminy,pmaxy);
    }
    else{
      cpgsch(0.5);
      cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
      cpgsch(3.5);
    }
  }
  else{
    cpgsubp(lmode,2);
    cpgpap(window_size*lmode,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx)))*2.0/(double)lmode);
    if(iDoLargeWedge) cpgsch(0.75);
    else cpgsch(0.5);
    //    cpgsch(0.5);
    for(i=0;i<lmode*2;i++) cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
    cpgsch(3.5);
  }
}

void end_graphics(int id){
  cpgslct(id);
  cpgclos();
}

void Pallet(int type,float contra,float bright){
  int r=9,h=5,hex=7,dozen=13,eighteen=19,twentyfour=25;

  float rl[]={ 0.0, 0.125,0.25, 0.375,0.5,  0.625,0.75,0.875,1.0};
  float rr[]={ 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
  float rg[]={ 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
  float rb[]={ 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};

  float hl[]={ 0.0, 0.2, 0.4, 0.6, 1.0};
  float hr[]={ 0.0, 0.5, 1.0, 1.0, 1.0};
  float hg[]={ 0.0, 0.0, 0.5, 1.0, 1.0};
  float hb[]={ 0.0, 0.0, 0.0, 0.3, 1.0};
  
  // RGB Green in the middle perfect hexagon section of the cube
  float hexl1[]={ 0.0, 1.0/6.0, 1.0/3.0,  0.5,2.0/3.0, 5.0/6.0, 1.0};
  float hexr1[]={ 1.0, 1.0,     0.5,      0.0,0.0,     0.5,     1.0};
  float hexg1[]={ 0.0, 0.5,     1.0,      1.0,0.5,     0.0,     0.0};
  float hexb1[]={ 0.5, 0.0,     0.0,      0.5,1.0,     1.0,     0.5};
  
    // RGB Blue in the middle perfect hexagon section of the cube
  float hexl2[]={ 0.0, 1.0/6.0, 1.0/3.0, 0.5, 2.0/3.0, 5.0/6.0, 1.0};
  float hexr2[]={ 1.0, 0.5,     0.0,     0.0, 0.5,     1.0,     1.0};
  float hexg2[]={ 0.5, 1.0,     1.0,     0.5, 0.0,     0.0,     0.5};
  float hexb2[]={ 0.0, 0.0,     0.5,     1.0, 1.0,     0.5,     0.0};
   
  // HSL Blue in the middle S=L=1.0
  float hexl3[]={ 0.0, 1.0/6.0, 1.0/3.0, 0.5, 2.0/3.0, 5.0/6.0, 1.0};
  float hexr3[]={ 1.0, 0.0,     0.0,     0.0, 1.0,     1.0,     1.0};
  float hexg3[]={ 1.0, 1.0,     1.0,     0.0, 0.0,     0.0,     1.0};
  float hexb3[]={ 0.0, 0.0,     1.0,     1.0, 1.0,     0.0,     0.0};
  /*
  // HSL Blue in the middle S=L=1.0
  hex=11;
  float hexl3[]={ 0.0, 1.0/8.0, 1.0/8.0, 2.0/8.0, 3.0/8.0, 0.5,  5.0/8.0, 6.0/8.0, 7.0/8.0, 7.0/8.0,1.0};
  float hexr3[]={ 0.0, 0.0,     1.0,     0.0,     0.0,     0.0,  1.0,     1.0,     1.0,     0.0,    0.0};
  float hexg3[]={ 0.0, 0.0,     1.0,     1.0,     1.0,     0.0,  0.0,     0.0,     1.0,     0.0,    0.0};
  float hexb3[]={ 0.0, 0.0,     0.0,     0.0,     1.0,     1.0,  1.0,     0.0,     0.0,     0.0,    0.0};
  */
  // HSL Blue @ the left edge S=L=1.0
  float hexr3_1[]={ 0.0, 1.0,     1.0,     1.0, 0.0,     0.0,     0.0};
  float hexg3_1[]={ 0.0, 0.0,     0.0,     1.0, 1.0,     1.0,     0.0};
  float hexb3_1[]={ 1.0, 1.0,     0.0,     0.0, 0.0,     1.0,     1.0};

  float dozenl[]={0.0, 1.0/12.0,1.0/6.0, 0.25,1.0/3.0, 5.0/12.0, 0.5, 7.0/12.0, 2.0/3.0, 0.75,5.0/6.0, 11.0/12.0, 1.0};
  float dozenr[]={0.0, 0.0,     0.5,     1.0, 1.0,     0.5,      0.0, 0.0,      0.5,     1.0,1.0,     0.5,       0.0};
  float dozeng[]={1.0, 0.5,     0.0,     0.0, 0.5,     1.0,      1.0, 0.5,      0.0,     0.0,0.5,     1.0,       1.0};
  float dozenb[]={0.5, 1.0,     1.0,     0.5, 0.0,     0.0,      0.5, 1.0,      1.0,     0.5,0.0,     0.0,       0.5};

  float eighteenl[]={ 0.0, 1.0/18.0,      1.0/9.0,  1.0/6.0,2.0/9.0, 5.0/18.0,     1.0/3.0, 7.0/18.0,      4.0/9.0,  0.5,5.0/9.0, 11.0/18.0,     2.0/3.0, 13.0/18.0,      7.0/9.0,  5.0/6.0,8.0/9.0, 17.0/18.0,     1.0};
  float eighteenr[]={ 1.0, 1.0,          0.5,      0.0,0.0,     0.5,         1.0, 1.0,          0.5,      0.0,0.0,     0.5,         1.0, 1.0,          0.5,      0.0,0.0,     0.5,         1.0};
  float eighteeng[]={ 0.0, 0.5,          1.0,      1.0,0.5,     0.0,         0.0, 0.5,          1.0,      1.0,0.5,     0.0,         0.0, 0.5,          1.0,      1.0,0.5,     0.0,         0.0};
  float eighteenb[]={ 0.5, 0.0,          0.0,      0.5,1.0,     1.0,         0.5, 0.0,          0.0,      0.5,1.0,     1.0,         0.5, 0.0,          0.0,      0.5,1.0,     1.0,         0.5};

  float twentyfourl[]={0.0, 1.0/24.0,1.0/12.0, 0.125,1.0/6.0, 5.0/24.0, 0.25, 7.0/24.0, 1.0/3.0, 0.375,5.0/12.0, 11.0/24.0, 0.5, 13.0/24.0,7.0/12.0, 5.0/8.0,2.0/3.0, 17.0/24.0, 3.0/4.0, 19.0/24.0, 5.0/6.0, 7.0/8.0,11.0/12.0, 23.0/24.0, 1.0};
  float twentyfourr[]={0.0, 0.0,     0.5,     1.0, 1.0,     0.5,      0.0, 0.0,      0.5,     1.0,1.0,     0.5,       0.0, 0.0,     0.5,     1.0, 1.0,     0.5,      0.0, 0.0,      0.5,     1.0,1.0,     0.5,       0.0};
  float twentyfourg[]={1.0, 0.5,     0.0,     0.0, 0.5,     1.0,      1.0, 0.5,      0.0,     0.0,0.5,     1.0,       1.0, 0.5,     0.0,     0.0, 0.5,     1.0,      1.0, 0.5,      0.0,     0.0,0.5,     1.0,       1.0};
  float twentyfourb[]={0.5, 1.0,     1.0,     0.5, 0.0,     0.0,      0.5, 1.0,      1.0,     0.5,0.0,     0.0,       0.5, 1.0,     1.0,     0.5, 0.0,     0.0,      0.5, 1.0,      1.0,     0.5,0.0,     0.0,       0.5};

  switch(type){
  case -2:
    cpgctab(rl, rr, rg, rb, r, contra, bright);
    break;
  case -1:
    cpgctab(hl, hr, hg, hb, h, contra, bright);
    break;
  case 0:
    cpgctab(hexl1, hexr1, hexg1, hexb1, hex, contra, bright);
    break;
  case 1:
    cpgctab(hexl2, hexr2, hexg2, hexb2, hex, contra, bright);
    break;
  case 2:
    //    printf("PALLET 2\n");
    cpgctab(hexl3, hexr3, hexg3, hexb3,hex, contra, bright);
    break;
  case 3:
    cpgctab(dozenl, dozenr, dozeng, dozenb, dozen, contra, bright);
    break;
  case 4:
    cpgctab(eighteenl, eighteenr, eighteeng, eighteenb, eighteen, contra, bright);
    break;
  case 5:
    cpgctab(twentyfourl, twentyfourr, twentyfourg, twentyfourb, twentyfour, contra, bright);
    break;
  case 2+PALLET_SHIFT_SHIFT:
    cpgctab(hexl3, hexr3_1, hexg3_1, hexb3_1, hex, contra, bright);
    break;
  case 2+PALLET_SPECIAL:
    cpgctab(hexl3, hexr3_1, hexg3_1, hexb3_1, hex, contra, bright);
    break;
  default:
    printf("Unknown pallet %i\n",type);
  }
}

void MakePallet(int type,float fS,float fPr,float fPo){
  float fContra=1.0;
  float fBright=0.5;
  float fMin,fMax;
  int iMin=0,iMax=5,iN=0;
  float fP;
  float *pfBuffer=NULL;
  int i,j;
  float *pfPalletL;
  float *pfPalletR;
  float *pfPalletG;
  float *pfPalletB;
  /*
  float pfPalletL[]={ 0.0, 1.0/6.0, 1.0/3.0, 0.5, 2.0/3.0, 5.0/6.0, 1.0};
  float pfPalletR[]={ 0.0, 1.0,     1.0,     1.0, 0.0,     0.0,     0.0};
  float pfPalletG[]={ 0.0, 0.0,     0.0,     1.0, 1.0,     1.0,     0.0};
  float pfPalletB[]={ 1.0, 1.0,     0.0,     0.0, 0.0,     1.0,     1.0};
  */
  fP=fS+fPr+fPo;
  fMin=fPr/fP;
  fMax=(fPr+fS)/fP;
  iMin=(int)floor((double)fMin*6.0);
  iMax=(int)floor((double)fMax*6.0);
  iN=iMax-iMin+2;
  if(!(pfBuffer=(float*)malloc(4*iN*sizeof(float)))){
    printf("GRA Cannot allocate for pfBuffer. No pallet change\n");
    return;
  }
  pfPalletL=pfBuffer;
  pfPalletR=pfBuffer+iN;
  pfPalletG=pfBuffer+2*iN;
  pfPalletB=pfBuffer+3*iN;
  pfPalletL[0]=0.0;
  pfPalletL[iN-1]=1.0;
  for(i=iMin+1,j=1;i<iMax+1;i++,j++)
    pfPalletL[j]=((float)i/6.0-fMin)/(fMax-fMin);
  if(type==PALLET_REDUCED){
    if(iMin==0){
      pfPalletR[0]=fMin*6.0;
      pfPalletG[0]=0.0;
      pfPalletB[0]=1.0;
    }
    if(iMin==1){
      pfPalletR[0]=1.0;
      pfPalletG[0]=0.0;
      pfPalletB[0]=1.0-(fMin*6.0-1.0);
    }
    if(iMin==2){
      pfPalletR[0]=1.0;
      pfPalletG[0]=fMin*6.0-2.0;
      pfPalletB[0]=0.0;
    }
    if(iMin==3){
      pfPalletR[0]=1.0-(fMin*6.0-3.0);
      pfPalletG[0]=1.0;
      pfPalletB[0]=0.0;
    }
    if(iMin==4){
      pfPalletR[0]=0.0;
      pfPalletG[0]=1.0;
      pfPalletB[0]=fMin*6.0-4.0;
    }
    if(iMin==5){
      pfPalletR[0]=0.0;
      pfPalletG[0]=1.0-(fMin*6.0-5.0);
      pfPalletB[0]=1.0;
    }
    if(iMax==0){
      pfPalletR[iN-1]=fMax*6.0;
      pfPalletG[iN-1]=0.0;
      pfPalletB[iN-1]=1.0;
    }
    if(iMax==1){
      pfPalletR[iN-1]=1.0;
      pfPalletG[iN-1]=0.0;
      pfPalletB[iN-1]=1.0-(fMax*6.0-1.0);
    }
    if(iMax==2){
      pfPalletR[iN-1]=1.0;
      pfPalletG[iN-1]=fMax*6.0-2.0;
      pfPalletB[iN-1]=0.0;
    }
    if(iMax==3){
      pfPalletR[iN-1]=1.0-(fMax*6.0-3.0);
      pfPalletG[iN-1]=1.0;
      pfPalletB[iN-1]=0.0;
    }
    if(iMax==4){
      pfPalletR[iN-1]=0.0;
      pfPalletG[iN-1]=1.0;
      pfPalletB[iN-1]=fMax*6.0-4.0;
    }
    if(iMax==5){
      pfPalletR[iN-1]=0.0;
      pfPalletG[iN-1]=1.0-(fMax*6.0-5.0);
      pfPalletB[iN-1]=1.0;
    }
    for(i=iMin+1,j=1;i<iMax+1;i++,j++){
      if(i==1){
	pfPalletR[j]=1.0;
	pfPalletG[j]=0.0;
	pfPalletB[j]=1.0;
      }
      if(i==2){
	pfPalletR[j]=1.0;
	pfPalletG[j]=0.0;
	pfPalletB[j]=0.0;
      }
      if(i==3){
	pfPalletR[j]=1.0;
	pfPalletG[j]=1.0;
	pfPalletB[j]=0.0;
      }
      if(i==4){
	pfPalletR[j]=0.0;
	pfPalletG[j]=1.0;
	pfPalletB[j]=0.0;
      }
      if(i==5){
	pfPalletR[j]=0.0;
	pfPalletG[j]=1.0;
	pfPalletB[j]=1.0;
      }
    }
  }

  cpgctab(pfPalletL, pfPalletR, pfPalletG, pfPalletB, iN, fContra, fBright);
  if(pfBuffer) free(pfBuffer);
}

#define BORDER_WRAP_NONE  0
#define BORDER_WRAP_RIGHT 1
#define BORDER_WRAP_LEFT  2
#define SMALL_FLOAT 1.0e-6

void MakeImagePallet(int type,float fS,float fPr,float fPo){
  float fContra=1.0;
  float fBright=0.5;
  float fMin,fMax;
  int iMin=0,iMax=5,iN=0;
  float fP;
  float *pfBuffer=NULL;
  int i,j;
  int iMinD=0,iMaxD=0;
  int iBorder=BORDER_WRAP_NONE;
  float *pfPalletL;
  float *pfPalletR;
  float *pfPalletG;
  float *pfPalletB;
  float fDum;
  int iNhex=7;
  float hexl[]={ 0.0, 1.0/6.0, 1.0/3.0, 0.5, 2.0/3.0, 5.0/6.0, 1.0};
  float hexr[]={ 1.0, 0.0,     0.0,     0.0, 1.0,     1.0,     1.0};
  float hexg[]={ 1.0, 1.0,     1.0,     0.0, 0.0,     0.0,     1.0};
  float hexb[]={ 0.0, 0.0,     1.0,     1.0, 1.0,     0.0,     0.0};
  float fR=0.0,fG=0.0,fB=0.0;

  fP=fS+fPr+fPo;
  fMin=0.5-fPo/fP;
  fMax=0.5+fPr/fP;
  if(fMin<0.0){ 
    iBorder|=BORDER_WRAP_LEFT;
    fMin+=1.0;
  }
  if(fMax>1.0){ 
    iBorder|=BORDER_WRAP_RIGHT;
    fMax-=1.0;
  }
  if(iBorder==3){
    printf("GRA Incorrect border\n");
    return;
  }
  if(fMin>fMax){
    fDum=fMin;
    fMin=fMax;
    fMax=fDum;
  }
  iN=iNhex+4;
  if(!(pfBuffer=(float*)malloc(4*iN*sizeof(float)))){
    printf("GRA Cannot allocate for pfBuffer. No pallet change\n");
    return;
  }
  pfPalletL=pfBuffer;
  pfPalletR=pfBuffer+iN;
  pfPalletG=pfBuffer+2*iN;
  pfPalletB=pfBuffer+3*iN;

  for(i=0,j=0;i<iNhex-1;i++,j++){
    pfPalletL[j]=hexl[i];
    pfPalletR[j]=hexr[i];
    pfPalletG[j]=hexg[i];
    pfPalletB[j]=hexb[i];
    
    if(fMin>= hexl[i] && fMin< hexl[i+1] && !iMinD){
      iMinD=1;
      j++;
      pfPalletL[j]=fMin;
      pfPalletL[j+1]=fMin+SMALL_FLOAT;
      pfPalletR[j+1]= pfPalletR[j]=fR;
      pfPalletG[j+1]= pfPalletG[j]=fG;
      pfPalletB[j+1]= pfPalletB[j]=fB;

      if(iBorder) iMin=j+1;
      else iMin=j;
      fDum=(iNhex-1)*(fMin-hexl[i]);
      pfPalletR[iMin]=hexr[i]+(hexr[i+1]-hexr[i])*fDum;
      pfPalletG[iMin]=hexg[i]+(hexg[i+1]-hexg[i])*fDum;
      pfPalletB[iMin]=hexb[i]+(hexb[i+1]-hexb[i])*fDum;
      j++;
    }
    if(fMax>= hexl[i] && fMax< hexl[i+1] && !iMaxD){
      iMaxD=1;
      j++;
      pfPalletL[j]=fMax-SMALL_FLOAT;
      pfPalletL[j+1]=fMax;
      pfPalletR[j+1]= pfPalletR[j]=fR;
      pfPalletG[j+1]= pfPalletG[j]=fG;
      pfPalletB[j+1]= pfPalletB[j]=fB;
      if(iBorder) iMax=j;
      else iMax=j+1;
      fDum=(iNhex-1)*(fMax-hexl[i]);
      pfPalletR[iMax]=hexr[i]+(hexr[i+1]-hexr[i])*fDum;
      pfPalletG[iMax]=hexg[i]+(hexg[i+1]-hexg[i])*fDum;
      pfPalletB[iMax]=hexb[i]+(hexb[i+1]-hexb[i])*fDum;
      j++;
    }
  }
  if(j==iN-1){
    pfPalletL[j]=hexl[iNhex-1];
    pfPalletR[j]=hexr[iNhex-1];
    pfPalletG[j]=hexg[iNhex-1];
    pfPalletB[j]=hexb[iNhex-1];
    if(iBorder){
      for(i=0;i<iMin;i++){
	pfPalletR[i]=fR;
	pfPalletG[i]=fG;
	pfPalletB[i]=fB;
      }
      for(i=iMax-1;i<iN;i++){
	pfPalletR[i]=fR;
	pfPalletG[i]=fG;
	pfPalletB[i]=fB;
      }      
    }
    else{
      for(i=iMin+1;i<iMax;i++){
	pfPalletR[i]=fR;
	pfPalletG[i]=fG;
	pfPalletB[i]=fB;
      }
    }
    cpgctab(pfPalletL, pfPalletR, pfPalletG, pfPalletB, iN, fContra, fBright);
  }
  else{
    printf("GRA Incorrect pallet %i(%i)\n",j,iN-1);
  }
  if(pfBuffer) free(pfBuffer);
}
