#include "iman1.h"

#define GET_STRING_MODE_ECHO_OFF 1

void MarkPoints(float fX,float fY,int iN,int iColor);

int GreenManager(int iVersion,FILE *fp,char *pc){
  int sid,pid;
  int i,j,k;
  int x,y,x0,y0,m;
  float dumf,xx,yy;
  char str[100],ch;
  FILEHEADER fh;
  SOFT_CHUNK SOFTChunk;
  unsigned long ulDataType;
  unsigned char *puc,*pucBuffer=NULL;
  unsigned short *pus,*pusBuffer=NULL;
  unsigned long ul,*pul,*pulBuffer=NULL;
  float *pf,*pfBuffer=NULL;
  int iReturn=0;
  int iCropROI=0;
  int iIgnoreFirstLine=0;
  //  int iOX=160,iOY=40;
  //  int iDX=480,iDY=920;
  int iOX=10,iOY=10;
  int iDX=10,iDY=10;
  int iXYShift=0;
  unsigned long ulXYDCrop;
  float fClipCoeffLo=0.1,fClipCoeffHi=-0.2;//-0.1;
  int iMarkPointsN=0;
  float *pfMarkPointsX=NULL;
  float *pfMarkPointsY=NULL;
  int iCount;

  //  dradius=DCRADIUS;
  radius=(int)floor(dradius);
  if(radius!=0){
    if(!(aindex_x=(int *)calloc((1+2*radius)*(1+2*radius),sizeof(int))) || 
       !(aindex_y=(int *)calloc((1+2*radius)*(1+2*radius),sizeof(int)))){
      printf("GRE Cannot allocate for aindex\n");
      return(6);
    }
    aindex_n=0;
    for(j=-radius;j<=radius;j++){
      for(k=-radius;k<=radius;k++){
	if(hypot((double)j,(double)k)<=dradius){ 
	  aindex_x[aindex_n]=k;
	  aindex_y[aindex_n++]=j;
	}
      }
    }
    printf("GRE aindex=%i\n",aindex_n);
  }
  rewind(fp);
  if(iVersion==VERSION_CHUNK){
    if((ul=FindChunkInFile(fp,"SOFT"))==~0UL){
      printf("GRE Cannot find SOFT chunk in file %s\n",pc);
      return(3);
    }
    if(ul+8!=sizeof(SOFT_CHUNK)){
      printf("GRE Not standard(%iu) SOFT chunk(%lu)\n",sizeof(SOFT_CHUNK),ul+8);
    }
    if(!ReadChunk((void *)&SOFTChunk,ul+8,fp)){
      printf("GRE Cannot read SOFT chunk from file %s\n",pc); 
      return(4); 
    }

    Xdim=SOFTChunk.XSize;
    Ydim=SOFTChunk.YSize;
    ulDataType=SOFTChunk.DataType;

    if(FindChunkInFileAndOffset(fp,"DATA",FALSE)==~0UL){
      printf("GRE Cannot find DATA chunk in file %s\n",pc);
      return(5);
    }
  }
  else{
    if(fread((void*)&fh,sizeof(FILEHEADER),1,fp) != 1){
      printf("GRE Cannot read fileheader from file %s\n",pc);
      return(3);
    }
    Xdim=fh.sumXsize;
    Ydim=fh.sumYsize;
    ulDataType=fh.sumLDataType;
  }

  XYdim=Xdim*Ydim;
  printf("GRE Dimensions: X %i Y %i Type %lu\n",Xdim,Ydim,ulDataType);

  if(!(pfBuffer=(float*)malloc(XYdim*sizeof(float)))){
    printf("GRE Cannot allocate for pfBuffer\n");
    return(1);
  }

  switch(ulDataType){
  case DATATYPE_UCHAR:
    if(!(pucBuffer=(unsigned char*)malloc(XYdim*sizeof(unsigned char)))){
      printf("GRE Cannot allocate for pucBuffer\n");
      iReturn=7;
      goto bailout;
    }
    if((i=fread((void*)pucBuffer,XYdim*sizeof(unsigned char),1,fp)) != 1){
      printf("GRE Cannot read frame from green file (%i)\n",i);
      iReturn=8;
      goto bailout;
    }
    for(ul=0,puc=pucBuffer,pf=pfBuffer;ul<XYdim;ul++) *(pf++)=(float)*(puc++);
    free(pucBuffer);
    break;
  case DATATYPE_USHORT:
    if(!(pusBuffer=(unsigned short*)malloc(XYdim*sizeof(unsigned short)))){
      printf("GRE Cannot allocate for pusBuffer\n");
      iReturn=9;
      goto bailout;
    }
    if((i=fread((void*)pusBuffer,XYdim*sizeof(unsigned short),1,fp)) != 1){
      printf("GRE Cannot read frame from green file (%i)\n",i);
      iReturn=10;
      goto bailout;
    }
    for(ul=0,pus=pusBuffer,pf=pfBuffer;ul<XYdim;ul++) *(pf++)=(float)*(pus++);
    free(pusBuffer);
    break;
  case DATATYPE_ULONG:
    if(!(pulBuffer=(unsigned long*)malloc(XYdim*sizeof(unsigned long)))){
      printf("GRE Cannot allocate for pulBuffer\n");
      iReturn=11;
      goto bailout;
    }
    if((i=fread((void*)pulBuffer,XYdim*sizeof(unsigned long),1,fp)) != 1){
      printf("GRE Cannot read frame from green file (%i)\n",i);
      iReturn=12;
      goto bailout;
    }
    for(ul=0,pul=pulBuffer,pf=pfBuffer;ul<XYdim;ul++) *(pf++)=(float)*(pul++);
    free(pulBuffer);
    break;
  case DATATYPE_FLOAT:
    if((i=fread((void*)pfBuffer,XYdim*sizeof(float),1,fp)) != 1){
      printf("GRE Cannot read frame from green file (%i)\n",i);
      iReturn=13;
      goto bailout;
    }
    break;
  default:
    printf("GRE Unknown data type %lu",ulDataType);
    iReturn=6;
    goto bailout;
  }

  fclose(fp);

  if(iCropROI){
    if(iOX+iDX>Xdim || iOY+iDY>Ydim){
      printf("GRE ROI out of the range\n");
      iReturn=7;
      goto bailout;
    }
    ulXYDCrop=(unsigned long)iDX*(unsigned long)iDY;
    iXYShift=iOX+iOY*Xdim;
    if(iXYShift){
      for(j=0;j<ulXYDCrop;j++){
	i=iXYShift+j%iDX+(j/iDX)*Xdim;
	pfBuffer[j]=pfBuffer[i];
      }
    }
    Xdim=iDX;
    Ydim=iDY;
    XYdim=Xdim*Ydim;
  }

  fg=0.0;
  bg=(float)(1<<16);
  if(iIgnoreFirstLine) i=Xdim;
  else i=0;
  for(;i<XYdim;i++){
    dumf=pfBuffer[i];
    if(fg<dumf){ 
      fg=dumf;
      ifor=i;
    }
    if(bg>dumf){ 
      bg=dumf;
      ibac=i;
    }
  }
  printf("GRE Not averaged green: Max=%f@%i Min=%f@%i\n",fg,ifor,bg,ibac);

  if(radius){
    printf("%i\n",aindex_n);
    fg=0.0;
    bg=(float)(1<<16);
    if(iIgnoreFirstLine) i=Xdim;
    else i=0;
    for(;i<XYdim;i++){
      dumf=0;
      m=0;
      x0=i%Xdim;
      y0=i/Xdim;
      for(j=0;j<aindex_n;j++){
	x=x0+aindex_x[j];
	y=y0+aindex_y[j];
	if(x>=0 && x<Xdim && y>=0 && y<Ydim){ 
	  dumf+=pfBuffer[x+y*Xdim];
	  m++;
	}
      }
      dumf/=(float)m;
      if(fg<dumf){ 
	fg=dumf;
	ifor=i;
      }
      if(bg>dumf){ 
	bg=dumf;
	ibac=i;
      }
    }
    printf("GRE Averaged green: Max=%f@%i Min=%f@%i\n",fg,ifor,bg,ibac);
  }

  dumf=fg-(fg-bg)*fClipCoeffHi;
  bg=bg+(fg-bg)*fClipCoeffLo;
  fg=dumf;

  pminx=-0.5;
  pmaxy=-0.5;
  pmaxx=(float)Xdim-0.5;
  pminy=(float)Ydim-0.5;
  //  if(iDrawAxesAndLabelOnGreen) 
  setenv("PGPLOT_ENVOPT","I",1);
  pid=cpgopen(device);
  cpgask(0);
  if(dradiusbig>0.0) cpgpap(dradiusbig,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
  else cpgpap(GREEN_IMAGE_SIZE,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
  if(iDrawAxesAndLabelOnGreen || 1){
    cpgsch(0.8);
    cpgenv(pminx,pmaxx,pminy,pmaxy,1,0);
    sprintf(str,"Green %s",pc);
    cpglab("","",str);
    cpgsch(1.0);
  }
  else{
    cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
  }
  cpggray(pfBuffer,Xdim,Ydim,1,Xdim,1,Ydim,fg*fContrast,bg,trans);

  if(do_intervals&2){
    if(finalframe>=0) iCount=finalframe;
    else iCount=1;

    while(1){
      cpgband(7,1,xx,yy,&xx,&yy,&ch); 
      if(do_verbose) printf("INT X %f Y %f C (%i)%c\n",xx,yy,ch,ch);

      if(ch=='q'){ break;}
      
      if(ch=='Q'){ exit(0);}

      if(ch=='m' || ch=='M'){
	if(ch=='M'){
	  if(!GetFloatPoint(&xx,&yy,7,3)) continue;
	}
	  
	if(!(pfMarkPointsX=(float*)realloc(pfMarkPointsX,(iMarkPointsN+1)*sizeof(float))) || 
	   !(pfMarkPointsY=(float*)realloc(pfMarkPointsY,(iMarkPointsN+1)*sizeof(float)))){
	  printf("INI Cannot realloc for pfMarkPointsX/pfMarkPointsY\n");
	  continue;
	}
	pfMarkPointsX[iMarkPointsN]=xx;
	pfMarkPointsY[iMarkPointsN]=yy;
	iMarkPointsN++;
	MarkPoints(xx,yy,iMarkPointsN,1);
      }

      x=(int)rint(xx);
      y=(int)rint(yy);
      if(initframe>0){
	x/=initframe;
	y/=initframe;
      }
      if(sacrifice_boundaries) printf("GRE X %i Y %i\n",x,y);
      else{
	GetString(str,GET_STRING_MODE_ECHO_OFF);
	printf("%i\t%i\t%s\t%i\n",x,y,str,iCount++);
      }
    }
  }

  printf("Save green image?(y/n): ");
  fgets(str,4,stdin);
  if(*str == 'y'){
    sprintf(str,"%s%s",pc,".ps/VCPS");
    printf("%s\n",str);
    sid=cpgopen("?");
    cpgpap(fSizeOfGreenSave,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));

    if(iDrawAxesAndLabelOnGreen){
      cpgsch(0.8);
      cpgenv(pminx,pmaxx,pminy,pmaxy,1,0);
      //      cpgeras();
      sprintf(str,"Green %s",pc);
      cpglab("","",str);
      cpgsch(1.0);
    }
    else{
      cpgsvp(0.0,1.0,0.0,1.0);
      cpgswin(pminx,pmaxx,pminy,pmaxy);
      //      cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
    }

    cpggray(pfBuffer,Xdim,Ydim,1,Xdim,1,Ydim,bg,fg*fContrast,trans);

    for(i=0;i<iMarkPointsN;i++) MarkPoints(pfMarkPointsX[i],pfMarkPointsY[i],i+1,0);

    cpgclos();
  }   
  cpgslct(pid);
  cpgclos();
  unsetenv("PGPLOT_ENVOPT");
 bailout:
  free(pfBuffer);
  if(pfMarkPointsX) free(pfMarkPointsX);
  if(pfMarkPointsY) free(pfMarkPointsY);
  return(iReturn);
}

char* GetString(char *str,int iMode){
  int i;
  char ch;
  float dumfx,dumfy;
  i=0;
  while(1){
    cpgband(0,0,0.0,0.0,&dumfx,&dumfy,&ch); 
    if(ch!=13 && ch!=32){
      if(ch==8){
	if(i){ 
	  i--;
	  str[i]='\0';
	  if(!(iMode&GET_STRING_MODE_ECHO_OFF)) printf("%c",8);
	}
      }
      else{
	str[i++]=ch;
	if(!(iMode&GET_STRING_MODE_ECHO_OFF)) printf("%c",ch);
      }
      fflush(stdout);
    }	
    else{ 
      if(ch=='q') return(NULL);
      str[i]='\0';
      if(!(iMode&GET_STRING_MODE_ECHO_OFF)){
	printf("%c",ch);
	fflush(stdout);
      }
      break;
    }
  }
  return(str);
}

#define MAX_GET_TEXT_STRING_LENGTH 128

int GetFloatPoint(float *pfX,float *pfY,int iMode,int iColor){
  char ch;
  char str[MAX_GET_TEXT_STRING_LENGTH];
  int iOldColor;
  int iReturn=1;
  cpgqci(&iOldColor);
  cpgsci(iColor);
  cpgband(iMode,0,*pfX,*pfY,pfX,pfY,&ch);
  if(ch=='q') iReturn=0;
  if(ch=='m'){
    if(GetString(str,0)){ 
      *pfX=atof(str);
      if(GetString(str,0)){
	*pfY=atof(str);
	iReturn=2;
      }
      else iReturn=0;
    }
    else iReturn=0;
  }
  cpgsci(iOldColor);
  return(iReturn);
}

int GetIntPoint(int *piX,int *piY,int iDX,int iDY,int iMode,int iColor){
  float fX,fY;
  int iReturn=0;
  if((iReturn=GetFloatPoint(&fX,&fY,iMode,iColor))){
    *piX=(int)rint(fX);
    if(*piX<0) *piX=0;
    if(*piX>iDX-1) *piX=iDX-1;
    *piY=(int)rint(fY);
    if(*piY<0) *piY=0;
    if(*piY>iDY-1) *piY=iDY-1;
  }
  return(iReturn);
}

#define TEXT_SHIFT_X_PCT 0.01
#define TEXT_SHIFT_Y_PCT 0.01

void MarkPoints(float fX,float fY,int iN,int iColor){
  int iOldColor;
  int iOldLineWidth;
  float fOldCharcterHeight;
  float fX1,fY1,fX2,fY2;
  float fXTextShift,fYTextShift;
  char str[128];
  int iMarkPointsEnumerate=1;
  int iMarkPointsSymbol=-3; //1-square, 3-*, 4-circle, 9-circle+dot
  //  int iMarkPointsColor=1;
  float fMarkPointsRadius=6.0; // if iMarkPointsSymbol=-3

  cpgqlw(&iOldLineWidth);
  cpgslw(2);
  cpgqci(&iOldColor);
  cpgsci(iColor);
  
  if(iMarkPointsSymbol!=-3){
    cpgpt1(fX*trans[1]+trans[0],fY*trans[5]+trans[3],iMarkPointsSymbol);
  }
  else{
    cpgsfs(2);
    cpgcirc(fX*trans[1]+trans[0],fY*trans[5]+trans[3],fMarkPointsRadius);
    cpgsfs(1);
  }

  if(iMarkPointsEnumerate){
    cpgqwin(&fX1,&fX2,&fY1,&fY2);
    fXTextShift=(fX2-fX1)*TEXT_SHIFT_X_PCT;
    fYTextShift=(fY2-fY1)*TEXT_SHIFT_Y_PCT;
    cpgqch(&fOldCharcterHeight);
    cpgsch(1.0);
    sprintf(str,"%i",iN++);
    cpgtext(fXTextShift+fX*trans[1]+trans[0],-fYTextShift+fY*trans[5]+trans[3],str);
    cpgsch(fOldCharcterHeight);
  }
  
  cpgsci(iOldColor);
  cpgslw(iOldLineWidth);
}
