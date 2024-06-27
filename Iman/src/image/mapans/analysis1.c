/* Map analysis                 */
/* Written by V.Kalatsky        */
/* 06Feb2001                    */
/* Modified 10May2001           */

#include "mapans1.h"

void Stencil(int xd,int yd,int n,POINTD *pDom,unsigned long *sn,unsigned long **s){
  static unsigned long *pul=(unsigned long *)0;
  unsigned long k;
  int j;
  unsigned long i,xy;
  double x,y;
  double angle;
  float scalex,scaley;
  double dVP,dSP;

  if(n==0){
    xy=(unsigned long)xd*(unsigned long)yd;
    pul=(unsigned long *)realloc(*s,xy*sizeof(unsigned long));
    *sn=xy;
    for(i=0;i<stenciln;i++) pul[i]=i;
    *s=pul;
    if(do_rescaling) memcpy(transs,trans,6*sizeof(float));
    return;
  }

  if(iAdditiveCropChanges){
    if(stenciln && stencil){
      xy=stenciln;
      pul=(unsigned long *)malloc(xy*sizeof(unsigned long));
    }
    else{
      printf("ANA EMPTY STENCIL\n");
      return;
    }
  }
  else{
    xy=(unsigned long)xd*(unsigned long)yd;
    pul=(unsigned long *)realloc(*s,xy*sizeof(unsigned long));
  }

  smaxx=smaxy=0;
  sminx=xd-1;
  sminy=yd-1;
  for(k=i=0;i<xy;i++){
    if(iAdditiveCropChanges){
      x=(double)(stencil[i]%xd);
      y=(double)(stencil[i]/xd);
   }
    else{
      x=(double)(i%xd);
      y=(double)(i/xd);
    }
    for(j=0,angle=0.0;j<n;j++){
      dVP=VECTPROD(pDom[j].x-x,pDom[j].y-y,pDom[j+1].x-x,pDom[j+1].y-y);
      dSP=DOTPROD(pDom[j].x-x,pDom[j].y-y,pDom[j+1].x-x,pDom[j+1].y-y);
      // if(dVP==0.0 && dSP==0.0) break;
      angle+=atan2(dVP,dSP);
    }
    //if(fabs(angle)>0.1) printf("ANA %li Angle %0.1f\n",i,angle);
    //    if(fabs(angle)>M_PI || j!=n){
    if(fabs(angle)>M_PI){
      if(iAdditiveCropChanges) pul[k]=stencil[i];
      else pul[k]=i;
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
  printf("ANA %lu pixels\n",k);
  if(do_rescaling){
    if(smaxx!=sminx) scalex=(float)(xd-1)/(float)(smaxx-sminx);
    else scalex=1.0;
    if(smaxy!=sminy) scaley=(float)(yd-1)/(float)(smaxy-sminy);
    else scaley=1.0;
    transs[1]=transs[5]= scalex<scaley ? scalex : scaley;
    transs[0]=0.0-transs[1]*sminx;
    transs[3]=0.0-transs[1]*sminy;
  }

  if(iFindMarkedInsideSelected){
    for(i=0;i<iMarkPointsN;i++){
      x=(double)(pMarkPoints[i].x);
      y=(double)(pMarkPoints[i].y);
      for(j=0,angle=0.0;j<n;j++){
	dVP=VECTPROD(pDom[j].x-x,pDom[j].y-y,pDom[j+1].x-x,pDom[j+1].y-y);
	dSP=DOTPROD(pDom[j].x-x,pDom[j].y-y,pDom[j+1].x-x,pDom[j+1].y-y);
	// if(dVP==0.0 && dSP==0.0) break;
	angle+=atan2(dVP,dSP);
      }
      if(fabs(angle)>M_PI){
	printf("%lu\n",i+1);
      }
    }
  }
  if(iSaveCropBorderAsLineSegments){
      if(!(pLineSegments=(PAIR*)realloc(pLineSegments,(iNLineSegments+n)*sizeof(PAIR)))){
	printf("DIS Cannot reallocate for pLineSegments\n");
	return;
      }
      for(j=0;j<n;j++){
	pLineSegments[iNLineSegments].n=iNLineSegments;
	pLineSegments[iNLineSegments].x1=pDom[j].x;
	pLineSegments[iNLineSegments].y1=pDom[j].y;
	pLineSegments[iNLineSegments].x2=pDom[j+1].x;
	pLineSegments[iNLineSegments].y2=pDom[j+1].y;
	iNLineSegments++; 
      }  
  }
}

int IsPointInsideDomain(double dX,int dY,int n,POINTD *pDom,RECT *pRect){
  double dAngle;
  int i;
  POINTD *pD;
  if(n==0 || pDom==NULL){
    if(!pRect){
      printf("ANA Domain is not specified\n");
      return(0);
    }
    if(dX>pRect->x1 && dX<pRect->x2 && dY>pRect->y1 && dY<pRect->y2) return(1);
    else return(0);
  }

  for(i=0,dAngle=0.0,pD=pDom;i<n;i++,pD++){
    dAngle+=atan2(VECTPROD(pD->x-dX,pD->y-dY,(pD+1)->x-dX,(pD+1)->y-dY),
		 DOTPROD(pD->x-dX,pD->y-dY,(pD+1)->x-dX,(pD+1)->y-dY));
  }
  if(fabs(dAngle)>M_PI) return(1);
  else return(0);
}

void StencilMap(int xd,int yd,float *pfX,float *pfY,int n,double *dx,double *dy,unsigned long *sn,unsigned long **s){
  unsigned long *pul=(unsigned long *)0;
  unsigned long k,l;
  int j;
  unsigned long i,xy;
  double x,y,dDum;
  double angle;
  double dVP,dSP,dThreshLo,dThreshHi;
  int iThresold=0;

  xy=(unsigned long)xd*(unsigned long)yd;
  pul=(unsigned long *)calloc(xy,sizeof(unsigned long));

  if(n==0){
    *sn=xy;
    for(i=0;i<xy;i++) pul[i]=i;
    if(*s) free(*s);
    *s=pul;
    return;
  }

  if(*sn==0 || *s==NULL){
    printf("WARNING: Empty stencil in StencilMap");
    free(pul);
    return;
  }

  if(n==1){
    // Threshold *dx - Lo (1) *dy - Hi (2) both (3)
    /*
    dThresh=*dx*(*dx);
    for(k=i=0;i<*sn;i++){
      l=(*s)[i];
      x=(double)(pfX[l]);
      y=(double)(pfY[l]);
      if(x*x+y*y>dThresh){
	pul[k++]=l;
      }
    }
    */
    iThresold=0;
    if(dx && *dx>=0.0){
      dThreshLo=*dx*(*dx);
      iThresold|=1;
    }
    if(dy && *dy>=0.0){
      dThreshHi=*dy*(*dy);
      iThresold|=2;
    }
    switch(iThresold){
    case 0:
      for(k=i=0;i<*sn;i++){
	l=(*s)[i];
	pul[k++]=l;
      }
      printf("ANA Nothing to do in StencilMap\n");
      break;
    case 1:
      for(k=i=0;i<*sn;i++){
	l=(*s)[i];
	x=(double)(pfX[l]);
	y=(double)(pfY[l]);
	if(x*x+y*y>=dThreshLo){
	  pul[k++]=l;
	}
      }
      break;
    case 2:
      for(k=i=0;i<*sn;i++){
	l=(*s)[i];
	x=(double)(pfX[l]);
	y=(double)(pfY[l]);
	if(x*x+y*y<=dThreshHi){
	  pul[k++]=l;
	}
      }
      break;
    case 3:
      for(k=i=0;i<*sn;i++){
	l=(*s)[i];
	x=(double)(pfX[l]);
	y=(double)(pfY[l]);
	dDum=x*x+y*y;
	if(dDum<=dThreshHi && dDum>=dThreshLo){
	  pul[k++]=l;
	}
      }
      break;
    }
  }
  else{
    // Domain
    for(k=i=0;i<*sn;i++){
      l=(*s)[i];
      x=(double)(pfX[l]);
      y=(double)(pfY[l]);
      for(j=0,angle=0.0;j<n;j++){
	dVP=VECTPROD(dx[j]-x,dy[j]-y,dx[j+1]-x,dy[j+1]-y);
	dSP=DOTPROD(dx[j]-x,dy[j]-y,dx[j+1]-x,dy[j+1]-y);
	// if(dVP==0.0 && dSP==0.0) break;
	angle+=atan2(dVP,dSP);
      }
      //    if(fabs(angle)>M_PI || j!=n){
      if(fabs(angle)>M_PI){
	pul[k++]=l;
      }
    }
  }
  *sn=k;
  pul=(unsigned long *)realloc(pul,k*sizeof(unsigned long));
  if(*s) free(*s);
  *s=pul;
}



int LoadStencilBitmap(int xd,int yd,unsigned long *sn,unsigned long **s,int iMode){
  char strFileName[256];
  FILE *pF;
  int i,iDum;
  int iType;
  int iBitsPerRecord;
  unsigned long ulBytes;
  char *pcRecords;
  int iError=0;
  unsigned long *pul=(unsigned long *)0;
  unsigned long j,k,l;
  unsigned long xy;

  xy=(unsigned long)xd*(unsigned long)yd;
  pul=(unsigned long *)calloc(xy,sizeof(unsigned long));

  fflush(stdin);
  printf("ANA Enter bitmap file name: ");
  fgets(strFileName,255,stdin);
  if(!strlen(strFileName)) printf("ANA Empty file name\n");
  if(strFileName[strlen(strFileName)-1]=='\n') strFileName[strlen(strFileName)-1]='\0';

  if((pF=fopen(strFileName,"r"))==NULL){
    printf("ANA Cannot open file <%s>\n",strFileName);
    return(1);
  }

  if(fread(&iType,sizeof(int),1,pF) != 1){
    printf("ANA Cannot read type from file %s\n",strFileName);
    iError=2;
  }

  if(fread(&i,sizeof(int),1,pF) != 1){
    printf("ANA Cannot read Xdim from file %s\n",strFileName);
    iError=3;
  }
  if(i!=xd){
    printf("ANA WARNING X dimension missmatch %i %i\n",xd,i);
    iError=4;
  }

  if(fread(&i,sizeof(int),1,pF) != 1){
    printf("ANA Cannot read Ydim from file %s\n",strFileName);
    iError=5;
  }
  if(i!=yd){
    printf("ANA WARNING Y dimension missmatch %i %i\n",yd,i);
    iError=6;
  }
  
  if(fread(&iBitsPerRecord,sizeof(int),1,pF) != 1){
    printf("ANA Cannot read iBitsPerRecord from file %s\n",strFileName);
    iError=7;
  }
  if(iBitsPerRecord!=1){
    printf("ANA Unsupported file format in %s\n",strFileName);
    iError=8;
  }

  if(iError){
    fclose(pF);
    return(iError);
  }

  ulBytes=xy/8;
  if(ulBytes%8) ulBytes++;

  if(!(pcRecords=malloc(ulBytes*sizeof(char)))){
    printf("ANA Cannot allocate for pcRecords\n");
    fclose(pF);
    return(9);
  }
  if(fread(pcRecords,ulBytes*sizeof(char),1,pF) != 1){
    printf("ANA Cannot read bitmap from file %s\n",strFileName);
    fclose(pF);
    return(10);
  }

  fclose(pF);

  // l - current map pixel
  // k - current pixels in the stencil
  // j - current bitmap byte
  // i - current bitmap byte bit mask
  switch(iMode){
  case MODE_LOAD_BITMAP_REPLACE:
    for(j=k=l=0,i=0;l<xy;l++){
      i=1<<(l%8);
      iDum=(int)pcRecords[l/8];
      if(iDum&i) pul[k++]=l;
    }
    break;
  default:
    printf("ANA Unknown mode\n");
    free(pcRecords);
    return(11);
    break;
  }
  free(pcRecords);
  *sn=k;
  pul=(unsigned long *)realloc(pul,k*sizeof(unsigned long));
  if(*s) free(*s);
  *s=pul;
  return(0);
}

int SaveStencilBitmap(int xd,int yd,unsigned long sn,unsigned long *s){
  char strFileName[256];
  FILE *pF;
  int i;
  int iType=0;
  unsigned long ulBytes;
  char *pcRecords,*pc;
  int iError=0;
  unsigned long k,l;
  unsigned long xy;
  static int iBitmapSaveCount=0;

  xy=(unsigned long)xd*(unsigned long)yd;
  ulBytes=xy/8;
  if(ulBytes%8) ulBytes++;

  if(!(pcRecords=calloc(ulBytes,sizeof(char)))){
    printf("ANA Cannot allocate for pcRecords\n");
    return(-1);
  }

  snprintf(strFileName,255,"%s_bitmap%i",filenames[0],iBitmapSaveCount++);
  printf("ANA Saving bitmap in %s\n",strFileName);

  if((pF=fopen(strFileName,"w"))==NULL){
    printf("ANA Cannot open file <%s>\n",strFileName);
    return(1);
  }

  if(fwrite(&iType,sizeof(int),1,pF) != 1){
    printf("ANA Cannot write type %i to file %s\n",i,strFileName);
    iError=2;
  }

  if(fwrite(&xd,sizeof(int),1,pF) != 1){
    printf("ANA Cannot write Xdim to file %s\n",strFileName);
    iError=3;
  }

  if(fwrite(&yd,sizeof(int),1,pF) != 1){
    printf("ANA Cannot write Ydim to file %s\n",strFileName);
    iError=5;
  }

  i=1;
  if(fwrite(&i,sizeof(int),1,pF) != 1){
    printf("ANA Cannot write iBitsPerRecord to file %s\n",strFileName);
    iError=7;
  }

  if(iError){
    free(pcRecords);
    fclose(pF);
    return(iError);
  }

  for(l=0;l<sn;l++){
    k=s[l];
    pc=&pcRecords[k/8];
    *pc|=(unsigned char)(1<<(k%8));
  }

  if(fwrite(pcRecords,ulBytes*sizeof(char),1,pF) != 1){
    printf("ANA Cannot write bitmap to file %s\n",strFileName);
    iError=11;
  }

  fclose(pF);
  free(pcRecords);
  return(iError);
}

int AddPixelToStencil(int x,int y,int xd,int yd,unsigned long *sn,unsigned long **s,int iMode){
  unsigned long i,ulPixel;
  if(x<0 || x>=xd) return(1);
  if(y<0 || y>=yd) return(2);
  ulPixel=(unsigned long)x+(unsigned long)y*(unsigned long)xd;
  for(i=0;i<*sn;i++) if(ulPixel==(*s)[i]) return(3);
  if(!(*s=(unsigned long *)realloc(*s,(*sn+1)*sizeof(unsigned long)))){
    printf("ANA Cannot realloc for stencil\n");
    return(-1);
  }
  else{
    (*s)[(*sn)++]=ulPixel;
    return(0);
  }
}

int AnalyzePW(int xd,int yd,double *mx,double *my){
  int i,j,k,m;
  int ix,iy,xpos,ypos,xrbits,yrbits,xbits,ybits,xangle,yangle;
  int p1,p2;
  double xmapx1,xmapy1,xmapx2,xmapy2;
  double ymapx1,ymapy1,ymapx2,ymapy2;
  double dumd,vp12,vp34,x,y;
  int c;
  int loop[4],loopx[4],loopy[4];
  FILE *pF=NULL;
  char strFileName[256];

  loop[0]=0;
  loop[1]=1;
  loop[2]=1+xd;
  loop[3]=xd;

  loopx[0]=0;
  loopx[1]=1;
  loopx[2]=1;
  loopx[3]=0;
  loopy[0]=0;
  loopy[1]=0;
  loopy[2]=1;
  loopy[3]=1;

  free(apinwheels);
  apinwheels_n=0;
  for(j=0;j<yd-1;j++){
    for(i=0;i<xd-1;i++){
      k=i+j*xd;

      for(xbits=ybits=ix=iy=m=0;m<4;m++){
	if(mx[k+loop[m]]>0.0){ 
	  ix++;
	  xbits|=1; 
	}
	xbits<<=1;
	if(my[k+loop[m]]>0.0){
	  iy++;
	  ybits|=1; 
	}
	ybits<<=1;
      }
      xbits>>=1;
      ybits>>=1;      

//      printf("%i %i %i %i (%i %i)\n",i,j,xbits,ybits,ix,iy);
      if((ix%4) && (iy%4)){

	if(!Rotator(xbits,&xrbits,&xpos,&xangle)) continue;
	p1=xpos;
	p2=(xpos+3)%4;
	dumd=mx[k+loop[p1]]/(mx[k+loop[p1]]-mx[k+loop[p2]]);
	xmapx1=loopx[p1]+dumd*(loopx[p2]-loopx[p1]);
	xmapy1=loopy[p1]+dumd*(loopy[p2]-loopy[p1]);
	p1=(xpos+(2-xrbits/3))%4;
	p2=(p1+1)%4;
	dumd=mx[k+loop[p1]]/(mx[k+loop[p1]]-mx[k+loop[p2]]);
	xmapx2=loopx[p1]+dumd*(loopx[p2]-loopx[p1]);
	xmapy2=loopy[p1]+dumd*(loopy[p2]-loopy[p1]);

	if(!Rotator(ybits,&yrbits,&ypos,&yangle)) continue;
	p1=ypos;
	p2=(ypos+3)%4;
	dumd=my[k+loop[p1]]/(my[k+loop[p1]]-my[k+loop[p2]]);
	ymapx1=loopx[p1]+dumd*(loopx[p2]-loopx[p1]);
	ymapy1=loopy[p1]+dumd*(loopy[p2]-loopy[p1]);
	p1=(ypos+(2-yrbits/3))%4;
	p2=(p1+1)%4;
	dumd=my[k+loop[p1]]/(my[k+loop[p1]]-my[k+loop[p2]]);
	ymapx2=loopx[p1]+dumd*(loopx[p2]-loopx[p1]);
	ymapy2=loopy[p1]+dumd*(loopy[p2]-loopy[p1]);
	
//	printf("(%f %f)(%f %f)(%f %f)(%f %f)\n",xmapx1,xmapy1,xmapx2,xmapy2,ymapx1,ymapy1,ymapx2,ymapy2);

	if(VECTPROD(xmapx2-xmapx1,xmapy2-xmapy1,ymapx1-xmapx1,ymapy1-xmapy1)*VECTPROD(xmapx2-xmapx1,xmapy2-xmapy1,ymapx2-xmapx1,ymapy2-xmapy1) <= 0.0){

	  if(!(apinwheels=(Pinwheel *)realloc(apinwheels,(++apinwheels_n)*sizeof(Pinwheel)))){
	    printf("Cannot reallocate for apinwheels\n");
	    return(1);
	  }
	  if((dumd=VECTPROD(xmapx2-xmapx1,xmapy2-xmapy1,ymapx2-ymapx1,ymapy2-ymapy1))>0.0) c=1.0;
	  else c=-1.0;

	  //	  c=Charge(xangle,yangle);
	  vp12=VECTPROD(xmapx1,xmapy1,xmapx2,xmapy2);
	  vp34=VECTPROD(ymapx1,ymapy1,ymapx2,ymapy2);
	  x=((ymapx1-ymapx2)*vp12-(xmapx1-xmapx2)*vp34)/dumd;
	  y=((ymapy1-ymapy2)*vp12-(xmapy1-xmapy2)*vp34)/dumd;
	  apinwheels[apinwheels_n-1].x=(x+(double)i);
	  apinwheels[apinwheels_n-1].y=(y+(double)j);

	  apinwheels[apinwheels_n-1].charge=(double)c;
//	  printf("PW: %2i %3i %3i x: %2i %i %i (%4.3f %4.3f)(%4.3f %4.3f) y: %2i %i %i (%4.3f %4.3f)(%4.3f %4.3f)\n",c,i,j,xbits,xrbits,xpos,xmapx1,xmapy1,xmapx2,xmapy2,ybits,yrbits,ypos,ymapx1,ymapy1,ymapx2,ymapy2);
//	  printf("%f %f %f %f %f\n",x,y,dumd,vp12,vp34);
	}
      }
    }
  }
  printf("NA=%i\n",apinwheels_n);
  for(i=0;i<apinwheels_n;i++){
    printf("%i (%f %f)\n",i,apinwheels[i].x,apinwheels[i].y);
    //    printf("%i (%f %f)\n",i,apinwheels[i].x,apinwheels[i].y);
  }
  if(iSaverPinwheels){
    sprintf(strFileName,"%s_PinWheels",filenames[0]);
    printf("ANA Saving pinwheels in %s\n",strFileName);
    if((pF=fopen(strFileName,"w"))==NULL){
      printf("ANA Cannot open file %s\n",strFileName);
    }
    else{
      if(fwrite((void*)&apinwheels_n,sizeof(int),1,pF) != 1){
	printf("INT Cannot write #Pinwheels %i to file %s\n",apinwheels_n,strFileName);
      }
      if(fwrite((void*)&xd,sizeof(int),1,pF) != 1){
	printf("INT Cannot write Xdim to file %s\n",strFileName);
      }
      if(fwrite((void*)&yd,sizeof(int),1,pF) != 1){
	printf("INT Cannot write Ydim to file %s\n",strFileName);
      }
      if(fwrite((void*)apinwheels,apinwheels_n*sizeof(Pinwheel),1,pF) != 1){
	printf("INT Cannot write apinwheels to file %s\n",strFileName);
      }
      fclose(pF);
    }
  }

  return(0);
}

int Rotator(int bits,int *rbits,int *position,int *angle){
  switch(bits){
  case 0:
    *rbits=0;
    *position=0;
    *angle=-1;
    printf("None bits\n");
    return(0);
    break;
  case 1:
    *rbits=1;
    *position=0;
    *angle=3;
    return(1);
    break;
  case 2:
    *rbits=1;
    *position=3;
    *angle=1;
    return(1);
    break;
  case 3:
    *rbits=3;
    *position=0;
    *angle=2;
    return(1);
    break;
  case 4:
    *rbits=1;
    *position=2;
    *angle=7;
    return(1);
    break;
  case 5:
    *rbits=5;
    *position=0;
    *angle=-1;
    //    printf("Bad bits: 5\n");
    return(0);
    break;
  case 6:
    *rbits=3;
    *position=3;
    *angle=0;
    return(1);
    break;
  case 7:
    *rbits=7;
    *position=0;
    *angle=1;
    return(1);
    break;
  case 8:
    *rbits=1;
    *position=1;
    *angle=5;
    return(1);
    break;
  case 9:
    *rbits=3;
    *position=1;
    *angle=4;
    return(1);
    break;
  case 10:
    *rbits=5;
    *position=1;
    *angle=-1;
    //    printf("Bad bits: 10\n");
    return(0);
    break;
  case 11:
    *rbits=7;
    *position=1;
    *angle=3;
    return(1);
    break;
  case 12:
    *rbits=3;
    *position=2;
    *angle=6;
    return(1);
    break;
  case 13:
    *rbits=7;
    *position=2;
    *angle=5;
    return(1);
    break;
  case 14:
    *rbits=7;
    *position=3;
    *angle=7;
    return(1);
    break;
  case 15:
    *rbits=0;
    *position=0;
    *angle=-1;
    printf("All bits\n");
    return(0);
    break;
  default:
    printf("Bad bits: %i\n",bits);
    return(0);
    break;    
  }
}

int Charge(int x,int y){
  int i;
  i=y-x;
  if(i<0) i+=8;
  if(!(i%4)) return(0);
  if(i/4) return(-1);
  else return(1);
}

int SmartContoursPolar(int xd,int yd,float *rho,float *phi,PAIR **sc,int *scn,float val,float threshold){
  unsigned long xy,ul,uli;
  int i,j,m,n;
  PAIR *psc;
  float to;
  unsigned long loop[4];
  int loopx[4],loopy[4];
  int positive,pattern;
  int position,rpattern,angle;
  int p1,p2;
  float dumf,dumf1,dumf2;

  xy=(unsigned long)xd*(unsigned long)yd;

  loop[0]=0;
  loop[1]=1;
  loop[2]=1+xd;
  loop[3]=xd;

  loopx[0]=0;
  loopx[1]=1;
  loopx[2]=1;
  loopx[3]=0;
  loopy[0]=0;
  loopy[1]=0;
  loopy[2]=1;
  loopy[3]=1;

  to=4.0*threshold;
  free(*sc);
  psc=(PAIR *)0;
  n=0;
  
  for(uli=0;uli<stenciln;uli++){
    ul=stencil[uli];
    i=ul%xd;
    j=ul/xd;
    if(i<xd-1 && j<yd-1){
      if(rho[ul]+rho[ul+1]+rho[ul+xd]+rho[ul+xd+1]>to){
	for(positive=pattern=m=0;m<4;m++){
	  if(phi[ul+loop[m]]>val){ 
	    positive++;
	    pattern|=1; 
	  }
	  pattern<<=1;
	}
	pattern>>=1;
	if((positive%4)){	  
	  if(Rotator(pattern,&rpattern,&position,&angle)){
	    n++;
	    if(!(psc=(PAIR *)realloc(psc,n*sizeof(PAIR)))){
	      printf("ANA Cannot reallocate for psc(1)\n");
	      return(1);
	    }
	    p1=position;
	    p2=(position+3)%4;
	    dumf1=(phi[ul+loop[p1]]-val)/(phi[ul+loop[p1]]-phi[ul+loop[p2]]);
	    psc[n-1].x1=(float)(i+loopx[p1])+dumf1*(loopx[p2]-loopx[p1]);
	    psc[n-1].y1=(float)(j+loopy[p1])+dumf1*(loopy[p2]-loopy[p1]);
	    p1=(position+(2-rpattern/3))%4;
	    p2=(p1+1)%4;
	    dumf2=(phi[ul+loop[p1]]-val)/(phi[ul+loop[p1]]-phi[ul+loop[p2]]);
	    psc[n-1].x2=(float)(i+loopx[p1])+dumf2*(loopx[p2]-loopx[p1]);
	    psc[n-1].y2=(float)(j+loopy[p1])+dumf2*(loopy[p2]-loopy[p1]);
	    psc[n-1].n=ul;
	    if(dumf1> 1.0 || dumf1<0.0 || dumf2> 1.0 || dumf2<0.0){
	      printf("ANA 1=%f 2=%f (pos=%i pat=%i ul=%li)\n",dumf1,dumf2,position,rpattern,ul);
	    }
	  }
	  else{
	    n+=2;
	    if(!(psc=(PAIR *)realloc(psc,n*sizeof(PAIR)))){
	      printf("ANA Cannot reallocate for psc(2)\n");
	      return(1);
	    }
	    printf("ANA Crossing @ %li\n",ul);
	    p1=position;
	    p2=(position+3)%4;
	    dumf=(phi[ul+loop[p1]]-val)/(phi[ul+loop[p1]]-phi[ul+loop[p2]]);
	    psc[n-2].x1=(float)(i+loopx[p1])+dumf*(loopx[p2]-loopx[p1]);
	    psc[n-2].y1=(float)(j+loopy[p1])+dumf*(loopy[p2]-loopy[p1]);
	    p1=(position+2)%4;
	    p2=(position+1)%4;
	    dumf=(phi[ul+loop[p1]]-val)/(phi[ul+loop[p1]]-phi[ul+loop[p2]]);
	    psc[n-2].x2=(float)(i+loopx[p1])+dumf*(loopx[p2]-loopx[p1]);
	    psc[n-2].y2=(float)(j+loopy[p1])+dumf*(loopy[p2]-loopy[p1]);
	    psc[n-2].n=ul;
	    p1=position;
	    p2=(position+1)%4;
	    dumf=(phi[ul+loop[p1]]-val)/(phi[ul+loop[p1]]-phi[ul+loop[p2]]);
	    psc[n-1].x1=(float)(i+loopx[p1])+dumf*(loopx[p2]-loopx[p1]);
	    psc[n-1].y1=(float)(j+loopy[p1])+dumf*(loopy[p2]-loopy[p1]);
	    p1=(position+2)%4;
	    p2=(position+3)%4;
	    dumf=(phi[ul+loop[p1]]-val)/(phi[ul+loop[p1]]-phi[ul+loop[p2]]);
	    psc[n-1].x2=(float)(i+loopx[p1])+dumf*(loopx[p2]-loopx[p1]);
	    psc[n-1].y2=(float)(j+loopy[p1])+dumf*(loopy[p2]-loopy[p1]);
	    psc[n-1].n=ul;
	  }
	}
      }
    }
  }
  
  *scn=n;
  *sc=psc;
  return(0);
}

int SmartContoursCartesian(int xd,int yd,float *x,float *y,PAIR **sc,int *scn,float val,float threshold){
  unsigned long xy,ul,uli;
  int i,j,m,n;
  PAIR *psc;
  float to;
  unsigned long loop[4];
  int loopx[4],loopy[4];
  int positive,pattern;
  int position,rpattern,angle;
  int p1,p2;
  float dumf,dumf1,dumf2,co,si;
  float phi[4];

  xy=(unsigned long)xd*(unsigned long)yd;

  loop[0]=0;
  loop[1]=1;
  loop[2]=1+xd;
  loop[3]=xd;

  loopx[0]=0;
  loopx[1]=1;
  loopx[2]=1;
  loopx[3]=0;
  loopy[0]=0;
  loopy[1]=0;
  loopy[2]=1;
  loopy[3]=1;

  co=(float)cos((double)val);
  si=(float)sin((double)val);

  if(do_verbose) printf("ANA co=%f si=%f\n",co,si);

  to=4.0*threshold*threshold;
  if(*sc) free(*sc);
  psc=(PAIR *)0;
  n=0;
  
  for(uli=0;uli<stenciln;uli++){
    ul=stencil[uli];
    i=ul%xd;
    j=ul/xd;
    if(i<xd-1 && j<yd-1){
      if(
	 HYPOT2(x[ul],y[ul])+HYPOT2(x[ul+1],y[ul+1])+HYPOT2(x[ul+xd],y[ul+xd])+HYPOT2(x[ul+xd+1],y[ul+xd+1])>to
	 && DOTPROD(x[ul],y[ul],co,si)>0.0 && DOTPROD(x[ul+1],y[ul+1],co,si)>0.0 
	 && DOTPROD(x[ul+xd],y[ul+xd],co,si)>0.0 && DOTPROD(x[ul+1+xd],y[ul+1+xd],co,si)>0.0
	 ){
	for(positive=pattern=m=0;m<4;m++){
	  phi[m]=VECTPROD(x[ul+loop[m]],y[ul+loop[m]],co,si);
	  if(phi[m]>0.0){ 
	    positive++;
	    pattern|=1; 
	  }
	  pattern<<=1;
	}
	pattern>>=1;
	if((positive%4)){	
	  if(Rotator(pattern,&rpattern,&position,&angle)){
	    n++;
	    if(!(psc=(PAIR *)realloc(psc,n*sizeof(PAIR)))){
	      printf("ANA Cannot reallocate for psc(1)\n");
	      return(1);
	    }
	    p1=position;
	    p2=(position+3)%4;
	    if((dumf1=phi[p1]-phi[p2])!=0.0) dumf1=phi[p1]/dumf1;
	    else{ 
	      dumf1=0.0;
	      //if(do_verbose) 
		printf("ANA Forced zero in 1 ul=%lu\n",ul);
	    }
	    psc[n-1].x1=(float)(i+loopx[p1])+dumf1*(loopx[p2]-loopx[p1]);
	    psc[n-1].y1=(float)(j+loopy[p1])+dumf1*(loopy[p2]-loopy[p1]);
	    p1=(position+(2-rpattern/3))%4;
	    p2=(p1+1)%4;
	    if((dumf2=phi[p1]-phi[p2])!=0.0) dumf2=phi[p1]/dumf2;
	    else{
	      dumf2=0.0;
	      //if(do_verbose) 
	      printf("ANA Forced zero in 2 ul=%lu\n",ul);
	    }
	    psc[n-1].x2=(float)(i+loopx[p1])+dumf2*(loopx[p2]-loopx[p1]);
	    psc[n-1].y2=(float)(j+loopy[p1])+dumf2*(loopy[p2]-loopy[p1]);
	    psc[n-1].n=ul;
	    if(do_verbose){ 
	      printf("ANA 1=%f 2=%f (pat=%i pos=%i rpat=%i ul=%li(%i %i))\n",dumf1,dumf2,pattern,position,rpattern,ul,i,j);
	      printf("ANA (%f %f) (%f %f)\n",psc[n-1].x1,psc[n-1].y1,psc[n-1].x2,psc[n-1].y2);
	    }
	    if(dumf1> 1.0 || dumf1<0.0 || dumf2> 1.0 || dumf2<0.0){
	      printf("ANA O 1=%f 2=%f (pos=%i pat=%i ul=%li)\n",dumf1,dumf2,position,rpattern,ul);
	    }
	  }
	  else{
	    n+=2;
	    if(!(psc=(PAIR *)realloc(psc,n*sizeof(PAIR)))){
	      printf("ANA Cannot reallocate for psc(2)\n");
	      return(1);
	    }
	    if(do_verbose) printf("ANA Crossing @ %li\n",ul);
	    p1=position;
	    p2=(position+3)%4;
	    dumf=phi[p1]/(phi[p1]-phi[p2]);
	    psc[n-2].x1=(float)(i+loopx[p1])+dumf*(loopx[p2]-loopx[p1]);
	    psc[n-2].y1=(float)(j+loopy[p1])+dumf*(loopy[p2]-loopy[p1]);
	    p1=(position+2)%4;
	    p2=(position+1)%4;
	    dumf=phi[p1]/(phi[p1]-phi[p2]);
	    psc[n-2].x2=(float)(i+loopx[p1])+dumf*(loopx[p2]-loopx[p1]);
	    psc[n-2].y2=(float)(j+loopy[p1])+dumf*(loopy[p2]-loopy[p1]);
	    psc[n-2].n=ul;
	    p1=position;
	    p2=(position+1)%4;
	    dumf=phi[p1]/(phi[p1]-phi[p2]);
	    psc[n-1].x1=(float)(i+loopx[p1])+dumf*(loopx[p2]-loopx[p1]);
	    psc[n-1].y1=(float)(j+loopy[p1])+dumf*(loopy[p2]-loopy[p1]);
	    p1=(position+2)%4;
	    p2=(position+3)%4;
	    dumf=phi[p1]/(phi[p1]-phi[p2]);
	    psc[n-1].x2=(float)(i+loopx[p1])+dumf*(loopx[p2]-loopx[p1]);
	    psc[n-1].y2=(float)(j+loopy[p1])+dumf*(loopy[p2]-loopy[p1]);
	    psc[n-1].n=ul;
	  }
	}
      }
    }
  }
  
  *scn=n;
  *sc=psc;
  return(0);
}

void Cartesian2Polar(unsigned long nn,float *x,float *y,float *rho,float *phi){
  unsigned long ul;
  float dumf;
  for(ul=0;ul<nn;ul++){
    dumf=(float)hypot((double)(y[ul]),(double)(x[ul]));
    phi[ul]=(float)atan2((double)(y[ul]),(double)(x[ul]));
    rho[ul]=dumf;
  }
  return;
}

void Polar2Cartesian(unsigned long nn,float *rho,float *phi,float *x,float *y){
  unsigned long ul;
  float dumf;
  for(ul=0;ul<nn;ul++){
    dumf =(float)((double)(rho[ul])*cos((double)(phi[ul])));
    y[ul]=(float)((double)(rho[ul])*sin((double)(phi[ul])));
    x[ul]=dumf;
  }
  return;
}

void Double2Float(unsigned long nn,double *pd,float *pf){
  unsigned long ul;
  for(ul=0;ul<nn;ul++){
    pf[ul]=(float)(pd[ul]);
  }
  return;  
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

int RectSegmentIntersect(POINT *pIniSegment,POINT *pFinSegment,RECT rect){
  POINT pSegment2[2];
  POINT pPoint[2];
  int iIntersection=0;
  int iX0=0,iX1=0,iY0=0,iY1=0;
  int i,iN=2,iPoint=0;

  /*  float fMin,fMax;

  if(pIniSegment[0].x>pIniSegment[1].x){
    fMin=pIniSegment[1].x;
    fMax=pIniSegment[0].x;
  }
  else{
    fMin=pIniSegment[0].x;
    fMax=pIniSegment[1].x;
  }
  */


  if(pIniSegment[0].x<rect.x1) iX0=-1;
  if(pIniSegment[0].x>rect.x2) iX0=1;
  if(pIniSegment[1].x<rect.x1) iX1=-1;
  if(pIniSegment[1].x>rect.x2) iX1=1;
  if(pIniSegment[0].y<rect.y1) iY0=-1;
  if(pIniSegment[0].y>rect.y2) iY0=1;
  if(pIniSegment[1].y<rect.y1) iY1=-1;
  if(pIniSegment[1].y>rect.y2) iY1=1;

  if(!iX0 && !iX1 && !iY0 && !iY1){
    memcpy(pFinSegment,pIniSegment,2*sizeof(POINT));
    return(1);
  }
  if(!iX0 && !iY0){
    memcpy(pFinSegment,pIniSegment,sizeof(POINT));
    iN--;
    iPoint=1;
  }
  if(!iX1 && !iY1){
    memcpy(pFinSegment+1,pIniSegment+1,sizeof(POINT));
    iN--;
    iPoint=0;
  }
  if(!iN){
    printf("RectSegmentIntersect Should not be here iN=0");
    return(1);
  }


  i=0;
  pSegment2[0].x=rect.x1;
  pSegment2[0].y=rect.y1;
  pSegment2[1].x=rect.x1;
  pSegment2[1].y=rect.y2;
  if(i!=iN && SegmentSegmentIntersect(pIniSegment,pSegment2,pPoint+i)){
    printf("Intersect 0\n");
    i++;
    iIntersection|=0x1;
  }
  pSegment2[0].x=rect.x1;
  pSegment2[0].y=rect.y2;
  pSegment2[1].x=rect.x2;
  pSegment2[1].y=rect.y2;
  if(i!=iN && SegmentSegmentIntersect(pIniSegment,pSegment2,pPoint+i)){
    printf("Intersect 1\n");
    i++;
    iIntersection|=0x2;
  }
  pSegment2[0].x=rect.x2;
  pSegment2[0].y=rect.y2;
  pSegment2[1].x=rect.x2;
  pSegment2[1].y=rect.y1;
  if(i!=iN && SegmentSegmentIntersect(pIniSegment,pSegment2,pPoint+i)){
    printf("Intersect 2\n");
    i++;
    iIntersection|=0x4;
  }
  pSegment2[0].x=rect.x2;
  pSegment2[0].y=rect.y1;
  pSegment2[1].x=rect.x1;
  pSegment2[1].y=rect.y1;
  if(i!=iN && SegmentSegmentIntersect(pIniSegment,pSegment2,pPoint+i)){
    printf("Intersect 3\n");
    i++;
    iIntersection|=0x8;
  }

  if(iIntersection){
    if(iN==1){
      pFinSegment[iPoint].x=pPoint[0].x;
      pFinSegment[iPoint].y=pPoint[0].y;
    }
    else{
      if(HYPOT2(pIniSegment[0].x-pPoint[0].x,pIniSegment[0].y-pPoint[0].y) <
	 HYPOT2(pIniSegment[0].x-pPoint[1].x,pIniSegment[0].y-pPoint[1].y)){
	pFinSegment[0].x=pPoint[0].x;
	pFinSegment[0].y=pPoint[0].y;
	pFinSegment[1].x=pPoint[1].x;
	pFinSegment[1].y=pPoint[1].y;
      }
      else{
	pFinSegment[0].x=pPoint[1].x;
	pFinSegment[0].y=pPoint[1].y;
	pFinSegment[1].x=pPoint[0].x;
	pFinSegment[1].y=pPoint[0].y;
      }
    }
    return(1);
  }
  else return(0);
}

int SegmentSegmentContinuationIntersect(POINT *pSegment1,POINT *pSegment2,POINT *pPoint){
  double vp,vp1;
  if((vp=VECTPROD(pSegment1[1].x-pSegment1[0].x,pSegment1[1].y-pSegment1[0].y,
	      pSegment2[1].x-pSegment2[0].x,pSegment2[1].y-pSegment2[0].y))==0.0){
    if(VECTPROD(pSegment1[1].x-pSegment1[0].x,pSegment1[1].y-pSegment1[0].y,
		pSegment2[0].x-pSegment1[0].x,pSegment2[0].y-pSegment1[0].y)==0.0){
      return(-1);
    }
    else{
      return(0);
    }
  }
  else{
    vp1=VECTPROD(pSegment2[0].x-pSegment1[0].x,pSegment2[0].y-pSegment1[0].y,
		 pSegment2[1].x-pSegment2[0].x,pSegment2[1].y-pSegment2[0].y);
    if(pPoint){
      pPoint->x=pSegment1[0].x+(pSegment1[1].x-pSegment1[0].x)*vp1/vp;
      pPoint->y=pSegment1[0].y+(pSegment1[1].y-pSegment1[0].y)*vp1/vp;
    }
    return(1);
  }
}

int SegmentSegmentIntersect(POINT *pSegment1,POINT *pSegment2,POINT *pPoint){
  double t3,t4,tt;
  int i;
  POINT point;

  if((i=SegmentSegmentContinuationIntersect(pSegment1,pSegment2,&point))){
    if(i==-1){
      if((pSegment1[1].x-pSegment1[0].x)!=0.0){
	t3=(pSegment2[0].x-pSegment1[0].x)/(pSegment1[1].x-pSegment1[0].x);
	t4=(pSegment2[1].x-pSegment1[0].x)/(pSegment1[1].x-pSegment1[0].x);
      }
      else{
	t3=(pSegment2[0].y-pSegment1[0].y)/(pSegment1[1].y-pSegment1[0].y);
	t4=(pSegment2[1].y-pSegment1[0].y)/(pSegment1[1].y-pSegment1[0].y);
      }
      if(t3>t4){
	tt=t3;
	t3=t4;
	t4=tt;
      }
      if(t3<=0.0){
	if(t4>=0.0 && t4<=1.0) tt=t4*0.5;
	else{
	  if(t4>=1.0) tt=0.5;
	  else return(0);
	}
      }
      else{
	if(t3<=1.0){
	  if(t4<=1.0) tt=(t3+t4)*0.5;
	  else{
	    if(t4>=1.0) tt=(t3+1.0)*0.5;
	    else return(0);
	  }
	}
	else return(0);
      }
      if(pPoint){
	pPoint->x=pSegment1[0].x+(pSegment1[1].x-pSegment1[0].x)*tt;
	pPoint->y=pSegment1[0].y+(pSegment1[1].y-pSegment1[0].y)*tt;
      }
      return(-1);
    }
    else{
      if(DOTPROD(pSegment1[0].x-point.x,pSegment1[0].y-point.y,
		 pSegment1[1].x-point.x,pSegment1[1].y-point.y)<=0.0 &&
	 DOTPROD(pSegment2[0].x-point.x,pSegment2[0].y-point.y,
		 pSegment2[1].x-point.x,pSegment2[1].y-point.y)<=0.0){
	if(pPoint) memcpy(pPoint,&point,sizeof(POINT));
	return(1);
      }
      else return(0);
    }
  }
  else return(0);
}

// Returns # of vortices inside circle (X-dX)**2+(Y-dY)**2=dRadius2
int IsUnitSquareInsideCircle(int iX,int iY,double dX,double dY,double dRadius2){
  int i=0,j;
  double dXX,dYY;
  for(j=0;j<4;j++){
    dXX=(double)(iX+j/2)-dX;
    dYY=(double)(iY+((j+1)%4)/2)-dY;
    if(dXX*dXX+dYY*dYY<=dRadius2) i++;
    else i--;
  }
  return(2+i/2);
}

// Returns non-zero if any part of a unit square is inside circle (X-dX)**2+(Y-dY)**2=dRadius2
int IsPartOfUnitSquareInsideCircle(int iX,int iY,double dX,double dY,double dRadius2){
  double dD[4];
  double dXX,dYY;
  int i,iSide;

  for(i=0;i<4;i++){
    dXX=(double)(iX+i/2)-dX;
    dYY=(double)(iY+((i+1)%4)/2)-dY;
    if((dD[i]=dXX*dXX+dYY*dYY-dRadius2)<0.0) return(1);
  }
  for(i=1,dYY=dD[0]+dD[1],iSide=0;i<4;i++){
    dXX=dD[i]+dD[(i+1)%4];
    if(dYY>dXX){
      dYY=dXX;
      iSide=i;
    }
  }
  
  if(iSide%2){
    dYY=(double)(iY+1-(iSide-1)/2)-dY;
    if(dYY*dYY<dRadius2) return(1);
    else return(0);
  }
  else{
    dXX=(double)(iX+iSide/2)-dX;
    if(dXX*dXX<dRadius2) return(1);
    else return(0);
  }
}


/*  Ran3 stuff */
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 0.000000001 //(1.0/MBIG)

/* Ran3 from Numerical Recipes */
double ran3(long *idum){
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  double dReturn;

  if (*idum < 0 || iff==0) {
    iff=1;
    mj = MSEED-(*idum<0 ? -*idum:*idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++){
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk<MZ) mk+=MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
	ma[i]-=ma[1+(i+30) % 55];
	if (ma[i]<MZ) ma[i]+=MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj<MZ) mj+=MBIG;
  ma[inext]=mj;
  dReturn=mj*FAC;
  if(dReturn<0.0 || dReturn>=1.0) printf("%f\n",dReturn);
  return mj*FAC;
}

