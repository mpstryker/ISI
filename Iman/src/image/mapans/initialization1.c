/* Map analysis                 */
/* Written by V.Kalatsky        */
/* 06Feb2001                    */
/* Modified 10May2001           */

#include <time.h>
#include "mapans1.h"
#include "fileheader.h"

//#define TEST

int LoadGreen(FILE *fp,char *cfile);
int GreenManager(int iVersion,FILE *fp,char *pc);

const char strISOI[]={"ISOI"};

int InitializeMaps(){
  FILE *fp;
  int i,dumi;
  unsigned short dumus;
  char *pcx=0,*pcy=0;
  double *binbuffer=0,*pd_s,*pd_d;
  double *pdTmpBuffer=NULL;
  double dDum;
  unsigned short xd,yd;
  unsigned long xy;
  int x,y,xb,yb;
  double rho,phi,*pd0,*pd1;
#ifdef TEST
  int j,k;
  double co,si;
  double xx,yy;
  double dSFCm=44.0,dScreenCm=40.0,dCoeff;
#endif

  lSeed=-(0xFFFFF&time(NULL));                                                          

  maps_orig=(double**)calloc(sizeof(double*),2);

  for(i=0;i<nfiles;i++){
    if((fp=fopen(filenames[i],"r"))==NULL){
      printf("INI Cannot open file %s\n",filenames[i]);
      return(1);
    }
    if(fread((void*)&dumi,sizeof(int),1,fp) != 1){
      printf("INI Cannot read type from file %s\n",filenames[i]);
      return(2);
    }
    if(memcmp(strISOI,&dumi,4)){
      printf("INI Data type %i\n",dumi);
    }
    else{
      printf("INI Green file\n");
      if(LoadGreen(fp,filenames[i])){
	fclose(fp);
	return(-2);
      }
      filestr=strdup(filenames[0]);
      type=SAVE_TYPE_BIT_GREEN;
      harmonic=1;
      harmonic_d=1.0;
      fclose(fp);
      break;
    }
    if(i==0){
      type=dumi;
    }
    else{
      if(dumi!=type){
	printf("INI Data type mismatch\n");
	return(4);
      }
    }
    if(i==0){
      if((type&SAVE_TYPE_BIT_XY_IN_FILES)){
	pcx=strstr(filenames[0],"mapraw");
	if(pcx) harmonic=(int)pcx[6]-48;
	filestr=strdup(filenames[0]);
	pcx=strstr(filestr,"mapraw");
	if(pcx) memmove(pcx,pcx+6,strlen(pcx+6)+1);
      }
      else{
	pcx=strstr(filenames[0],"mapxraw");
	pcy=strstr(filenames[0],"mapyraw");
	if(pcx) harmonic=(int)pcx[7]-48;
	if(pcy){ 
	  harmonic=(int)pcy[7]-48;
	  printf("INI WARNING: Y component is input before X\n");
	}
	filestr=strdup(filenames[0]);
	pcx=strstr(filestr,"mapxraw");
	pcy=strstr(filestr,"mapyraw");
	if(pcx) memmove(pcx,pcx+7,strlen(pcx+7)+1);
	if(pcy) memmove(pcy,pcy+7,strlen(pcy+7)+1);
      }
    }
    if(fread((void*)&dumus,sizeof(unsigned short),1,fp) != 1){
      printf("INI Cannot read Xdim from file %s\n",filenames[i]);
      return(2);
    }
    if(i==0){
      Xdim=dumus;
    }
    else{
      if(dumus!=Xdim){
	printf("INI Xdim mismatch %i(%i)\n",Xdim,dumus);
	return(4);
      }
    }
    if(fread((void*)&dumus,sizeof(unsigned short),1,fp) != 1){
      printf("INI Cannot read Ydim from file %s\n",filenames[i]);
      return(2);
    }
    if(i==0){
      Ydim=dumus;
      printf("INI Type=%i Xdim=%i Ydim=%i\n",type,Xdim,Ydim);
    }
    else{
      if(dumus!=Ydim){
	printf("INI Ydim mismatch %i(%i)\n",Ydim,dumus);
	return(4);
      }
    }
    if((type&SAVE_TYPE_BIT_XY_IN_FILES)){
      for(dumi=0;dumi<2;dumi++){
	if(!(maps_orig[dumi]=(double *)calloc((unsigned long)Xdim*(unsigned long)Ydim,sizeof(double)))){
	  printf("INI Cannot allocate for maps_orig[%i]\n",dumi);
	  return(3);
	}
	if(fread((void*)maps_orig[dumi],(unsigned long)Xdim*(unsigned long)Ydim*sizeof(double),1,fp) != 1){
	  printf("INI Cannot read maps_orig[%i] from file %s\n",dumi,filenames[i]);
	  return(2);
	}
      }
    }
    else{
      if(!(maps_orig[i]=(double *)calloc((unsigned long)Xdim*(unsigned long)Ydim,sizeof(double)))){
	printf("INI Cannot allocate for maps_orig[%i]\n",i);
	return(3);
      }
      if(fread((void*)maps_orig[i],(unsigned long)Xdim*(unsigned long)Ydim*sizeof(double),1,fp) != 1){
	printf("INI Cannot read maps_orig[%i] from file %s\n",i,filenames[i]);
	return(2);
      }
    }
    if((type&SAVE_TYPE_BIT_ADD_DOUBLE_AT_END)){
      if(fread((void*)&harmonic_d,sizeof(double),1,fp) != 1){
	printf("INI Cannot read harmonic from file %s\n",filenames[i]);
	harmonic_d=1.0;
      }
    }
    fclose(fp);
   }

  if(iBlankSource){
    memset(maps_orig[0],0,(unsigned long)Xdim*(unsigned long)Ydim*sizeof(double));
    memset(maps_orig[1],0,(unsigned long)Xdim*(unsigned long)Ydim*sizeof(double));
  }

  if(iRotateMaps){
    if(iRotateMaps%2){
      xd=Ydim;
      yd=Xdim;
    }
    else{
      xd=Xdim;
      yd=Ydim;
    }

    xy=(unsigned long)xd*(unsigned long)yd;
    if((pd_d=pdTmpBuffer=(double*)calloc(xy,sizeof(double)))){
      for(i=0;i<2;i++){
	pd_s=maps_orig[i];
	switch(iRotateMaps){
	case 1:
	  for(y=0;y<Ydim;y++){
	    for(x=0;x<Xdim;x++) pd_d[y+(Xdim-1-x)*Ydim]=pd_s[x+y*Xdim];
	  }
	  break;
	case 2:
	  for(y=0;y<Ydim;y++){
	    for(x=0;x<Xdim;x++) pd_d[(Xdim-1-x)+(Ydim-1-y)*Xdim]=pd_s[x+y*Xdim];
	  }
	  break;
	case 3:
	  for(y=0;y<Ydim;y++){
	    for(x=0;x<Xdim;x++) pd_d[(Ydim-1-y)+x*Ydim]=pd_s[x+y*Xdim];
	  }
	  break;
	default:
	  printf("INI Bad rotation angle %i\n",iRotateMaps);
	  memcpy(pd_d,pd_s,xy*sizeof(double));
	  break;
	}
	memcpy(maps_orig[i],pdTmpBuffer,xy*sizeof(double));
      }
      Xdim=xd;
      Ydim=yd;
      free(pdTmpBuffer);
    }
    else{
      printf("INI Cannot allocate for pdTmpBuffer. Will go on without rotation.\n");
    }    
  }

  if(iFlipMaps){
    for(i=0;i<2;i++){
      pd_s=maps_orig[i];
      if(iFlipMaps&INVERSION_X){
	for(y=0;y<Ydim;y++){
	  for(x=0;x<Xdim/2;x++){ 
	    dDum=pd_s[x+y*Xdim];
	    pd_s[x+y*Xdim]=pd_s[(Xdim-1-x)+y*Xdim];
	    pd_s[(Xdim-1-x)+y*Xdim]=dDum;
	  }	
	}
      }
      if(iFlipMaps&INVERSION_Y){
	for(x=0;x<Xdim;x++){ 
	  for(y=0;y<Ydim/2;y++){
	    dDum=pd_s[x+y*Xdim];
	    pd_s[x+y*Xdim]=pd_s[x+(Ydim-1-y)*Xdim];
	    pd_s[x+(Ydim-1-y)*Xdim]=dDum;
	  }	
	}
      }
    }
  }


  if(spatial_binning>1){
    printf("INI Running %ix%i binning\n",spatial_binning,spatial_binning);
    xd=Xdim/spatial_binning;
    yd=Ydim/spatial_binning;
    xy=(unsigned long)xd*(unsigned long)yd;
    if((binbuffer=(double*)malloc(xy*sizeof(double)))){
      for(i=0;i<2;i++){
	memset(binbuffer,0,xy*sizeof(double));
	for(y=0;y<yd;y++){
	  for(x=0;x<xd;x++){
	    pd_s=maps_orig[i]+(x+y*Xdim)*spatial_binning;
	    pd_d=binbuffer+(x+y*xd);
	    for(yb=0;yb<spatial_binning;yb++){
	      for(xb=0;xb<spatial_binning;xb++) *pd_d += *(pd_s+xb+yb*Xdim);
	    }
	  }
	}
	maps_orig[i]=(double*)realloc(maps_orig[i],xy*sizeof(double));
	memcpy(maps_orig[i],binbuffer,xy*sizeof(double));
      }
      Xdim=xd;
      Ydim=yd;
      free(binbuffer);
    }
    else{
      printf("INI Cannot allocate for binbuffer. Will go on without binning.\n");
    }
  }

  if(iDoPhaseScale && dPhaseScale<0.0){
    iDoPhaseScale=0;
    dPhaseScale=-dPhaseScale;
    xy=(unsigned long)Xdim*(unsigned long)Ydim;
    pd0=maps_orig[0];
    pd1=maps_orig[1];
    for(i=0;i<xy;i++,pd0++,pd1++){
      rho=hypot(*pd1,*pd0);
      /*
      phi=atan2(*pd1,*pd0)-M_PI;
      if(phi<=-M_PI) phi+=2.0*M_PI;
      phi*=dPhaseScale;
      if(phi<-M_PI) phi=-M_PI;
      if(phi>M_PI) phi=M_PI;
      phi+=M_PI;
      */
      phi=atan2(*pd1,*pd0);
      if(phi<0.0) phi+=2*M_PI;
      phi*=dPhaseScale;
      //      phi=dPhaseScale*atan2(*pd1,*pd0);
      *pd0=rho*cos(phi);
      *pd1=rho*sin(phi);
    }
  }

  if(iDoAmplitudeNorm){
    xy=(unsigned long)Xdim*(unsigned long)Ydim;
    pd0=maps_orig[0];
    pd1=maps_orig[1];
    for(i=0;i<xy;i++,pd0++,pd1++){
      *pd0/=dAmplitudeNormCoeff;
      *pd1/=dAmplitudeNormCoeff;
    }
 }
  //extern double dLensRadius;
  //extern int iLensCorrection;

  if(iLensCorrection){
    xy=(unsigned long)Xdim*(unsigned long)Ydim;
    pd0=maps_orig[0];
    pd1=maps_orig[1];
    for(y=0;y<Ydim;y++){
      for(x=0;x<Xdim;x++){
	dDum=1.0/sqrt(1.0-((Xdim*0.5-x)*(Xdim*0.5-x)+(Ydim*0.5-y)*(Ydim*0.5-y))/(dLensRadius*dLensRadius));
	pd0[y*Xdim+x]*=dDum;
	pd1[y*Xdim+x]*=dDum;
      }
    }
  }

  if(iBlendBorder && Xdim>iBlendBorderPixel && Ydim>iBlendBorderPixel){
    xy=(unsigned long)Xdim*(unsigned long)Ydim;
    pd0=maps_orig[0];
    pd1=maps_orig[1];
    for(y=0;y<Ydim;y++){
      for(x=0;x<Xdim;x++){
	if(y<iBlendBorderPixel){
	  pd0[y*Xdim+x]*=(float)y/(float)iBlendBorderPixel;
	  pd1[y*Xdim+x]*=(float)y/(float)iBlendBorderPixel;
	}
	if(y>Ydim-iBlendBorderPixel-1){
	  pd0[y*Xdim+x]*=(float)(Ydim-y-1)/(float)iBlendBorderPixel;
	  pd1[y*Xdim+x]*=(float)(Ydim-y-1)/(float)iBlendBorderPixel;
	}

	if(x<iBlendBorderPixel){
	  pd0[y*Xdim+x]*=(float)x/(float)iBlendBorderPixel;
	  pd1[y*Xdim+x]*=(float)x/(float)iBlendBorderPixel;
	}
	if(x>Xdim-iBlendBorderPixel-1){
	  pd0[y*Xdim+x]*=(float)(Xdim-x-1)/(float)iBlendBorderPixel;
	  pd1[y*Xdim+x]*=(float)(Xdim-x-1)/(float)iBlendBorderPixel;
	}
      }
    }
  }


#ifdef TEST

  //Xdim=100;
  //Ydim=100;
  Xdim=400;
  Ydim=400;
  //    Xdim=128;
  //    Ydim=96;
  //  Xdim=612;
  //  Ydim=783;
  //  Xdim=612*2;
  //  Ydim=783*2;
  for(i=0;i<2;i++){
    if(!(maps_orig[i]=(double *)realloc(maps_orig[i],(unsigned long)Xdim*(unsigned long)Ydim*sizeof(double)))){
      printf("INI Cannot reallocate for maps_orig[%i]\n",i);
      return(4);
    }
  }

  //  Xdim=48;
  //  Ydim=48;
  //Wedge dimensions
  // Fat wedge
  //    Xdim=256;
  //    Ydim=16;
  // Slim wedge
    //  Xdim=256;
    //  Ydim=8;

  harmonic_d=1.0;

  dCoeff=4*M_PI*dScreenCm/(dSFCm*Xdim);

   co=cos(M_PI/4.0); si=sin(M_PI/4.0);
  //  co=cos(M_PI/2.0); si=sin(M_PI/2.0);
  //  phi=2*M_PI*(80.0/66.0)/(double)Xdim;
  phi=2*M_PI/(0.5*(Xdim+Ydim));
  for(k=0;k<Ydim;k++){
    for(j=0;j<Xdim;j++){
      yy=(double)(k-Ydim/2);
      xx=(double)(j-Xdim/2);

      //Polar coord
      rho=hypot(xx,yy);

      //counterclockwise, blue up
      //       phi=-2.0*atan2(yy,xx)+M_PI;
       phi=atan2(yy,xx);
      // Radially modulated
      //      xx=rho*cos(phi);
      //      yy=rho*sin(phi);
      // Flat
       //      xx=cos(phi);
       //      yy=sin(phi);
      
      //Circular diagram
      
      //Square fill
      //maps_orig[0][k*Xdim+j]=co*xx-si*yy;
      //maps_orig[1][k*Xdim+j]=si*xx+co*yy;
      
      //Circle fill 
      
       //if(rho<Xdim*0.4999){
	  if(rho<Xdim*0.45 && rho>Xdim*0.35){ // Ring fill
	 // maps_orig[0][k*Xdim+j]=co*xx-si*(-yy);
	 // maps_orig[1][k*Xdim+j]=si*xx+co*(-yy);
	 maps_orig[0][k*Xdim+j]=cos(-phi);
	 maps_orig[1][k*Xdim+j]=sin(-phi);
       }
       else{
	 maps_orig[0][k*Xdim+j]=0.0;
	 maps_orig[1][k*Xdim+j]=0.0;
       }
      
      //Rect diagram

      //Horizontal stripes - yy
      //Vertical stripes - xx
      // modulated along x
      //      maps_orig[0][k*Xdim+j]=(xx+Xdim/2)*cos(yy*phi);
      //      maps_orig[1][k*Xdim+j]=-(xx+Xdim/2)*sin(yy*phi);
      // flat
      // y gradient
       // maps_orig[0][k*Xdim+j]=cos(yy*dCoeff);
       // maps_orig[1][k*Xdim+j]=-sin(yy*dCoeff);
      // x gradient
       //maps_orig[0][k*Xdim+j]=cos(xx*dCoeff);
       //maps_orig[1][k*Xdim+j]=sin(xx*dCoeff);
      
      
      //maps_orig[0][k*Xdim+j]=(double)(k-Ydim/2);
      //maps_orig[1][k*Xdim+j]=-(double)(j-Xdim/2);

      // Step
      //      if(yy>=0.0) maps_orig[0][k*Xdim+j]=1.0;
      //      else maps_orig[0][k*Xdim+j]=-1.0;
      //      maps_orig[1][k*Xdim+j]=0;

      //Random
      //maps_orig[0][k*Xdim+j]=ran3(&lSeed)-0.5;
      //maps_orig[1][k*Xdim+j]=ran3(&lSeed)-0.5;

      //Flat Phi, Gradient Rho
       //      maps_orig[0][k*Xdim+j]=0;
       //maps_orig[1][k*Xdim+j]=-j-k-Xdim*0.5+2*Xdim;

      //Gradient Phi x, Gradient Rho y
      //maps_orig[0][k*Xdim+j]=k*cos(j*2*M_PI/Xdim)/Ydim;
      //maps_orig[1][k*Xdim+j]=k*sin(j*2*M_PI/Xdim)/Ydim;

      // Pinwheel Spiral
      /*    
      if(rho<Xdim/2){
	//maps_orig[0][k*Xdim+j]=cos(2*phi)*(Xdim/4-rho)*(Xdim/4+rho);
	//maps_orig[1][k*Xdim+j]=sin(2*phi)*(Xdim/4-rho)*(Xdim/4+rho);
	//	if(rho<Xdim/6){
	// maps_orig[0][k*Xdim+j]=cos(2*phi)*cos(3.0*M_PI*rho/Xdim);
	// maps_orig[1][k*Xdim+j]=sin(2*phi)*cos(3.0*M_PI*rho/Xdim);
	  //}
	  //else{
	  maps_orig[0][k*Xdim+j]=cos(2*phi+M_PI*rho/Xdim)*cos(3.0*M_PI*rho/Xdim);
	  maps_orig[1][k*Xdim+j]=sin(2*phi+M_PI*rho/Xdim)*cos(3.0*M_PI*rho/Xdim);
	  //}
	//maps_orig[0][k*Xdim+j]=cos(2*phi+2.0*M_PI*rho/Xdim)*cos(M_PI*rho/Xdim);
	//maps_orig[1][k*Xdim+j]=sin(2*phi+2.0*M_PI*rho/Xdim)*cos(M_PI*rho/Xdim);
      }
      else{
	maps_orig[0][k*Xdim+j]=0;
	maps_orig[1][k*Xdim+j]=0;
      }
      */
    } 
  }
  //  maps_orig[0][(Xdim*Ydim)/2]=0.1;
  //  maps_orig[1][(Xdim*Ydim)/2]=0.01;
  /*
  for(k=0;k<Ydim;k++){
    for(j=0;j<Xdim;j++){
      printf("INI %i %i %i %f %f\n",k*Xdim+j,j,k,maps_orig[0][k*Xdim+j],maps_orig[1][k*Xdim+j]);
    } 
  }
  */
#endif

  mainRect.x1=0.0;
  mainRect.y1=0.0;
  mainRect.x2=(float)(Xdim-1);
  mainRect.y2=(float)(Ydim-1);
  
  XYdim=(unsigned long)Xdim*(unsigned long)Ydim;
  pminx=-0.5;
  pmaxx=(float)Xdim-0.5;
  pmaxy=-0.5;
  pminy=(float)Ydim-0.5;
  //pminy=-0.5;
  //pmaxy=(float)Ydim-0.5;
  
  printf("INI Maps initialized\n");
  return(0);
}

int InitializeBuffers(){
  int i;
  unsigned long ul;

  maps=(double**)calloc(sizeof(double*),2);
  for(i=0;i<2;i++){
    if(!(maps[i]=(double *)calloc(XYdim,sizeof(double)))){
      printf("INI Cannot allocate for maps[%i]\n",i);
      return(4);
    }
    memcpy(maps[i],maps_orig[i],XYdim*sizeof(double));
  }

  fmaps_orig=(float**)calloc(sizeof(float*),2);
  for(i=0;i<2;i++){
    if(!(fmaps_orig[i]=(float *)calloc(XYdim,sizeof(float)))){
      printf("INI Cannot allocate for fmaps_orig[%i]\n",i);
      return(4);
    }
    Double2Float(XYdim,maps[i],fmaps_orig[i]);
  }

  if(!(fbufferphi=(float*)malloc(XYdim*sizeof(float)))){
    printf("INI Cannot allocate for fbufferphi\n");
    return(1);
  }
  if(!(fbufferrho=(float*)malloc(XYdim*sizeof(float)))){
    printf("INI Cannot allocate for fbufferrho\n");
    return(1);
  }
  if(!(fbufferx=(float *)calloc(XYdim,sizeof(float)))){
    printf("INI Cannot allocate for fbufferx\n");
    return(2);
  }
  if(!(fbuffery=(float *)calloc(XYdim,sizeof(float)))){
    printf("INI Cannot allocate for fbuffery\n");
    return(3);
  }
  if(!(pdBufferX=(double *)calloc(XYdim,sizeof(double)))){
    printf("INI Cannot allocate for pdBufferX\n");
    return(2);
  }
  if(!(pdBufferY=(double *)calloc(XYdim,sizeof(double)))){
    printf("INI Cannot allocate for pdBufferY\n");
    return(3);
  }

  if(!(fbufferi=(float*)malloc(XYdim*sizeof(float)))){
    printf("INI Cannot allocate for fbufferi\n");
    return(1);
  }

  if(!(fbufferdiv=(float*)malloc((Xdim-1)*(Ydim-1)*sizeof(float)))){
    printf("INI Cannot allocate for fbufferdiv\n");
    return(1);
  }
  if(!(fbufferrot=(float*)malloc((Xdim-1)*(Ydim-1)*sizeof(float)))){
    printf("INI Cannot allocate for fbufferrot\n");
    return(1);
  }

  if(!(fbufferdivrot_phi=(float*)malloc((Xdim-1)*(Ydim-1)*sizeof(float)))){
    printf("INI Cannot allocate for fbufferdivrot_phi\n");
    return(1);
  }
  if(!(fbufferdivrot_rho=(float*)malloc((Xdim-1)*(Ydim-1)*sizeof(float)))){
    printf("INI Cannot allocate for fbufferdivrot_rho\n");
    return(1);
  }

  dbuffers=(double**)calloc(sizeof(double*),2);
  for(i=0;i<2;i++){
    if(!(dbuffers[i]=(double *)calloc(XYdim,sizeof(double)))){
      printf("INI Cannot allocate for dbuffers[%i]\n",i);
      return(4);
    }
  }

  smoothfbuffers=(float**)calloc(sizeof(float*),4);
  for(i=0;i<4;i++){
    if(!(smoothfbuffers[i]=(float *)calloc(XYdim,sizeof(float)))){
      printf("INI Cannot allocate for smoothfbuffers[%i]\n",i);
      return(4);
    }
  }

  if(!(chopfbuffer=(float *)calloc(XYdim,sizeof(float)))){
    printf("INI Cannot allocate for chopfbuffer\n");
    return(3);
  }

  stenciln=XYdim;
  if(!(stencil=(unsigned long *)malloc(stenciln*sizeof(unsigned long)))){
    printf("INI Cannot allocate for stencil\n");
    return(1);      
  }
  for(ul=0;ul<stenciln;ul++) stencil[ul]=ul;

  if(!(contours=(float*)malloc((depth_lines+1)*sizeof(float)))){
    printf("INI Cannot allocate for contours\n");    
  }

  smaxx=Xdim-1;
  sminx=0;
  smaxy=Ydim-1;
  sminy=0;

  radius=(int)floor(dradius);
  radiusbig=(int)floor(dradiusbig);
  cradius=(int)floor(dcradius);
  sradius=(int)floor(dsradius);

  if(radiusbig==0){ 
    if(display_mode==DISPLAY_MODE6) display_mode=DISPLAY_MODE4;
    if(save_mode==DISPLAY_MODE6) save_mode=DISPLAY_MODE4;
  }

  printf("INI Buffers initialized\n");
  return(0);
}


int LoadGreen(FILE *fp,char *cfile){
  DATA_CHUNK DATAChunk;
  SOFT_CHUNK SOFTChunk;
  FILEHEADER fileheader;
  char strTag[5];
  unsigned long ul;
  unsigned int nFileHeaderSize=0;
  unsigned int nFrameHeaderSize=0;
  int iVersion=0;
  int imagetype;

  rewind(fp);
  if(!ReadChunkHeader((void*)&DATAChunk,fp)){
    printf("INI Cannot read chunk header from file %s\n",cfile);
    return(1);
  }
  strTag[4]='\0';

  if(HasName((void*)&DATAChunk,"ISOI")){
    printf("INI Got ISOI fileheader from file %s\n",cfile);
      if((ul=FindChunkInFile(fp,"SOFT"))!=~0UL){
	if(ul+CHUNK_HEAD_SIZE!=sizeof(SOFT_CHUNK)){
	  printf("INI Not standard(%iu) SOFT chunk(%lu)\n",sizeof(SOFT_CHUNK),ul+CHUNK_HEAD_SIZE);
	}
	if(!ReadChunk((void *)&SOFTChunk,ul+CHUNK_HEAD_SIZE,fp)){
	  printf("INI Cannot read SOFT chunk from file %s\n",cfile);  
	  return(4); 
	}
	memcpy(strTag,SOFTChunk.Tag,CHUNK_TAG_SIZE);
	iVersion=VERSION_CHUNK;
	nFrameHeaderSize=SOFTChunk.FrameHeaderSize;
	printf("INI SubjectID %s\n",SOFTChunk.SubjectID);
      }
      else{
	printf("INI Cannot find SOFT chunk in file %s\n",cfile);
	return(3);
      }
      if(FindChunkInFileAndOffset(fp,"DATA",FALSE)!=~0UL){
	nFileHeaderSize=ftell(fp);
      }
      else{
	printf("INI Cannot find DATA chunk in file %s\n",cfile);
	return(3);
      }
    }
    else{
      if(fread((void*)&fileheader,sizeof(FILEHEADER),1,fp) != 1){
	printf("INI Cannot read fileheader from file %s\n",cfile);
	return(3);
      }
      iVersion=VERSION_HEADER;
      memcpy(strTag,fileheader.sumTag,CHUNK_TAG_SIZE);
      nFileHeaderSize=sizeof(FILEHEADER);
      nFrameHeaderSize=sizeof(FRAMEHEADER);
    }
    printf("INI Got header from %s TAG=%s\n",cfile,strTag);

    switch(*strTag){
    case 'A':
      imagetype = FILE_TYPE_ANALYSIS;
      printf("INI Analysis file. Use other routines to view image\n");
      return(1);
      break;
    case 'C':
      imagetype = FILE_TYPE_COMPRESSED;
      printf("INI Comressed raw data file\n");
      return(2);
      break;
    case 'G':
    case 'E':
      imagetype = FILE_TYPE_GREEN;
      return(GreenManager(iVersion,fp,cfile));
      break;
    case 'T':
      imagetype = FILE_TYPE_STREAM;
      return(3);
      break;
    default:
      imagetype = FILE_TYPE_UNKNOWN;
      printf("INI Unknown file type. Bailing out\n");
      return(-1);
    }
    return(0);
}

int GreenManager(int iVersion,FILE *fp,char *pc){
  int i;
  FILEHEADER fh;
  SOFT_CHUNK SOFTChunk;
  unsigned long ulDataType;
  unsigned char *puc,*pucBuffer=NULL;
  unsigned short *pus,*pusBuffer=NULL;
  unsigned long ul,*pul,*pulBuffer=NULL;
  float *pf,*pfBuffer=NULL;
  int iReturn=0;
  unsigned long ulXY;
  double *pd;

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

  ulXY=(unsigned long)Xdim*(unsigned long)Ydim;
  printf("GRE Dimensions: X %i Y %i Type %lu\n",Xdim,Ydim,ulDataType);

  for(i=0;i<2;i++){
    if(!(maps_orig[i]=(double *)calloc(ulXY,sizeof(double)))){
      printf("GRE Cannot allocate for maps_orig[%i]\n",i);
      return(3);
    }
  }

  switch(ulDataType){
  case DATATYPE_UCHAR:
    if(!(pucBuffer=(unsigned char*)malloc(ulXY*sizeof(unsigned char)))){
      printf("GRE Cannot allocate for pucBuffer\n");
      iReturn=7;
      goto bailout;
    }
    if((i=fread((void*)pucBuffer,ulXY*sizeof(unsigned char),1,fp)) != 1){
      printf("GRE Cannot read frame from green file (%i)\n",i);
      iReturn=8;
      goto bailout;
    }
    for(ul=0,puc=pucBuffer,pd=maps_orig[0];ul<ulXY;ul++) *(pd++)=(float)*(puc++);
    free(pucBuffer);
    break;
  case DATATYPE_USHORT:
    if(!(pusBuffer=(unsigned short*)malloc(ulXY*sizeof(unsigned short)))){
      printf("GRE Cannot allocate for pusBuffer\n");
      iReturn=9;
      goto bailout;
    }
    if((i=fread((void*)pusBuffer,ulXY*sizeof(unsigned short),1,fp)) != 1){
      printf("GRE Cannot read frame from green file (%i)\n",i);
      iReturn=10;
      goto bailout;
    }
    for(ul=0,pus=pusBuffer,pd=maps_orig[0];ul<ulXY;ul++) *(pd++)=(float)*(pus++);
    free(pusBuffer);
    break;
  case DATATYPE_ULONG:
    if(!(pulBuffer=(unsigned long*)malloc(ulXY*sizeof(unsigned long)))){
      printf("GRE Cannot allocate for pulBuffer\n");
      iReturn=11;
      goto bailout;
    }
    if((i=fread((void*)pulBuffer,ulXY*sizeof(unsigned long),1,fp)) != 1){
      printf("GRE Cannot read frame from green file (%i)\n",i);
      iReturn=12;
      goto bailout;
    }
    for(ul=0,pul=pulBuffer,pd=maps_orig[0];ul<ulXY;ul++) *(pd++)=(float)*(pul++);
    free(pulBuffer);
    break;
  case DATATYPE_FLOAT:
    if(!(pfBuffer=(float*)malloc(ulXY*sizeof(float)))){
      printf("GRE Cannot allocate for pfBuffer\n");
      iReturn=11;
      goto bailout;
    }
    if((i=fread((void*)pfBuffer,ulXY*sizeof(float),1,fp)) != 1){
      printf("GRE Cannot read frame from green file (%i)\n",i);
      iReturn=13;
      goto bailout;
    }
    for(ul=0,pf=pfBuffer,pd=maps_orig[0];ul<ulXY;ul++) *(pd++)=(float)*(pf++);
    free(pfBuffer);
    break;
  default:
    printf("GRE Unknown data type %lu",ulDataType);
    iReturn=6;
    goto bailout;
  }

 bailout:
  fclose(fp);
  return(iReturn);
}

int SetFilter(FILTER *pFilter,int xd,int yd){
  int iRadius;
  int iNP;
  int i,j;
  double dDum;
  double dCoef;

  iRadius=(int)pFilter->dRadius;
  if(iRadius> ((xd>yd) ? xd/2 : yd/2)){
    printf("INI Filter radius is too large (%i)\n",iRadius);
    return(1);
  }

  if(iRadius>0){
    SAFE_FREE(pFilter->piX);
    SAFE_FREE(pFilter->piY);
    SAFE_FREE(pFilter->pdValues);
    iNP=(1+2*iRadius)*(1+2*iRadius);
    if(!(pFilter->piX=(int *)calloc(iNP,sizeof(int))) || 
       !(pFilter->piY=(int *)calloc(iNP,sizeof(int))) || 
       !(pFilter->pdValues=(double *)calloc(iNP,sizeof(double)))){
      printf("AVE Cannot allocate for pFilter->piX/pFilter->piY/pFilter->pdValues\n");
      return(2);
    }
    dCoef=M_PI/pFilter->dRadius;
    pFilter->iNPoints=0;
    for(j=-iRadius;j<=iRadius;j++){
      for(i=-iRadius;i<=iRadius;i++){
	if((dDum=hypot((double)j,(double)i))<=pFilter->dRadius){ 
	  pFilter->piX[pFilter->iNPoints]=i;
	  pFilter->piY[pFilter->iNPoints]=j;
	  switch(pFilter->iKind){
	  case FILTER_KIND_COS:
	    pFilter->pdValues[pFilter->iNPoints]=pFilter->dDepth*0.5*(1.0+cos(dDum*dCoef));
	    break;
	  }
	  pFilter->iNPoints++;
	}
      }
    }
  }
  return(0);
}
