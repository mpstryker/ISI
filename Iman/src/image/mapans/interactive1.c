/* Map analysis                 */
/* Written by V.Kalatsky        */
/* 06Feb2001                    */
/* Modified 10May2001           */

#include "mapans1.h"
#include <termios.h> 

#define BIAS_WINDOW_SIZE 8.0
#define BIAS_PANEL_N 201
#define BIAS_RING_N 5
#define BIAS_RAY_N 8
#define BIAS_PANEL_POLAR 0
#define BIAS_PANEL_CARTESIAN 1
#define BIAS_PANEL_DEFAULT BIAS_PANEL_CARTESIAN

#define CORRELATION_TEST_NONE   0
#define CORRELATION_TEST_SHIFT  1
#define CORRELATION_TEST_RANDOM 2

#define DRAW_DOMAIN_FIRST  0
#define DRAW_DOMAIN_LAST   1
#define DRAW_DOMAIN_FINISH 2
#define DRAW_DOMAIN_CLOSED 3
#define DRAW_DOMAIN_OPEN   4

#define SWITCH_PANEL_MODE_REFRESH 0
#define SWITCH_PANEL_MODE_SWITCH  1
#define SWITCH_PANEL_MODE_STENCIL 2

void Biaspanel(int id,int type,float fmax);
void Bullseye(float rmax,int rdiv,int phidiv);
void Mesh(float xmin,float xmax,float ymin,float ymax,int xdiv,int ydiv);
void SwitchPanel(int xd,int yd,float *buf,double **mp,int *mode,float *pfMin,float *pfMax,int iMode);
void RestorePanel(int xd,int yd,float *buf,int *mode,float fmin,float fmax);
char* GetString(char *str);
int GetFloatPoint(float *pfX,float *pfY,int iMode,int iColor);
int GetIntPoint(int *piX,int *piY,int iDX,int iDY,int iMode,int iColor);
int GetStringFromSTDIN(char *str,int iNChars,char *strPrompt);
int ReadDomain(int *piN,POINTD **ppDom);
int WriteDomain(int iN,POINTD *pDom);
int DrawDomain(int iN,POINTD *pDom,int iMode);

void Interactive(int xd,int yd,double **mp,int *mode){
  float fIMin=0.0,fIMax=0.0,x=0.0,y=0.0;
  int interactiveid=0,biasid=0,callerid;
  int xi,yi,x0,y0;
  char ch,str[256],strDum[256],*pc;
  int i,k,l;
  double phi,rho,phip=0.0,rhop=0.0,dumd;
  unsigned long j,xy;
  unsigned short dumus;
  int bias_panel_type=BIAS_PANEL_DEFAULT;
  FILE *fp=NULL,*pF=NULL;
  double co,si;
  POINT tmpSegment[2];
  int iCropZerosOnSave=1;
  int iXDCrop,iYDCrop;
  unsigned long ulXYDCrop;
  int iXMinCrop,iYMinCrop;
  int iXMaxCrop,iYMaxCrop;
  int iXYShift;
  int iError=0;
  float fX1,fY1,fX2,fY2;
  int iXL,iYT,iXR,iYB;
  double dX,dY,dZ,dN1,dN2;
  double *pdX,*pdY;
  int iReturn;
  int iBatchOutput=0;
  RECORD *pR;
  int iDum;
  int iCorrelationTest=0;
  int iCorrelationTestShiftX=0,iCorrelationTestShiftY=0;
  //  int iCorrelationTestRandomX=0,iCorrelationTestRandomY=100;
  //  int iCorrelationTestRandomXDim=500,iCorrelationTestRandomYDim=350;
  RECT stRect;
  double dDum;
  int iAddPixelsRadius=1;
  int iPanelMode=0;
  int iShuffle=0;

  if(*mode!=INT_MODE_RHO && *mode!=INT_MODE_PHI) return;

  cpgqid(&callerid);

  xy=(unsigned long)xd*(unsigned long)yd;

  stRect.x1=0;stRect.x2=xd-1;stRect.y1=0;stRect.y2=yd-1;

  setenv("PGPLOT_ENVOPT","I",1); 
  interactiveid=cpgopen("/xw");
  cpgask(0);
  cpgpap(fInteractiveWindowWidth,(float)yd/(float)xd);
  cpgsch(0.8);
  cpgenv(pminx,pmaxx,pminy,pmaxy,1,0);
  unsetenv("PGPLOT_ENVOPT");

  SwitchPanel(xd,yd,fbufferi,mp,mode,&fIMin,&fIMax,iPanelMode|SWITCH_PANEL_MODE_SWITCH);
  
  rhop=hypot(mp[0][0],mp[1][0]);
  phip=180.0/M_PI*atan2(mp[1][0],mp[0][0]);

  while(1){
    cpgslct(interactiveid);
    cpgband(7,1,x,y,&x,&y,&ch); 
    if(do_verbose) printf("INT X %f Y %f C (%i)%c\n",x,y,ch,ch);
    //printf("INT qC (%i)%c\n",ch,ch);
    // Keys in use: 0,A,B,D,F,G,H,I,J,K,M,O,P,Q,R,S,T,W,Z,a,b,c,d,e,g,h,i,l,m,n,o,p,q,r,s,t,u,v

    if(ch=='q'){ break;}

    if(ch=='Q'){ exit(0);}

    if(ch=='m'){
      if(GetString(str)) x=atof(str);
      else continue;
      if(GetString(str)) y=atof(str);
      else continue;
      printf("INT %f %f ",x,y);
    }

    if(ch=='c'){
      iNDomainPoints=0;
      SAFE_FREE(pDomain);
      cpgsci(1);
      
    label1:
      RestorePanel(xd,yd,fbufferi,mode,fIMin,fIMax);
      cpgband(0,0,x,y,&x,&y,&ch);
      if(ch=='q'){
	continue;
      }
      if(ch=='c'){
	Stencil(xd,yd,0,pDomain,&stenciln,&stencil);	          
	DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
	RestorePanel(xd,yd,fbufferi,mode,fIMin,fIMax);
	continue;
      }
      if(ch=='r'){
	if(GetString(str)) fX1=atof(str);
	else continue;
	if(GetString(str)) fY1=atof(str);
	else continue;
	if(GetString(str)) fX2=atof(str);
	else continue;
	if(GetString(str)) fY2=atof(str);
	else continue;

	if(iShiftDomainForIntegerPoints){
	  fX1-=0.5;
	  fX2-=0.5;
	  fY1-=0.5;
	  fY2-=0.5;
	}

	iNDomainPoints=4;
	pDomain=(POINTD *)malloc((iNDomainPoints+1)*sizeof(POINTD));
	pDomain[0].x=pDomain[iNDomainPoints].x=(double)fX1;
	pDomain[0].y=pDomain[iNDomainPoints].y=(double)fY1;
	pDomain[1].x=(double)fX2;
	pDomain[1].y=(double)fY1;
	pDomain[2].x=(double)fX2;
	pDomain[2].y=(double)fY2;
	pDomain[3].x=(double)fX1;
	pDomain[3].y=(double)fY2;

	DrawDomain(iNDomainPoints,pDomain,DRAW_DOMAIN_CLOSED);

	// Dum hack for auto-correlation, one pixel shift
	for(iXR=iYB=-1,iXL=xd+1,iYT=yd+1,i=0;i<4;i++){
	  if(iXL>(int)pDomain[i].x) iXL=(int)(pDomain[i].x);
	  if(iYT>(int)pDomain[i].y) iYT=(int)(pDomain[i].y);
	  if(iXR<(int)pDomain[i].x) iXR=(int)(pDomain[i].x);
	  if(iYB<(int)pDomain[i].y) iYB=(int)(pDomain[i].y);
	}
	if(iXL<0) iXL=0;
	if(iYT<0) iYT=0;
	if(iXR>=xd) iXR=xd-1;
	if(iYB>=yd) iYB=yd-1;
	//	printf("\nINT X L %i R %i Y T %i B %i\n",iXL,iXR,iYT,iYB);
	if(iXR<xd-1 && iYT<yd-1){
	  pdX=mp[0];
	  pdY=mp[1];
	  dX=dY=dN1=dN2=0.0;
	  for(j=iYT;j<iYB;j++){
	    for(i=iXL;i<iXR;i++){
	      l=i+xd*j;
	      k=l+xd+1;
	      dX+=pdX[l]*pdX[k]+pdY[l]*pdY[k];
	      dY+=pdX[l]*pdY[k]-pdY[l]*pdX[k];
	      dN1+=pdX[l]*pdX[l]+pdY[l]*pdY[l];
	      dN2+=pdX[k]*pdX[k]+pdY[k]*pdY[k];
	    }
	  }
	  printf("\nINT Corr %f\n",sqrt((dX*dX+dY*dY)/(dN1*dN2)));
	}

	Stencil(xd,yd,iNDomainPoints,pDomain,&stenciln,&stencil);	    
	DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
	if(do_fields){ 
	  DisplayFields(field_mode,Xdim,Ydim,maps[0],maps[1],dbuffers[0],dbuffers[1]);
	}
	continue;
      }
      if(ch=='m'){
	if(GetString(str)) x=atof(str);
	else continue;
	if(GetString(str)) y=atof(str);
	else continue;
      }
      iNDomainPoints=1;
      pDomain=(POINTD *)malloc((iNDomainPoints)*sizeof(POINTD));
      pDomain[iNDomainPoints-1].x=(double)x;
      pDomain[iNDomainPoints-1].y=(double)y;

      DrawDomain(iNDomainPoints,pDomain,DRAW_DOMAIN_FIRST);

      printf("%i %f %f\n",iNDomainPoints,x,y);
      while(1){
	cpgband(1,0,pDomain[iNDomainPoints-1].x,pDomain[iNDomainPoints-1].y,&x,&y,&ch);
	if(ch=='q'){
	  iNDomainPoints=0;
	  SAFE_FREE(pDomain);
	  RestorePanel(xd,yd,fbufferi,mode,fIMin,fIMax);
	  break;
	}
	if(ch=='A' || ch=='m'){
	  if(ch=='m'){
	    if(GetString(str)) x=atof(str);
	    else continue;
	    if(GetString(str)) y=atof(str);
	    else continue;
	  }
	  iNDomainPoints++;
	  pDomain=(POINTD *)realloc((void *)pDomain,iNDomainPoints*sizeof(POINTD));
	  pDomain[iNDomainPoints-1].x=(double)x;
	  pDomain[iNDomainPoints-1].y=(double)y;
	  DrawDomain(iNDomainPoints,pDomain,DRAW_DOMAIN_LAST);
	  printf("%i %f %f\n",iNDomainPoints,x,y);
	}
	if(ch=='D'){
	  iNDomainPoints=0;
	  SAFE_FREE(pDomain);
	  goto label1;
	}
	if(ch=='X'){
	  if(iNDomainPoints>2){ 
	    pDomain=(POINTD *)realloc((void *)pDomain,(iNDomainPoints+1)*sizeof(POINTD));
	    pDomain[iNDomainPoints].x=pDomain[0].x;
	    pDomain[iNDomainPoints].y=pDomain[0].y;

	    DrawDomain(iNDomainPoints,pDomain,DRAW_DOMAIN_FINISH);

	    Stencil(xd,yd,iNDomainPoints,pDomain,&stenciln,&stencil);	    
	    DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
	    if(do_fields){ 
	      DisplayFields(field_mode,Xdim,Ydim,maps[0],maps[1],dbuffers[0],dbuffers[1]);
	    }
	    break;
	  }
	  else printf("INT Needs three points to complete selection\n");
	}
	if(ch=='p'){
	  SwitchPanel(xd,yd,fbufferi,mp,mode,&fIMin,&fIMax,iPanelMode|SWITCH_PANEL_MODE_SWITCH);
	  DrawDomain(iNDomainPoints,pDomain,DRAW_DOMAIN_OPEN);
	  continue;
	}
     }
      continue;
    }

    if(ch=='p'){
      SwitchPanel(xd,yd,fbufferi,mp,mode,&fIMin,&fIMax,iPanelMode|SWITCH_PANEL_MODE_SWITCH);
      DrawDomain(iNDomainPoints,pDomain,DRAW_DOMAIN_CLOSED);
      continue;
    }

    if(ch=='i'){
      inversion=-inversion;
      printf("INT Inversion=%i\n",(int)inversion);
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      continue;
    }

    if(ch=='a'){
      if(GetString(str)) alpha=atof(str);
      else continue;
      printf("INT Alpha=%f\n",alpha);
      alpha*=M_PI/180.0;
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      continue;
    }

    if(ch=='d'){
      if(GetString(str)) depth_lines=atoi(str);
      else continue;
      if(depth_lines<0){
	depth_lines=0;
      }
      if(depth_lines){
	if(!(contours=(float*)realloc(contours,(depth_lines+1)*sizeof(float)))){
	  printf("INT Cannot reallocate for contours\n");
	  contours=(float*)0;
	  depth_lines=0;
	}
	else{
	  printf("INT Depth_lines=%i\n",depth_lines);
	  DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
	}
      }
      continue;
    }

    if(ch=='s'){
      if(GetString(str)) contour_scheme=atoi(str);
      else continue;
      if(contour_scheme>10 || contour_scheme<0 || contour_scheme==3 || contour_scheme==7){
	printf("INT Bad contour scheme %i, using default %i\n",contour_scheme,CONTOUR_DEFAULT);
	contour_scheme=CONTOUR_DEFAULT;
      }
      contour_scheme_phi=contour_scheme & 3;
      contour_scheme_rho=contour_scheme >> 2;
      
      printf("INT Contour_scheme=%i(%i %i)\n",contour_scheme,contour_scheme_phi,contour_scheme_rho);
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      continue;
    }

    if(ch=='r'){
      if(do_rescaling){
	do_rescaling=0;
	memcpy(transs,trans,6*sizeof(float));
      }
      else{
	do_rescaling=1;
      }
      continue;
    }

    if(ch=='v'){
      do_verbose=1-do_verbose;
      continue;
    }

    if(ch=='b'){
      /*
      fmin=BIG_FLOAT;
      fmax=-fmin;
      for(i=0;i<xy;i++){
	dumf=(float)hypot(mp[0][i],mp[1][i]);
	if(i>=xd){ // Do not take into account first screwed row of each frame
	  if(dumf>fmax){ 
	    fmax=dumf;
	  }
	  if(dumf<fmin){
	    fmin=dumf;
	  }
	}
      }
      */

      biasid=cpgopen("/xw");
      cpgask(0);
      Biaspanel(biasid,bias_panel_type,fIMax);
      while(1){
	cpgslct(biasid);
	cpgband(7,1,x,y,&x,&y,&ch);

	// Keys in use: q,Q,0,p,h,m

	if(ch=='Q'){
	  exit(1);
	}

	if(ch=='q'){
	  cpgslct(biasid);
	  cpgclos();
	  biasid=0;
	  break;
	}

	if(ch=='0'){
	  x=xshift=0.0;
	  y=yshift=0.0;
	}

	if(ch=='p'){
	  Biaspanel(biasid,(bias_panel_type=1-bias_panel_type),fIMax);
	  continue;
	}

	if(ch=='m'){
	  if(GetString(str)) x=atof(str);
	  else continue;
	  if(GetString(str)) y=atof(str);
	  else continue;
	}

	if(ch=='h'){
	  printf("q - quit interactive dialog\n");
	  printf("Q - quit program\n");
	  printf("0 - null shift\n");
	  printf("p - switch, polar <-> cartesian\n");
	  printf("h - print this help\n");
	  continue;
	}

	if(bias_panel_type==BIAS_PANEL_POLAR){
	  xshift=x;
	  yshift=y;
	}
	if(bias_panel_type==BIAS_PANEL_CARTESIAN){
	  xshift=(float)(y*cos((double)x*M_PI/180.0));
	  yshift=(float)(y*sin((double)x*M_PI/180.0));
	  if(y<0.0){ 
	    y=-y;
	    x+=180.0;
	  }
	  if(x>180.0) x-=360.0;
	  if(x<-180.0) x+=360.0;
	}
	printf("INT Bias X=%f Y=%f (R=%f P=%f)\n",xshift,yshift,hypot((double)yshift,(double)xshift),180.0/M_PI*atan2((double)yshift,(double)xshift));
	DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      }
      continue;
    }

    if(ch=='D'){
      do_distribution=1;
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);      
    }

    if(ch=='0'){
      inversion=1;
      xshift=0.0;
      yshift=0.0;
      alpha=0.0;
      amplitude_chop_l=0.0;
      amplitude_chop_u=1.0;
      smartcontours_threshold=DEFAULT_SMART_CONTOURS_THRESHOLD;
      Stencil(xd,yd,0,pDomain,&stenciln,&stencil);	          
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);      
      RestorePanel(xd,yd,fbufferi,mode,fIMin,fIMax);
      continue;
    }

    if(ch=='n'){
      xshift=yshift=0.0;
      for(j=0;j<stenciln;j++){
	xshift-=mp[0][stencil[j]];
	yshift-=mp[1][stencil[j]];
      }
      xshift/=(double)stenciln;
      yshift/=(double)stenciln;
      printf("INT Bias X=%f Y=%f (R=%f P=%f)\n",xshift,yshift,hypot((double)yshift,(double)xshift),180.0/M_PI*atan2((double)yshift,(double)xshift));
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);      
      continue;
    }
    if(ch=='S'){
      if(iCropZerosOnSave){
	iXMinCrop=iXMaxCrop=stencil[0]%xd;
	iYMinCrop=iYMaxCrop=stencil[0]/xd;
	for(j=1;j<stenciln;j++){
	  i=stencil[j];
	  xi=i%xd;
	  yi=i/xd;
	  if(iXMinCrop>xi) iXMinCrop=xi;
	  if(iXMaxCrop<xi) iXMaxCrop=xi;
	  if(iYMinCrop>yi) iYMinCrop=yi;
	  if(iYMaxCrop<yi) iYMaxCrop=yi;
	}
	iXDCrop=1+iXMaxCrop-iXMinCrop;
	iYDCrop=1+iYMaxCrop-iYMinCrop;
	printf("INT Crop box Dim (%i %i) Origin (%i %i)\n",iXDCrop,iYDCrop,iXMinCrop,iYMinCrop);
      }
      else{
	iXDCrop=xd;
	iYDCrop=yd;
	iXMinCrop=iYMinCrop=0;
	iXMaxCrop=xd-1;
	iYMaxCrop=yd-1;
      }
      ulXYDCrop=(unsigned long)iXDCrop*(unsigned long)iYDCrop;
      iXYShift=iXMinCrop+iYMinCrop*xd;

      memset((void *)dbuffers[0],0,(size_t)xd*(size_t)yd*sizeof(double));
      memset((void *)dbuffers[1],0,(size_t)xd*(size_t)yd*sizeof(double));
      co=cos(alpha);
      si=sin(alpha);
      for(j=0;j<stenciln;j++){
	i=stencil[j];
	dbuffers[0][i]=(mp[0][i]+xshift)*co-inversion*(mp[1][i]+yshift)*si;
	dbuffers[1][i]=inversion*(mp[1][i]+yshift)*co+(mp[0][i]+xshift)*si;
      }
      //      if(iXYShift){
      //	printf("INT XYShift=%i\n",iXYShift);
	for(j=0;j<ulXYDCrop;j++){
	  i=iXYShift+j%iXDCrop+(j/iXDCrop)*xd;
	  dbuffers[0][j]=dbuffers[0][i];
	  dbuffers[1][j]=dbuffers[1][i];
	}
	//      }
      if((type&SAVE_TYPE_BIT_XY_IN_FILES)==SAVE_XY_IN_ONE_FILE*SAVE_TYPE_BIT_XY_IN_FILES){
	if(pcSaveFileSuffix && *pcSaveFileSuffix) sprintf(str,"%s%s",filenames[0],pcSaveFileSuffix);
	else sprintf(str,"%s_a",filenames[0]);
	printf("INT Saving curent maps in %s\n",str);
	if((fp=fopen(str,"w"))==NULL){
	  printf("INT Cannot open file %s\n",str);
	}
	else{
	  i=type&0xE;
	  if(fwrite((void*)&i,sizeof(int),1,fp) != 1){
	    printf("INT Cannot write type %i to file %s\n",i,str);
	  }
	  dumus=(unsigned short)iXDCrop;
	  if(fwrite((void*)&dumus,sizeof(unsigned short),1,fp) != 1){
	    printf("INT Cannot write Xdim=%i to file %s\n",dumus,str);
	  }
	  dumus=(unsigned short)iYDCrop;
	  if(fwrite((void*)&dumus,sizeof(unsigned short),1,fp) != 1){
	    printf("INT Cannot write Ydim=%i to file %s\n",dumus,str);
	  }
	  if(fwrite((void*)dbuffers[0],ulXYDCrop*sizeof(double),1,fp) != 1){
	    printf("INT Cannot write dbuffer[0] to file %s\n",str);
	  }
	  if(fwrite((void*)dbuffers[1],ulXYDCrop*sizeof(double),1,fp) != 1){
	    printf("INT Cannot write dbuffer[1] to file %s\n",str);
	  }
	  if((SAVE_TYPE_BIT_ADD_DOUBLE_AT_END&type)){
	    if(fwrite((void*)&harmonic_d,sizeof(double),1,fp) != 1){
	      printf("INT Cannot write harmonic %f to file %s\n",harmonic_d,str);
	    }
	  }
	  fclose(fp);
	}
      }
      else{
	for(k=0;k<2;k++){
	  if(pcSaveFileSuffix && *pcSaveFileSuffix) sprintf(str,"%s%s",filenames[k],pcSaveFileSuffix);
	  else sprintf(str,"%s_a",filenames[k]);
	  printf("INT Saving curent %cmap in %s\n",88+k,str);
	  if((fp=fopen(str,"w"))==NULL){
	    printf("INT Cannot open file %s\n",str);
	  }
	  else{
	    i=type&0xE;
	    if(fwrite((void*)&i,sizeof(int),1,fp) != 1){
	      printf("INT Cannot write type %i to file %s\n",i,str);
	    }
	    dumus=(unsigned short)iXDCrop;
	    if(fwrite((void*)&dumus,sizeof(unsigned short),1,fp) != 1){
	      printf("INT Cannot write Xdim=%i to file %s\n",dumus,str);
	    }
	    dumus=(unsigned short)iYDCrop;
	    if(fwrite((void*)&dumus,sizeof(unsigned short),1,fp) != 1){
	      printf("INT Cannot write Ydim=%i to file %s\n",dumus,str);
	    }
	    if(fwrite((void*)dbuffers[k],ulXYDCrop*sizeof(double),1,fp) != 1){
	      printf("INT Cannot write dbuffer[%i] to file %s\n",k,str);
	    }
	    if(type>1){
	      if(fwrite((void*)&harmonic_d,sizeof(double),1,fp) != 1){
		printf("INT Cannot write harmonic %f to file %s\n",harmonic_d,str);
	      }
	    }
	    fclose(fp);
	  }
	}
      }
      continue;
    }

    if(ch=='e'){
      switch(pallet){
      case 2: 
	pallet=3;
	break;
      case 3:
	pallet=4;
	break;
      case 4:
	pallet=5;
	break;
      case 5:
	pallet=2;
	break;
      }
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);      
    }

    if(ch=='l'){
      if(GetString(str)) amplitude_chop_l=atof(str)/100.0;
      else continue;
      if(amplitude_chop_l<0.0 || amplitude_chop_l>amplitude_chop_u){
	printf("INT Fixing lower chop=%f\n",amplitude_chop_l);	
	amplitude_chop_l=0.0;
      }
      printf("INT Lower chop=%f\n",amplitude_chop_l);
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      continue;      
    }

    if(ch=='u'){
      if(GetString(str)) amplitude_chop_u=atof(str)/100.0;
      else continue;
      if(amplitude_chop_u>1.0 || amplitude_chop_u<amplitude_chop_l){
	printf("INT Fixing lower chop=%f\n",amplitude_chop_u);	
	//Val: quick hack
	//	amplitude_chop_u=1.0;
      }
      printf("INT Upper chop=%f\n",amplitude_chop_u);
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      continue;      
    }

    if(ch=='t'){
      if(do_smartcontours){
	if(GetString(str)) smartcontours_threshold=atof(str)/100.0;
	else continue;
	DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
	RestorePanel(xd,yd,fbufferi,mode,fIMin,fIMax);
      }
      continue;      
    }

    if(ch=='O'){
      color_scheme=(color_scheme+1)%(COLOR_N);
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      continue;      
    }

    if(ch=='g'){
      cpgsci(1);
      
      RestorePanel(xd,yd,fbufferi,mode,fIMin,fIMax);
      cpgband(0,0,x,y,&(tmpSegment[0].x),&(tmpSegment[0].y),&ch);
      if(ch=='q'){
	continue;
      }

      if(ch=='m'){
	printf("INT B ");
	if(GetString(str)) tmpSegment[0].x=atof(str);
	else continue;
	printf(" X %f ",tmpSegment[0].x);
	if(GetString(str)) tmpSegment[0].y=atof(str);
	else continue;
	printf(" Y %f\n",tmpSegment[0].y);
      }

      cpgmove(tmpSegment[0].x,tmpSegment[0].y);
      cpgsci(2);
      cpgpt1(tmpSegment[0].x,tmpSegment[0].y,2);
      cpgsci(1);
      do{
	cpgband(1,0,tmpSegment[0].x,tmpSegment[0].y,&(tmpSegment[1].x),&(tmpSegment[1].y),&ch);
	if(ch=='q'){
	  RestorePanel(xd,yd,fbufferi,mode,fIMin,fIMax);
	  goto looplabel;
	}
	if(ch=='m'){
	  printf("INT E ");
	  if(GetString(str)) tmpSegment[1].x=atof(str);
	  else continue;
	  printf(" X %f ",tmpSegment[1].x);
	  if(GetString(str)) tmpSegment[1].y=atof(str);
	  else continue;
	  printf(" Y %f\n",tmpSegment[1].y);
	}
      }while(!RectSegmentIntersect(tmpSegment,section,mainRect));
      
      cpgdraw(tmpSegment[1].x,tmpSegment[1].y);
      cpgsci(2);
      cpgpt1(tmpSegment[1].x,tmpSegment[1].y,2);
      cpgsci(2);
      cpgpt1(section[0].x,section[0].y,2);
      cpgsci(3);
      cpgpt1(section[1].x,section[1].y,2);
      cpgsci(1);
      do_sectionplot=1;
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      RestorePanel(xd,yd,fbufferi,mode,fIMin,fIMax);
      goto looplabel;
    }


    if(ch=='R'){
      if(GetStringFromSTDIN(str,255,"INT Enter contour file: ")) continue;

      if((fp=fopen(str,"r"))==NULL){
	printf("INT Cannot open file <%s>\n",str);
	continue;
      }
      if(fread(&i,sizeof(int),1,fp) != 1){
	printf("INT Cannot read type %i from file %s\n",i,str);
      }
      if(fread(&i,sizeof(int),1,fp) != 1){
	printf("INT Cannot read Xdim from file %s\n",str);
      }
      if(i!=xd) printf("INT WARNING X dimension missmatch %i %i\n",xd,i);
      if(fread(&i,sizeof(int),1,fp) != 1){
	printf("INT Cannot read Ydim from file %s\n",str);
      }
      if(i!=yd) printf("INT WARNING Y dimension missmatch %i %i\n",yd,i);

      if(piSmartContoursNExternal){
	free(piSmartContoursNExternal);
	piSmartContoursNExternal=NULL;
      }
      if(ppPAIRSmartContoursExternal){
	for(i=0;i<iNDepthLinesExternal;i++){
	  if(ppPAIRSmartContoursExternal[i]) free(ppPAIRSmartContoursExternal[i]);
	}
	if(ppPAIRSmartContoursExternal){
	  free(ppPAIRSmartContoursExternal);
	  ppPAIRSmartContoursExternal=NULL;
	}
      }

      if(fread(&iNDepthLinesExternal,sizeof(int),1,fp) != 1){
	printf("INT Cannot read number of contours from file %s\n",str);
      }
      if(!(piSmartContoursNExternal=(int*)calloc(iNDepthLinesExternal,sizeof(int)))){
	printf("INT Cannot allocate for piSmartContoursNExternal\n");
	fclose(fp);
	continue;
      }
      if(!(ppPAIRSmartContoursExternal=(PAIR**)calloc(iNDepthLinesExternal,sizeof(PAIR*)))){
	printf("INT Cannot allocate for ppPAIRSmartContoursExternal\n");
	fclose(fp);
	continue;
      }

      if(fread((void*)piSmartContoursNExternal,iNDepthLinesExternal*sizeof(int),1,fp) != 1){
	printf("INT Cannot read numbers of contour segments from file %s\n",str);
      }
      for(i=0;i<iNDepthLinesExternal;i++){
	if(piSmartContoursNExternal[i]){
	  if(!(ppPAIRSmartContoursExternal[i]=(PAIR*)calloc(piSmartContoursNExternal[i],sizeof(PAIR)))){
	    printf("INT Cannot allocate for ppPAIRSmartContoursExternal[%i]\n",i);
	    break;
	  }
	  if(fread(ppPAIRSmartContoursExternal[i],piSmartContoursNExternal[i]*sizeof(PAIR),1,fp) != 1){
	    printf("INT Cannot read contour segments(%i) from file %s\n",i,str);
	  }
	}
      }
      fclose(fp);

      printf("INT read %i contour lines\n",iNDepthLinesExternal);

      if(iPlotExternalContours) DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      continue;
    }

    if(ch=='W'){
      if(!depth_lines){
	printf("INT No contour lines to save (depth_lines=0)\n");
	continue;
      }
      if(!piSmartContoursN || !ppPAIRSmartContours){
	printf("INT No contour lines to save\n");
	continue;
      }
      for(i=0;i<depth_lines;i++) if(piSmartContoursN[i]) break;
      if(i==depth_lines){
	printf("INT Empty contour lines\n");
	continue;
      }
      if(contours_iniangle!=0.0)
	sprintf(str,"%s_contours%i_%i_%i_%i",filenames[0],depth_lines,
		(int)(smartcontours_threshold*100.0),(int)(dsradius),(int)(contours_iniangle));
      else
	sprintf(str,"%s_contours%i_%i_%i",filenames[0],depth_lines,
		(int)(smartcontours_threshold*100.0),(int)(dsradius));
      printf("INT Saving curent contours in %s\n",str);
      if((fp=fopen(str,"w"))==NULL){
	printf("INT Cannot open file %s\n",str);
      }
      i=0;
      if(fwrite((void*)&i,sizeof(int),1,fp) != 1){
	printf("INT Cannot write type %i to file %s\n",i,str);
      }
      if(fwrite((void*)&xd,sizeof(int),1,fp) != 1){
	printf("INT Cannot write Xdim=%i to file %s\n",xd,str);
      }
      if(fwrite((void*)&yd,sizeof(int),1,fp) != 1){
	printf("INT Cannot write Ydim=%i to file %s\n",yd,str);
      }
      if(fwrite((void*)&depth_lines,sizeof(int),1,fp) != 1){
	printf("INT Cannot write number of contours to file %s\n",str);
      }
      if(fwrite((void*)piSmartContoursN,depth_lines*sizeof(int),1,fp) != 1){
	printf("INT Cannot write numbers of contour segments to file %s\n",str);
      }
      for(i=0;i<depth_lines;i++){
	if(piSmartContoursN[i]){
	  if((iError=fwrite(ppPAIRSmartContours[i],piSmartContoursN[i]*sizeof(PAIR),1,fp)) != 1){
	    printf("INT Cannot write contour segments(%i) to file %s (%i)\n",i,str,iError);
	    break;
	  }
	}
      }
      fclose(fp);
      continue;
    }

    if(ch=='T'){
      if(!iNDepthLinesExternal){
	printf("INT No contour lines to transfer (iNDepthLinesExternal=0)\n");
	continue;
      }
      if(!piSmartContoursNExternal || !ppPAIRSmartContoursExternal){
	printf("INT No contour lines to transfer\n");
	continue;
      }
      for(i=0;i<iNDepthLinesExternal;i++) if(piSmartContoursNExternal[i]) break;
      if(i==iNDepthLinesExternal){
	printf("INT Empty external contour lines\n");
	continue;
      }
      
      if(piSmartContoursN){
	free(piSmartContoursN);
	piSmartContoursN=NULL;
      }
      if(ppPAIRSmartContours){
	for(i=0;i<depth_lines;i++){
	  if(ppPAIRSmartContours[i]) free(ppPAIRSmartContours[i]);
	}
	free(ppPAIRSmartContours);
	ppPAIRSmartContours=NULL;
      }
      
      depth_lines=iNDepthLinesExternal;

      if(!(contours=(float*)realloc(contours,(depth_lines+1)*sizeof(float)))){
	printf("INT Cannot reallocate for contours\n");
	contours=(float*)0;
	depth_lines=0;
	continue;
      }

      if(!(piSmartContoursN=(int*)calloc(depth_lines,sizeof(int)))){
	printf("INT Cannot allocate for piSmartContoursN\n");
	continue;
      }
      if(!(ppPAIRSmartContours=(PAIR**)calloc(depth_lines,sizeof(PAIR*)))){
	printf("INT Cannot allocate for ppPAIRSmartContours\n");
	continue;
      }
      
      memcpy(piSmartContoursN,piSmartContoursNExternal,iNDepthLinesExternal*sizeof(int));
      
      for(i=0;i<depth_lines;i++){
	if(piSmartContoursN[i]){
	  if(!(ppPAIRSmartContours[i]=(PAIR*)calloc(piSmartContoursN[i],sizeof(PAIR)))){
	    printf("INT Cannot allocate for ppPAIRSmartContours[%i]\n",i);
	    break;
	  }
	  memcpy(ppPAIRSmartContours[i],ppPAIRSmartContoursExternal[i],piSmartContoursN[i]*sizeof(PAIR));
	}
      }
      
      printf("INT transferred %i contour lines\n",depth_lines);
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
    }

    if(ch=='B'){
      if(GetStringFromSTDIN(str,255,"INT Enter batch file: ")) continue;

      if((fp=fopen(str,"r"))==NULL){
	printf("INT Cannot open file <%s>\n",str);
	continue;
      }

      i=0;
      while(fgets(strDum,255,fp)){
	if(*strDum=='#') continue;
	pc=strtok(strDum," \t");
	if(pc) xi=atoi(pc);
	else{
	  break;
	}
	if(xi<0 || xi>xd-1){
	  printf("INT Batch: x(%i) out of range(0-%i)\n",xi,xd-1);
	  continue;
	}
	pc=strtok(NULL," \t");
	if(pc) yi=atoi(pc);
	else{
	  break;
	}
	if(yi<0 || yi>yd-1){
	  printf("INT Batch: y(%i) out of range(0-%i)\n",yi,yd-1);
	  continue;
	}
	pc=strtok(NULL," \t");
	i++;
	//	printf("INT X %i Y %i\n",xi,yi);
	l=xi+yi*xd;
	rho=hypot(mp[0][l],mp[1][l]);
	phi=180.0/M_PI*atan2(mp[1][l],mp[0][l]);

	if(iDoWedgeAnnotation){
	  dumd=phi/360.0;
	  if(dumd<0.0) dumd+=1.0;

	  if(iPrePostStimTimeRescale){
	    dumd=(dumd-dPreStimTime)/(1.0-dPreStimTime-dPostStimTime);
	  }

	  dumd=dumd*(fWedgeMax-fWedgeMin)+fWedgeMin;
	  dumd=exp(dumd*log(2.0));
	  //	  if(pc) printf("INI X %i Y %i F %.3f %s\n",xi,yi,dumd,pc);
	  //	  else printf("INI X %i Y %i F %.3f\n",xi,yi,dumd);
	  //	  if(pc) printf("#%i\n%.3f %s\n",i,dumd,pc);
	  //if(pc) printf("%.3f\n",dumd);
	  //else printf("%.3f\n",dumd);
	  switch(iBatchOutput){
	  case 1:
	    printf("%.3f\t%.4f\n",dumd,rho);
	    break;
	  case 0:
	  default:
	    printf("%.3f\n",dumd);
	    break;
	  }
	}
	else{
	  printf("X %i Y %i Rho %f Phi %f\n",xi,yi,rho,phi);
	}
      }
      printf("INT Processed %i items\n",i);
      fclose(fp);
      continue;
    }

    if(ch=='M'){
      iMarkPoints=1;
      if((iReturn=GetFloatPoint(&x,&y,7,3))>0){
	if(!(pMarkPoints=(MARKRECORD*)realloc(pMarkPoints,(iMarkPointsN+1)*sizeof(MARKRECORD)))){
	  printf("INI Cannot realloc for pMarkPoints\n");
	  continue;
	}
	pMarkPoints[iMarkPointsN].x=x;
	pMarkPoints[iMarkPointsN++].y=y;
	printf("INT Mark %0.1f %0.1f\n",x,y);
	if(!iMarkPointsHide)  DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
	else continue;
      }
      else{
	if(!iReturn) continue;
	if(-iReturn=='M'){
	  if(GetStringFromSTDIN(str,255,"INT Enter mark batch file: ")) continue;
	  
	  if((fp=fopen(str,"r"))==NULL){
	    printf("INT Cannot open file <%s>\n",str);
	    continue;
	  }
	  iMarkPointsLoadN++;
	  
	  i=0;
	  while(fgets(strDum,255,fp)){
	    if(*strDum=='#') continue;
	    pc=strtok(strDum," \t");
	    if(pc) x=atof(pc);
	    else{
	      break;
	    }
	    if(x<0 || x>xd-1){
	      printf("INT Batch: x(%0.1f) out of range(0-%i)\n",x,xd-1);
	      continue;
	    }
	    pc=strtok(NULL," \t");
	    if(pc) y=atof(pc);
	    else{
	      break;
	    }
	    if(y<0 || y>yd-1){
	      printf("INT Batch: y(%0.1f) out of range(0-%i)\n",y,yd-1);
	      continue;
	    }
	    pc=strtok(NULL," \t");
	    i++;
	    //	printf("INT X %i Y %i\n",xi,yi);
	    
	    if(!(pMarkPoints=(MARKRECORD*)realloc(pMarkPoints,(iMarkPointsN+1)*sizeof(MARKRECORD)))){
	      printf("INI Cannot realloc for pMarkPoints\n");
	      continue;
	    }
	    pMarkPoints[iMarkPointsN].iIndex=iMarkPointsN+1;
	    pMarkPoints[iMarkPointsN].iSymbol=iMarkPointsLoadN;
	    pMarkPoints[iMarkPointsN].iSelect=0;
	    pMarkPoints[iMarkPointsN].x  =x;
	    pMarkPoints[iMarkPointsN++].y=y;
	  }
	  printf("INT Processed %i items\n",i);
	  fclose(fp);
	  if(!iMarkPointsHide)  DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
	  else continue;
	}
	if(-iReturn=='S'){
	  if(GetStringFromSTDIN(str,255,"INT Enter select marked file: ")) continue;
	  
	  if((fp=fopen(str,"r"))==NULL){
	    printf("INT Cannot open file <%s>\n",str);
	    continue;
	  }
	  iMarkPointsSelectN++;
	  
	  i=0;
	  while(fgets(strDum,255,fp)){
	    if(*strDum=='#') continue;
	    if((pc=strtok(strDum," \t"))) xi=atoi(pc);
	    else break;

	    if(xi<1 || xi>iMarkPointsN){
	      printf("INT Select: %i out of range(1-%i)\n",xi,iMarkPointsN);
	      continue;
	    }
	    i++;
	    //	printf("INT X %i Y %i\n",xi,yi);
	    if(pMarkPoints[xi-1].iSelect) printf("INT WARNING mark point %i has been selected already\n",xi);
	    switch(iMarkPointsSelectMode){
	    case MARK_SELECTION_MODE_EXCLUDE:
	      pMarkPoints[xi-1].iSelect=-1;
	      break;
	    case MARK_SELECTION_MODE_SPECIAL:
	      pMarkPoints[xi-1].iSelect=iMarkPointsSelectN;
	      break;
	    case MARK_SELECTION_MODE_INCLUDE:
	    default:
	      pMarkPoints[xi-1].iSelect=0;
	      break;
	    }
	  }
	  printf("INT Processed %i items\n",i);
	  fclose(fp);
	  if(!iMarkPointsHide)  DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
	  else continue;
	}
      }
    }

    if(ch=='Z'){
      if(GetStringFromSTDIN(str,255,"INT Enter correlation file[s]: ")) continue;

      pc=strtok(str," ");
      if((pc=strtok(NULL," "))){
	k=2;
	if((pF=fopen(pc,"r"))==NULL){
	  printf("INT Cannot open file 2 <%s>\n",pc);
	  continue;
	}
      }
      else k=1;
      if((fp=fopen(str,"r"))==NULL){
	printf("INT Cannot open file 1 <%s>\n",str);
	continue;
      }
      
      iNCorrRecords=0;
      SAFE_FREE(pCorrRecords);
      while(!GetStringFromFile(fp,strDum,255,1)){
	pc=strtok(strDum," \t");
	if(pc) xi=(int)(dX=atof(pc));
	else{
	  printf("INT Incomplete X record in %s @ %i\n",str,i);
	  break;
	}
	if(xi<0 || xi>xd-1){
	  printf("INT Correlation: x(%i) out of range(0-%i)\n",xi,xd-1);
	  continue;
	}

	pc=strtok(NULL," \t");
	if(pc) yi=(int)(dY=atof(pc));
	else{
	  printf("INT Incomplete Y record in %s @ %i\n",str,i);
	  break;
	}
	if(yi<0 || yi>yd-1){
	  printf("INT Correlation: y(%i) out of range(0-%i)\n",yi,yd-1);
	  continue;
	}
	if(k==1){
	  pc=strtok(NULL," \t");
	  if(pc) dZ=atof(pc);
	  else{
	    printf("INT Incomplete Z record in %s @ %i\n",str,i);
	    break;
	  }
	  pc=strtok(NULL," \t");
	}
	else{
	  if(!GetStringFromFile(pF,strDum,255,1)){
	    pc=strtok(strDum," \t");
	    if(pc) dZ=atof(pc);
	    else{
	      printf("INT File 2 is shorter than 1\n");
	      break;
	    }
	    pc=NULL;
	  }
	  else{
	    printf("INT File 2 is shorter than 1\n");
	    break;    
	  }
	}
	
	if(!(pCorrRecords=(RECORD*)realloc(pCorrRecords,(iNCorrRecords+1)*sizeof(RECORD)))){
	  printf("INI Cannot realloc for pCorrRecords\n");
	  continue;
	}

	pR=&(pCorrRecords[iNCorrRecords++]);
	if(pc) pR->iIndex=atoi(pc);
	else pR->iIndex=iNCorrRecords;

	switch(iCorrelationTest){
	case CORRELATION_TEST_NONE:
	  break;
	case CORRELATION_TEST_SHIFT:
	  dX+=iCorrelationTestShiftX;
	  if(dX<0.0){
	    printf("INT Clipping X=%f below %i\n",dX,iNCorrRecords);
	    dX=0.0;
	  }
	  if(dX>=xd){
	    printf("INT Clipping X=%f above %i\n",dX,iNCorrRecords);
	    dX=xd-1;
	  }
	  dY+=iCorrelationTestShiftY;
	  if(dY<0.0){
	    printf("INT Clipping Y=%f below %i\n",dY,iNCorrRecords);
	    dY=0.0;
	  }
	  if(dY>=yd){
	    printf("INT Clipping Y=%f above %i\n",dY,iNCorrRecords);
	    dY=yd-1;
	  }
	  break;
	case CORRELATION_TEST_RANDOM:
	  //dX=iCorrelationTestRandomX+iCorrelationTestRandomXDim*ran3(&lSeed);
	  //dY=iCorrelationTestRandomY+iCorrelationTestRandomYDim*ran3(&lSeed);
	  do{
	    dX=(xd-1)*ran3(&lSeed);
	    dY=(yd-1)*ran3(&lSeed);
	  }while(!IsPointInsideDomain(dX,dY,iNDomainPoints,pDomain,&stRect));

	  break;
	}

	pR->dX=dX;
	pR->dY=dY;
	pR->dZ=dZ;
	pR->iHide=0;
      }
      printf("INT Processed %i records\n",iNCorrRecords);
      if(fp){
	fclose(fp);
	fp=NULL;
      }
      if(pF){
	fclose(pF);
	pF=NULL;
      }

      if(iMarkCorrelationPoints){
	iMarkPoints=1;
	if(!iMarkPointsHide)  DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      }
      CorrelationToExternalSource(iNCorrRecords,pCorrRecords,xd,yd,mp,0);
      continue;
    }
      
    if(ch=='I'){
      if(GetStringFromSTDIN(str,255,"INT Enter pinwheel file: ")) continue;

      if((fp=fopen(str,"r"))==NULL){
	printf("INT Cannot open file <%s>\n",str);
	continue;
      }
      iError=0;
      if(fread(&iDum,sizeof(int),1,fp) != 1){
	printf("INT Cannot read # pinwheels from file %s\n",str);
	iError|=1;
      }
      if(fread(&i,sizeof(int),1,fp) != 1){
	printf("INT Cannot read Xdim from file %s\n",str);
	iError|=2;
      }
      if(i!=xd){
	printf("INT WARNING X dimension missmatch %i %i\n",xd,i);
	iError|=4;
      }
      if(fread(&i,sizeof(int),1,fp) != 1){
	printf("INT Cannot read Ydim from file %s\n",str);
	iError|=8;
      }
      if(i!=yd){
	printf("INT WARNING Y dimension missmatch %i %i\n",yd,i);
	iError|=16;
      }
      if(iError){
	fclose(fp);
	continue;
      }
      else{
	apinwheels_n=iDum;
	SAFE_FREE(apinwheels);
	if(!(apinwheels=(Pinwheel*)calloc(apinwheels_n,sizeof(Pinwheel)))){
	printf("INT Cannot allocate for apinwheels\n");
	apinwheels_n=0;
	fclose(fp);
	continue;
      }

      if(fread(apinwheels,sizeof(Pinwheel)*apinwheels_n,1,fp) != 1){
	printf("INT Cannot read pinwheels from file %s\n",str);
	fclose(fp);
	continue;
      }
      
      printf("INT read pinwheels %i\n",apinwheels_n);
      fclose(fp);

      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      continue;
      }
    }
  
    if(ch=='J'){
      if(!LoadStencilBitmap(xd,yd,&stenciln,&stencil,iLoadBitmapMode)){
	DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      }
      continue;
    }

    if(ch=='K'){
      SaveStencilBitmap(xd,yd,stenciln,stencil);
      continue;
    }

    if(ch=='G'){
      if(!ReadDomain(&iNDomainPoints,&pDomain)){
	Stencil(xd,yd,iNDomainPoints,pDomain,&stenciln,&stencil);	    
	DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
	if(do_fields){ 
	  DisplayFields(field_mode,Xdim,Ydim,maps[0],maps[1],dbuffers[0],dbuffers[1]);
	}
      }
      continue;
    }

    if(ch=='H'){
      WriteDomain(iNDomainPoints,pDomain);
      continue;
    }

    if(ch=='F'){
      stFilter.iKind=FILTER_KIND_COS;
      stFilter.dRadius=25.5;
      stFilter.dDepth=0.2;
      if(SetFilter(&stFilter,xd,yd)){
	printf("INT Cannot SetFilter\n");
	continue;
      }

      while(1){
	cpgband(0,1,x,y,&x,&y,&ch);
	if(ch=='q' || ch=='Q') break;
	  dDum=0.05;//noise level
 	  x0=(int)x;
	  y0=(int)y;
	  for(j=0;j<stFilter.iNPoints;j++){
	    xi=x0+stFilter.piX[j];
	    yi=y0+stFilter.piY[j];
	    if(xi>=0 && xi<xd && yi>=0 && yi<yd){ 
	      mp[0][xi+yi*xd]*=(1-stFilter.pdValues[j]*ran3(&lSeed));
	      mp[1][xi+yi*xd]*=(1-stFilter.pdValues[j]*ran3(&lSeed));
	      //Noisify
	      mp[0][xi+yi*xd]+=dDum*(ran3(&lSeed)-0.5);
	      mp[1][xi+yi*xd]+=dDum*(ran3(&lSeed)-0.5);
	    }
	  }
	  //Shuffle
	  iShuffle=0;
	  if(iShuffle){
	    for(j=0;j<stFilter.iNPoints;j++){
	      i=(int)floor(stFilter.iNPoints*ran3(&lSeed));
	      k=(int)floor(stFilter.iNPoints*ran3(&lSeed));
	      if(INSIDE(x0+stFilter.piX[i],0,xd) && INSIDE(y0+stFilter.piY[i],0,yd) &&
		 INSIDE(x0+stFilter.piX[k],0,xd) && INSIDE(y0+stFilter.piY[k],0,yd)){
		dDum=mp[0][x0+stFilter.piX[i]+(y0+stFilter.piY[i])*xd];
		mp[0][x0+stFilter.piX[i]+(y0+stFilter.piY[i])*xd]=
		  mp[0][x0+stFilter.piX[k]+(y0+stFilter.piY[k])*xd];
		mp[0][x0+stFilter.piX[k]+(y0+stFilter.piY[k])*xd]=dDum;
		dDum=mp[1][x0+stFilter.piX[i]+(y0+stFilter.piY[i])*xd];
		mp[1][x0+stFilter.piX[i]+(y0+stFilter.piY[i])*xd]=
		  mp[1][x0+stFilter.piX[k]+(y0+stFilter.piY[k])*xd];
		mp[1][x0+stFilter.piX[k]+(y0+stFilter.piY[k])*xd]=dDum;
	      }
	    }
	  }
	  SwitchPanel(xd,yd,fbufferi,mp,mode,&fIMin,&fIMax,iPanelMode);
      }
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      continue;
    }
    
    if(ch=='P'){
      xi=(int)rint(x);
      yi=(int)rint(y);
      for(i=xi-iAddPixelsRadius;i<=xi+iAddPixelsRadius;i++){
	for(k=yi-iAddPixelsRadius;k<=yi+iAddPixelsRadius;k++){
	  AddPixelToStencil(i,k,xd,yd,&stenciln,&stencil,0);
	}
      }
      SwitchPanel(xd,yd,fbufferi,mp,mode,&fIMin,&fIMax,iPanelMode);
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      continue;
    }

    if(ch=='o'){

      // Keys in use: A,B,C,E,G,H,M,N,P,R,S,W,Z,a,b,c,f,h,l,m,n,p,r,s,w,z
      cpgband(7,1,x,y,&x,&y,&ch); 

      if(ch=='h'){
	printf("  A - use RhoMin on section amplitude plot %i\n",iUseRhoMinOnSectionAmplitudePlot);
	printf("  B - batch output mode %i\n",iBatchOutput);
	printf("  C - contour type switch, current value %i\n",do_smartcontours);
	printf("  E - enumarate marked switch, current value %i\n",iMarkPointsEnumerate);
	printf("  G - batch select marked points mode %i\n",iMarkPointsSelectMode);
	printf("  H - hide marked external points %i\n",iMarkPointsHide);
	printf("  M - mark external points %i\n",iMarkPoints);
	printf("  N - save noise level, current value %i\n",iGetNoiseLevel);
	printf("  P - do not plot phase, current value %i\n",iDoNotPlotPhase);
	printf("  R - do shift for correlation, current value %i\n",iCorrelationModeShift);
	printf("  S - save crop border as lines switch, current value %i\n",iSaveCropBorderAsLineSegments);
	printf("  W - plot reduced color wedge %i\n",iReducedColorWedge);
	printf("  Z - correlation Rho threshold %f\n",dCorrelationThreshold);
	printf("  a - additive crop switch, current value %i\n",iAdditiveCropChanges);
	printf("  b - plot black-white wedge switch, current value %i\n",iDrawBWWedge);
	printf("  c - plot external contours switch, current value %i\n",iPlotExternalContours);
	printf("  e - emphasize first contour switch, current value %i\n",emphasize_contours_iniangle);
	printf("  f - find marked inside selected region switch, current value %i\n",iFindMarkedInsideSelected);
	printf("  h - print this help\n");
	printf("  l - plot large werge, current value %i\n",iDoLargeWedge);
	printf("  m - mark min/max locations %i\n",iMarkMinMax);
	printf("  n - null blank responses %i\n",iNullOutsiders);
	printf("  p - panel stencil switch, current value %i\n",iPanelMode);
	printf("  r - replot contours switch, current value %i\n",iReplotSavedContours);
	printf("  s - use single pixel value for section plot, current value %i\n",iSectionSinglePixel);
	printf("  w - plot color wedge switch, current value %i\n",iDrawColorWedge);
	printf("  z - crop Zeros on saves, current value %i\n",iCropZerosOnSave);
	continue;
      }

      if(ch=='A'){ 
	iUseRhoMinOnSectionAmplitudePlot=1-iUseRhoMinOnSectionAmplitudePlot;
	continue;
      }

      if(ch=='B'){
	if(GetString(str)) iBatchOutput=atoi(str);
	continue;
      }

      if(ch=='C'){
	do_smartcontours=1-do_smartcontours;
	//      smartcontours_threshold=DEFAULT_SMART_CONTOURS_THRESHOLD;
	if(depth_lines) DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
	continue;      
      }
      
      if(ch=='E') iMarkPointsEnumerate=1-iMarkPointsEnumerate;

      if(ch=='G'){
	iMarkPointsSelectMode=(iMarkPointsSelectMode+1)%MARK_SELECTION_MODE_N;
	continue;
      }

      if(ch=='H') iMarkPointsHide=1-iMarkPointsHide;

      if(ch=='M'){ 
	if(iMarkPoints){
	  iMarkPoints=0;
	  iMarkPointsHide=0;
	  SAFE_FREE(pMarkPoints);
	  iMarkPointsN=0;
	}
	else{
	  iMarkPoints=1;
	  iMarkPointsHide=0;
	}
      }

      if(ch=='N') iGetNoiseLevel=1;
      if(ch=='P') iDoNotPlotPhase=1-iDoNotPlotPhase;
      if(ch=='R'){
	iCorrelationModeShift=1-iCorrelationModeShift;
	continue;
      }
       if(ch=='S'){
	iSaveCropBorderAsLineSegments=1-iSaveCropBorderAsLineSegments;
	continue;
      }
     if(ch=='W') iReducedColorWedge=1-iReducedColorWedge;

      if(ch=='Z'){
	if(GetString(str) && atof(str)>=0.0) dCorrelationThreshold=atof(str);
	CorrelationToExternalSource(iNCorrRecords,pCorrRecords,xd,yd,mp,0);
	continue;
      }
     
      if(ch=='a'){ 
	iAdditiveCropChanges=1-iAdditiveCropChanges;
	continue;
      }
      if(ch=='b') iDrawBWWedge=1-iDrawBWWedge;
      if(ch=='c') iPlotExternalContours=1-iPlotExternalContours;
      if(ch=='e') emphasize_contours_iniangle=1-emphasize_contours_iniangle;
      if(ch=='f'){ 
	iFindMarkedInsideSelected=1-iFindMarkedInsideSelected;
	continue;
      }
      if(ch=='l'){ 
	iDoLargeWedge=1-iDoLargeWedge;
	//	if(!iDrawColorWedge) continue;
      }
      if(ch=='m'){ 
	iMarkMinMax=1-iMarkMinMax;
      }
      if(ch=='n'){ 
	iNullOutsiders=1-iNullOutsiders;
      }
      if(ch=='p'){ 
	if(iPanelMode&SWITCH_PANEL_MODE_STENCIL) iPanelMode&=~SWITCH_PANEL_MODE_STENCIL;
	else iPanelMode|=SWITCH_PANEL_MODE_STENCIL;
	SwitchPanel(xd,yd,fbufferi,mp,mode,&fIMin,&fIMax,iPanelMode);
	continue;      
      }
      if(ch=='r'){
	iReplotSavedContours=1-iReplotSavedContours;
	continue;      
      }
      if(ch=='s'){
	iSectionSinglePixel=1-iSectionSinglePixel;
	continue;      
      }

      if(ch=='w') iDrawColorWedge=1-iDrawColorWedge;

      if(ch=='z'){ 
	iCropZerosOnSave=1-iCropZerosOnSave;
	continue;
      }

      
      DisplayMap(display_mode,xd,yd,mp[0],mp[1],dbuffers[0],dbuffers[1]);
      continue;      
    }
    
    if(ch=='h'){
      printf("0 - set original image\n");
      printf("B - read in batch file and output stuff out\n");
      printf("D - plot distribution\n");
      printf("G - read domain segments from file\n");
      printf("H - write domain segments to file\n");
      printf("O - color scheme switch\n");
      printf("Q - quit program\n");
      printf("R - read contours from file\n");
      printf("S - save current maps\n");
      printf("T - transfer contours from external to internal source\n");
      printf("W - write contours to file\n");
      printf("Z - read in correlation file\n");
      printf("a - <float> angle shift (default 0.0)\n");
      printf("b - spawn bias interactive panel (default off)\n");
      printf("c - <mouse input> outline domain of interest (default all image)\n");
      printf("d - <int> number of contour lines (default 0)\n");
      printf("g - <mouse input(2 points)> 2D plot along a section\n");
      printf("h - print this help\n");
      printf("i - switch, invert images (x <-> -x) (default x)\n");
      printf("l - <float>pcnt lower chop for amplitude (default 0.0)\n");
      printf("m - <int> <int> manual input of coordinates\n");
      printf("n - normalize on selected domain\n");
      printf("o - other options\n");
      printf("p - switch, phi <-> rho image on interactive panel (defined by options)\n");
      printf("q - quit interactive dialog\n");
      printf("r - switch, rescale <-> do not rescale selected area (default not)\n");
      printf("s - <int> contour line scheme (default 0)\n");
      printf("t - <float>pcnt threshold for smart contours (default %f)\n",DEFAULT_SMART_CONTOURS_THRESHOLD);
      printf("u - <float>pcnt upper chop for amplitude (default 100.0)\n");
      printf("v - switch, verbose mode (default off)\n");
     continue;
    }

    xi=(int)rint(x);
    if(xi<0) xi=0;
    if(xi>xd-1) xi=xd-1;
    yi=(int)rint(y);
    if(yi<0) yi=0;
    if(yi>yd-1) yi=yd-1;
    printf("INT X %i Y %i\n",xi,yi);
    
    l=xi+yi*xd;
    rho=hypot(mp[0][l],mp[1][l]);
    phi=180.0/M_PI*atan2(mp[1][l],mp[0][l]);

    printf("INT Rho=%f(%f %f) Phi=%f(%f)\n",rho,rho-rhop,rho/rhop,phi,phi-phip);
    if(iDoWedgeAnnotation){
      dumd=phi/360.0;
      if(dumd<0.0) dumd+=1.0;
      if(iPrePostStimTimeRescale){
	dumd=(dumd-dPreStimTime)/(1.0-dPreStimTime-dPostStimTime);
      }
      dumd=dumd*(fWedgeMax-fWedgeMin)+fWedgeMin;
      dumd=exp(dumd*log(2.0));
      printf("INI Frequency %.3f\n",dumd);
    }

    if(do_phase_transform){
      printf("INT Transformed Phi=%f\n",phi*phase_transform*M_PI/180.0);
    }

    rhop=rho;
    phip=phi;
    
  looplabel:
    ;
  }
  cpgslct(interactiveid);
  cpgclos();
  cpgslct(callerid);
}

void Biaspanel(int id,int type,float fmax){
  int bias_paneln;
  float *bias_panel,bias_trans[6];
  int i,j,l;
  double dumd;

  printf("Panel type %i\n",type);
  setenv("PGPLOT_ENVOPT","I",1); 
  cpgslct(id);
  cpgsci(1);

  if(type==BIAS_PANEL_POLAR){
    cpgpap(BIAS_WINDOW_SIZE,1.0);
    cpgsch(0.8);
    cpgenv(-fmax,fmax,-fmax,fmax,1,0);
    bias_paneln=BIAS_PANEL_N;
    if((bias_panel=(float*)malloc((size_t)bias_paneln*(size_t)bias_paneln*sizeof(float)))){
      bias_trans[0]=-fmax*(1.0+1.0/(float)(bias_paneln/2));
      bias_trans[1]=fmax/(float)(bias_paneln/2);
      bias_trans[2]=0.0;
      bias_trans[3]=-fmax*(1.0+1.0/(float)(bias_paneln/2));
      bias_trans[4]=0.0;
      bias_trans[5]=fmax/(float)(bias_paneln/2);
      dumd=(double)fmax/(double)(bias_paneln/2);
      for(j=0;j<bias_paneln;j++){
	for(i=0;i<bias_paneln;i++){
	  l=i+j*bias_paneln;
	  bias_panel[l]=(float)atan2(dumd*(j-bias_paneln/2),dumd*(i-bias_paneln/2));
	}
      }
      Pallet(DEFAULT_PALLET,1.0,0.5);
      cpgimag(bias_panel,bias_paneln,bias_paneln,1,bias_paneln,1,bias_paneln,-M_PI,M_PI,bias_trans);
      free(bias_panel);
    }
    else{
      printf("INT Cannot allocate for bias_panel\n");
    }
    Bullseye(fmax,BIAS_RING_N,BIAS_RAY_N);
  }
  if(type==BIAS_PANEL_CARTESIAN){
    cpgpap(BIAS_WINDOW_SIZE,0.5*(sqrt(5.0)-1.0));
    cpgsch(0.8);
    cpgenv(-180.0,180.0,0.0,fmax,0,0);
    bias_paneln=BIAS_PANEL_N;
    if((bias_panel=(float*)malloc((size_t)bias_paneln*sizeof(float)))){
      bias_trans[0]=-180.0*(1.0+1.0/(float)(bias_paneln/2));
      bias_trans[1]=180.0/(float)(bias_paneln/2);
      bias_trans[2]=0.0;
      bias_trans[3]=-0.5*fmax;
      bias_trans[4]=0.0;
      bias_trans[5]=fmax;

      for(i=0;i<bias_paneln;i++){
	bias_panel[i]=180.0*(i-bias_paneln/2)/(float)(bias_paneln/2);
      }
      Pallet(DEFAULT_PALLET,2.0,0.5);
      cpgimag(bias_panel,bias_paneln,1,1,bias_paneln,1,1,-180.0,180.0,bias_trans);
      free(bias_panel);
    }
    else{
      printf("INT Cannot allocate for bias_panel\n");
    }
    Mesh(-180.0,180.0,0.0,fmax,BIAS_RAY_N,BIAS_RING_N);    
  }

  unsetenv("PGPLOT_ENVOPT");  
}

void Bullseye(float rmax,int rdiv,int phidiv){
  int i;
  float x,y,r;
  double a;
 
  cpgsfs(2);
  cpgsci(0);
  cpgsls(2);
  r=rmax/rdiv;
  //  r=1.0;
  for(i=0;i<rdiv;i++){
    cpgcirc(0.0,0.0,r*(i+1));
  }
  a=M_PI/phidiv;
  r=rmax*sqrt(2.0);
  for(i=0;i<phidiv;i++){
    x=r*(float)cos(a*i);
    y=r*(float)sin(a*i);
    cpgmove(x,y);
    cpgdraw(-x,-y);
  }
  cpgsls(1);
}

void Mesh(float xmin,float xmax,float ymin,float ymax,int xdiv,int ydiv){
  int i;
  float x,y;
 
  cpgsfs(2);
  cpgsci(0);
  cpgsls(2);
  x=0.5*(xmax-xmin)/xdiv;
  for(i=0;i<xdiv;i++){
    cpgmove(x*i,ymin);
    cpgdraw(x*i,ymax);
    cpgmove(-x*i,ymin);
    cpgdraw(-x*i,ymax);
  }
  y=(ymax-ymin)/ydiv;
  for(i=1;i<ydiv;i++){
    cpgmove(xmin,ymin+y*i);
    cpgdraw(xmax,ymin+y*i);
  }
  cpgsls(1);

}

void SwitchPanel(int xd,int yd,float *buf,double **mp,int *mode,float *pfMin,float *pfMax,int iMode){
  unsigned int xy,i,j,lmax,lmin;
  float dumf;
  double *pdX,*pdY;

  xy=(unsigned int)xd*(unsigned int)yd;
  pdX=mp[0];
  pdY=mp[1];

  if(iMode&SWITCH_PANEL_MODE_SWITCH){
    switch(*mode){
    case INT_MODE_PHI:
      *mode=INT_MODE_RHO;
      break;
    case INT_MODE_RHO:
    default:
      *mode=INT_MODE_PHI;
    }
  }

  switch(*mode){
  case INT_MODE_RHO:
    *pfMax=-(*pfMin=BIG_FLOAT);
    if(iMode&SWITCH_PANEL_MODE_STENCIL){
      memset(buf,0,sizeof(float)*xy);
      for(j=0;j<stenciln;j++){
	i=stencil[j];
	buf[i]=dumf=(float)hypot(pdX[i],pdY[i]);
	if(dumf>*pfMax){ 
	  *pfMax=dumf;
	  lmax=i;
	}
	if(dumf<*pfMin){
	  *pfMin=dumf;
	  lmin=i;
	}
      }
    }
    else{
      for(i=0;i<xy;i++){
	buf[i]=dumf=(float)hypot(pdX[i],pdY[i]);
	if(i>=xd){ // Do not take into account first screwed row of each frame
	  if(dumf>*pfMax){ 
	    *pfMax=dumf;
	    lmax=i;
	  }
	  if(dumf<*pfMin){
	    *pfMin=dumf;
	    lmin=i;
	  }
	}
      }
    }

    if(*pfMin==*pfMax){
      *pfMin-=1.0;
      *pfMax+=1.0;
    }
    break;
  case INT_MODE_PHI:
  default:
    if(iMode&SWITCH_PANEL_MODE_STENCIL){
      memset(buf,0,sizeof(float)*xy);
      for(j=0;j<stenciln;j++){
	i=stencil[j];
	buf[i]=(float)atan2(pdY[i],pdX[i]);
      }
    }
    else{
      for(i=0;i<xy;i++){
	buf[i]=(float)atan2(pdY[i],pdX[i]);
      }
    }
    Pallet(DEFAULT_PALLET,1.0,0.5);
  }
  RestorePanel(xd,yd,buf,mode,*pfMin,*pfMax);
}


void RestorePanel(int xd,int yd,float *buf,int *mode,float fmin,float fmax){
  int i,j;
  PAIR *pPAIR;
  float fminc;

  if(*mode==INT_MODE_RHO){
    if(fmin==fmax){
      fmin-=1.0;
      fmax+=1.0;
    }
    cpggray(buf,xd,yd,1,xd,1,yd,fmax,fmin,trans);

    fminc=fmin+smartcontours_threshold*(fmax-fmin);
    for(i=0;i<depth_lines;i++){ 
      contours[i]=fminc+(float)(i)*(fmax-fminc)/(float)(depth_lines);
    }
    cpgsci(1);

    cpgcont(buf,xd,yd,1,xd,1,yd,contours,depth_lines,trans);

    /*
    for(i=0;i<depth_lines;i++){
      pPAIR=ppPAIRSmartContours[i];
      if(emphasize_contours_iniangle && !i) cpgsci(2);
      else cpgsci(3);
      for(j=0;j<piSmartContoursN[i];j++){
	cpgmove(pPAIR[j].x1,pPAIR[j].y1);
	cpgdraw(pPAIR[j].x2,pPAIR[j].y2);
      }
    }
    */
  }
  if(*mode==INT_MODE_PHI){
    cpgimag(buf,xd,yd,1,xd,1,yd,-M_PI,M_PI,trans);

    for(i=0;i<depth_lines;i++){
      pPAIR=ppPAIRSmartContours[i];
      if(emphasize_contours_iniangle && !i) cpgsci(1);
      else cpgsci(0);
      for(j=0;j<piSmartContoursN[i];j++){
	cpgmove(pPAIR[j].x1,pPAIR[j].y1);
	cpgdraw(pPAIR[j].x2,pPAIR[j].y2);
      }
    }

    if(iPlotExternalContours){
      for(i=0;i<iNDepthLinesExternal;i++){
	if(!(pPAIR=ppPAIRSmartContoursExternal[i])) continue;
	if(emphasize_contours_iniangle){
	  if(!i){
	    cpgsci(0);
	    cpgsls(1);
	  }
	  else{ 
	    cpgsci(1);
	    cpgsls(1);
	  }
	}
	else{
	  cpgsci(1);
	  cpgsls(1);
	}
	for(j=0;j<piSmartContoursNExternal[i];j++){
	  cpgmove(pPAIR[j].x1,pPAIR[j].y1);
	  cpgdraw(pPAIR[j].x2,pPAIR[j].y2);
	}
      }
    }

  }
  if(iPlotLineSegments) PlotLineSegments(0,1,iPlotLineSegmentsStyle);
  if(iMarkPoints && !iMarkPointsHide){
    if(*mode==INT_MODE_RHO) MarkPoints(xd,iMarkPointsColor+1,INT_MODE_RHO);
    else MarkPoints(xd,iMarkPointsColor,INT_MODE_PHI);
  }
}

char* GetString(char *str){
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
	  printf("%c",8);
	}
      }
      else{
	str[i++]=ch;
	printf("%c",ch);
      }
      fflush(stdout);
    }	
    else{ 
      if(ch=='q') return(NULL);
      str[i]='\0';
      printf("%c",ch);
      fflush(stdout);
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
  if(ch=='M') iReturn=-(int)'M';
  if(ch=='S') iReturn=-(int)'S';
  if(ch=='m'){
    if(GetString(str)){ 
      *pfX=atof(str);
      if(GetString(str)){
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

//iMode=0 - any, 1 - uncommented only, no empty strings
int GetStringFromFile(FILE *pF,char *str,int iNChars,int iMode){
  while(fgets(str,iNChars,pF)){
    if((*str=='#' || *str=='\n') && iMode==1) continue;
    else return(0);
  }
  return(1);
}

int GetStringFromSTDIN(char *str,int iNChars,char *strPrompt){
  tcflush(0,TCIOFLUSH);
  printf(strPrompt);
  fflush(stdout);
  fgets(str,iNChars,stdin);
  if(str[strlen(str)-1]=='\n') str[strlen(str)-1]='\0';
  if(!strlen(str)){ 
    printf("INT Empty string\n");
    return(1);
  }
  return(0);
}

int DrawDomain(int iN,POINTD *pDom,int iMode){
  int i;
  switch(iMode){
  case DRAW_DOMAIN_FIRST:
    cpgsci(2);
    cpgpt1((float)(pDom[0].x),(float)(pDom[0].y),2);
    cpgsci(1);
    break;
  case DRAW_DOMAIN_LAST:
    if(iN<2) return(1);
    cpgsci(1);
    cpgmove((float)pDom[iN-2].x,(float)pDom[iN-2].y);
    cpgdraw((float)pDom[iN-1].x,(float)pDom[iN-1].y);
    cpgsci(2);
    cpgpt1((float)pDom[iN-1].x,(float)pDom[iN-1].y,2);
    cpgsci(1);
    break;
  case DRAW_DOMAIN_FINISH:
    if(iN<1) return(2);
    cpgsci(1);
    cpgdraw((float)pDom[iN].x,(float)pDom[iN].y);
    break;
  case DRAW_DOMAIN_CLOSED:
    if(iN<1) return(3);
    cpgmove((float)pDom[0].x,(float)pDom[0].y);
    cpgsci(1);
    for(i=1;i<iN+1;i++) cpgdraw((float)(pDom[i].x),(float)(pDom[i].y));
    cpgsci(2);
    for(i=0;i<iN;i++) cpgpt1((float)(pDom[i].x),(float)(pDom[i].y),2);
    break;
  case DRAW_DOMAIN_OPEN:
    if(iN<1) return(4);
    cpgmove((float)pDom[0].x,(float)pDom[0].y);
    cpgsci(1);
    for(i=1;i<iN;i++) cpgdraw((float)(pDom[i].x),(float)(pDom[i].y));
    cpgsci(2);
    for(i=0;i<iN;i++) cpgpt1((float)(pDom[i].x),(float)(pDom[i].y),2);
    cpgsci(1);
   break;
  }
  return(0);
}

int ReadDomain(int *piN,POINTD **ppDom){
  FILE *pF=NULL;
  int iNRecords;
  //static int iDomainCount=0;
  char str[256],strDum[256],*pc;
  int iError=0;
  double dX,dY;

  if(GetStringFromSTDIN(str,255,"INT Enter domain file: ")) return(1);

  if((pF=fopen(str,"r"))==NULL){
    printf("INT Cannot open file <%s>\n",str);
    return(2);
  }

  iNRecords=0;
  SAFE_FREE(*ppDom);
  while(!GetStringFromFile(pF,strDum,255,1)){
    pc=strtok(strDum," \t");
    if(pc) dX=atof(pc);
    else{
      printf("INT Incomplete X record in %s @ %i\n",str,iNRecords);
      break;
    }
    
    pc=strtok(NULL," \t");
    if(pc) dY=atof(pc);
    else{
      printf("INT Incomplete Y record in %s @ %i\n",str,iNRecords);
      break;
    }
	
    if(!(*ppDom=(POINTD*)realloc(*ppDom,(iNRecords+2)*sizeof(POINTD)))){
      printf("INI Cannot realloc for *ppDom\n");
      iNRecords=0;
      iError=3;
      goto label_bailout;
    }
    
    (*ppDom)[iNRecords].x=dX;
    (*ppDom)[iNRecords].y=dY;

    if(iNRecords)
      cpgdraw((float)dX,(float)dY);	    
    else
      cpgmove((float)dX,(float)dY);	    
      
    iNRecords++;
 }

  (*ppDom)[iNRecords].x=(*ppDom)[0].x;
  (*ppDom)[iNRecords].y=(*ppDom)[0].y;

  cpgdraw((float)(*ppDom)[iNRecords].x,(float)(*ppDom)[iNRecords].y);	    

  printf("INT Processed %i records\n",iNRecords);
  
 label_bailout:

  if(iNRecords<3){
    printf("INT Incomplete domain\n");
    SAFE_FREE(*ppDom);
    *piN=0;
  }
  else{
    *piN=iNRecords;
  }
  
  if(pF){
    fclose(pF);
    pF=NULL;
  }
  
  return(iError);
}


int WriteDomain(int iN,POINTD *pDom){
  FILE *pF=NULL;
  int i;
  static int iDomainCount=0;
  char str[256];
  if(iN==0 || pDom==NULL) return(1);
  snprintf(str,255,"%sDomain%i",filenames[0],iDomainCount);
  if((pF=fopen(str,"w"))==NULL){
    printf("INT Cannot open file %s\n",str);
    return(1);
  }
  printf("INT Saving curent domain in %s\n",str);
  for(i=0;i<iN;i++){
    if(fprintf(pF,"%f\t%f\n",pDom[i].x,pDom[i].y)<0){
      printf("INT Cannot write to <%s>\n",str);
      break;
    }
  }
  printf("INT Wrote %i domain points into %s\n",i--,str);
  iDomainCount++;
  if(pF) fclose(pF);
  return(0);
}

/*
int GetSelection(){
  int iColorSave;
  int iReturn=0;
  double *pdDomainX=(double*)0,*pdDomainY =(double*)0;
  int iDomainN=0;
  float fX1,fY1,fX2,fY2;
  char ch,str[256];

  cpgqci(&iColorSave);

  cpgsci(1);
  
 label_start:
  RestorePanel(xd,yd,fbufferi,mode,fmin,fmax);
  cpgband(0,0,x,y,&x,&y,&ch);
  if(ch=='q'){
    iReturn=0;
    goto label_exit;
  }
  if(ch=='r'){
    if(GetString(str)) fX1=atof(str);
    else continue;
    if(GetString(str)) fY1=atof(str);
    else continue;
    if(GetString(str)) fX2=atof(str);
    else continue;
    if(GetString(str)) fY2=atof(str);
    else continue;
    
    iDomainN=4;
    pdDomainX=(double *)malloc((iDomainN+1)*sizeof(double));
    pdDomainY=(double *)malloc((iDomainN+1)*sizeof(double));
    pdDomainX[0]=pdDomainX[iDomainN]=(double)fX1;
    pdDomainY[0]=pdDomainY[iDomainN]=(double)fY1;
    pdDomainX[1]=(double)fX2;
    pdDomainY[1]=(double)fY1;
    pdDomainX[2]=(double)fX2;
    pdDomainY[2]=(double)fY2;
    pdDomainX[3]=(double)fX1;
    pdDomainY[3]=(double)fY2;
    cpgmove(fX1,fY1);
    cpgsci(1);
    for(i=1;i<5;i++) cpgdraw((float)(pdDomainX[i]),(float)(pdDomainY[i]));
    cpgsci(2);
    for(i=0;i<4;i++) cpgpt1((float)(pdDomainX[i]),(float)(pdDomainY[i]),2);

    Stencil(xd,yd,iDomainN,pdDomainX,pdDomainY,&stenciln,&stencil);	    
    
    iReturn=1;
    goto label_exit;
  }
  if(ch=='m'){
    if(GetString(str)) x=atof(str);
    else continue;
    if(GetString(str)) y=atof(str);
    else continue;
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
      RestorePanel(xd,yd,fbufferi,mode,fmin,fmax);
      iReturn=1;
      goto label_exit;
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
      goto label_start;
    }
    if(ch=='X'){
      if(iDomainN>2){ 
	pdDomainX=(double *)realloc((void *)pdDomainX,(iDomainN+1)*sizeof(double));
	pdDomainY=(double *)realloc((void *)pdDomainY,(iDomainN+1)*sizeof(double));
	pdDomainX[iDomainN]=pdDomainX[0];
	pdDomainY[iDomainN]=pdDomainY[0];
	cpgdraw((float)pdDomainX[iDomainN],(float)pdDomainY[iDomainN]);	    
	Stencil(xd,yd,iDomainN,pdDomainX,pdDomainY,&stenciln,&stencil);	    
	iReturn=1;
	break;
      }
      else printf("INT Needs three points to complete selection\n");
    }
  }

 label_exit:
  iDomainN=0;
  if(pdDomainX) free(pdDomainX);
  if(pdDomainY) free(pdDomainY);

  cpgsci(iColorSave);
  return(iReturn);
}

*/
