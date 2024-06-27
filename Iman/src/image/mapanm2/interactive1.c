#include "mapanm21.h"

#define INTERACTIVE_WINDOW_SIZEX 7.0

void RestorePanel(int xd,int yd,double **dbuf,float *buf,int mode);
char* GetString(char *str);

void Interactive(int panel,int xd,int yd,double **dmp){
  static float x=0.0,y=0.0;
  int interactiveid=0,callerid=0;
  int xi,yi;
  char ch,str[128];
  int i,j,l;
  double phi,rho,phip,rhop,dDum;
  unsigned long xy;
  unsigned short dumus;
  FILE *fp;
  int panelmode=0;
  int iX,iY;

  cpgqid(&callerid);
  xy=(unsigned long)xd*(unsigned long)yd;

  setenv("PGPLOT_ENVOPT","I",1); 
  interactiveid=cpgopen("/xw");
  cpgask(0);
  cpgpap(INTERACTIVE_WINDOW_SIZEX,(float)yd/(float)xd);
  cpgsch(0.8);
  cpgenv(0.5,(float)xd+0.5,(float)yd+0.5,0.5,1,0);
  cpglab("","",modestr[panel-1]);
  unsetenv("PGPLOT_ENVOPT");
  
  rhop=hypot(dmp[0][0],dmp[1][0]);
  phip=180.0/M_PI*atan2(dmp[1][0],dmp[0][0]);
  
  RestorePanel(xd,yd,dmp,fbufferi,panelmode);

  while(1){
    cpgslct(interactiveid);
    cpgband(7,1,x,y,&x,&y,&ch); 
    //    printf("INT X %f Y %f C (%i)%c\n",x,y,ch,ch);

    // Keys in use: q,D,Q,m,S,p,h,X,T,Z

    if(ch=='q' || ch=='D'){ break;}

    if(ch=='Q'){ exit(0);}

    if(ch=='m'){
      cpgslct(interactiveid);
      x=atof(GetString(str));
      printf(" %s(%f) ",str,x);
      y=atof(GetString(str));
      printf("INT %s(%f) ",str,y);
    }

    if(ch=='S'){
      sprintf(str,"%s.%s",savefile,modestr[panel-1]);
      if(pcSaveFileSuffix) strcat(str,pcSaveFileSuffix);
      printf("INT Saving curent map in %s\n",str);
      if((fp=fopen(str,"w"))){
	i=(SAVE_TYPE_BIT_DOUBLE_FLOAT&type) |
	  (SAVE_TYPE_BIT_ADD_DOUBLE_AT_END&type) |
	  (SAVE_TYPE_BIT_XY_IN_FILES*SAVE_XY_IN_ONE_FILE);       
	if(fwrite((void*)&i,sizeof(int),1,fp) != 1){
	  printf("INT Cannot write type %i to file %s\n",i,str);
	}
	dumus=(unsigned short)xd;
	if(fwrite((void*)&dumus,sizeof(unsigned short),1,fp) != 1){
	  printf("INT Cannot write Xdim %i to file %s\n",dumus,str);
	}
	dumus=(unsigned short)yd;
	if(fwrite((void*)&dumus,sizeof(unsigned short),1,fp) != 1){
	  printf("INT Cannot write Ydim %i to file %s\n",dumus,str);
	}
	if(fwrite((void*)(dmp[0]),XYdim*sizeof(double),1,fp) != 1){
	  printf("INT Cannot write dmp[0] to file %s\n",str);
	}
	if(fwrite((void*)(dmp[1]),XYdim*sizeof(double),1,fp) != 1){
	  printf("INT Cannot write dmp[1] to file %s\n",str);
	}
	if((type&SAVE_TYPE_BIT_ADD_DOUBLE_AT_END)){
	  if(fwrite((void*)&harmonic_d,sizeof(double),1,fp) != 1){
	    printf("INT Cannot write harmonic %f to file %s\n",harmonic_d,str);
	  }
	}
	fclose(fp);
      }
      else{
	printf("INT Cannot open file %s\n",str);
      }
      continue;
    }

    if(ch=='p'){
      RestorePanel(xd,yd,dmp,fbufferi,(panelmode=1-panelmode));
    }

    if(ch=='h'){
      printf("q or D - quit interactive dialog\n");
      printf("Q - quit program\n");
      printf("m - <int> <int> manual input of coordinates\n");
      printf("p - switch, phi <-> rho image on interactive panel (defined by options)\n");
      printf("S - save current maps\n");
      printf("h - print this help\n");
      printf("X - flip phase by Pi in a region (size %i)\n",iFlipPhaseSize);
      continue;
    }

    if(ch=='T'){
      cpgslct(interactiveid);
      dFlipPhaseThreshold=atof(GetString(str));
      printf("INT Phase Flip Threshold %f\n",dFlipPhaseThreshold);
      continue;
    }

    if(ch=='Z'){
      cpgslct(interactiveid);
      iFlipPhaseSize=atoi(GetString(str));
      printf("INT Phase Flip Threshold %i\n",iFlipPhaseSize);
      continue;
    }

    xi=(int)rint(x);
    if(xi<0) xi=0;
    if(xi>xd-1) xi=xd-1;
    yi=(int)rint(y);
    if(yi<0) yi=0;
    if(yi>yd-1) yi=yd-1;
    printf("INT X %i Y %i\n",xi,yi);    


    if(ch=='X'){
      //      printf("INT Flipping phase at %i %i\n",xi,yi);
      l=xi+yi*xd;
      phi=atan2(dmp[1][l],dmp[0][l]);
      if(phi<0.0) phi+=2*M_PI;
      
      phip=phi;
      for(i=1;i<=iFlipPhaseSize;i++){
	iX=xi+i;
	if(iX<0 || iX>xd-1) continue;
	l=iX+yi*xd;
	dDum=atan2(dmp[1][l],dmp[0][l]);
	if(dDum<0.0) dDum+=2*M_PI;
	dDum=fabs(phip-dDum);
	if(dDum>M_PI*dFlipPhaseThreshold && 
	   dDum<M_PI*(2.0-dFlipPhaseThreshold)){
	  dmp[0][l]=-dmp[0][l];
	  dmp[1][l]=-dmp[1][l];
	}	
	phip=atan2(dmp[1][l],dmp[0][l]);
	if(phip<0.0) phip+=2*M_PI;
      }

      phip=phi;
      for(i=-1;i>=-iFlipPhaseSize;i--){
	iX=xi+i;
	if(iX<0 || iX>xd-1) continue;
	l=iX+yi*xd;
	dDum=atan2(dmp[1][l],dmp[0][l]);
	if(dDum<0.0) dDum+=2*M_PI;
	dDum=fabs(phip-dDum);
	if(dDum>M_PI*dFlipPhaseThreshold && 
	   dDum<M_PI*(2.0-dFlipPhaseThreshold)){
	  dmp[0][l]=-dmp[0][l];
	  dmp[1][l]=-dmp[1][l];
	}	
	phip=atan2(dmp[1][l],dmp[0][l]);
	if(phip<0.0) phip+=2*M_PI;
      }

      for(i=-iFlipPhaseSize;i<=iFlipPhaseSize;i++){
	iX=xi+i;
	if(iX<0 || iX>xd-1) continue;

	l=iX+yi*xd;
	phi=atan2(dmp[1][l],dmp[0][l]);
	if(phi<0.0) phi+=2*M_PI;
	
	phip=phi;
	for(j=1;j<=iFlipPhaseSize;j++){
	iY=yi+j;
	if(iY<0 || iY>yd-1) continue;
	  l=iX+iY*xd;
	  dDum=atan2(dmp[1][l],dmp[0][l]);
	  if(dDum<0.0) dDum+=2*M_PI;
	  dDum=fabs(phip-dDum);
	  if(dDum>M_PI*dFlipPhaseThreshold && 
	     dDum<M_PI*(2.0-dFlipPhaseThreshold)){
	    dmp[0][l]=-dmp[0][l];
	    dmp[1][l]=-dmp[1][l];
	  }	  
	  phip=atan2(dmp[1][l],dmp[0][l]);
	  if(phip<0.0) phip+=2*M_PI;
	}

	phip=phi;
	for(j=-1;j>=-iFlipPhaseSize;j--){
	iY=yi+j;
	if(iY<0 || iY>yd-1) continue;
	  l=iX+iY*xd;
	  dDum=atan2(dmp[1][l],dmp[0][l]);
	  if(dDum<0.0) dDum+=2*M_PI;
	  dDum=fabs(phip-dDum);
	  if(dDum>M_PI*dFlipPhaseThreshold && 
	     dDum<M_PI*(2.0-dFlipPhaseThreshold)){
	    dmp[0][l]=-dmp[0][l];
	    dmp[1][l]=-dmp[1][l];
	  }	  
	  phip=atan2(dmp[1][l],dmp[0][l]);
	  if(phip<0.0) phip+=2*M_PI;
	}

      }
      
      /*
      for(j=-iFlipPhaseSize;j<=iFlipPhaseSize;j++){
	iY=yi+j;
	if(iY<0 || iY>yd-1) continue;
	for(i=-iFlipPhaseSize;i<=iFlipPhaseSize;i++){
	  iX=xi+i;
	  if(iX<0 || iX>xd-1) continue;
	  l=iX+iY*xd;
	  dDum=atan2(dmp[1][l],dmp[0][l]);
	  if(dDum<0.0) dDum+=2*M_PI;
	  dDum=fabs(phi-dDum);
	  if(dDum>M_PI*dFlipPhaseThreshold && 
	     dDum<M_PI*(2.0-dFlipPhaseThreshold)){
	    dmp[0][l]=-dmp[0][l];
	    dmp[1][l]=-dmp[1][l];
	  }	  
	}
      }
      */
      RestorePanel(xd,yd,dmp,fbufferi,panelmode);
    }

    l=xi+yi*xd;

    rho=hypot(dmp[0][l],dmp[1][l]);
    phi=180.0/M_PI*atan2(dmp[1][l],dmp[0][l]);
    printf("INT Rho=%f(%f %f) Phi=%f(%f) (%f)\n",rho,rho-rhop,rho/rhop,phi,phi-phip,180.0/M_PI*fbufferi[l]);
    rhop=rho;
    phip=phi;
  }

  cpgslct(interactiveid);
  cpgclos();
  cpgslct(callerid);
}

void RestorePanel(int xd,int yd,double **dbuf,float *buf,int mode){
  int i;
  float fmin,fmax;

  if(mode==0){
    for(i=0;i<xd*yd;i++){
      buf[i]=(float)atan2(dbuf[1][i],dbuf[0][i]);
    }
    fmin=-M_PI;
    fmax=M_PI;
    Pallet(2,1.0,0.5);
    cpgimag(buf,xd,yd,1,xd,1,yd,fmin,fmax,trans);
  }
  else{
    fmax=-(fmin=1.0e10);
    for(i=0;i<xd*yd;i++){
      buf[i]=(float)hypot(dbuf[1][i],dbuf[0][i]);
      if(buf[i]>fmax) fmax=buf[i];
      if(buf[i]<fmin) fmin=buf[i];
    }
    cpggray(buf,xd,yd,1,xd,1,yd,fmax,fmin,trans);  
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
      str[i]='\0';
      printf("%c",ch);
      fflush(stdout);
      break;
    }
  }
  return(str);
}
