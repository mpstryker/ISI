/* Bin analysis                 */
/* Written by V.Kalatsky        */
/* 20Jul2001                    */

#include "binan1.h"

#define INTERACTIVE_WINDOW_SIZEX 8.0
//#define SAVE_WINDOW_SIZE 4.0
#define SAVE_WINDOW_SIZE 2.0

#define TUNING_WINDOW_SIZE 7.0
#define TUNING_WINDOW_RATIO 0.6
#define TUNING_CURVE_N 5

void RestorePanel(int xd,int yd,float *buf,int mode,float fmin,float fmax);
char* GetString(char *str);
void DrawTuningCurve(int id,int n,float **mp,unsigned long pixel,int iShift,int iSaveCurve);

void Interactive(int xd,int yd,float **mp,int *mode){
  float fmin=0.0,fmax=0.0,x=0.0,y=0.0;
  int interactiveid=0,callerid,tuningid=0,dumid=0,saveid=0;
  int xi,yi;
  char ch,str[128],savefile[128];
  int i,l;
  unsigned long lmax=0,lmin=0;
  unsigned long xy;
  int iShift=0;

  if(*mode<0 && *mode>=nbins) return;

  cpgqid(&callerid);

  xy=(unsigned long)xd*(unsigned long)yd;

  setenv("PGPLOT_ENVOPT","I",1); 
  interactiveid=cpgopen("/xw");
  cpgask(0);
  cpgpap(INTERACTIVE_WINDOW_SIZEX,(float)yd/(float)xd);
  cpgsch(0.8);
  cpgenv(pminx,pmaxx,pminy,pmaxy,1,0);
  unsetenv("PGPLOT_ENVOPT");
  
  FindExtremaEx(xy,(unsigned long)xd,fmaps[*mode],&fmin,&fmax,&lmin,&lmax);
  RestorePanel(xd,yd,fmaps[*mode],*mode,fmin,fmax);

  tuningid=cpgopen("/xw");
  cpgask(0);
  cpgpap(TUNING_WINDOW_SIZE,TUNING_WINDOW_RATIO);
  cpgsch(0.8);
  cpgenv(0.0,(float)nbins,-1.0,1.0,0,1);

  l=0;
  while(1){
    cpgslct(interactiveid);
    cpgband(7,1,x,y,&x,&y,&ch); 
    if(do_verbose) printf("INT X %f Y %f C (%i)%c\n",x,y,ch,ch);

    // Keys in use: q,Q,m,+/-,v, ,S,h,s,p

    if(ch=='q'){ break;}

    if(ch=='Q'){ exit(0);}

    if(ch=='m'){
      x=atof(GetString(str));
      printf(" %s(%f) ",str,x);
      y=atof(GetString(str));
      printf("INT %s(%f) ",str,y);
    }

    if(ch=='-' || ch=='+'){
      if(ch=='-'){
	*mode=(*mode+nbins-1)%nbins;
	iShift=(iShift+nbins-1)%nbins;
      }
      else{
	*mode=(*mode+1)%nbins;
	iShift=(iShift+1)%nbins;
      }
      FindExtremaEx(xy,(unsigned long)xd,fmaps[*mode],&fmin,&fmax,&lmin,&lmax);
      RestorePanel(xd,yd,fmaps[*mode],*mode,fmin,fmax);
      continue;
    }

    if(ch=='v'){
      do_verbose=1-do_verbose;
      continue;
    }

    if(ch==' '){
      DrawTuningCurve(tuningid,nbins,(float**)0,0L,0,0);
      continue;      
    }

    if(ch=='S'){
      cpgqid(&dumid);
      for(i=0;i<nbins;i++){
	//	sprintf(savefile,"%s_%03i%s",filenames[0],i,".gif/GIF");
	sprintf(savefile,"%s%02i%s","mov",i,".gif/GIF");
	printf("Saving %3i in %s\n",i,savefile);
	saveid=cpgopen(savefile);
	cpgpap(SAVE_WINDOW_SIZE,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
	cpgsch(0.1);
	cpgenv(pminx,pmaxx,pminy,pmaxy,1,-1);
	cpgask(0);

	FindExtremaEx(xy,(unsigned long)xd,fmaps[*mode],&fmin,&fmax,&lmin,&lmax);
	cpggray(fmaps[i],xd,yd,1,xd,1,yd,fmax,fmin,trans);
	//	getchar();
	cpgslct(saveid);
	cpgclos();
      }
      cpgslct(dumid);
      continue;      
    }

    if(ch=='s'){ 
      DrawTuningCurve(tuningid,nbins,mp,l,iShift,1);
    }

    if(ch=='p'){
      iPlotFComponents=(iPlotFComponents+1)%3;
      continue;
    }

    if(ch=='h'){
      printf("q - quit interactive dialog\n");
      printf("Q - quit program\n");
      printf("m - <int> <int> manual input of coordinates\n");
      printf("-/+ - switch, bin image on interactive panel (default 0)\n");
      printf("v - switch, verbose mode (default off)\n");
      printf(" (Space) - clear accumulated tuning curves\n");
      printf("S - save as PS files\n");
      printf("s - save curve\n");
      printf("p - plot Fourier components switch\n");
      printf("h - print this help\n");
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

    DrawTuningCurve(tuningid,nbins,mp,l,iShift,0);
  }
  cpgslct(interactiveid);
  cpgclos();
  cpgslct(callerid);
}

void RestorePanel(int xd,int yd,float *buf,int mode,float fmin,float fmax){
  static char str[16]={"\0"};
  int i;
  cpgqci(&i);
  cpggray(buf,xd,yd,1,xd,1,yd,fmax,fmin,trans);
  cpgsci(0);
  cpglab("","",str);
  sprintf(str,"%i",mode);  
  cpgsci(7);
  cpglab("","",str);
  cpgsci(i);
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

#define SIN_PLOT_N_EXTRA_POINTS 5

void DrawTuningCurve(int id,int n,float **mp,unsigned long pixel,int iShift,int iSaveCurve){
  int callerid;
  int i,j,dumi;
  unsigned long ul;
  static int n_saved_curves=0;
  static int current_curve=-1;
  static float *tcurve[TUNING_CURVE_N];
  static float *ft_tcurve[TUNING_CURVE_N];
  static float fmin[TUNING_CURVE_N],fmax[TUNING_CURVE_N];
  static unsigned long lmin[TUNING_CURVE_N],lmax[TUNING_CURVE_N];
  static unsigned long tpixel[TUNING_CURVE_N];
  static float globalmin=BIG_FLOAT,globalmax=-BIG_FLOAT;
  static char str[64]={"\0"};
  char strSaveFile[256];
  FILE *fp;
  double dPhi,dRho,dOmega,dDum;

  cpgqid(&callerid);
  cpgslct(id);

  if(!mp){
    n_saved_curves=0;
    current_curve=-1;
    for(i=0;i<TUNING_CURVE_N;i++) free(tcurve[i]);
    for(i=0;i<TUNING_CURVE_N;i++) free(ft_tcurve[i]);
    cpgsci(1);
    cpgenv(0.0,(float)n,globalmin,globalmax,0,1);
    if(callerid) cpgslct(callerid);
    return;
  }

  if(!n_saved_curves){
    for(i=0;i<TUNING_CURVE_N;i++) tcurve[i]=(float*)malloc(n*sizeof(float));
    for(i=0;i<TUNING_CURVE_N;i++) ft_tcurve[i]=(float*)malloc(n*sizeof(float));
  }
  if(n_saved_curves<TUNING_CURVE_N) n_saved_curves++;
  current_curve=(current_curve+1)%TUNING_CURVE_N;
  tpixel[current_curve]=pixel;
  for(i=0;i<n;i++){
    tcurve[current_curve][i]=mp[(i+iShift)%n][pixel];
  }
  FindExtremaEx((unsigned long)n,0L,tcurve[current_curve],fmin+current_curve,fmax+current_curve,lmin+current_curve,lmax+current_curve);
  FindMin((unsigned long)n_saved_curves,fmin,&globalmin,&ul);
  FindMax((unsigned long)n_saved_curves,fmax,&globalmax,&ul);

  //Fourier stuff
  for(i=iMinFComponent;i<iMinFComponent+iNFComponents;i++){
    FindFTComponents(n,tcurve[current_curve],i,
		     ft_tcurve[current_curve]+(2*i),
		     ft_tcurve[current_curve]+(2*i+1));
  }

  cpgsci(1);
  cpgenv(0.0,(float)n,globalmin,globalmax,0,1);
  for(j=0;j<n_saved_curves;j++){
    sprintf(str,"%lu-%lu ",tpixel[j]%Xdim,tpixel[j]/Xdim);
    cpgsci(2+j);
    if(j==current_curve){
      cpgqcf(&dumi);
      cpgscf(2+dumi);
      cpgtext((j+1)*((float)n)*0.09,globalmax+(globalmax-globalmin)*0.04,str);
      cpgscf(dumi);
    }
    else{
      cpgtext((j+1)*((float)n)*0.09,globalmax+(globalmax-globalmin)*0.04,str);
    }
    cpgmove(0.0,tcurve[j][0]);
    for(i=0;i<n;i++){
      cpgdraw((float)i,tcurve[j][i]);
    }
    cpgdraw((float)n,tcurve[j][0]);
  }

  if(iPlotFComponents && iNFComponents){
    // Separate components
    if(iPlotFComponents&1){
      for(i=iMinFComponent;i<iMinFComponent+iNFComponents;i++){
	cpgsls(2+2*(i-iMinFComponent));
	dPhi=atan2(ft_tcurve[current_curve][2*i+1],ft_tcurve[current_curve][2*i]);
	dRho=hypot(ft_tcurve[current_curve][2*i+1],ft_tcurve[current_curve][2*i]);
	dOmega=2.0*M_PI*(double)i/(double)n;
	cpgmove(0.0,(float)(dRho*cos(dPhi)));      
	for(j=1;j<SIN_PLOT_N_EXTRA_POINTS*n;j++){
	  cpgdraw((float)j/SIN_PLOT_N_EXTRA_POINTS,(float)(dRho*cos(-dPhi+j*dOmega/SIN_PLOT_N_EXTRA_POINTS)));
	}
      }
    }
    
    // Reconstruct curve
    if(iPlotFComponents&2){
      cpgsls(2);
      cpgmove(0.0,FourierSum(n,ft_tcurve[current_curve],iNFComponents,iMinFComponent,0.0));
      for(j=1;j<SIN_PLOT_N_EXTRA_POINTS*n;j++){
	cpgdraw((float)j/SIN_PLOT_N_EXTRA_POINTS,FourierSum(n,ft_tcurve[current_curve],iNFComponents,iMinFComponent,(float)j/SIN_PLOT_N_EXTRA_POINTS));
      }
      dDum=FindExtremaFComponents(n,ft_tcurve[current_curve],iNFComponents,iMinFComponent,1);
      printf("INT Max %f(%f)\n",dDum,FourierSum(n,ft_tcurve[current_curve],iNFComponents,iMinFComponent,dDum));
      cpgslw(4);
      cpgpt1((float)dDum,(float)FourierSum(n,ft_tcurve[current_curve],iNFComponents,iMinFComponent,dDum),0);
      cpgslw(1);
    }

    cpgsls(1);
  }

  if(callerid) cpgslct(callerid);

  if(iSaveCurve){
    sprintf(strSaveFile,"%s_%li-%licurve",filenames[0],tpixel[current_curve]%Xdim,tpixel[current_curve]/Xdim);
    printf("Saving curve in %s\n",strSaveFile);
    if((fp=fopen(strSaveFile,"w"))==NULL){
      printf("INI Cannot open file %s\n",strSaveFile);
      return;
    }
    
    for(i=0;i<n;i++){
      fprintf(fp,"%f %f\n",fBinSizeSec*(float)i,tcurve[current_curve][i]);
    }
    fclose(fp);
  }
}
