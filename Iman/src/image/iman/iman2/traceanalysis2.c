#include "iman1.h"

//#define TEST
//#define PLOT_PHASE
//#define PLOT_DIFFERENCE

float TestFunction(float x,float p1,float p2,float p3,float p4){
  double c;

  //  c=(x-p1-p2*( 0.5+floor((double)((x-p1)/p2)) ))/p3;
  //  c=exp(-c*c*0.5)-p4;
  //  c=c<0.0 ? 0.0 : c;

  //  c=cos(x*p1);

  //   c=-(x-(initframe-1))*(x-(finalframe+1));

  //  c=c>0.5 ? 0.25*c : 0.0;
  
  c=(x-initframe)/(finalframe-initframe);
  if(p1>0.5) c*=c;
  if(p1>1.5) c*=c;
  if(p1>2.5) c*=c;
  return((float)(c));
}

#define BIG_UNSIGNEDSHORT 0xFFFE
#define BIG_FLOAT 10e10
// used for paper
#define SAVE_PLOT_WINDOW_SIZEX 12.0
//#define SAVE_PLOT_WINDOW_SIZEX 15.0
#define PLOT_WINDOW_RATIO 0.8
// used for paper
#define SAVE_PLOT_WINDOW_RATIO 1.1
//#define SAVE_PLOT_WINDOW_RATIO 0.7
//#define CHARACTERSIZE1 1.1
//#define CHARACTERSIZE2 1.5
#define CHARACTERSIZE1 1.6
#define CHARACTERSIZE2 2

#define NTRACES 4
#define NPLOTS 3

#define SAVE_MODE1 1
#define SAVE_MODE2 2
#define SAVE_MODE3 4

#define POWER_SPECTRUM_LOG_PLOT

#define FOURIER_THRESHOLD -2.0
#define FOURIER_THRESHOLD_BIT_DEPTH 16

void realft(float*,unsigned long,int);

int TraceAnalysis(TRAIN *tr,int ifr,int ffr,int a){
  int interactiveid,plotid,saveid;
  float x,y,dumf,dumfx,dumfy;
  char ch,str[128];
  float *trace[NTRACES],tracemax[NTRACES],tracemin[NTRACES];
  float *ftrace[NTRACES];
  float *ftracem[NTRACES],ftracemaxm[NTRACES],ftraceminm[NTRACES];
  float *fatrace=NULL,*fatracem=NULL;
  float *atrace[NTRACES];
  double resp=0.0,imsp=0.0,phil;
  float pssp;
  int tracepixel[NTRACES];
  unsigned long tracefttn,nframes;
  float tracemaxglobal=0.0,traceminglobal=0.0;
  float atracemaxglobal,atraceminglobal;
  float ftracemaxglobal=0.0,ftraceminglobal=0.0;
  double tracemean,traceshift,tracemeanx;
  double s,sx,sxx,delta;
  double fit_a=0.0,fit_b=0.0;
  int ntraces=0,ntraces_count=0,current_trace=0,plot_trace;
  float iwminx,iwmaxx,iwminy,iwmaxy, pwminx,pwmaxx,pwminy,pwmaxy;
  float accumulated_timeaverage,timeinterval,timeinterval_local;
  int xi=0,yi=0,l=0;
  int i=0,j,k;
  int do_shift=0,do_fft=0,plot_fit=0;
  CAR *pcar;
  double co,si,coN,siN,CO,SI,CON,SIN,re,im,dumd;
  int cut_offll,cut_offuu;
  int trace_poly_fitn;
  double *poly_fit=(double*)0;
  // used 2 for paper?
  int save_lw=1;
  int save_mode=SAVE_MODE1 | SAVE_MODE3;
  //     int save_mode=SAVE_MODE1;
  FILE *fp;
  unsigned short us;
  int iDoFrequencyConvert=1;
  int iDoTimeConvert=1;
  int iDoLabelOnSavePlot=0;
  float fInterFrameTime=0.133332;
  float fTimeCoeff=1.0;
  float fFrequencyCoeff=1.0;
  float fFourierThreshold=FOURIER_THRESHOLD;
  float fTraceCoeff=0.0001;
  int iPlotAverageTraceOnSavePlot=1;
  double dNFrames2;
  int iIniAFrame,iFinAFrame;
  int iPlotHPFourierTrace=1;
  float *pfDum;
  int iOldTraces=0;
  int iDoSavePlotBeginTimeOffset=1;
  int iSavePlotBeginTimeOffset=0;
  int iDoSavePlotPoints=0;
  int iDoSavePlotColor=0;
  double dFrameTotal;

#ifdef TEST
  float omegadum,alphadum;
  float p1,p2,p3,p4;
#endif
#ifdef PLOT_PHASE
  int phaseid;
  float *ftracep[NTRACES];
#endif

  if(ifr<0 || ifr>=tr->max_n){
    printf("TRA Initial frame (%i) is out of range(%i-%i)\n",ifr,0,tr->max_n-1);
    return(1);
  }
  if(ffr<0 || ffr>=tr->max_n){
    printf("TRA Final frame (%i) is out of range(%i-%i)\n",ffr,0,tr->max_n-1);
    return(2);
  }
  if(ffr<ifr){
    printf("TRA Negative interval Initial frame %i Final frame %i\n",ifr,ffr);
    return(3);
  }

  nframes=ffr-ifr+1;
  printf("TRA Tracing frames: %i-%i (%li)\n",ifr,ffr,nframes);

  tracefttn=1;
  while(tracefttn<nframes) tracefttn<<=1;
  tracefttn<<=1;

  if(!synch_phi){
    cyclesexact=(double)(n_cycles)*harmonic;
  }
  
  if(poly_fitn){
    trace_poly_fitn=poly_fitn;
  }
  else{
    trace_poly_fitn=(int)ceil((double)n_cycles*POLYNOM_COEFF)+POLYNOM_ADD;
  }
  printf("TRA Will fit %i smart polies\n",trace_poly_fitn);

  if(!(poly_fit=(double *)calloc(trace_poly_fitn,sizeof(double)))){
    printf("TRA Cannot allocate for poly_fit(%iB)\n",trace_poly_fitn*sizeof(double));
  }

  for(i=0;i<NTRACES;i++){
    if(!(trace[i]=(float *)calloc(nframes,sizeof(float)))){
      printf("TRA Cannot allocate for trace[%i]\n",i);
      return(3);
    }
    if(!(ftrace[i]=(float *)calloc(tracefttn,sizeof(float)))){
      printf("TRA Cannot allocate for ftrace[%i]\n",i);
      return(3);
    }
    if(!(ftracem[i]=(float *)calloc(nframes/2+1,sizeof(float)))){
      printf("TRA Cannot allocate for ftracem[%i]\n",i);
      return(3);
    }

#ifdef PLOT_PHASE
    if(!(ftracep[i]=(float *)calloc(nframes/2+1,sizeof(float)))){
      printf("TRA Cannot allocate for ftracep[%i]\n",i);
      return(3);
    }
#endif
    // Using only one array out of NTRACES, 
    // will use array 0 instead of current_trace for ploting and 
    // array 1 for difference between trace and atrace
    // array 2 for difference between trace and fit
    if(!(atrace[i]=(float *)calloc(nframes,sizeof(float)))){
      printf("TRA Cannot allocate for atrace[%i]\n",i);
      return(3);
    }
  }

  if(iPlotHPFourierTrace){
    if(!(fatrace=(float *)calloc(tracefttn,sizeof(float)))){
      printf("TRA Cannot allocate for fatrace. Will go on without HPFourier Plot\n");
      iPlotHPFourierTrace=0;
    }
    else{
      if(!(fatracem=(float *)calloc(nframes/2+1,sizeof(float)))){
	printf("TRA Cannot allocate for fatracem. Will go on without HPFourier Plot\n");
	free(fatrace);
	iPlotHPFourierTrace=0;
      }
    }
  }

  timeinterval=(float)(2*timeradius+1);

  if(iDoFrequencyConvert) fFrequencyCoeff=1.0/(fInterFrameTime*(1+ffr-ifr));
  else fFrequencyCoeff=1.0;
  if(iDoTimeConvert) fTimeCoeff=fInterFrameTime;
  else fTimeCoeff=1.0;

  iwminx=-0.5;
  iwmaxx=(float)(Xdim)-0.5;
  iwminy=(float)(Ydim)-0.5;
  iwmaxy=-0.5;
  pwminx=(float)ifr;
  pwmaxx=(float)ffr;
  pwminy=0.0;
  pwmaxy=(float)BIG_UNSIGNEDSHORT;
  
  if(!(pcar=AddFrame(tr,ifr,a))){
    printf("TRA AddFrame return 0 for initframe %i\n",ifr);     
    return(2);
  }
  x=BIG_FLOAT;
  y=0.0;
  dFrameTotal=0.0;
  for(i=0;i<XYdim;i++){
    fbuffer[i]=dumf=(float)pcar->pimage[i];
    dFrameTotal+=dumf;
    if(i>=Xdim){ // Do not take into account first screwed row of each frame
      if(dumf>y) y=dumf;
      if(dumf<x) x=dumf;
    }
  }
  dFrameTotal/=XYdim;

  if(iUseExternalFourierThreshold){
    fFourierThreshold=fExternalFourierThreshold;
  }
  else{
    fFourierThreshold=FOURIER_THRESHOLD+rint(log10(dFrameTotal/50000.0));
  }


  // Quick hack to print first frame
  /*
  cpgopen("watch1.ps/vps");
  setenv("PGPLOT_ENVOPT","I",1); 
  cpgask(0);
  cpgpap(5.5,1.0*(float)fabs((double)((iwmaxy-iwminy)/(iwmaxx-iwminx))));
  cpgsch(0.8);
  cpgenv(iwminx,iwmaxx,iwminy,iwmaxy,1,1);
  unsetenv("PGPLOT_ENVOPT");   
  cpggray(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,x,y,trans);
  cpgclos();
  */
  interactiveid=cpgopen(device);
  setenv("PGPLOT_ENVOPT","I",1); 
  cpgask(0);
  cpgpap(fSizeOfTraceInteractiveWindow,1.0*(float)fabs((double)((iwmaxy-iwminy)/(iwmaxx-iwminx))));
  cpgsch(0.8);
  cpgenv(iwminx,iwmaxx,iwminy,iwmaxy,1,1);
  unsetenv("PGPLOT_ENVOPT");
   
  cpggray(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,y,x,trans);

#ifdef PLOT_PHASE
  phaseid=cpgopen(device);
  cpgsubp(1,1);
  cpgpap(fSizeOfTracePlotWindow,PLOT_WINDOW_RATIO/3.0);
  cpgask(0);
#endif

  plotid=cpgopen(device);
  cpgsubp(1,3);
  cpgpap(fSizeOfTracePlotWindow,PLOT_WINDOW_RATIO);
  cpgask(0);
  x=y=0.0;
  ntraces=ntraces_count=current_trace=0;
  while(1){
    cpgslct(interactiveid);
    cpgband(7,0,x,y,&x,&y,&ch); 
    //    printf("TRA X %f Y %f C (%i)%c\n",x,y,ch,ch);

    // Keys in use s,S,q,Q,m,0
    if(ch=='s'){
      if(ntraces){
	if(iDoSavePlotBeginTimeOffset) iSavePlotBeginTimeOffset=ifr;
	else iSavePlotBeginTimeOffset=0;
	if(save_mode){
	  save_mode%=8;
	  k=(save_mode&1)+((save_mode&2)>>1)+save_mode/4;

	  plot_trace=(ntraces_count-1)%NTRACES;
	
	  sprintf(str,"%s%s%i%s%i%s%i_%ix%i-%i%s",tr->filename,".",n_cycles,"_",initframe,"_",finalframe,xi,yi,save_mode,".trace.ps/CPS");	
	  printf("%s\n",str);
	  saveid=cpgopen("?");
	  cpgsubp(1,k);
	  cpgpap(SAVE_PLOT_WINDOW_SIZEX,k*SAVE_PLOT_WINDOW_RATIO/3.0);

	  pwminx=(float)(ifr-iSavePlotBeginTimeOffset)*fTimeCoeff;
	  pwmaxx=(float)(ffr-iSavePlotBeginTimeOffset)*fTimeCoeff;
	  
	  if(save_mode&1){
	    pwminy=tracemin[plot_trace]*fTraceCoeff;
	    pwmaxy=tracemax[plot_trace]*fTraceCoeff;
	    cpgsci(1);
	    cpgsch(CHARACTERSIZE1);
	    cpgslw(save_lw+1);

	    cpgenv(pwminx,pwmaxx,pwminy,pwmaxy,0,2);
	    if(iDoLabelOnSavePlot){
	      sprintf(str,"X %i Y %i",xi,yi);
	      cpgsch(CHARACTERSIZE2);
	      cpglab("","",str);
	    }
	    if(iDoSavePlotColor) cpgsci(3);
	    cpgmove((float)(ifr-iSavePlotBeginTimeOffset)*fTimeCoeff,trace[plot_trace][0]*fTraceCoeff);
	    for(i=ifr+1;i<=ffr;i++){
	      cpgdraw((float)(i-iSavePlotBeginTimeOffset)*fTimeCoeff,trace[plot_trace][i-ifr]*fTraceCoeff);	
	    } 
	    if(iDoSavePlotColor){
	      cpgsci(2);
	      cpgslw(save_lw+1);
	      for(i=ifr;i<=ffr;i++){
		cpgpt1((float)(i-iSavePlotBeginTimeOffset)*fTimeCoeff,trace[plot_trace][i-ifr]*fTraceCoeff,-1);
	      } 
	    }
	    
	    // Plot time average
	    /*
	    cpgsci(2);
	    cpgslw(2);
	    for(i=ifr;i<=ffr;i++){
	      cpgpt1((float)i,trace[current_trace][i-ifr],-1);
	    } 
	    */
	    if(iPlotAverageTraceOnSavePlot){
	      cpgslw(save_lw+1);
	      cpgsci(1);
	      cpgmove((float)(ifr-iSavePlotBeginTimeOffset)*fTimeCoeff,atrace[0/*current_trace*/][0]*fTraceCoeff);
	      for(i=ifr+1;i<=ffr;i++){
		cpgdraw((float)(i-iSavePlotBeginTimeOffset)*fTimeCoeff,atrace[0/*current_trace*/][i-ifr]*fTraceCoeff);	
	      } 
	    }
	    
	    // Plot fit
	    if(plot_fit){
	      /*
		// Linear fit
	      cpgsci(7);
	      cpgmove((float)ifr,(float)fit_a);
	      for(i=ifr+1;i<=ffr;i++){
		cpgdraw((float)i,(float)(fit_a+fit_b*(i-ifr)));	
	      } 
	      */
	      if(poly_fit){
		if(iDoSavePlotColor) cpgsci(4);
		PolynomialFitM(nframes,trace_poly_fitn,0.0,&dumd,poly_fit,ASSEMBLE);
		cpgmove((float)(ifr-iSavePlotBeginTimeOffset)*fTimeCoeff,(float)dumd*fTraceCoeff);
		for(i=ifr+1;i<=ffr;i++){
		  PolynomialFitM(nframes,trace_poly_fitn,(double)(i-ifr),&dumd,poly_fit,ASSEMBLE);
		  cpgdraw((float)(i-iSavePlotBeginTimeOffset)*fTimeCoeff,(float)(dumd)*fTraceCoeff);	
		}
	      } 
	    }
	  }

	  if(save_mode&2){
	    pwminy=traceminglobal;
	    pwmaxy=tracemaxglobal;
	    cpgsci(1);
	    cpgsch(CHARACTERSIZE1);
	    cpgslw(save_lw+1);

	    cpgenv(pwminx,pwmaxx,pwminy,pwmaxy,0,2);
	    cpgsch(CHARACTERSIZE2);
	    for(j=0;j<ntraces;j++){
	      cpgsci(3+j);
	      sprintf(str,"%i-%i ",tracepixel[j]%Xdim,tracepixel[j]/Xdim);
	      cpgtext(pwminx+(j+1)*(pwmaxx-pwminx)*0.08,pwmaxy+(pwmaxy-pwminy)*0.04,str);
	      cpgmove((float)ifr*fTimeCoeff,trace[j][0]);
	      for(i=ifr+1;i<=ffr;i++){
		cpgdraw((float)i*fTimeCoeff,trace[j][i-ifr]);	
	      } 
	    }
	    if(synch_trace){
	      cpgsci(1);
	      cpgmove((float)ifr*fTimeCoeff,pwminy+0.2*(pwmaxy-pwminy)*synch_trace[0]);
	      for(i=ifr+1;i<=ffr;i++){
		cpgdraw((float)i*fTimeCoeff,pwminy+0.2*(pwmaxy-pwminy)*synch_trace[i-ifr]);	
	      } 
	    }	    
	  }

	  //Save Fourie plot
	  if(save_mode&4){

	    pwminy=ftraceminglobal;
	    pwmaxy=ftracemaxglobal;
	    if(pwminy<fFourierThreshold) pwminy=fFourierThreshold;
	    pwminx=(float)cut_offl*fFrequencyCoeff;
	    pwmaxx=(float)cut_offu*fFrequencyCoeff;
	    cpgsci(1);
	    cpgsch(CHARACTERSIZE1);
	    
	    cpgslw(save_lw+1);	    	    

#ifdef POWER_SPECTRUM_LOG_PLOT    
	    cpgenv((float)log10((double)pwminx),(float)log10((double)pwmaxx),
		   pwminy,1.01*pwmaxy,0,30);
	    cpgsch(CHARACTERSIZE2);


	    // Plot fourie traces 
	    cpgslw(save_lw+1);	    	    
	    for(j=0;j<ntraces;j++){
	      if(iDoSavePlotColor) cpgsci(3+j);
	      cpgmove((float)log10((double)cut_offl*fFrequencyCoeff),ftracem[j][cut_offl]);
	      for(i=cut_offl+1;i<=cut_offu;i++){
		cpgdraw((float)log10((double)i*fFrequencyCoeff),ftracem[j][i]);	
	      } 
	    }
	    if(iDoSavePlotColor){
	      cpgsci(2);
	      cpgslw(save_lw+1);
	      for(i=cut_offl;i<=cut_offu;i++){
		cpgpt1((float)log10((double)i*fFrequencyCoeff),ftracem[plot_trace][i],-1);
	      } 
	    }

	    if(iPlotHPFourierTrace){
	      cpgslw(save_lw+1);	    	    
	      if(iDoSavePlotColor) cpgsci(4);
	      cpgmove((float)log10((double)cut_offl*fFrequencyCoeff),fatracem[cut_offl]);
	      for(i=cut_offl+1;i<=cut_offu;i++){
		cpgdraw((float)log10((double)i*fFrequencyCoeff),fatracem[i]);	
	      } 
	      if(iDoSavePlotColor){
		cpgsci(1);
		cpgslw(save_lw+1);
		for(i=cut_offl;i<=cut_offu;i++){
		  cpgpt1((float)log10((double)i*fFrequencyCoeff),fatracem[i],-1);
		}
	      }
	      /*
	      // Emphasize raw trace points 
	      cpgsci(2);
	      cpgslw(save_lw+1);
	      for(i=cut_offl;i<=cut_offu;i++){
		cpgpt1((float)log10((double)i*fFrequencyCoeff),ftracem[plot_trace][i],-1);
	      } 
	      */
	      cpgslw(save_lw+1);
	    }
#else
	    cpgenv(pwminx,pwmaxx,pwminy,pwmaxy,0,2);
	    cpgsch(CHARACTERSIZE2);
	    for(j=0;j<ntraces;j++){
	      cpgsci(3+j);
	      cpgmove((float)cut_offl*fFrequencyCoeff,ftracem[j][cut_offl]);
	      for(i=cut_offl+1;i<=cut_offu;i++){
		cpgdraw((float)i*fFrequencyCoeff,ftracem[j][i]);	
	      } 
	    }
	    cpgsci(2);
	    cpgslw(2);
	    for(i=cut_offl;i<=cut_offu;i++){
	      cpgpt1((float)i*fFrequencyCoeff,ftracem[plot_trace][i],-1);
	    } 
#endif

	  }
	  cpgclos();
	}
	else{
	  printf("TRA Save mode is not defined\n");
	}
      }
      else{
	printf("TRA Nothing to save\n");
      }
      continue;
    }

    if(ch=='S'){
      if(ntraces){
	plot_trace=(ntraces_count-1)%NTRACES;	
	
	sprintf(str,"%s.%i_%i_%i_%ix%i%s",tr->filename,n_cycles,initframe,finalframe,xi,yi,".trace.dat");	
	printf("TRA Saving raw trace in %s\n",str);
	
	if(!(fp=fopen(str,"w"))){
	  printf("TRA Cannot open raw trace file %s\n",str);
	}
	else{
	  for(i=ifr;i<=ffr;i++){
	    fprintf(fp,"%i\n",(int)(trace[plot_trace][i-ifr]));
	  } 
	}
	fclose(fp);
      }
      else{
	printf("TRA Nothing to save\n");
      }
      continue;
    }
    
    if(ch=='q'){ break;}
    if(ch=='Q'){ exit(0);}
    if(ch=='m'){
      i=0;
      while(1){
	cpgband(0,0,0.0,0.0,&dumfx,&dumfy,&ch); 
	if(ch!=13 && ch!=32){ 
	  str[i++]=ch;
	  printf("TRA %c",ch);
	  fflush(stdout);
	}	
	else{ 
	  str[i]='\0';
	  printf("%c",ch);
	  fflush(stdout);
	  break;
	}
      }
      x=atof(str);
      printf(" %s(%f) ",str,x);
      i=0;
      while(1){
	cpgband(0,0,0.0,0.0,&dumfx,&dumfy,&ch); 
	if(ch!=13 && ch!=32){
	  str[i++]=ch;
	  printf("TRA%c",ch);
	  fflush(stdout);
	}
	else{ 
	  str[i]='\0';
	  printf("%c",ch);
	  fflush(stdout);
	  break;
	}
      }
      printf("\n");
      y=atof(str);
      printf("TRA %s(%f) ",str,y);
      //	printf("TRA Input X Y: ");
      //	scanf("%f %f",&x,&y);
    }

    if(ch=='0'){
      ntraces=ntraces_count=current_trace=0;
      continue;
    }
    
    if(ch=='M'){
      if(GetString(str,1)) save_mode=atoi(str);
      continue;
    }

    if(ch=='F'){
      iPlotHPFourierTrace=1-iPlotHPFourierTrace;
      iOldTraces=1;
      goto plot;
    }

    
    xi=(int)rint((double)x);
    if(xi<0) xi=0;
    if(xi>Xdim-1) xi=Xdim-1;
    yi=(int)rint((double)y);
    if(yi<0) yi=0;
    if(yi>Ydim-1) yi=Ydim-1;
    printf("TRA X %i Y %i\n",xi,yi);
    
    l=xi+yi*Xdim;
    tracepixel[current_trace]=l;
    if(ntraces<NTRACES) ntraces++;

  plot:

    if(ntraces<=0) continue;

#ifdef TEST
    omegadum=2.0*M_PI/20.;//27837626523;
    alphadum=omegadum/50.0;
    //    dumf=TestFunction(ifr*omegadum,alphadum,0.0,0.0,0.0);
    //    p1=omegadum;
    p1=(float)current_trace;
    printf("%i %f\n",current_trace,p1);
    //    p1=6000.0;
    p2=5000.0;
    p3=50.0;
    p4=0.5;
    dumf=TestFunction((float)ifr,p1,p2,p3,p4);
#else
    
    if(synch_phi) phil=synch_phi[0]*harmonic;
    else phil=0.0;
    /*
    if(!(pcar=AddFrame(tr,ifr,a))){
      printf("TRA AddFrame return 0 for initframe %i\n",ifr);     
      return(2);
    }
    dumf=(float)(dumd=(double)pcar->pimage[l]);
    */


    // Accumulate time averaged
    //Val: Using  iIniAFrame and iFinAFrame instead of ifr and ffr (8-25-2002)
    //        iIniAFrame=ifr;
    //        iFinAFrame=ffr;
    iIniAFrame=initframe_ta;
    iFinAFrame=finalframe_ta;

    if(GetRecordUS(tr,iIniAFrame,l,&us)){
      printf("TRA GetRecordUS error for initframe %i\n",iIniAFrame);     
      return(2);
    }
    dumf=(float)(dumd=(double)us);
    accumulated_timeaverage=dumf;
    timeinterval_local=1.0;
    for(i=iIniAFrame+1;i<ifr+timeradius;i++){
      if(i<iFinAFrame){
	timeinterval_local+=1.0;
	/*
	if(!(pcar=AddFrame(tr,i,a))){
	  printf("TRA AddFrame return 0 at AT accumulation for frame %i\n",i);
	  return(2);
	}
	accumulated_timeaverage +=(float)(pcar->pimage[l]);
	*/
	if(GetRecordUS(tr,i,l,&us)){
	  printf("TRA GetRecordUS error at AT accumulation for frame %i\n",i);
	  return(2);
	}
	accumulated_timeaverage +=(float)us;
      }
    }
    
    //    timeinterval_local=(float)(timeradius);
#endif
   
    atracemaxglobal=-(atraceminglobal=BIG_FLOAT);
    trace[current_trace][0]=dumf;
    tracemean=dumd;
    tracemeanx=0.0;
    if(poly_fit){
      memset(poly_fit,0,trace_poly_fitn*sizeof(double));
      PolynomialFitM(nframes,trace_poly_fitn,0.0,&dumd,poly_fit,DISASSEMBLE);
    }
    tracemin[current_trace]=tracemax[current_trace]=dumf;

    resp=cos(phil)*dumd;
    imsp=sin(phil)*dumd;

    pwminy=dumf-2.0*(float)sqrt(dumd);
    pwmaxy=-pwminy+2.0*dumf;
    if(pwmaxy==pwminy){ 
      pwmaxy=1.0;
      pwminy=-0.1;
    }

    pwminx=(float)ifr;
    pwmaxx=(float)ffr;

    cpgslct(plotid);
    cpgpanl(1,3);
    cpgsci(1);
    cpgsch(CHARACTERSIZE1);
    cpgenv(pwminx,pwmaxx,pwminy,pwmaxy,0,2);
    sprintf(str,"X %i Y %i",xi,yi);
    cpgsch(CHARACTERSIZE2);
    cpglab("","",str);
    cpgsci(3);
    
    cpgmove((float)ifr,trace[current_trace][0]);
    
    for(i=ifr+1;i<=ffr;i++){

#ifdef TEST
      //      dumf=TestFunction(i*omegadum,alphadum,0.0,0.0,0.0);
      dumf=TestFunction((float)i,p1,p2,p3,p4);
#else
      /*      
      if(!(pcar=AddFrame(tr,i,a))){
	printf("TRA AddFrame return 0 for frame %i\n",i);     
	return(2);
      }
      dumf=(float)(dumd=(double)pcar->pimage[l]);
      */
      if(GetRecordUS(tr,i,l,&us)){
	printf("TRA GetRecordUS error for frame %i\n",i);     
	return(2);
      }
      dumf=(float)(dumd=(double)us);
      
      if(synch_phi) phil=synch_phi[i-ifr]*harmonic;
      else phil=harmonic*2.0*M_PI*n_cycles*(i-ifr)/(double)nframes;
#endif

      trace[current_trace][i-ifr]=dumf;
      tracemean+=dumd;
      tracemeanx+=dumd*(double)(i-ifr);      
      if(poly_fit) PolynomialFitM(nframes,trace_poly_fitn,(double)(i-ifr),&dumd,poly_fit,DISASSEMBLE);

      resp+=cos(phil)*dumd;
      imsp+=sin(phil)*dumd;

      if(dumf>tracemax[current_trace]) tracemax[current_trace]=dumf;
      if(dumf<tracemin[current_trace]) tracemin[current_trace]=dumf;
      cpgdraw((float)i,dumf);	

      if(i<iFinAFrame-timeradius){
	/*
	if(!(pcar=AddFrame(tr,i+timeradius,a))){
	  printf("TRA AddFrame return 0 at time+ for %i\n",i+timeradius);     
	  return(4);
	}
	accumulated_timeaverage+=(float)pcar->pimage[l];
	*/
	if(GetRecordUS(tr,i+timeradius,l,&us)){
	  printf("TRA GetRecordUS error at time+ for frame %i\n",i+timeradius);
	  return(4);
	}
	accumulated_timeaverage +=(float)us;
	timeinterval_local+=1.0;
      }
      if(i>iIniAFrame+timeradius){
	/*
	if(!(pcar=AddFrame(tr,i-timeradius,a))){
	  printf("TRA AddFrame return 0 at time- for %i\n",i-timeradius);     
	  return(5);
	}
	accumulated_timeaverage-=(float)pcar->pimage[l];
	*/
	if(GetRecordUS(tr,i-timeradius,l,&us)){
	  printf("TRA GetRecordUS error at time- for frame %i\n",i-timeradius);
	  return(5);
	}
	accumulated_timeaverage -=(float)us;
	timeinterval_local-=1.0;
      }

      dumf=atrace[0][i-ifr]=accumulated_timeaverage/timeinterval_local;
      // Val: Using atrace[1] for high-pass filtered signal (bad hack)
      dumf=atrace[1][i-ifr]=trace[current_trace][i-ifr]-atrace[0][i-ifr];
      if(dumf>atracemaxglobal) atracemaxglobal=dumf;
      if(dumf<atraceminglobal) atraceminglobal=dumf;

    }

    traceshift=(float)tracemean/(float)nframes;

    s=(double)nframes;
    sx=0.5*s*(double)(nframes-1);
    sxx=sx*(double)(2*nframes-1)/3.0;
    delta=sx*sx-s*sxx;
    fit_a=(sx*tracemeanx-sxx*tracemean)/delta;
    fit_b=(sx*tracemean-s*tracemeanx)/delta;
    // See note Val(7-26-02)
    pssp=10.0*(float)((resp*resp+imsp*imsp)/((double)nframes*(double)nframes));

    printf("TRA Trace %i min=%f max=%f mean=%f\n",current_trace,tracemin[current_trace],tracemax[current_trace],traceshift);
    printf("TRA %f Re=%f Im=%f P=%f A=%f\n",cyclesexact,resp/sqrt((double)nframes),imsp/sqrt((double)nframes),sqrt(pssp),atan2(imsp,resp));

    pwminy=tracemin[current_trace];
    pwmaxy=tracemax[current_trace];
    if(pwmaxy==pwminy){ 
      pwmaxy=1.0;
      pwminy=-0.1;
    }

    cpgsci(1);
    cpgpanl(1,3);
    cpgsch(CHARACTERSIZE1);
    cpgenv(pwminx,pwmaxx,pwminy,pwmaxy,0,2);
    sprintf(str,"X %i Y %i",xi,yi);
    cpgsch(CHARACTERSIZE2);
    cpglab("","",str);

    cpgsci(3);
    cpgmove((float)ifr,trace[current_trace][0]);
    for(i=ifr+1;i<=ffr;i++){
      cpgdraw((float)i,trace[current_trace][i-ifr]);	
    } 

    // Plot time average
    cpgsci(2);
    cpgslw(8);
    for(i=ifr;i<=ffr;i++){
      cpgpt1((float)i,trace[current_trace][i-ifr],-1);
    } 
    cpgslw(1);
    cpgsci(1);
    cpgmove((float)ifr,atrace[0/*current_trace*/][0]);
    for(i=ifr+1;i<=ffr;i++){
      cpgdraw((float)i,atrace[0/*current_trace*/][i-ifr]);	
    } 

    // Plot fit
    if(plot_fit){
      cpgsci(7);
      cpgmove((float)ifr,(float)fit_a);
      for(i=ifr+1;i<=ffr;i++){
	cpgdraw((float)i,(float)(fit_a+fit_b*(i-ifr)));	
      } 
      if(poly_fit){
	cpgsci(4);
	PolynomialFitM(nframes,trace_poly_fitn,0.0,&dumd,poly_fit,ASSEMBLE);
	cpgmove((float)ifr,(float)dumd);
	for(i=ifr+1;i<=ffr;i++){
	  PolynomialFitM(nframes,trace_poly_fitn,(double)(i-ifr),&dumd,poly_fit,ASSEMBLE);
	  cpgdraw((float)i,(float)(dumd));
	}
      } 
    }

    if(do_shift){
      tracemin[current_trace]-=traceshift;
      tracemax[current_trace]-=traceshift;
      for(i=ifr;i<=ffr;i++) trace[current_trace][i-ifr]-=traceshift;
    }

    // Fourier transforms

    if(cut_offl<0) cut_offl=0;
    if(cut_offu>nframes/2) cut_offu=nframes/2;
    if(cut_offl>cut_offu){
      printf("TRA Screwed up cut-offs: lower=%i upper=%i\n",cut_offl,cut_offu);
      cut_offl=1;
      cut_offu=nframes/2;
      printf("TRA Setting to max range: lower=%i upper=%i\n",cut_offl,cut_offu);
    }
    cut_offll=cut_offl;
    if(cut_offl==0)  cut_offll++;
    
    cut_offuu=cut_offu;
    if(cut_offu==nframes/2) cut_offuu--;

    if(do_fft){
      memcpy(ftrace[current_trace],trace[current_trace],nframes*sizeof(float));
      memset(ftrace[current_trace]+nframes,0,(tracefttn-nframes)*sizeof(float));
      realft(ftrace[current_trace]-1,tracefttn,1);
      if(iPlotHPFourierTrace){
	memcpy(fatrace,atrace[1],nframes*sizeof(float));
	memset(fatrace+nframes,0,(tracefttn-nframes)*sizeof(float));
	realft(fatrace-1,tracefttn,1);
      }
    }
    else{
      pfDum=trace[current_trace];
      re=im=0.0;
      for(i=0;i<nframes/2;i++){
	re+=(double)pfDum[2*i]+(double)pfDum[2*i+1];
	im+=(double)pfDum[2*i]-(double)pfDum[2*i+1];
      }
      ftrace[current_trace][0]=(float)re;
      ftrace[current_trace][1]=(float)im;

      CO=cos(2.0*M_PI/(double)nframes);
      SI=sin(2.0*M_PI/(double)nframes);
      CON=1.0;
      SIN=0.0;
      for(k=1;k<=cut_offuu;k++){
	co=dumd=CON*CO-SIN*SI;
	si=SIN=SIN*CO+CON*SI;
	CON=dumd;
	//	co=cos(2.0*M_PI*(double)k/(double)nframes);
	//	si=sin(2.0*M_PI*(double)k/(double)nframes);
	coN=1.0; 
	siN=0.0; 
	re=(double)pfDum[0];
	im=0.0;
	for(i=1;i<nframes;i++){
	  dumd=coN*co-siN*si;
	  siN=siN*co+coN*si;
	  coN=dumd;

	  dumd=(double)pfDum[i];
	  re+=dumd*coN;
	  im+=dumd*siN;
	}
	ftrace[current_trace][k*2]=(float)re;
	ftrace[current_trace][k*2+1]=(float)im;
      }
      if(iPlotHPFourierTrace){
	pfDum=atrace[1];
	re=im=0.0;
	for(i=0;i<nframes/2;i++){
	  re+=(double)pfDum[2*i]+(double)pfDum[2*i+1];
	  im+=(double)pfDum[2*i]-(double)pfDum[2*i+1];
	}
	fatrace[0]=(float)re;
	fatrace[1]=(float)im;
	
	CO=cos(2.0*M_PI/(double)nframes);
	SI=sin(2.0*M_PI/(double)nframes);
	CON=1.0;
	SIN=0.0;
	for(k=1;k<=cut_offuu;k++){
	  co=dumd=CON*CO-SIN*SI;
	  si=SIN=SIN*CO+CON*SI;
	  CON=dumd;
	  coN=1.0; 
	  siN=0.0; 
	  re=(double)pfDum[0];
	  im=0.0;
	  for(i=1;i<nframes;i++){
	    dumd=coN*co-siN*si;
	    siN=siN*co+coN*si;
	    coN=dumd;
	    
	    dumd=(double)pfDum[i];
	    re+=dumd*coN;
	    im+=dumd*siN;
	  }
	  fatrace[k*2]=(float)re;
	  fatrace[k*2+1]=(float)im;
	}
      }
    }

    /*
    for(k=0;k<9;k++){
      co=(double)ftrace[current_trace][k*2]/(double)nframes;
      si=(double)ftrace[current_trace][k*2+1]/(double)nframes;
      printf("TRA %i %f %f %f\n",k,co,si,hypot(co,si));
    }
    */

    traceminglobal=tracemin[0];
    tracemaxglobal=tracemax[0];
    for(i=1;i<ntraces;i++){
      if(tracemaxglobal<tracemax[i]) tracemaxglobal=tracemax[i];
      if(traceminglobal>tracemin[i]) traceminglobal=tracemin[i];
    }
    pwminy=traceminglobal;
    pwmaxy=tracemaxglobal;
    if(pwmaxy==pwminy){ 
      pwmaxy=1.0;
      pwminy=-0.1;
    }

    cpgpanl(1,1);
    cpgsch(CHARACTERSIZE1);
    cpgsci(1);
    cpgenv(pwminx,pwmaxx,pwminy,pwmaxy,0,2);
    cpgsch(CHARACTERSIZE2);
    for(j=0;j<ntraces;j++){
      cpgsci(4+j);
      sprintf(str,"%i-%i ",tracepixel[j]%Xdim,tracepixel[j]/Xdim);
      cpgtext(pwminx+(j+1)*(pwmaxx-pwminx)*0.08,pwmaxy+(pwmaxy-pwminy)*0.04,str);
      cpgmove((float)ifr,trace[j][0]);
      for(i=ifr+1;i<=ffr;i++){
	cpgdraw((float)i,trace[j][i-ifr]);	
      } 
    }
    if(synch_trace){
      cpgsci(1);
      cpgmove((float)ifr,pwminy+0.2*(pwmaxy-pwminy)*synch_trace[0]);
      for(i=ifr+1;i<=ffr;i++){
	cpgdraw((float)i,pwminy+0.2*(pwmaxy-pwminy)*synch_trace[i-ifr]);	
      } 
    }
   /*
    if(synch_phi){
      cpgmove((float)ifr,pwminy+0.2*(pwmaxy-pwminy)*(cos(synch_phi[0])+1.0));
      for(i=ifr+1;i<=ffr;i++){
	cpgdraw((float)i,pwminy+0.2*(pwmaxy-pwminy)*(cos(synch_phi[i-ifr])+1.0));	
      } 
    }
    */

#ifdef PLOT_PHASE
    ftracep[current_trace][0]=ftracep[current_trace][nframes/2]=0.0;
#endif
    //Val(7-26-02) Will normalize on nframes^2
    //    dNFrames2=1.0/(double)nframes;
    dNFrames2=1.0/((double)nframes*(double)nframes);

    pfDum=ftrace[current_trace];
    if(cut_offl==0){
      ftracem[current_trace][0]=(float)log10(((double)pfDum[0])*((double)pfDum[0])*dNFrames2);
    }
    if(cut_offu==nframes/2){
      ftracem[current_trace][nframes/2]=(float)log10(((double)pfDum[1])*((double)pfDum[1])*dNFrames2);
    }

    for(i=cut_offll;i<=cut_offuu;i++){ 
      ftracem[current_trace][i]=(float)log10((pfDum[2*i]*pfDum[2*i]+pfDum[2*i+1]*pfDum[2*i+1])*dNFrames2);

#ifdef PLOT_PHASE
      ftracep[current_trace][i]=(float)atan2(pfDum[2*i+1],pfDum[2*i]);
#endif
    }

    if(iPlotHPFourierTrace){
      pfDum=fatrace;
      if(cut_offl==0){
	fatracem[0]=(float)log10(((double)pfDum[0])*((double)pfDum[0])*dNFrames2);
      }
      if(cut_offu==nframes/2){
	fatracem[nframes/2]=(float)log10(((double)pfDum[1])*((double)pfDum[1])*dNFrames2);
      }
      
      for(i=cut_offll;i<=cut_offuu;i++){ 
	fatracem[i]=(float)log10((pfDum[2*i]*pfDum[2*i]+pfDum[2*i+1]*pfDum[2*i+1])*dNFrames2);
      }
    }

    k=(int)(cycleslower);
    if(k<=cut_offu && k>=cut_offl) printf("TRA %i Re=%f Im=%f P=%f\n",k,ftrace[current_trace][2*k]/sqrt((double)nframes),ftrace[current_trace][2*k+1]/sqrt((double)nframes),exp(0.5*log(10.0)*(double)ftracem[current_trace][k]));
    k=(int)(cyclesupper);
    if(k<=cut_offu && k>=cut_offl) printf("TRA %i Re=%f Im=%f P=%f\n",k,ftrace[current_trace][2*k]/sqrt((double)nframes),ftrace[current_trace][2*k+1]/sqrt((double)nframes),exp(0.5*log(10.0)*(double)ftracem[current_trace][k]));

    ftracemaxm[current_trace]=ftraceminm[current_trace]=ftracem[current_trace][cut_offl];
    for(i=cut_offl+1;i<=cut_offu;i++){ 
      dumf=ftracem[current_trace][i];
      //      printf("TRA %i %f\n",i,dumf);
      if(dumf>ftracemaxm[current_trace]) ftracemaxm[current_trace]=dumf;
      if(dumf<ftraceminm[current_trace]) ftraceminm[current_trace]=dumf;
    }

    ftraceminglobal=ftraceminm[0];
    ftracemaxglobal=ftracemaxm[0];
    //    ftracemaxglobal=log10(pssp);
    for(i=0;i<ntraces;i++){
      if(ftracemaxglobal<ftracemaxm[i]) ftracemaxglobal=ftracemaxm[i];
      if(ftraceminglobal>ftraceminm[i]) ftraceminglobal=ftraceminm[i];
    }
    pwminy=ftraceminglobal;
    pwmaxy=ftracemaxglobal;
    if(pwminy<fFourierThreshold) pwminy=fFourierThreshold;
    if(pwmaxy==pwminy){ 
      pwmaxy=1.0;
      pwminy=-0.1;
    }

    pwminx=(float)cut_offl;
    pwmaxx=(float)(cut_offu-1);

    cpgpanl(1,2);
    cpgsci(1);
    cpgsch(CHARACTERSIZE1);

#ifdef POWER_SPECTRUM_LOG_PLOT    
    cpgenv((float)log10((double)pwminx),(float)log10((double)pwmaxx),pwminy,1.01*pwmaxy,0,30);
    cpgsch(CHARACTERSIZE2);

    // Plot spike @ exact frequency
    if(pssp!=0.0){
      cpgsci(1);
      cpgmove((float)log10(cyclesexact),pwminy);
      cpgdraw((float)log10(cyclesexact),log10(pssp));
      cpgsci(2);
      cpgslw(14);
      cpgpt1((float)log10(cyclesexact),log10(pssp),-1);
      cpgslw(1);
    }

    // Plot fourie traces 
    // PLot HP filtered first 
    if(iPlotHPFourierTrace){
      cpgslw(1);
      cpgsci(3);
      cpgmove((float)log10((double)cut_offl),fatracem[cut_offl]);
      for(i=cut_offl+1;i<=cut_offu;i++){
	cpgdraw((float)log10((double)i),fatracem[i]);	
      } 
      cpgsci(1);
      cpgslw(8);
      for(i=cut_offl;i<=cut_offu;i++){
	cpgpt1((float)log10((double)i),fatracem[i],-1);
      } 
      cpgslw(1);
    }

    for(j=0;j<ntraces;j++){
      cpgsci(4+j);
      cpgmove((float)log10((double)cut_offl),ftracem[j][cut_offl]);
      for(i=cut_offl+1;i<=cut_offu;i++){
	cpgdraw((float)log10((double)i),ftracem[j][i]);	
      } 
    }
    cpgsci(2);
    cpgslw(8);
    for(i=cut_offl;i<=cut_offu;i++){
      cpgpt1((float)log10((double)i),ftracem[current_trace][i],-1);
    } 

    cpgslw(1);
    cpgsci(2);
    cpgmove((float)log10((double)cut_offl),(float)log10(tracemean));
    cpgdraw((float)log10((double)cut_offu),(float)log10(tracemean));

    //    cpgsci(7);
    //    cpgmove((float)log10((double)cut_offl),1);
    //    cpgdraw((float)log10((double)cut_offu),1);

#else
    cpgenv(pwminx,pwmaxx,pwminy,1.01*pwmaxy,0,2);
    cpgsch(CHARACTERSIZE2);

    // Plot spike @ exact frequency
    if(pssp!=0.0){
      cpgsci(1);
      cpgmove((float)cyclesexact,pwminy);
      cpgdraw((float)cyclesexact,log10(pssp));
      cpgsci(2);
      cpgslw(14);
      cpgpt1((float)cyclesexact,log10(pssp),-1);
      cpgslw(1);
    }

    // Plot fourie traces 
    for(j=0;j<ntraces;j++){
      cpgsci(4+j);
      cpgmove((float)cut_offl,ftracem[j][cut_offl]);
      for(i=cut_offl+1;i<=cut_offu;i++){
	cpgdraw((float)i,ftracem[j][i]);	
      } 
    }
    cpgsci(2);
    cpgslw(8);
    for(i=cut_offl;i<=cut_offu;i++){
      cpgpt1((float)i,ftracem[current_trace][i],-1);
    } 

    if(iPlotHPFourierTrace){
      cpgsci(3);
      cpgmove((float)cut_offl,fatracem[cut_offl]);
      for(i=cut_offl+1;i<=cut_offu;i++){
	cpgdraw((float)i,fatracem[i]);	
      } 
      cpgsci(1);
      cpgslw(8);
      for(i=cut_offl;i<=cut_offu;i++){
	cpgpt1((float)i,fatracem[i],-1);
      } 
    }

    cpgslw(1);
    cpgsci(2);
    cpgmove((float)cut_offl,(float)log10(tracemean));
    cpgdraw((float)cut_offu,(float)log10(tracemean));

#endif


#ifdef PLOT_PHASE
    cpgslct(phaseid); 
#ifdef PLOT_DIFFERENCE

    //Quick hack, it's supposed to plot phase
    //Plot difference of trace and time averaged trace

    cpgsci(1);
    cpgsch(CHARACTERSIZE1);
    if(atraceminglobal==atracemaxglobal) atraceminglobal=-(atracemaxglobal=1.0);
    cpgenv((float)ifr,(float)ffr,atraceminglobal,atracemaxglobal,0,2);
    printf("%i %f %f\n",current_trace,atraceminglobal,atracemaxglobal);
    cpgsci(3);
    cpgmove((float)ifr,atrace[1][0]);
    for(i=ifr+1;i<=ffr;i++){
      cpgdraw((float)i,atrace[1][i-ifr]);	
    } 
    if(poly_fit){
      cpgsci(4);
      PolynomialFitM(nframes,trace_poly_fitn,0.0,&dumd,poly_fit,ASSEMBLE);
      cpgmove((float)ifr,trace[current_trace][0]-(float)dumd);
      for(i=ifr+1;i<=ffr;i++){
	PolynomialFitM(nframes,trace_poly_fitn,(double)(i-ifr),&dumd,poly_fit,ASSEMBLE);
	cpgdraw((float)i,trace[current_trace][i-ifr]-(float)(dumd)); 
      }
    }
#else
    cpgsci(1);
    cpgsch(CHARACTERSIZE1);
    cpgenv(pwminx,pwmaxx,-M_PI,M_PI,0,2);
    cpgsch(CHARACTERSIZE2);
    for(j=0;j<ntraces;j++){
      cpgsci(4+j);
      cpgmove((float)cut_offl,ftracem[j][cut_offl]);
      for(i=cut_offl+1;i<=cut_offu;i++){
	cpgdraw((float)i,ftracep[j][i]);	
      } 
    }
#endif 
#endif    
    if(!iOldTraces){
      ntraces_count++;
      current_trace=ntraces_count%NTRACES;
    }
    else iOldTraces=0;
  }

  cpgslct(plotid);
  cpgclos();

#ifdef PLOT_PHASE
  cpgslct(phaseid);
  cpgclos();
#endif

  cpgslct(interactiveid);
  cpgclos();

  if(iPlotHPFourierTrace){
    if(fatrace) free(fatrace);
    if(fatracem) free(fatracem);
  }
  free(poly_fit);
  if(synch_trace) free(synch_trace);
  for(i=0;i<NTRACES;i++){
    free(trace[i]);
    free(atrace[i]);
    free(ftrace[i]);
    free(ftracem[i]);
#ifdef PLOT_PHASE
    free(ftracep[i]);
#endif
  }
  return(0);
}

