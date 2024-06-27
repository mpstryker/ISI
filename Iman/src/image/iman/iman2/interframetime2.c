#include "iman1.h"

#define BIG_LONG 4294967296.0
#define PLOT_WINDOW_SIZEX 12.0
#define PLOT_WINDOW_RATIO 0.8
#define CHARACTERSIZE1 1.1
#define CHARACTERSIZE2 1.5

#define SEEK_CHUNK 0.1

double ConvertI64(I64 *pi64){
  return((double)pi64->LowPart+BIG_LONG*(double)pi64->HighPart);
}

double ConvertULs(unsigned long ulLowPart,unsigned long ulHighPart){
  return((double)ulLowPart+BIG_LONG*(double)ulHighPart);
}

int InterframeTime(TRAIN *tr,int ifr,int ffr,int *iifr,int *fffr,int a){
  int i=0,k;
  int nframes;
  double *frame_intervals;
  double average_frame_interval,nframes_cycle,residue,residue_min;
  double frequency=0.0,tc=0.0,tp=0.0,dumd,average_ti,average_ti2,sigma;
  int interframeid;
  CAR *pcar;
  float pwminx,pwmaxx,pwminy,pwmaxy,dumf;
  char str[64];
  int n_cycles_min,n_cycles_best;
  int error_status=0;
  int orig_ifr,orig_ffr;
  int iVersion;

  iVersion=tr->iGlobalVersion;

  nframes=ffr-ifr+1;
  if(!(frame_intervals=(double *)calloc(nframes-1,sizeof(double)))){
    printf("INT Cannot allocate for frame_intervals (%i)\n",nframes-1);
    return(1);
  }
  if(!(pcar=AddFrame(tr,ifr,a))){
    printf("INT AddFrame return 0 (%i)\n",ifr);     
    return(2);
  }

  orig_ifr=ifr;
  orig_ffr=ffr;

  if(iVersion==VERSION_CHUNK){
    tp=1.0e-6*ConvertULs(((FRAM_CHUNK*)pcar->pFrameHeader)->TimeArrivalUsecLo,
			 ((FRAM_CHUNK*)pcar->pFrameHeader)->TimeArrivalUsecHi);
  }
  else{
    frequency=ConvertI64(&(((FRAMEHEADER*)pcar->pFrameHeader)->perfreq));
    tp=ConvertI64(&(((FRAMEHEADER*)pcar->pFrameHeader)->perfcount))/frequency;
  }

  printf("INT |");
  for(i=ifr+1,k=1;i<=ffr;i++){      
    if(!(pcar=AddFrame(tr,i,a))){
      printf("INT AddFrame return 0 (%i)\n",i);     
      return(2);
    }
    if(iVersion==VERSION_CHUNK){
      tc=1.0e-6*ConvertULs(((FRAM_CHUNK*)pcar->pFrameHeader)->TimeArrivalUsecLo,
			   ((FRAM_CHUNK*)pcar->pFrameHeader)->TimeArrivalUsecHi);
    }
    else{
      tc=ConvertI64(&(((FRAMEHEADER*)pcar->pFrameHeader)->perfcount))/frequency;
    }
    frame_intervals[i-ifr-1]=tc-tp;
    tp=tc;

    k++;
    if(!(k%100)){ 
      if(!(k%1000)){ 
	printf("|");
	fflush(stdout);
      }
      else printf(".");
      fflush(stdout);
    }
  }
  printf("\n");

  pwminx=(float)ifr;
  pwmaxx=(float)(ffr-1);
  pwminy=pwmaxy=(float)(average_ti=frame_intervals[0]);
  average_ti2=average_ti*average_ti;
  for(i=1;i<nframes-1;i++){
    average_ti+=(dumd=frame_intervals[i]);
    average_ti2+=dumd*dumd;
    dumf=(float)dumd;
    if(dumf>pwmaxy) pwmaxy=dumf;
    if(dumf<pwminy) pwminy=dumf;      
  }
  average_frame_interval=average_ti/(double)(nframes-1);
  sigma=sqrt(average_ti2/(double)(nframes-1)-average_frame_interval*average_frame_interval);
  printf("INT Interframe intervals %e(%e)\n",average_frame_interval,sigma);
  nframes_cycle=stimulus_period/average_frame_interval;
  cycles_per_frame=harmonic/nframes_cycle;

  if(do_intervals%2){
    if(synch_phi){
      frame_omega=2.0*M_PI/(nframes_cycle*n_cycles);
      printf("INT b=%e(%e) f=%e\n",1.0/nframes_cycle,sigma/stimulus_period,frame_omega/(2.0*M_PI));
      printf("INT NF=%f\n",2.0*M_PI/frame_omega);
    }
    else{
      n_cycles=(int)floor(nframes/nframes_cycle);
      if(n_cycles<1){
	printf("INT Cannot grab cycles. Too few frames\n");
	return(3);
      }
      nframes=(int)rint(n_cycles*nframes_cycle);
      residue=rint(n_cycles*nframes_cycle)-n_cycles*nframes_cycle;
      printf("INT      NC=%i NF=%i Res=%f\n",n_cycles,nframes,residue);

      if(find_best_cycle){
	residue_min=residue;
	n_cycles_best=n_cycles;
	
	n_cycles_min=(int)((1.0-SEEK_CHUNK)*(double)n_cycles);
	for(i=n_cycles_min;i<n_cycles;i++){
	  dumd=(double)i*nframes_cycle;
	  printf("INT NC %i Res %f\n",i,rint(dumd)-dumd);
	  if(fabs(rint(dumd)-dumd)<fabs(residue_min)){
	    n_cycles_best=i;
	    residue_min=rint(dumd)-dumd;
	  }
	}
	n_cycles=n_cycles_best;
	nframes=(int)rint(n_cycles*nframes_cycle);
	residue=residue_min;
	printf("INT Best NC=%i NF=%i Res=%f\n",n_cycles,nframes,residue);
      }

      if(do_precise){
	frame_omega=2.0*M_PI/(nframes_cycle*n_cycles);
      }
      else{
	frame_omega=2.0*M_PI/(double)(nframes);
      }

      printf("INT Fixing final frame from %i to %i\n",ffr,nframes+ifr-1);
      ffr=nframes+ifr-1;
      *fffr=ffr;
      initframe_ta=*iifr;
      finalframe_ta=*fffr;

      if(do_timeaveraging){
	timeradius=(int)rint(timeaverage_cycles*nframes_cycle);
	if(InitializeTimeAveraging(tr)){ 
	  printf("INT Time averaging initialization failed. Will try curve fitting\n");
	}
	else{
	  if(sacrifice_boundaries){
	    if(n_cycles>2*(int)(timeaverage_cycles)){
	      printf("INT Will chop %i cycles(%i frames) from each side\n",(int)(timeaverage_cycles),timeradius);
	      n_cycles-=2*(int)(timeaverage_cycles);
	      (*iifr)+=timeradius;
	      nframes=(int)rint((double)n_cycles*nframes_cycle);
	      *fffr=nframes+(*iifr)-1;
	      
	      if(do_precise){
		frame_omega=2.0*M_PI/(nframes_cycle*n_cycles);
	      }
	      else{
		frame_omega=2.0*M_PI/(double)(nframes);
	      }
	      if(synch_phi){
		memmove(synch_phi,synch_phi+timeradius,nframes*sizeof(double));
		synch_phi=realloc((void *)synch_phi,nframes*sizeof(double));	  
	      }
	      printf("INT NFrames=%i NCycles=%i\n",nframes,n_cycles);
	    }
	    else{
	      sacrifice_boundaries=0;
	      printf("INT Too few cycles. No boundary chopping.\n");
	    }
	  }
	}
      }

      if(do_remove_curve || do_remove_linear){
	if(InitializeFitting(tr)){
	  printf("INT Poly fit initialization failed.\n");
	  error_status=1;
	}
      }
     
      if(!(synch_phi=(double *)calloc((nframes),sizeof(double)))){
	printf("INT Cannot allocate for synch_phi (%luB)\n",(unsigned long)nframes*sizeof(double));
	error_status|=2;
      }
      else{
	dumd=2.0*M_PI/nframes_cycle;
	for(i=0;i<nframes;i++){
	  synch_phi[i]=dumd*(double)i; 
	}
      }
    }
  }

  if((do_intervals&2)){
    if((interframeid=cpgopen("/xw"))<1){
      printf("INT Cannot open device /xw\n");
      error_status|=4;
    }
    else{
      cpgsubp(1,1);
      cpgpap(PLOT_WINDOW_SIZEX,PLOT_WINDOW_RATIO);
      cpgsci(1);
      cpgsch(CHARACTERSIZE1);
      cpgenv(pwminx,pwmaxx,pwminy,pwmaxy,0,2);
      sprintf(str,"Interframe intervals");
      cpgsch(CHARACTERSIZE2);
      cpglab("","",str);
      
      cpgsci(3);
      cpgmove((float)orig_ifr,frame_intervals[0]);
      for(i=1;i<orig_ffr-orig_ifr;i++){
	cpgdraw((float)(i+orig_ifr),frame_intervals[i]);	
      } 
      cpgsci(2);
      cpgslw(8);
      for(i=0;i<orig_ffr-orig_ifr;i++){
	cpgpt1((float)(i+orig_ifr),frame_intervals[i],-1);
      } 
      cpgslw(1);
      cpgask(0);
      cpgband(0,0,0.0,0.0,&dumf,&dumf,str); 
      cpgclos();
    }
  }
  free(frame_intervals);
  
  return(error_status);
}
