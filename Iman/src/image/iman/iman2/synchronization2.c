#include "iman1.h"

//#define DEBUG_SYNCHRONIZATION

//#define UPPERLOWERMAPS

#define THRESHOLD 0.5
#define SEEK_CHUNK 0.1

#define TRY_SALVAGE_OUTER_FRAMES_NOT 0
#define TRY_SALVAGE_OUTER_FRAMES_PARTIAL 1
#define TRY_SALVAGE_OUTER_FRAMES_COMPLETE 2

//#define TRY_SALVAGE_COMPLETE

int Synchronization(TRAIN *tr,int ifr,int ffr,int *iifr,int *fffr,int a){
  int i,j,k,iDum;
  int fi,ff=0;
  int nframes;
  int *synch_tracei=NULL;
  CAR *pcar;
  int threshold;
  double S,Sx,Sy,Sxx,Sxy,DD,aa,bb,aad,bbd,aau,bbu;
  double sigma,sigmaa,sigmab,sigmaad,sigmabd,sigmaau,sigmabu;
  double dNframesynch,residue,residued,residueu;
  int iReturn=0;
  double dSynchMax=0.0;
  ULONG ulSynchMax=0,ul;
  int iCyclesOrig,iCyclesLeft,iCyclesRight;
  int iTryToSalvageOuterFrames=TRY_SALVAGE_OUTER_FRAMES_COMPLETE;
#ifdef UPPERLOWERMAPS
  double rel,iml,reu,imu,norm;
  double *maplx=NULL,*maply=NULL,*mapux=NULL,*mapuy=NULL;
  double cou,siu,col,sil,coL,siL,coU,siU;
  double dumd;
#endif
  FRAM_COST_CHUNK *pFCC;
  ULONG ulSynchChannelOffset;
  char str[16];

#ifdef DEBUG_SYNCHRONIZATION
  printf("SYN DEBUG ifr=%i ffr=%i\n",ifr,ffr);
#endif

  nframes=ffr-ifr+1;
  if(!(synch_tracei=(int *)calloc(nframes,sizeof(int)))){
    printf("SYN Cannot allocate for synch_tracei (%luB)\n",(unsigned long)nframes*sizeof(int));
    return(1);
  }
  if(!(synch_trace=(float *)calloc(nframes,sizeof(float)))){
    printf("SYN Cannot allocate for synch_trace (%luB)\n",(unsigned long)nframes*sizeof(float));
    iReturn=3;
    goto bailout;
  }


  if(tr->iGlobalVersion==VERSION_CHUNK){
    if(tr->iGlobalExperiment==EXPERIMENT_CONTINUOUS_MODE_CONTINUOUS){

      if(tr->pulSynchChannelMax[iSynchChannel]) dSynchMax=tr->pulSynchChannelMax[iSynchChannel];
      else dSynchMax=4294967296.0;

      if(!(pcar=AddFrame(tr,ifr,a))){
	printf("SYN Ini AddFrame return 0 (%i)\n",ifr);     
	iReturn=4;
	goto bailout;
      }

      pFCC=(FRAM_COST_CHUNK*)(pcar->pFrameHeader+CHUNK_HEAD_SIZE+*(ULONG*)(pcar->pFrameHeader+CHUNK_ID_SIZE));
      ulSynchChannelOffset=(void*)&(pFCC->SynchChannel[iSynchChannel]) - pcar->pFrameHeader;
      printf("SYN Synch channel offset %lu\n",ulSynchChannelOffset);

      printf("SYN |");
      for(i=ifr,k=1;i<=ffr;i++,k++){

	/*
	if(!(pcar=AddFrame(tr,i,a))){
	  printf("SYN AddFrame return 0 (%i)\n",i);     
	  iReturn=5;
	  goto bailout;
	}
	synch_trace[i-ifr]=(synch_tracei[i-ifr]=ul=*(ULONG*)(pcar->pFrameHeader+ulSynchChannelOffset))/dSynchMax;
	*/

	if((j=GetRecordUL(tr,i,(int)ulSynchChannelOffset,&ul))){
	  printf("SYN GetRecordUL return %i\n",j);     
	  iReturn=5;
	  goto bailout;
	}
	synch_tracei[i-ifr]=ul;
	if(ulSynchMax<ul) ulSynchMax=ul;

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
    }
    if(tr->iGlobalExperiment==EXPERIMENT_CONTINUOUS_MODE_EPISODIC){
      printf("SYN EPST experiment uses other means of synchronization\n");
      printf("SYN Use option -N with proper FRAME_RULE\n");
      iReturn=-1;
      goto bailout;
    }
  }
  else{
    dSynchMax=DEFAULT_SYNCH_MAX;
    printf("SYN |");
    for(i=ifr,k=1;i<=ffr;i++,k++){
      if(!(pcar=AddFrame(tr,i,a))){
	printf("SYN AddFrame return 0 (%i)\n",i);     
	iReturn=5;
	goto bailout;
      }
      ul=synch_tracei[i-ifr]=((FRAMEHEADER*)pcar->pFrameHeader)->synch_in>>2;
      if(ulSynchMax<ul) ulSynchMax=ul;

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
  }
  /*
  if(ulSynchMax+1!=(unsigned long)dSynchMax){
    printf("SYN WARNING Synch Max=%lu, expected %lu\n",ulSynchMax,(unsigned long)(dSynchMax)-1);
  }
  */
  if(ulSynchMax+1>(unsigned long)dSynchMax){
    printf("SYN WARNING Synch Max=%lu is larger than absolute max %lu\n",ulSynchMax,(unsigned long)(dSynchMax)-1);
    printf("SYN Continue (y/n)?");
    fgets(str,4,stdin);
    if(*str != 'y'){
      iReturn=55;
      goto bailout;
    }
  }
  
  threshold=(int)(THRESHOLD*dSynchMax);
  j=synch_tracei[0];
  k=synch_tracei[1];
  i=1;
  while(j==k){
    if(i==nframes-1){
      printf("SYN Synchronization failed, too few frames (1)\n");
      iReturn=7;
      goto bailout;
    }
    k=synch_tracei[++i];
  }

  //Val: got rid of the direction business (8Aug2002)
  if(abs(k-j)>threshold){
    if(k>j){
      iDum=(int)(dSynchMax)-1;
      for(i=0;i<nframes;i++){
	synch_tracei[i]=iDum-synch_tracei[i];
      }
    }
  }
  else{
    if(k<j){
      iDum=(int)(dSynchMax)-1;
      for(i=0;i<nframes;i++){
	synch_tracei[i]=iDum-synch_tracei[i];
      }
    }
  }
  
  for(i=0;i<nframes;i++){
    synch_trace[i]=synch_tracei[i]/dSynchMax;
  }

  j=synch_tracei[0];
  k=synch_tracei[1];
  i=1;
  while(j<=k){
    if(i==nframes-1){
      printf("SYN Synchronization failed, too few frames (2)\n");
      iReturn=9;
      goto bailout;
    }
    j=k;
    k=synch_tracei[++i];
  }
  fi=i;
  
  for(k=0,i=fi+1;i<nframes;i++){
    if(synch_tracei[i]<synch_tracei[i-1]){ 
      ff=i-1;
      k++;
    }
  }
  if(k==0){
    printf("SYN Synchronization failed, too few frames (3)\n");
    iReturn=11;
    goto bailout;
  }

  n_cycles=k;
  ff+=additionalframes;
  nframes=ff-fi+1;

  printf("SYN Grabbed %i cycles IF=%i(%i) FF=%i(%i) NF=%i\n",n_cycles,fi+ifr,fi,ff+ifr,ff,nframes);

  if(!(synch_phi=(double *)calloc((nframes+2),sizeof(double)))){
    printf("SYN Cannot allocate for synch_phi (%luB)\n",(unsigned long)(nframes+2)*sizeof(double));
    iReturn=13;
    goto bailout;
  }
 
  synch_phi[0]=(double)(synch_tracei[fi])/dSynchMax;
  for(i=1,k=0;i<nframes;i++){
    if(synch_tracei[fi+i]<synch_tracei[fi+i-1]){
      k++;
    }
    synch_phi[i]=(double)(synch_tracei[fi+i]+k*(int)(dSynchMax))/dSynchMax;
  }

  sigma=0.5/dSynchMax;
  S=(double)nframes;
  Sx=0.5*S*(S-1.0);
  Sxx=Sx*(2.0*S-1.0)/3.0;
  Sy=Sxy=0.0;
  for(i=0;i<nframes;i++){
    Sy+=synch_phi[i];
    Sxy+=synch_phi[i]*(double)i;
  }
  DD=S*Sxx-Sx*Sx;
  aad=(Sxx*Sy-Sx*Sxy)/DD;
  bbd=(S*Sxy-Sx*Sy)/DD;
  sigmaad=sigma*sqrt(Sxx/DD);
  sigmabd=sigma*sqrt(S/DD);
  residued=(double)nframes-((double)((int)n_cycles))/bbd;

  S+=1.0;
  Sx+=(double)i;
  Sxx+=(double)i*(double)i;
  Sy+=synch_phi[i];
  Sxy+=synch_phi[i]*(double)i;
  DD=S*Sxx-Sx*Sx;
  aau=(Sxx*Sy-Sx*Sxy)/DD;
  bbu=(S*Sxy-Sx*Sy)/DD;
  sigmaau=sigma*sqrt(Sxx/DD);
  sigmabu=sigma*sqrt(S/DD);
  residueu=(double)(nframes+1)-((double)((int)n_cycles))/bbu;

  printf("SYN Residue D=%f(b=%e) U=%f(b=%e)\n",residued,bbd,residueu,bbu);

  j=0;
//  nframes++;
  if(fabs(residued)<fabs(residueu)){
    bb=bbd;
    aa=aad;
    sigmaa=sigmaad;
    sigmab=sigmabd;
  }
  else{
    bb=bbu;
    aa=aau;
    sigmaa=sigmaau;
    sigmab=sigmabu;
    nframes++;
    j++;
  }

  cycles_per_frame=fabs(bb)*harmonic;
  cyclesexact=nframes*(bb)*harmonic;
  cycleslower=floor(cyclesexact);
  cyclesupper=cycleslower+1.0;
  dNframesynch=(double)n_cycles/bb;
  residue=(double)nframes-dNframesynch;

  if(do_precise){
    synch_omega=2.0*M_PI*bb/(double)n_cycles;
    synch_phase_shift=2.0*M_PI*aa;
  }
  else{
    synch_omega=2.0*M_PI/(double)nframes;
    synch_phase_shift=0.0;
  }

  printf("SYN a=%e(%e) b=%e(%e)\n",aa,sigmaa,bb,sigmab);
  printf("SYN f=%e Exact: Cycles=%f Frames/cycle=%f\n",synch_omega/(2.0*M_PI),cyclesexact,1.0/fabs(bb));
  printf("SYN NF=%i NFC=%f Residue=%f(%i)\n",nframes,dNframesynch,residue,j);
  for(i=0;i<nframes;i++){
    synch_phi[i]=2.0*M_PI*(fabs(bb)*(double)i); // Phase shift (aa) removed
    //    synch_phi[i]=2.0*M_PI*(aa+bb*(double)i);
  }

  memmove(synch_trace,synch_trace+fi,nframes*sizeof(float));
  synch_trace=realloc((void *)synch_trace,nframes*sizeof(float));

  fi+=ifr;
  ff+=ifr;
  initframe_ta=*iifr=fi;
  finalframe_ta=*fffr=ff;
 
  printf("SYN IFr=%i FFr=%i FCycles=%i NFramesC=%f\n",*iifr,*fffr,n_cycles,((double)((int)n_cycles))/fabs(bb));

  if(do_timeaveraging){
    timeradius=(int)rint(timeaverage_cycles/fabs(bb));
    if(InitializeTimeAveraging(tr)){ 
      printf("SYN Time averaging initialization failed. Will try curve fitting\n");
    }
    else{
      if(sacrifice_boundaries){
	if(n_cycles>2*(int)(ceil(timeaverage_cycles))){
	  iCyclesOrig=iCyclesRight=iCyclesLeft=(int)ceil(timeaverage_cycles);
	  
	  //Left border
	  switch(iTryToSalvageOuterFrames){
#ifdef TRY_SALVAGE_COMPLETE
	  case TRY_SALVAGE_OUTER_FRAMES_COMPLETE:
	    if(fi>timeradius+1){
	      printf("SYN LEFT COMPLETE\n");
	      iCyclesLeft=0;
	      break;
	    }
#endif
	  case TRY_SALVAGE_OUTER_FRAMES_PARTIAL:
	    if((int)ceil(timeaverage_cycles)>(int)floor(timeaverage_cycles)){
	      if(fi>(int)rint((timeaverage_cycles-floor(timeaverage_cycles))/fabs(bb))+1){
		printf("SYN LEFT PARTIAL\n");
		iCyclesLeft--;
		break;
	      }
	    }
	  case TRY_SALVAGE_OUTER_FRAMES_NOT:
	    printf("SYN LEFT NOT\n");
	    break;
	  default:
	    break;
	  }
	  
	  //Right border
	  switch(iTryToSalvageOuterFrames){
#ifdef TRY_SALVAGE_COMPLETE
	  case TRY_SALVAGE_OUTER_FRAMES_COMPLETE:
	    if(tr->max_n-ff>timeradius+1){
	      printf("SYN RIGHT COMPLETE\n");
	      iCyclesRight=0;
	      break;
	    }
#endif
	  case TRY_SALVAGE_OUTER_FRAMES_PARTIAL:
	    if((int)ceil(timeaverage_cycles)>(int)floor(timeaverage_cycles)){
	      if(tr->max_n-ff>(int)rint((timeaverage_cycles-floor(timeaverage_cycles))/fabs(bb))+1){
		printf("SYN RIGHT PARTIAL\n");
		iCyclesRight--;
		break;
	      }
	    }
	  case TRY_SALVAGE_OUTER_FRAMES_NOT:
	    printf("SYN RIGHT NOT\n");
	    break;
	    default:
	    break;
	  }
	  
	  printf("SYN Cycles: Orig %i Chopping Left %i, Right %i (frames %i)\n",
		 iCyclesOrig,iCyclesLeft,iCyclesRight,timeradius);
	  //	  n_cycles -= (int)ceil(2.0*timeaverage_cycles);
	  if(n_cycles<=iCyclesLeft+iCyclesRight){
	    sacrifice_boundaries=0;
	    printf("SYN WARNING: n_cycles=%i\n",n_cycles);
	      goto bailout;
	  }
	  n_cycles -= iCyclesLeft+iCyclesRight;
	  //	  fi+=timeradius;
	  fi+=(int)rint(iCyclesLeft/fabs(bb));
	  initframe_ta=fi-timeradius; //New
	  //Check ones again
	  if(initframe_ta<1) printf("SYN WARNING: initframe_ta=%i\n",initframe_ta);
	  dNframesynch=(double)n_cycles/fabs(bb);
	  nframes=(int)rint(dNframesynch);
	  ff=nframes+fi-1;
	  finalframe_ta=ff+timeradius; //New
	  //Check ones again
	  if(finalframe_ta>=tr->max_n){
	    printf("SYN WARNING: finalframe_ta=%i initframe_ta=%i\n",finalframe_ta,initframe_ta);
	      printf("SYN WARNING: fi=%i ff=%i timeradius=%i\n",fi,ff,timeradius);
	  }
	  cyclesexact=(double)nframes*fabs(bb)*harmonic;
	  cycleslower=floor(cyclesexact);
	  cyclesupper=cycleslower+1.0;
	  residue=(double)nframes-dNframesynch;
	    
	  if(do_precise){
	    synch_omega=2.0*M_PI*fabs(bb)/(double)n_cycles;
	    synch_phase_shift=2.0*M_PI*(aa+timeradius*bb);
	  }
	  else{
	    synch_omega=2.0*M_PI/(double)nframes;
	    synch_phase_shift=0.0;
	  }
	  memmove(synch_trace,synch_trace+timeradius,nframes*sizeof(float));
	  synch_trace=realloc((void *)synch_trace,nframes*sizeof(float));
	  memmove(synch_phi,synch_phi+timeradius,nframes*sizeof(double));
	  synch_phi=realloc((void *)synch_phi,nframes*sizeof(double));	  
	  
	  printf("SYN f=%e Exact: Cycles=%f Frames/cycle=%f\n",synch_omega/(2.0*M_PI),cyclesexact,1.0/bb);
	  printf("SYN NF=%i NFC=%f Residue=%f\n",nframes,dNframesynch,residue);
	  *iifr=fi;
	  *fffr=ff;
	  
	  printf("SYN IFr=%i FFr=%i FCycles=%i NFramesC=%f\n",*iifr,*fffr,n_cycles,((double)((int)n_cycles))/fabs(bb));
	}
	else{
	  sacrifice_boundaries=0;
	  printf("SYN Too few cycles. No boundary chopping.\n");
	}
      }
      else{
	//Keep the boundaries and do the best on frame salvaging
	if(fi+timeradius>ffr || ff-timeradius<ifr){
	  printf("SYN WARNING timeradius=%i is too large\n",timeradius);
	  goto bailout;
	}
	initframe_ta =fi-timeradius>ifr ? fi-timeradius : ifr;
	finalframe_ta=ff+timeradius<ffr ? ff+timeradius : ffr;
	printf("SYN Time averaging boundaries %i-%i\n",initframe_ta,finalframe_ta);
      }
    }
  }
  
  if(do_remove_curve || do_remove_linear){
    if(InitializeFitting(tr)){
      printf("SYN Poly fit initialization failed.\n");
      iReturn=15;
      goto bailout;
    }
  }

#ifdef UPPERLOWERMAPS

  if(!(maplx=(double *)calloc(XYdim,sizeof(double)))){
    printf("SYN Cannot allocate for maplx (OK)\n");
    goto UPPERLOWERMAPSbailout;
  }
  if(!(maply=(double *)calloc(XYdim,sizeof(double)))){
    printf("SYN Cannot allocate for maply (OK)\n");
    goto UPPERLOWERMAPSbailout;
  }
  if(!(mapux=(double *)calloc(XYdim,sizeof(double)))){
    printf("SYN Cannot allocate for mapux (OK)\n");
    goto UPPERLOWERMAPSbailout;
  }
  if(!(mapuy=(double *)calloc(XYdim,sizeof(double)))){
    printf("SYN Cannot allocate for mapuy (OK)\n");
    goto UPPERLOWERMAPSbailout;
  }

  norm=1.0/(double)nframes;

  coL=cos(2.0*M_PI*cycleslower*norm);
  siL=sin(2.0*M_PI*cycleslower*norm);
  coU=cos(2.0*M_PI*cyclesupper*norm);
  siU=sin(2.0*M_PI*cyclesupper*norm);
 
  col=1.0; 
  sil=0.0; 
  cou=1.0; 
  siu=0.0; 

  if(!(pcar=AddFrame(tr,fi,a))){
    printf("SYN AddFrame return 0 (%i) (OK)\n",fi);
    goto UPPERLOWERMAPSbailout;
  }
  for(i=0;i<XYdim;i++){
    dumd=(double)pcar->pimage[i];
    maplx[i]=mapux[i]=dumd;
    maply[i]=mapuy[i]=0.0;    
  }

  printf("SYN |");
  for(j=1;j<nframes;j++){
    dumd=col*coL-sil*siL;
    sil=sil*coL+col*siL;
    col=dumd;
    dumd=cou*coU-siu*siU;
    siu=siu*coU+cou*siU;
    cou=dumd;
    
    if(!(pcar=AddFrame(tr,fi+j,a))){
      printf("SYN AddFrame return 0 in (%i) (OK)\n",fi+j);
      goto UPPERLOWERMAPSbailout;
    }
    for(i=0;i<XYdim;i++){
      dumd=(double)pcar->pimage[i];
      maplx[i]+=dumd*col;
      maply[i]+=dumd*sil;
      mapux[i]+=dumd*cou;
      mapuy[i]+=dumd*siu;
    }
   
    if(!(j%100)){ 
      if(!(j%1000)){ 
	printf("|");
	fflush(stdout);
      }
      else printf(".");
      fflush(stdout);
    }
  }
  printf("\n");

  for(i=0;i<XYdim;i++){
    maplx[i]*=norm;
    maply[i]*=norm;
    mapux[i]*=norm;
    mapuy[i]*=norm;
  }
  
  DisplayMap(tr,maplx,maply);
  DisplayMap(tr,mapux,mapuy);

 UPPERLOWERMAPSbailout:

  if(maplx) free(maplx);  
  if(maply) free(maply);  
  if(mapux) free(mapux);  
  if(mapuy) free(mapuy);  

#endif

 bailout:

  if(synch_tracei) free(synch_tracei);  
  return(iReturn);
}
