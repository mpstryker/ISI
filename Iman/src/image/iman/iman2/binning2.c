#include "iman1.h"

void AverageMaps(int xd,int yd,double **dm,float **fm,int nn,double dr,double drb,int mode);

// All Frame Rules return bin number

// Collapses all frames onto one bined cycle
// respects cycle boundary (corrected 5 April 2002)
int FrameRule0(TRAIN *tr,int fr,CAR *pcar){
  double dPhase;
  dPhase=cycles_per_frame*(double)(fr-initframe);
  dPhase-=floor(dPhase);
  return( (int)(dPhase*(double)n_bins) );
}

#define N_JUNK_FRAMES 4

// Bins frames by condition 
// does condition sub-binning if requested
int FrameRule1(TRAIN *tr,int fr,CAR *pcar){
  static int iJunkFrames=0;
  FRAM_EPST_CHUNK *pepstChunk=NULL;  
  int iCond,iFrameType;

  if(tr->iGlobalVersion==VERSION_CHUNK){
    if(!(pepstChunk=(FRAM_EPST_CHUNK*)FindChunkInBuffer(pcar->pFrameHeader,tr->nFrameHeaderSize,"epst"))){
      printf("BIN FrameRule1 Cannot find epst chunk in frame %i\n",fr);     
      return(-1);
    }
    if(pepstChunk->FrameType==EPST_FRAME_STIM){
      iJunkFrames++;
      if(iJunkFrames>N_JUNK_FRAMES) return(pepstChunk->Condition);
      else return(-1);
    }
    else{
      iJunkFrames=0;
      return(-1);
    }
  }
  else{
    iCond=((FRAMEHEADER*)(pcar->pFrameHeader))->condition;
    iFrameType=((FRAMEHEADER*)(pcar->pFrameHeader))->frame_type;
    //    printf("BIN Cond %2i Type %i\n",iCond,iFrameType);
    if(iFrameType){
      iJunkFrames++;
      if(iJunkFrames>1) return(iCond);
      else return(-1);
    }
    else{
      iJunkFrames=0;
      return(-1);
    }    
  }
}

#define NCONDITIONS_L 8
#define NCONDITIONS_R 8
#define NCONDITIONS_B 4

int BinFrames(TRAIN *tr,int a){
  int fr,nfr;
  int i;
  int ok=1;
  unsigned long ul;
  double dumd,*pd;
  char str[128],small_str[64],*pc;
  FILE *fp;
  double *pdBlank=NULL;
  double *pdCos=NULL,*pdSin=NULL;
  double dCos,dSin;
  int iDoCocktailBlank=0,iNormailizeByBlank=0;
  int iDoHPFiltering=0;
  double dHPFilteringRadius=60.0;

  if(n_bins>0){
    if(tr->iGlobalExperiment!=EXPERIMENT_CONTINUOUS_MODE_EPISODIC && 
       n_bins>(int)(1.0/cycles_per_frame)){
      printf("BIN WARNING: There are more bins than frames in a cycle\n");
    }
    if(!(frame_bins=(double**)malloc(n_bins*sizeof(double*))) || !(frames_in_bin=(int*)calloc(n_bins,sizeof(int)))){
      printf("BIN Cannot allocate for frame_bins/frames_in_bin\n");
      ok=0;
    }
    else{
      for(i=0;i<n_bins;i++){
	if(!(frame_bins[i]=(double*)calloc(tr->XYdim,sizeof(double)))){
	  printf("BIN Cannot allocate for frame_bins[%i]\n",i);
	  ok=0;
	  break;
	}
      }
    }
    if(!ok){
      printf("BIN Time binning initialization failed\nBIN Enter 'q' to quit: ");
      fgets(str,4,stdin);
      if(*str == 'q') return(1);
      n_bins=0;
    }
  }
  else{
    printf("BIN Bad number of bins %i\n",n_bins);
    return(2);
  }

  printf("BIN ");
  for(nfr=0,fr=initframe;fr<=finalframe;fr+=increment1,nfr++){
    if(AddFrameToBins(tr,fr,a)){
      printf("BIN Cannot add frame %i to bins\n",fr);
      return(3);
    }
    if(!(nfr%displayimageN)){
      LoadFrame(tr,fr,fbuffer,a);
      //Double2Float(tr->XYdim,timeaverage_buffer,fbuffer);
      DisplayImage(tr,fr,fbuffer);
    }
    if(!(nfr%100)){ 
      if(!(nfr%1000)){ 
	printf("|"); fflush(stdout);
      }
      else printf(".");
      fflush(stdout);
    }
  }
  printf("\n");


  for(i=0;i<n_bins;i++){
    if(frames_in_bin[i]){
      dumd=(double)(frames_in_bin[i]);
      pd=frame_bins[i];
      for(ul=0;ul<tr->XYdim;ul++){
	pd[ul]/=dumd;
      }
    }
    else{
      printf("BIN Bin %i is empty\n",i);
    }
    printf("BIN B=%3i F=%3i\n",i,frames_in_bin[i]);
  }

  if(iDoDivisionByAverage){
    if(iDivisionByAverageCount!=1+finalframe-initframe) 
      printf("BIN WARNING Division Frame Count missmatch %i %i\n",iDivisionByAverageCount,1+finalframe-initframe);
    dumd=1.0/(dDivisionByAverageFactor*(double)iDivisionByAverageCount);

    for(i=0;i<tr->XYdim;i++) map0[i]*=dumd;

    if(do_timeaveraging || do_remove_curve || do_remove_linear){
      for(i=0;i<n_bins;i++){
	pd=frame_bins[i];
	for(ul=0;ul<tr->XYdim;ul++){
	  pd[ul]/=map0[ul];
	}
      }
    }
    else{
      for(i=0;i<n_bins;i++){
	pd=frame_bins[i];
	for(ul=0;ul<tr->XYdim;ul++){
	  pd[ul]=pd[ul]/map0[ul]-dDivisionByAverageFactor;
	}
      }
    }
  }

  DisplaySeries(tr,n_bins,frame_bins);
  
  printf("BIN Save time bins ?(y/n): ");
  fgets(str,4,stdin);
  if(*str == 'y'){
    pc=small_str;
    if(iDoDivisionByAverage){ 
      *pc='d';
      pc++;
    }
    if(do_timeaveraging) sprintf(pc,"t%i",timeaverage_cycles_int);
    else
      if(do_remove_curve) sprintf(pc,"c%i",poly_fitn);
      else
	if(do_remove_linear) sprintf(pc,"l");
	else *pc='\0';

    i=DOUBLE;
    if(radius!=0 || radiusbig!=0){
      sprintf(str,"%s%s%i%s%i%s%i%s%i%s%i%s%i%s",tr->filename,".bins",(int)harmonic,"_",n_bins,"_",radius,"_",radiusbig,"_",increment1,"_",initframe,small_str);
    }
    else{
      if(increment1==1) sprintf(str,"%s%s%i%s%i%s%i%s%i%s%i%s",tr->filename,".bins",(int)harmonic,"_",n_bins,"_",n_cycles,"_",initframe,"_",finalframe,small_str);
      else sprintf(str,"%s%s%i%s%i%s%i%s%i%s%i%s%i%s",tr->filename,".bins",(int)harmonic,"_",n_bins,"_",n_cycles,"_",initframe,"_",finalframe,"_",increment1,small_str);
    }
    //      printf("BIN %s\n",str);
    if((fp=fopen(str,"w"))==NULL){
      printf("BIN Cannot open file %s\n",str);
    }
    //Val: a quick hack to save episodic bins
    //    if(tr->iGlobalVersion==VERSION_CHUNK){
      if(fwrite((void*)&i,sizeof(int),1,fp) != 1){
	printf("BIN Cannot write type to file %s\n",str);
      }
      if(fwrite((void*)&Xdim,sizeof(unsigned short),1,fp) != 1){
	printf("BIN Cannot write Xdim to file %s\n",str);
      }
      if(fwrite((void*)&Ydim,sizeof(unsigned short),1,fp) != 1){
	printf("BIN Cannot write Ydim to file %s\n",str);
      }
      if(fwrite((void*)&n_bins,sizeof(int),1,fp) != 1){
	printf("BIN Cannot write number of bins to file %s\n",str);
      }
      if(fwrite((void*)&harmonic,sizeof(double),1,fp) != 1){
	printf("BIN Cannot write harmonic to file %s\n",str);
      }
      /*
    }
    else{
      if(fwrite(tr->files[0].pFileHeader,sizeof(FILEHEADER),1,fp) != 1){
	printf("BIN Cannot write FileHeader to file %s\n",str);
      }
    }
      */
    for(i=0;i<n_bins;i++){
      if(fwrite((void*)(frame_bins[i]),XYdim*sizeof(double),1,fp) != 1){
	printf("BIN Cannot write frame_bins[%i] to file %s\n",i,str);
      }
    }
    fclose(fp);
  }
  if(tr->iGlobalExperiment==EXPERIMENT_CONTINUOUS_MODE_EPISODIC){
    printf("BIN Make maps from bins ?(y/n): ");
    fgets(str,4,stdin);
    if(*str == 'y'){
      if(NCONDITIONS_L+NCONDITIONS_R+NCONDITIONS_B!=n_bins){
	printf("BIN Fix NCONDITIONS_x\n");
	return(0);
      }
      if(NCONDITIONS_L!=NCONDITIONS_R){
	printf("BIN Fix NCONDITIONS_L!=NCONDITIONS_R\n");
	return(0);
      }
      if(!(pdCos=(double *)calloc(NCONDITIONS_L,sizeof(double)))){
	printf("BIN Cannot allocate for pdCos\n");
	return(12);
      }
      if(!(pdSin=(double *)calloc(NCONDITIONS_L,sizeof(double)))){
	printf("BIN Cannot allocate for pdSin\n");
	return(13);
      }
      for(i=0;i<NCONDITIONS_L;i++){
	pdCos[i]=cos(i*2.0*M_PI/NCONDITIONS_L);
	pdSin[i]=sin(i*2.0*M_PI/NCONDITIONS_L);
      }
      /*
      if(iDoHPFiltering){
	AverageMaps(tr->Xdim,tr->Ydim,frame_bins,NULL,n_bins-NCONDITIONS_B,0.0,dHPFilteringRadius,DOUBLE);
      }
      */

      if(!do_timeaveraging && iNormailizeByBlank){
	if(!(pdBlank=(double *)calloc(tr->XYdim,sizeof(double)))){
	  printf("BIN Cannot allocate for dBlank\n");
	  return(11);
	}
	if(iDoCocktailBlank){
	  for(i=0;i<NCONDITIONS_L;i++){
	    pd=frame_bins[i];
	    for(ul=0;ul<tr->XYdim;ul++){
	      pdBlank[ul]+=pd[ul];
	    }
	  }
	  dumd=1.0/(dDivisionByAverageFactor*NCONDITIONS_L);
	  for(ul=0;ul<tr->XYdim;ul++){
	    pdBlank[ul]*=dumd;
	  }
	  for(i=0;i<NCONDITIONS_L;i++){
	    pd=frame_bins[i];
	    for(ul=0;ul<tr->XYdim;ul++){
	      pd[ul]/=pdBlank[ul];
	    }
	  }

	  memset(pdBlank,0,tr->XYdim*sizeof(double));
	  for(i=NCONDITIONS_L;i<NCONDITIONS_L+NCONDITIONS_R;i++){
	    pd=frame_bins[i];
	    for(ul=0;ul<tr->XYdim;ul++){
	      pdBlank[ul]+=pd[ul];
	    }
	  }
	  dumd=1.0/(dDivisionByAverageFactor*NCONDITIONS_R);
	  for(ul=0;ul<tr->XYdim;ul++){
	    pdBlank[ul]*=dumd;
	  }
	  for(i=NCONDITIONS_L;i<NCONDITIONS_L+NCONDITIONS_R;i++){
	    pd=frame_bins[i];
	    for(ul=0;ul<tr->XYdim;ul++){
	      pd[ul]/=pdBlank[ul];
	    }
	  }
	}
	else{
	  for(i=n_bins-NCONDITIONS_B;i<n_bins;i++){
	    pd=frame_bins[i];
	    for(ul=0;ul<tr->XYdim;ul++){
	      pdBlank[ul]+=pd[ul];
	    }
	  }
	  dumd=1.0/(dDivisionByAverageFactor*NCONDITIONS_B);
	  for(ul=0;ul<tr->XYdim;ul++){
	    pdBlank[ul]*=dumd;
	  }
	  for(i=0;i<NCONDITIONS_L+NCONDITIONS_R;i++){
	    pd=frame_bins[i];
	    for(ul=0;ul<tr->XYdim;ul++){
	      pd[ul]/=pdBlank[ul];
	    }
	  }
	}
	free(pdBlank);
      }

      if(iDoHPFiltering){
	AverageMaps(tr->Xdim,tr->Ydim,frame_bins,NULL,n_bins-NCONDITIONS_B,-1,dHPFilteringRadius,DOUBLE);
      }
      DisplaySeries(tr,n_bins,frame_bins);
      
      memset(mapx,0,tr->XYdim*sizeof(double));
      memset(mapy,0,tr->XYdim*sizeof(double));

      for(i=0;i<NCONDITIONS_L;i++){
	pd=frame_bins[i];
	dCos=pdCos[i];
	dSin=pdSin[i];
	for(ul=0;ul<tr->XYdim;ul++){
	  dumd=pd[ul];
	  mapx[ul]+=dumd*dCos;
	  mapy[ul]+=dumd*dSin;
	}
      }

      DisplayMap(tr,mapx,mapy);

      memset(mapx,0,tr->XYdim*sizeof(double));
      memset(mapy,0,tr->XYdim*sizeof(double));

      for(i=0;i<NCONDITIONS_R;i++){
	pd=frame_bins[i+NCONDITIONS_L];
	dCos=pdCos[i];
	dSin=pdSin[i];
	for(ul=0;ul<tr->XYdim;ul++){
	  dumd=pd[ul];
	  mapx[ul]+=dumd*dCos;
	  mapy[ul]+=dumd*dSin;
	}
      }

      DisplayMap(tr,mapx,mapy);

      free(pdCos);
      free(pdSin);
     }
  }

  return(0);
}


int AddFrameToBins(TRAIN *tr,int fr,int a){
  int bin;
  unsigned long i;
  double *pd;
  CAR *pcar;

  if(!(pcar=AddFrame(tr,fr,a))){
    printf("BIN AddFrame return 0 for frame %i\n",fr);
    return(1);
  }
  bin=pfunFrameRule[frame_rule](tr,fr,pcar);
  //  if(bin>0) printf("BIN bin = %i %i\n",bin,fr);

  if(iDoDivisionByAverage){
    iDivisionByAverageCount++;
    for(i=0;i<tr->XYdim;i++) map0[i]+=(double)pcar->pimage[i];
  }

  if(bin>=0){
    frames_in_bin[bin]++;
    pd=frame_bins[bin];
  }
  else pd=NULL;
  if(do_timeaveraging){    
    if(AdvanceTimeAverage(tr,fr,a,&timeaverage_buffer_frames_in,&timeaverage_buffer)){
      printf("BIN Cannot advance time average for frame %i\n",fr);
      return(2);
    }
    if(pd){
      for(i=0;i<tr->XYdim;i++){
	pd[i]+=(double)pcar->pimage[i]-timeaverage_buffer[i]/timeaverage_buffer_frames_in;
      }
    }
    return(0);
  }
  if(pd){
    for(i=0;i<tr->XYdim;i++){
      pd[i]+=(double)pcar->pimage[i];
    }
  }
  return(0);
}


int Double2Float(unsigned long xy ,double *dbuf,float *fbuf){
  unsigned long i;
  float dumf;

  fg=0.0;
  bg=(float)(1<<16);
  for(i=0;i<xy;i++){
    dumf=fbuf[i]=(float)dbuf[i];
    if(i>=Xdim){ // Do not take into account first screwed row of each frame
      if(fg<dumf){ 
	fg=dumf;
	ifor=i;
      }
      if(bg>dumf){ 
	bg=dumf;
	ibac=i;
      }
    }
  }
  return(0);
}

int EPSTParser(TRAIN *tr,int ifr,int ffr,int *iifr,int *fffr){
  int iNConditions=0;
  int iNRepetitions=0;
  int iNFramesITI=0,iFramesITICount=0;
  int iNFramesStim=0,iFramesStimCount=0;
  int iNFramesPre=0,iFramesPreCount=0;
  int iNFramesPost=0,iFramesPostCount=0;
  int iFramesPauseCount=0;
  int *piConditionsCount=NULL;
  int i,k;
  void *pv;
  int iepstChunkOffset;
  CAR *pcar;
  EPST_CHUNK *pEPSTChunk;
  FRAM_EPST_CHUNK *pepstChunk;
  char str[128];
  int iExitCode=0;
  int iFirstUsableFrame=-1,iLastUsableFrame=-1;

  printf("BIN IFr=%i FFr=%i\n",ifr,ffr);

  if(!(pEPSTChunk=FindChunkInBuffer(tr->files->pFileHeader,tr->files->nFileHeaderSize,"EPST"))){
    printf("BIN Cannot find EPST chunk in file 0\n");     
    return(1);
  }

  iNConditions=pEPSTChunk->NConditions;
  iNRepetitions=pEPSTChunk->NRepetitions;
  iNFramesITI=pEPSTChunk->NFramesITI;
  iNFramesStim=pEPSTChunk->NFramesStim;
  iNFramesPre=pEPSTChunk->NFramesBlankPre;
  iNFramesPost=pEPSTChunk->NFramesBlankPost;

  printf("BIN Conditions=%i Repetitions=%i\n",iNConditions,iNRepetitions);
  printf("BIN Frames ITI=%i Stim=%i Pre=%i Post=%i\n",iNFramesITI,iNFramesStim,iNFramesPre,iNFramesPost);
  printf("BIN Frames Total (not paused) = %i\n",iNConditions*iNRepetitions*(iNFramesITI+iNFramesStim)+iNFramesPre+iNFramesPost);

  if(!(pcar=AddFrame(tr,ifr,AVERAGE_NOT))){
    printf("BIN AddFrame return NULL (%i)\n",ifr);     
    return(2);
  }
  if(!(pv=FindChunkInBuffer(pcar->pFrameHeader,tr->nFrameHeaderSize,"epst"))){
    printf("BIN Cannot find epst chunk in frame %i\n",ifr);     
    return(3);
  }

  iepstChunkOffset=pv - pcar->pFrameHeader;
  //  printf("BIN iepstChunkOffset %i\n",iepstChunkOffset);

  if(!(piConditionsCount=(int*)calloc(iNConditions,sizeof(int)))){
    printf("BIN Cannot allocate for piConditionsCount\n");     
    return(4);
  }

  printf("BIN |");
  for(i=ifr,k=1;i<=ffr;i++,k++){
    if(!(pcar=AddFrame(tr,i,AVERAGE_NOT))){
      printf("BIN AddFrame return NULL (%i)\n",i);     
      iExitCode=1;
      goto bailout;
    }
    pepstChunk=(FRAM_EPST_CHUNK *)(pcar->pFrameHeader+iepstChunkOffset);
    if(memcmp(pepstChunk,"epst",4)){
      sprintf(str,"BIN Frame header %i does not have valid epst chunk",i);
      if(AskToQuit(str)){
	iExitCode=5;
	goto bailout;
      }
      continue;
    }

    //    printf("%i FSN1 = %lu FSN2 = %lu FT = %lu \n",i,((FRAM_CHUNK*)pcar->pFrameHeader)->FrameSeqNumber,pepstChunk->SeqNumber,pepstChunk->FrameType);

    switch(pepstChunk->FrameType){
    case EPST_FRAME_ITI:
      iFramesITICount++;
      if(iFirstUsableFrame==-1) iFirstUsableFrame=i;
      iLastUsableFrame=i;
      break;

    case EPST_FRAME_STIM:
      iFramesStimCount++;
      if(iFirstUsableFrame==-1) iFirstUsableFrame=i;
      iLastUsableFrame=i;
      if(pepstChunk->Condition<iNConditions){
	piConditionsCount[pepstChunk->Condition]++;
      }
      else{
	sprintf(str,"BIN Condition %lu out of bound in frame %i\n",pepstChunk->Condition,i);
	if(AskToQuit(str)){
	  iExitCode=6;
	  goto bailout;
	}
      }
      break;

    case EPST_FRAME_PRE_STIM:
      iFramesPreCount++;
      break;

    case EPST_FRAME_POST_STIM:
      iFramesPostCount++;
      break;

    case EPST_FRAME_PAUSE:
      iFramesPauseCount++;
      break;

    default:
      printf("BIN Unknown EPST frame type %lu\n",pepstChunk->FrameType);
      //      return(11);
    }
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

  for(i=0;i<iNConditions;i++){
    if(piConditionsCount[i]){
      printf("BIN Condition %2i Frames %3i\n",i,piConditionsCount[i]);
    }
    else{
      printf("BIN Condition %2i Frames   0 <--- NOT PRESENT\n",i);
    }
  }

  printf("BIN Frame Count: ITI=%i Stim=%i Pre=%i Post=%i Pause=%i\n",iFramesITICount,iFramesStimCount,iFramesPreCount,iFramesPostCount,iFramesPauseCount);
  printf("BIN Frames Total = %i\n",iFramesITICount+iFramesStimCount+iFramesPreCount+iFramesPostCount+iFramesPauseCount);
  printf("BIN FirstUsableFrame %i LastUsableFrame %i\n",iFirstUsableFrame,iLastUsableFrame);

  if(iFirstUsableFrame<0){
    printf("BIN No usable frames found\n");
    iExitCode=11;
    goto bailout;
  }

  if(iFirstUsableFrame-ifr>iFramesPreCount){
    printf("BIN There are %i NON PRE_STIM (pause?) frames in the header\n",iFirstUsableFrame-ifr-iFramesPreCount);
  }
  else{
    if(iFirstUsableFrame-ifr<iFramesPreCount){
      printf("BIN There are %i PRE_STIM frames outside the header\n",iFramesPreCount-(iFirstUsableFrame-ifr));
    }
  }

  if(ffr-iLastUsableFrame>iFramesPostCount){
    printf("BIN There are %i NON POST_STIM (pause?) frames in the footer\n",ffr-iLastUsableFrame-iFramesPreCount);
  }
  else{
    if(ffr-iLastUsableFrame<iFramesPostCount){
      printf("BIN There are %i POST_STIM frames outside the footer\n",iFramesPostCount-(ffr-iLastUsableFrame));
    }
  }

  if(do_timeaveraging){
    timeradius=(int)rint(timeaverage_blocks*(iNFramesITI+iNFramesStim));
    if(InitializeTimeAveraging(tr)){ 
      printf("BIN Time averaging initialization failed\n");
      iExitCode=13;
      goto bailout;
    }

    // Using cycles_per_frame for condition_blocks_per_frame
    cycles_per_frame=1.0/(double)(iNFramesITI+iNFramesStim);

    i=ifr;
    if(iFirstUsableFrame-ifr>=timeradius){
      ifr=iFirstUsableFrame;
      initframe_ta=ifr-timeradius;
      printf("BIN Lower boundary OK\n");
    }
    else{
      initframe_ta=ifr;
      if(sacrifice_boundaries){
	if(ifr+timeradius<iLastUsableFrame) ifr+=timeradius;
	else{
	  ifr=iFirstUsableFrame;
	  sacrifice_boundaries=0;
	  printf("BIN Cannot sacrifice boundaries. Too few frames\n");
	}
      }
      else{
	ifr=iFirstUsableFrame;
      }
    }
    *iifr=ifr;
    
    if(i!=ifr) printf("BIN Fixed init frame from %i to %i\n",i,ifr);

    i=ffr;
    if(ffr-iLastUsableFrame>=timeradius){
      ffr=iLastUsableFrame;
      finalframe_ta=ffr+timeradius;
      printf("BIN Upper boundary OK\n");
    }
    else{
      finalframe_ta=ffr;
      if(sacrifice_boundaries){
	if(ffr-timeradius>ifr) ffr-=timeradius;
	else{
	  ffr=iLastUsableFrame;
	  sacrifice_boundaries=0;
	  printf("BIN Cannot sacrifice upper boundaries. Too few frames\n");
	}
      }
      else{
	ffr=iLastUsableFrame;
      }
    }
    *fffr=ffr;
    
    if(i!=ffr) printf("BIN Fixed final frame from %i to %i\n",i,ffr);
  }

  if(n_bins>1){
    printf("BIN Sub binning is not implemented\n");
  }
  n_bins=iNConditions;
  printf("BIN N bins = %i\n",n_bins);
  

 bailout:

  if(piConditionsCount) free(piConditionsCount);
  return(iExitCode);
}

void AverageMaps(int xd,int yd,double **dm,float **fm,int nn,double dr,double drb,int mode){
  int i,j,k,m,n,x,y,x0,y0;
  int r,rb;
  unsigned long xy;
  double dumd;
  double *localdbuffer=NULL,*pd;
  float dumf;
  float *localfbuffer=NULL,*pf;
  int *piAIndex_x=NULL,*piAIndexbig_x=NULL,*piAIndex_y=NULL,*piAIndexbig_y=NULL;
  int iAIndex_n=0,iAIndexbig_n=0;

  if(dr<=0.0) r=0;
  else r=(int)dr;
  if(r> ((xd>yd) ? xd/2 : yd/2)){
    printf("AVE Low-pass radius is too large (%i)\n",r);
    return;
  }

  xy=(unsigned long)xd*(unsigned long)yd;
  if(r>0){
    if(!(piAIndex_x=(int *)calloc((1+2*r)*(1+2*r),sizeof(int))) || 
       !(piAIndex_y=(int *)calloc((1+2*r)*(1+2*r),sizeof(int)))){
      printf("AVE Cannot allocate for aindex\n");
      return;
    }
    iAIndex_n=0;
    for(j=-r;j<=r;j++){
      for(k=-r;k<=r;k++){
	if(hypot((double)j,(double)k)<=dr){ 
	  piAIndex_x[iAIndex_n]=k;
	  piAIndex_y[iAIndex_n++]=j;
	}
      }
    }
    printf("AVE Low-pass averaging over %i pixels\n",iAIndex_n);
    switch(mode){
    case DOUBLE:
      if(!(localdbuffer=(double *)calloc(xy,sizeof(double)))){
	printf("INI Cannot allocate for localdbuffer\n");
	return;
      }
      for(n=0;n<nn;n++){
	pd=dm[n];
	for(i=0;i<xy;i++){
	  dumd=0.0;
	  m=0;
	  x0=i%xd;
	  y0=i/xd;
	  for(j=0;j<iAIndex_n;j++){
	    x=x0+piAIndex_x[j];
	    y=y0+piAIndex_y[j];
	    if(x>=0 && x<xd && y>=0 && y<yd){ 
	      dumd+=pd[x+y*xd];
	      m++;
	    }
	  }
	  localdbuffer[i]=dumd/(double)m;
	}
	memcpy(pd,localdbuffer,xy*sizeof(double));    
      }    
      free(localdbuffer);
      break;
    case FLOAT:
      if(!(localfbuffer=(float *)calloc(xy,sizeof(float)))){
	printf("INI Cannot allocate for localfbuffer\n");
	return;
      }
      for(n=0;n<nn;n++){
	pf=fm[n];
	for(i=0;i<xy;i++){
	  dumf=0.0;
	  m=0;
	  x0=i%xd;
	  y0=i/xd;
	  for(j=0;j<iAIndex_n;j++){
	    x=x0+piAIndex_x[j];
	    y=y0+piAIndex_y[j];
	    if(x>=0 && x<xd && y>=0 && y<yd){ 
	      dumf+=pf[x+y*xd];
	      m++;
	    }
	  }
	  localfbuffer[i]=dumf/(float)m;
	}
	memcpy(pf,localfbuffer,xy*sizeof(float));    
      }
      free(localfbuffer);      
      break;
    default:
      printf("AVE Unknown mode %i\n",mode);
      return;
    }
    free(piAIndex_x);
    free(piAIndex_y);
    printf("AVE Maps low-pass averaged\n");
  }

  if(drb<=0.0) rb=0;
  else rb=(int)drb;
  if(rb> ((xd>yd) ? xd/2 : yd/2)){
    printf("AVE High-pass radius is too large (%i)\n",rb);
    return;
  }

  if(rb>0){
    if(!(piAIndexbig_x=(int *)calloc((1+2*rb)*(1+2*rb),sizeof(int))) ||
       !(piAIndexbig_y=(int *)calloc((1+2*rb)*(1+2*rb),sizeof(int)))){
      printf("AVE Cannot allocate for aindexbig\n");
      return;
    }
    iAIndexbig_n=0;
    for(j=-rb;j<=rb;j++){
      for(k=-rb;k<=rb;k++){
	if(hypot((double)j,(double)k)<=drb){
	  piAIndexbig_x[iAIndexbig_n]=k;
	  piAIndexbig_y[iAIndexbig_n++]=j;
	}
      }
    }
    printf("AVE High-pass averaging %i maps over %i pixels\n",nn,iAIndexbig_n);

    switch(mode){
    case DOUBLE:
      if(!(localdbuffer=(double *)calloc(xy,sizeof(double)))){
	printf("INI Cannot allocate for localdbuffer\n");
	return;
      }
      for(n=0;n<nn;n++){
	pd=dm[n];
	for(i=0;i<xy;i++){
	  dumd=0;
	  m=0;
	  x0=i%xd;
	  y0=i/xd;
	  for(j=0;j<iAIndexbig_n;j++){
	    x=x0+piAIndexbig_x[j];
	    y=y0+piAIndexbig_y[j];
	    if(x>=0 && x<xd && y>=0 && y<yd){ 
	      dumd+=pd[x+y*xd];
	      m++;
	    }
	  }
	  localdbuffer[i]=dumd/(double)m;
	}
	for(i=0;i<xy;i++) pd[i]-=localdbuffer[i];
	//	memcpy(pd,localdbuffer,xy*sizeof(double));
	printf(".");
	fflush(stdout);
      }
      printf("\n");
      free(localdbuffer);
      break;
    case FLOAT:
      if(!(localfbuffer=(float *)calloc(xy,sizeof(float)))){
	printf("INI Cannot allocate for localfbuffer\n");
	return;
      }
      for(n=0;n<nn;n++){
	pf=fm[n];
	for(i=0;i<xy;i++){
	  dumf=0;
	  m=0;
	  x0=i%xd;
	  y0=i/xd;
	  for(j=0;j<iAIndexbig_n;j++){
	    x=x0+piAIndexbig_x[j];
	    y=y0+piAIndexbig_y[j];
	    if(x>=0 && x<xd && y>=0 && y<yd){ 
	      dumf+=pf[x+y*xd];
	      m++;
	    }
	  }
	  localfbuffer[i]=dumf/(float)m;
	}
	for(i=0;i<xy;i++) pf[i]-=localfbuffer[i];
	printf(".");
	fflush(stdout);
      }
      printf("\n");
      free(localfbuffer);
      break;
    default:
      printf("AVE Unknown mode %i\n",mode);
      return;
    }
    free(piAIndexbig_x);
    free(piAIndexbig_y);
    printf("AVE Maps high-pass averaged\n");
  }
}
