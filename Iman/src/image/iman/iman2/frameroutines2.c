#include "iman1.h"

int ReadIn(TRAIN *tr,int fr,CAR *pcar,int a);

CAR  *AddFrame(TRAIN *tr,int fr,int a){
  CAR *pcar;

  if(fr<0 || fr>=tr->max_n){
    printf("FRA(AF) Frame number(%i) is out of range(%i-%i)\n",fr,0,tr->max_n-1);
    return((CAR*)0);
  }
  if(tr->frameread[fr]==(char)FRAME_IN){
    //    printf("Frame %i has been added already\n",fr);
    return(tr->framemap[fr]);
    for(pcar=tr->tail;pcar!=tr->head;pcar=pcar->next) if(pcar->n==fr) break;
    return(pcar);
  }

  if(tr->n >= max_frames_in_cue){
    return(ReplaceFrame(tr,tr->tail->n,fr,replace_scheme,a));
  }

  if(ReadIn(tr,fr,tr->head,a)){
    printf("FRA ReadIn failed\n");
    return((CAR*)0);
  }

  if(!(tr->head->next=(CAR*)malloc(sizeof(CAR)))){
    printf("\nFRA Cannot allocate for car (%luB)\n",(unsigned long)sizeof(CAR));
    return(ShrinkCue(tr,fr,a));
  } 
  if(!(tr->head->next->pFrame=malloc(tr->nFrameSize))){
    printf("\nFRA Cannot allocate for tr->head->pFrame (%luB)\n",(unsigned long)tr->nFrameSize);
    return(ShrinkCue(tr,fr,a));
  }

  tr->head->next->iLock=0;
  tr->head->next->pFrameHeader = tr->head->next->pFrame;
  tr->head->next->pFrameImage  = tr->head->next->pFrame + tr->nFrameHeaderSize;
  tr->head->next->pimage = (unsigned short *)(tr->head->next->pFrameImage);

  tr->newest=tr->head;
  tr->head->n=fr;
  tr->frameread[fr]=(char)FRAME_IN;
  tr->framemap[fr]=tr->head;
  tr->n++;  
  tr->head->next->prev=tr->head;
  tr->head=tr->head->next;
  tr->head->next=(CAR*)0;
  
  if(do_verbose) PrintFrameheader(tr,fr);
  return(tr->newest);
}

CAR *ShrinkCue(TRAIN *tr,int fr,int a){
  int i;
  printf("\nFRA Cannot allocate for car (%luB)\n",(unsigned long)sizeof(CAR));
  printf("FRA Frames in cue %i Fixing max_frames_in_cue=%i -> %i\n",tr->n,max_frames_in_cue,tr->n-NFRAMES_TO_FREE);
  for(i=0;i<NFRAMES_TO_FREE;i++){
    if(RemoveFrame(tr,tr->tail->n)){}
  }
  max_frames_in_cue=tr->n;
  printf("FRA Fixed max_frames_in_cue=%i \n",max_frames_in_cue);    
  return(ReplaceFrame(tr,tr->tail->n,fr,replace_scheme,a));
}

int RemoveFrame(TRAIN *tr,int fr){
  CAR *pcar;
  if(tr->frameread[fr] == (char)FRAME_OUT){
    printf("FRA Frame %i is not in train\n",fr);
    return(1);
  }
  if(fr<0 || fr>=tr->max_n){
    printf("FRA Frame number(%i) is out of range(%i-%i)\n",fr,0,tr->max_n-1);
    return(2);
  }
  if(tr->n > 0){
    for(pcar=tr->tail;pcar!=tr->head;pcar=pcar->next) if(pcar->n==fr) break;
    free(pcar->pFrame);
    if(pcar == tr->newest) tr->newest=(CAR*)0;
    if(pcar == tr->replace){
      if(pcar->next != tr->head) tr->replace=pcar->next;
      else tr->replace=tr->tail;
    }
    if(pcar == tr->tail){
      tr->tail=tr->tail->next;
      free(tr->tail->prev);
    }
    else{
      pcar->next->prev=pcar->prev;
      pcar->prev->next=pcar->next;
      free(pcar);
    }
  }
  tr->frameread[fr] = (char)FRAME_OUT;
  tr->framemap[fr]=NULL;
  tr->n--;

  return(0);
}

CAR *ReplaceFrame(TRAIN *tr,int frold,int frnew,int mode,int a){
  CAR *pcar;
  int i;
  if(frnew<0 || frnew>=tr->max_n){
    printf("FRA New frame number(%i) is out of range(%i-%i)\n",frnew,0,tr->max_n-1);
    return((CAR*)0);
  }
  if(tr->frameread[frnew]==(char)FRAME_IN){
    //    printf("FRA Frame %i has been added already\n",frnew);
    for(pcar=tr->tail;pcar!=tr->head;pcar=pcar->next) if(pcar->n==frnew) break;
    return(pcar);
  }

  switch(mode){

  case REPLACE_OLDEST:
    pcar=tr->replace;
    //    frold=tr->replace->n;
    if(pcar->next != tr->head) tr->replace=pcar->next;
    else tr->replace=tr->tail;    
    break;

  case REPLACE_REQUESTED:
    if(tr->frameread[frold]==(char)FRAME_OUT){
      printf("FRA Frame %i is not loaded\n",frold);
      return((CAR*)0);
    }
    if(frold<0 || frold>=tr->max_n){
      printf("FRA Old frame number(%i) is out of range(%i-%i)\n",frnew,0,tr->max_n-1);
      return((CAR*)0);
    }
    
    for(pcar=tr->tail;pcar!=tr->head;pcar=pcar->next) if(pcar->n==frold) break;
    if(pcar->iLock){
      printf("FRA Old frame (%i) is locked\n",frold);
      return(ReplaceFrame(tr,frold,frnew,REPLACE_RANDOM,a));
    }
    if(pcar == tr->replace){
      if(pcar->next != tr->head) tr->replace=pcar->next;
      else tr->replace=tr->tail;
    }
    break;

  case REPLACE_RANDOM:
    frold=rand()%tr->n;
    for(i=0,pcar=tr->tail;i<frold;i++,pcar=pcar->next);
    for(i=0;i<tr->n;i++){
      if(!pcar->iLock){
	frold=pcar->n;
	break;
      }
      if(pcar!=tr->head) pcar=pcar->next;
      else pcar=tr->tail;
    }
    if(i==tr->max_n){
      printf("FRA Cannot find loaded unlocked frame\n");
      return((CAR*)0);      
    }
    if(pcar == tr->replace){
      if(pcar->next != tr->head) tr->replace=pcar->next;
      else tr->replace=tr->tail;
    }
    break;

  default:
    printf("FRA Unknown replace scheme %i\n",mode);
    return((CAR*)0);    
  }

  if(ReadIn(tr,frnew,pcar,a)){
    printf("FRA ReadIn failed\n");
    return((CAR*)0);
  }

  pcar->n=frnew;
  tr->frameread[frold]=(char)FRAME_OUT;
  tr->frameread[frnew]=(char)FRAME_IN;
  tr->framemap[frold]=NULL;
  tr->framemap[frnew]=pcar;
  tr->newest=pcar;
   
  return(pcar);
}

int ReadIn(TRAIN *tr,int fr,CAR *pcar,int a){
  int i,j,m,fn;
  int x,y,x0,y0;
  unsigned long dl;

  for(fn=0;fn<tr->n_files && fr >= tr->frame_file_n[fn];fn++);
  tr->fp=tr->files[fn].fp;
  if(tr->files[fn].state_of_file!=FILE_OPEN){
    printf("FRA File %s is not open. How come?\n",tr->files[fn].filename);
    return(1);
  }

  if(fseek(tr->fp,tr->size_fileheader + (fr-tr->files[fn].offset_frame)*tr->nFrameSize,SEEK_SET)){
    printf("FRA Cannot rewind file %s %i %i(%i %i)\n",tr->files[fn].filename,(fr-tr->files[fn].offset_frame)*tr->nFrameSize,errno,EBADF,EINVAL); 
    return(2);
  }
  
  if(fread(pcar->pFrame,tr->nFrameSize,1,tr->fp)!=1){
    printf("FRA Cannot read frame from file %s\n",tr->files[fn].filename);
    return(3);
  }
  
  if(a==AVERAGE && radius!=0){
    memcpy((void*)dumbuffer,pcar->pFrameImage,tr->nFrameImageSize);
      /* Average the frame */
    for(i=0;i<XYdim;i++){
      dl=0;
      m=0;
      x0=i%Xdim;
      y0=i/Xdim;
      for(j=0;j<aindex_n;j++){
	x=x0+aindex_x[j];
	y=y0+aindex_y[j];
	if(x>=0 && x<Xdim && y>=0 && y<Ydim){ 
	  dl+=dumbuffer[x+y*Xdim];
	  m++;
	}
      }
      pcar->pimage[i]=(unsigned short)(dl/(float)m);
    }
  }
  return(0);
}

int LoadFrame(TRAIN *tr,int fr,float *fbuf,int a){
  int i;
  float dumf;
  CAR *pcar;

  if(!(pcar=AddFrame(tr,fr,a))){
    printf("FRA AddFrame return 0 fr=%i\n",fr);
    return(1);
  }
  fg=0.0;
  bg=(float)(1<<16);
  for(i=0;i<tr->XYdim;i++){
    dumf=fbuf[i]=(float)pcar->pimage[i];
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

int LockFrame(TRAIN *tr,CAR *pcar,int iLock){
  CAR *pcar1;

  for(pcar1=tr->tail;pcar1!=tr->head && pcar1!=pcar;pcar1=pcar1->next);
  if(pcar1!=pcar){
    printf("FRA Cannot lock: frame is not loaded\n");
    return(1);
  }
  if(iLock) pcar->iLock=1;
  else pcar->iLock=0;
  return(0);
}

int GetRecordUS(TRAIN *tr,int fr,int iRecord,unsigned short *pus){
  int fn;
  FILE *pF;

  if(fr<0 || fr>=tr->max_n){
    printf("FRA(GetRecordUL) Frame number(%i) is out of range(%i-%i)\n",fr,0,tr->max_n-1);
    return(1);
  }
  if(tr->frameread[fr]==(char)FRAME_IN){
    *pus=tr->framemap[fr]->pimage[iRecord];
  }
  else{
    for(fn=0;fn<tr->n_files && fr >= tr->frame_file_n[fn];fn++);
    pF=tr->files[fn].fp;
    if(fseek(pF,tr->size_fileheader + (fr-tr->files[fn].offset_frame)*tr->nFrameSize + tr->nFrameHeaderSize + iRecord*sizeof(unsigned short),SEEK_SET)){
      printf("FRA Cannot rewind file %s %i (%i %i)\n",tr->files[fn].filename,errno,EBADF,EINVAL); 
      return(2);
    }
    if(fread(pus,2,1,pF)!=1){
      printf("FRA Cannot read record from file %s\n",tr->files[fn].filename);
      return(3);
    }
  }

  return(0);
}

int GetRecordUL(TRAIN *tr,int fr,int iByteOffset,unsigned long *pul){
  int fn;
  FILE *pF;

  if(fr<0 || fr>=tr->max_n){
    printf("FRA(GetRecordUL) Frame number(%i) is out of range(%i-%i)\n",fr,0,tr->max_n-1);
    return(1);
  }
  if(tr->frameread[fr]==(char)FRAME_IN){
    *pul=*((unsigned long *)(tr->framemap[fr]->pFrameHeader+iByteOffset));
  }
  else{
    for(fn=0;fn<tr->n_files && fr >= tr->frame_file_n[fn];fn++);
    pF=tr->files[fn].fp;
    if(fseek(pF,tr->size_fileheader + (fr-tr->files[fn].offset_frame)*tr->nFrameSize + iByteOffset,SEEK_SET)){
      printf("FRA Cannot rewind file %s %i %i(%i %i)\n",tr->files[fn].filename,(fr-tr->files[fn].offset_frame)*tr->nFrameSize,errno,EBADF,EINVAL); 
      return(2);
    }
    if(fread(pul,4,1,pF)!=1){
      printf("FRA Cannot read record from file %s\n",tr->files[fn].filename);
      return(3);
    }
  }

  return(0);
}

int GetFrameAverage(TRAIN *tr,int fr,unsigned long *pul,int a){
  unsigned long ul;
  CAR *pcar;
  if(fr<0 || fr>=tr->max_n){
    printf("FRA(GetRecordUL) Frame number(%i) is out of range(%i-%i)\n",fr,0,tr->max_n-1);
    return(1);
  }
  if(!(pcar=AddFrame(tr,fr,a))){
    printf("TRA AddFrame return 0 for initframe %i\n",fr);     
    return(2);
  }
  
  *pul=0;
  for(*pul=ul=0;ul<XYdim;ul++){
    (*pul)+=(unsigned long)pcar->pimage[ul];
  }
  *pul/=XYdim;
  return(0);
}

