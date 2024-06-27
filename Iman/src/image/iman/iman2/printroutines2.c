#include "iman1.h"

#define FILE_NAME_SIZE 16

void PrintFileheader(TRAIN *tr,int fi){
  char str[3*(FILE_NAME_SIZE+1)],strTag[5];
  unsigned long i,x,y,fr;
  unsigned long ulDataBytes,ulFrameHeadSize;
  FILEHEADER *pFILEHEADER=NULL;
  DATA_CHUNK *pDATAChunk=NULL;
  SOFT_CHUNK *pSOFTChunk=NULL;
  short int siNconds=0;
  short int siNreps=0;
  short int siNframeiti=0;
  short int siNframestim=0;

  if(tr->files[fi].iVersion==VERSION_CHUNK){
    if(!(pDATAChunk = (DATA_CHUNK*)FindChunkInBuffer(tr->files[fi].pFileHeader,tr->files[fi].nFileHeaderSize,"DATA"))){
      printf("PRI Cannot find DATA chunk in file header\n");
      return;
    }
    ulDataBytes=pDATAChunk->Size;

    if(!(pSOFTChunk = (SOFT_CHUNK*)FindChunkInBuffer(tr->files[fi].pFileHeader,tr->files[fi].nFileHeaderSize,"SOFT"))){
      printf("PRI Cannot find SOFT chunk in file header\n");
      return;
    }

    strncpy(strTag,pSOFTChunk->Tag,4);
    i=pSOFTChunk->SizeOfDataType;
    x=pSOFTChunk->XSize;
    y=pSOFTChunk->YSize;
    fr=pSOFTChunk->NFramesThisFile;
    ulFrameHeadSize=pSOFTChunk->FrameHeaderSize;
    strncpy(str,pSOFTChunk->PrevFilename,FILE_NAME_SIZE);
    strncpy(str+FILE_NAME_SIZE+1,pSOFTChunk->ThisFilename,FILE_NAME_SIZE);
    strncpy(str+2*(FILE_NAME_SIZE+1),pSOFTChunk->NextFilename,FILE_NAME_SIZE);
  }
  else{
    pFILEHEADER=(FILEHEADER*)(tr->files[fi].pFileHeader);

    ulDataBytes=pFILEHEADER->sumBytes-sizeof(FILEHEADER)+4;
    strncpy(strTag,pFILEHEADER->sumTag,4);
    i=pFILEHEADER->sumLSizeOf;
    x=(long)pFILEHEADER->sumXsize;
    y=(long)pFILEHEADER->sumYsize;
    fr=0;
    ulFrameHeadSize=pFILEHEADER->sumFrameHeaderSize;
    strncpy(str,pFILEHEADER->sumPrevFilename,FILE_NAME_SIZE);
    strncpy(str+FILE_NAME_SIZE+1,pFILEHEADER->sumOrigFilename,FILE_NAME_SIZE);
    strncpy(str+2*(FILE_NAME_SIZE+1),pFILEHEADER->sumNextFilename,FILE_NAME_SIZE);
    if(n_bins){
      siNconds=pFILEHEADER->sumNconds;
      siNreps=pFILEHEADER->sumNreps;
      siNframeiti=pFILEHEADER->sumNframeiti;
      siNframestim=pFILEHEADER->sumNframestim;
    }
  }
  strTag[4]='\0';
  str[FILE_NAME_SIZE]='\0';
  str[2*(FILE_NAME_SIZE+1)-1]='\0';
  str[3*(FILE_NAME_SIZE+1)-1]='\0';

  printf("PRI Tag:           %s           File# %i Version %i\n",strTag,fi,tr->files[fi].iVersion);
  printf("PRI DataBytes:     %lu (frames: %lu %lu)\n",ulDataBytes,fr,ulDataBytes/(x*y*i+ulFrameHeadSize));
  printf("PRI Xsize:         %lu\n",x);
  printf("PRI Ysize:         %lu\n",y);
  printf("PRI Bin X Y T:     %lu %lu %lu\n",ulBinX,ulBinY,ulBinT);
  printf("PRI PrevFilename:  %s\n",str);
  printf("PRI OrigFilename:  %s\n",str+FILE_NAME_SIZE+1);
  printf("PRI NextFilename:  %s\n",str+2*(FILE_NAME_SIZE+1));
  if(n_bins){
    printf("PRI Nconds %i Nreps %i\n",siNconds,siNreps);
    printf("PRI Nframes: iti %i stim %i\n",siNframeiti,siNframestim);
    printf("PRI Total Nframes %li\n",((siNframeiti+siNframestim)/ulBinT)*siNconds*siNreps);
  }
}

void PrintFrameheader(TRAIN *tr,int fr){
  CAR *pcar;
  FRAMEHEADER *pFH;
  FRAM_CHUNK *pFC;
  FRAM_COST_CHUNK *pFCC;
  FRAM_EPST_CHUNK *pFEC;

  if(tr->frameread[fr] == (char)FRAME_OUT){
    printf("PRI Frame %i is not in train\n",fr);
    return;
  }

  for(pcar=tr->tail;pcar!=tr->head;pcar=pcar->next) if(pcar->n==fr) break;

  if(tr->iGlobalVersion==VERSION_CHUNK){
    pFC=(FRAM_CHUNK*)pcar->pFrameHeader;
    printf("\nPRI Tag                %c%c%c%c\n",pFC->Tag[0],pFC->Tag[1],pFC->Tag[2],pFC->Tag[3]);
    printf("PRI FrameSeqNumber     %lu\n",pFC->FrameSeqNumber);
    printf("PRI FrameRingNumber    %lu\n",pFC->FrameRingNumber);
    printf("PRI TimeArrivalUsecLo  %lu\n",pFC->TimeArrivalUsecLo);
    printf("PRI TimeArrivalUsecHi  %lu\n",pFC->TimeArrivalUsecHi);
    printf("PRI TimeDelayUsec      %lu\n",pFC->TimeDelayUsec);
    printf("PRI PotentiallyBad     %lu\n",pFC->PotentiallyBad);
    printf("PRI Locked             %lu\n",pFC->Locked);
    printf("PRI FrameSeqCount      %lu\n",pFC->FrameSeqCount);
    printf("PRI CallbackResult     %lu\n",pFC->CallbackResult);
    if(tr->iGlobalExperiment==EXPERIMENT_CONTINUOUS_MODE_CONTINUOUS){
      pFCC=(FRAM_COST_CHUNK*)(pcar->pFrameHeader+8+*(ULONG*)(pcar->pFrameHeader+4));

      printf("PRI cost Tag                  %c%c%c%c\n",pFCC->Tag[0],pFCC->Tag[1],pFCC->Tag[2],pFCC->Tag[3]);
      printf("PRI cost SynchChannel[0]      %lu\n",pFCC->SynchChannel[0]);
      printf("PRI cost SynchChannelDelay[0] %lu\n",pFCC->SynchChannelDelay[0]);
    }
    if(tr->iGlobalExperiment==EXPERIMENT_CONTINUOUS_MODE_EPISODIC){
      pFEC=(FRAM_EPST_CHUNK*)(pcar->pFrameHeader+8+*(ULONG*)(pcar->pFrameHeader+4));
      printf("PRI epst Tag               %c%c%c%c\n",pFEC->Tag[0],pFEC->Tag[1],pFEC->Tag[2],pFEC->Tag[3]);
      printf("PRI epst SeqNumber         %lu\n",pFEC->SeqNumber);
      printf("PRI epst Repetition        %lu\n",pFEC->Repetition);
      printf("PRI epst Trial             %lu\n",pFEC->Trial);
      printf("PRI epst Condition         %lu\n",pFEC->Condition);
      printf("PRI epst FrameOfCondition  %lu\n",pFEC->FrameOfCondition);
      printf("PRI epst FramePaused       %lu\n",pFEC->FramePaused);
      printf("PRI epst FrameType         %lu\n",pFEC->FrameType);
    }
  }
  else{
    pFH=(FRAMEHEADER*)pcar->pFrameHeader;
    printf("PRI ID            %lu\n",pFH->frameheadID);
    printf("PRI Length        %lu\n",pFH->frameheadlength);
    printf("PRI SeqNum        %lu\n",pFH->seqnum);
    printf("PRI Time          %lu\n",pFH->time);
    printf("PRI PCount        %lu %lu\n",pFH->perfcount.HighPart,pFH->perfcount.LowPart);
    printf("PRI PFreq         %lu %lu\n",pFH->perfreq.HighPart,pFH->perfreq.LowPart);
    printf("PRI Rep           %li\n",pFH->rep);
    printf("PRI Trial         %li\n",pFH->trial);
    printf("PRI Condition     %li\n",pFH->condition);
    printf("PRI Frame of cond %li\n",pFH->frame_of_cond);
    printf("PRI Frame paused  %li\n",pFH->frame_paused);
    printf("PRI Frame type    %li\n",pFH->frame_type);
    printf("PRI seqnum_all    %li\n",pFH->seqnum_all);
    printf("PRI synch_in      %li\n",pFH->synch_in);
  }
}
