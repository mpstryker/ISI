#include "iman1.h"

#include <time.h>

int InitializeTrain(TRAIN *tr,char *filen){
  struct stat stat_buf;
  unsigned long f_s,d_s,fr_n;
  FILE *fp;
  char str[128],dir[256],*pc;
  char cfile[256],pfile[256],rfile[256],prfile[256],rfileU[256],rfileL[256];
  int delim='/';
  FILESTRUC *pFS=NULL;
  FILEHEADER fileheader,*pFILEHEADER=NULL;
  DATA_CHUNK DATAChunk,*pDATAChunk=NULL;
  SOFT_CHUNK SOFTChunk,*pSOFTChunk=NULL;
  int iVersion=0;
  char strTag[5];
  unsigned int nFileHeaderSize=0;
  unsigned int nFrameHeaderSize=0;
  int i;
  void *pv,*pvExp=NULL;
  unsigned long ul;
  int iMemoryMultiplier=1;

  strTag[4]='\0';

  if(stat(filen,&stat_buf)){
    printf("INI Cannot stat file %s\n",filen);
    return(-1);
  }

  if(!(fp=fopen("/proc/meminfo","r"))){
    printf("INI Cannot open file /proc/meminfo\nWill allow %i frames in cue\n",max_frames_in_cue);
    memory_size=0;
  }
  else{
    //fgets(str,256,fp);
    fgets(str,256,fp);
    fclose(fp);
    strtok(str," ");
    memory_size=atoi(strtok(NULL," "));
    if((pc=strtok(NULL," "))){
      if(memcmp(pc,"kB",2)){
	iMemoryMultiplier=1;
	strcpy(str,"B");
      }
      else{
	iMemoryMultiplier=1000;
 	strcpy(str,"kB");
     }
    }
    printf("INI This box has %lu %s of RAM\n",memory_size,str);
  }

  if((pc=strrchr(filen,delim))!=NULL){  
    strncpy(dir,filen,(pc-filen)+1);
    dir[(pc-filen)+1]='\0';
    printf("INI Search dir: %s\n",dir);
    strcpy(rfile,pc+1);
  }
  else{
    dir[0]='\0';
    strcpy(rfile,filen);    
  }
  printf("INI Local file: %s\n",rfile);

  while(*rfile){
    strcpy(pfile,cfile);
    strcpy(cfile,dir);
    strcat(cfile,rfile);
    if((fp=fopen(cfile,"r"))==NULL){
      printf("INI Cannot open file %s\nTrying lower case... ",cfile);
      for(i=0;i<strlen(rfile);i++) rfileL[i]=tolower(rfile[i]);
      rfileL[strlen(rfile)]='\0';
      strcpy(cfile,dir);
      strcat(cfile,rfileL);
      if((fp=fopen(cfile,"r"))==NULL){
	printf("\nINI Cannot open file %s(%i)\nTrying upper case... ",cfile,strlen(cfile));
	for(i=0;i<strlen(rfile);i++) rfileU[i]=toupper(rfile[i]);
	rfileU[strlen(rfile)]='\0';
	strcpy(cfile,dir);
	strcat(cfile,rfileU);
	if((fp=fopen(cfile,"r"))==NULL){
	  printf("\nINI Cannot open file %s\n Exhausted\n",cfile);
	  printf("INI Head data loss should be assumed\nEnter 'q' to quit: ");
	  fgets(str,4,stdin);
	  if(*str == 'q'){
	    return(2);
	  }
	  else{
	    strcpy(cfile,pfile);
	    break;
	  }
	}
	else{
	  printf("Succeeded!\n");
	  strcpy(rfile,rfileU);
	}
      }
      else{
	printf("Succeeded!\n");
	strcpy(rfile,rfileL);
      }
    }

    if(!ReadChunkHeader((void*)&DATAChunk,fp)){
      printf("INI Cannot read chunk header from file %s\n",cfile);
      return(3);
    }

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
    printf("INI Got header from %s TAG=%s\n",rfile,strTag);

    switch(*strTag){
    case 'A':
      imagetype = FILE_TYPE_ANALYSIS;
      printf("INI Analysis file. Use other routines to view image\n");
      return(1);
      break;
    case 'C':
      imagetype = FILE_TYPE_COMPRESSED;
      if(compress!=DECOMPRESS){
	printf("INI Use -Q option to decompress this file\n");
	return(2);
      }
      //      compress=DECOMPRESS;
      break;
    case 'G':
    case 'E':
      imagetype = FILE_TYPE_GREEN;
      GreenManager(iVersion,fp,rfile);
      exit(0);
      break;
    case 'T':
      imagetype = FILE_TYPE_STREAM;
      break;
    default:
      imagetype = FILE_TYPE_UNKNOWN;
      printf("INI Unknown file type. Bailing out\n");
      return(-1);
    }
    fclose(fp);
    strcpy(prfile,rfile);

    if(iVersion==VERSION_CHUNK) strncpy(rfile,SOFTChunk.PrevFilename,16);
    else  strcpy(rfile,fileheader.sumPrevFilename);
  }
  printf("INI Reached head of series: %s\n",prfile);
  printf("INI Traversing forward\n");
  strcpy(rfile,prfile);

  tr->iGlobalVersion=iVersion;
  tr->n_files=0;
  tr->files=(FILESTRUC*)0;
  tr->size_fileheader=nFileHeaderSize;
  tr->nFrameHeaderSize=nFrameHeaderSize;

  while(*rfile){
    strcpy(pfile,cfile);
    strcpy(cfile,dir);
    strcat(cfile,rfile);
    if(stat(cfile,&stat_buf)){
      printf("INI Cannot stat file |%s|\nINI Trying lower case... ",cfile);
      for(i=0;i<strlen(rfile);i++) rfileL[i]=tolower(rfile[i]);
      strcpy(cfile,dir);
      strcat(cfile,rfileL);
      if(stat(cfile,&stat_buf)){
	printf("\nINI Cannot stat file |%s|\nINI Trying upper case... ",cfile);
	for(i=0;i<strlen(rfile);i++) rfileL[i]=toupper(rfile[i]);
	strcpy(cfile,dir);
	strcat(cfile,rfileU);
	if(stat(cfile,&stat_buf)){
	  printf("\nINI Cannot stat file |%s|\nINI Exhausted\n",cfile);
	  printf("INI Tail data loss should be assumed\nINI Enter 'q' to quit: ");
	  fgets(str,4,stdin);
	  if(*str == 'q') return(2);
	  else break;
	}
	else{
	  printf(" Succeeded!\n");
	  strcpy(rfile,rfileU);
	}
      }
      else{
	printf(" Succeeded!\n");
	strcpy(rfile,rfileL);
      }
    }

    if((fp=fopen(cfile,"r"))==NULL){
      printf("INI Cannot open file %s\n",cfile);
      printf("INI Tail data loss should be assumed\nINI Enter 'q' to quit: ");
      fgets(str,4,stdin);
      if(*str == 'q') return(2);
      else break;
    }

    if(!(tr->files=(FILESTRUC*)realloc(tr->files,(tr->n_files+1)*sizeof(FILESTRUC)))){
      printf("INI Cannot reallocate for files struct\n");
      return(1);
    }

    pFS=tr->files+tr->n_files;
    pFS->fp=NULL;
    pFS->pFileFooter=NULL;
    pFS->pFileHeader=NULL;
    pFS->fullfilename=NULL;
    pFS->filename=NULL;

    pFS->iVersion=iVersion;

    if(iVersion==VERSION_CHUNK){
      if(FindChunkInFileAndOffset(fp,"DATA",FALSE)!=~0UL){
	pFS->nFileHeaderSize=ftell(fp);
	rewind(fp);
      }
      else{
	printf("INI Cannot find DATA chunk in file %s\n",cfile);
	return(3);
      }
    }
    else pFS->nFileHeaderSize=sizeof(FILEHEADER);

    if(!(pv=pFS->pFileHeader=calloc(1,pFS->nFileHeaderSize))){
      printf("INI Cannot allocate for pFileHeader\n");
      return(1);
    }

    pFS->fullfilename=strdup(cfile);
    pFS->filename=strdup(rfile);
    pFS->fp=fp;
    pFS->state_of_file=FILE_OPEN;

    if(fread(pv,pFS->nFileHeaderSize,1,fp) != 1){
      printf("INI Cannot read fileheader from file %s\n",cfile);
      return(3);
    }

    if(iVersion==VERSION_CHUNK){
      if(!(pSOFTChunk = (SOFT_CHUNK*)FindChunkInBuffer(pv,pFS->nFileHeaderSize,"SOFT"))){
	printf("INI Cannot find SOFT chunk in file header\n");
	return(5);
      }

      // Experiment type
      for(pc = (char*)strListOfKnownExperimentChunks ; *pc && !(pvExp = FindChunkInBuffer(pv, pFS->nFileHeaderSize, pc)); pc+=CHUNK_ID_SIZE+1); 
      
      if(!pvExp){
	printf("INI Cannot find experiment chunk in file header\n");
	return(6);
      }
      
      pFS->iExperiment=1+(pc-strListOfKnownExperimentChunks)/(CHUNK_ID_SIZE+1);
      strncpy(str,pc,CHUNK_ID_SIZE);
      str[CHUNK_ID_SIZE]='\0';
      printf("INI Found %s experiment chunk(%i)\n",str,pFS->iExperiment);

      if(pFS->iExperiment==EXPERIMENT_CONTINUOUS_MODE_CONTINUOUS){
	printf("INI NSynchChannels = %lu NStimulusChanels = %lu\n",((COST_CHUNK*)pvExp)->NSynchChannels,((COST_CHUNK*)pvExp)->NStimulusChanels);

	for(i=0;i<((COST_CHUNK*)pvExp)->NSynchChannels;i++)
	  printf("INI SynchChannelMax[%i] = %lu\n",i,((COST_CHUNK*)pvExp)->SynchChannelMax[i]);

	for(i=0;i<((COST_CHUNK*)pvExp)->NStimulusChanels;i++)
	  printf("INI StimulusPeriod[%i] = %lu\n",i,((COST_CHUNK*)pvExp)->StimulusPeriod[i]);
      }
    }
    else{
      if(frame_rule==1) pFS->iExperiment=EXPERIMENT_CONTINUOUS_MODE_EPISODIC;
      else pFS->iExperiment=EXPERIMENT_CONTINUOUS_MODE_CONTINUOUS;
      pFILEHEADER=(FILEHEADER*)pv;
    }

    if(!tr->n_files){
      if(iVersion==VERSION_CHUNK){
	ulBinX=pSOFTChunk->SpatialBinningX;
	ulBinY=pSOFTChunk->SpatialBinningY;
	ulBinT=pSOFTChunk->TemporalBinning;
	tr->Xdim=Xdim=pSOFTChunk->XSize;
	tr->Ydim=Ydim=pSOFTChunk->YSize;
	tr->XYdim=XYdim=(unsigned long)Xdim*(unsigned long)Ydim;
	tr->nFrameImageSize=pSOFTChunk->SizeOfDataType*XYdim;
	tr->ulDataType=pSOFTChunk->DataType;

	tr->ulNSynchChannels=((COST_CHUNK*)pvExp)->NSynchChannels;
	if(tr->ulNSynchChannels){
	  if(iSynchChannel>=(int)tr->ulNSynchChannels){
	    printf("INI Incorrect Synch channel requested %i Max %lu\n",iSynchChannel,tr->ulNSynchChannels-1);
	    return(11);
	  }
	  
	  if(!(tr->pulSynchChannelMax=(ULONG*)calloc(tr->ulNSynchChannels,sizeof(ULONG)))){
	    printf("INI Cannot allocate for pulSynchChannelMax\n");
	    return(21);
	  }
	  for(i=0;i<tr->ulNSynchChannels;i++){
	    tr->pulSynchChannelMax[i]=((COST_CHUNK*)pvExp)->SynchChannelMax[i];
	  }
	}
	else{
	  if(iSynchChannel>=0){
	    printf("INI No Synch channels available\n");
	    return(31);
	  }
	}
      }
      else{
	ulBinX=ulBinY=pFILEHEADER->sumBinning;
	ulBinT=pFILEHEADER->sumTempBinning;
	tr->Xdim=Xdim=pFILEHEADER->sumXsize;
	tr->Ydim=Ydim=pFILEHEADER->sumYsize;
	tr->XYdim=XYdim=(unsigned long)Xdim*(unsigned long)Ydim;
	tr->nFrameImageSize=pFILEHEADER->sumLSizeOf*XYdim;
	tr->ulDataType=pFILEHEADER->sumLDataType;
      }

      tr->nFrameSize=tr->nFrameHeaderSize+tr->nFrameImageSize;

      printf("INI Dimensions: X %i Y %i Type %lu\n",Xdim,Ydim,tr->ulDataType);
      if(memory_size){
	max_frames_in_cue=(unsigned int)(MAX_MEMORY_USAGE*iMemoryMultiplier*(double)memory_size/(double)tr->nFrameSize);
	printf("INI Will allow %i frames in cue\n",max_frames_in_cue);
      }
    }

    printf("INI Filled header from %s\n",rfile);

    pFS->ini_frame=0;

    if(iVersion==VERSION_CHUNK){
      f_s=((ISOI_CHUNK*)pv)->Size+CHUNK_HEAD_SIZE;

      fr_n=pSOFTChunk->NFramesThisFile;

      if(f_s!=stat_buf.st_size){
	printf("INI Bogus File Size %li(%lu) field in %s\n",f_s,stat_buf.st_size,pFS->filename);
       
	if(f_s>stat_buf.st_size){
	  f_s=stat_buf.st_size;
	  fr_n=(f_s-pFS->nFileHeaderSize)/tr->nFrameSize;
	  printf("INI Fixing NFrames from %li to %lu\n",pSOFTChunk->NFramesThisFile,fr_n);
	}
	else{
	  f_s=stat_buf.st_size;
	  fr_n=(f_s-pFS->nFileHeaderSize)/tr->nFrameSize;
	  printf("INI File size is larger than expected. Keeping orig frame number %lu\n",fr_n);
	}
      }
      if(!(pDATAChunk = (DATA_CHUNK*)FindChunkInBuffer(pv,pFS->nFileHeaderSize,"DATA"))){
	printf("INI Cannot find DATA chunk in file header\n");
	return(5);
      }
      d_s=pDATAChunk->Size;
      //fr_n=pSOFTChunk->NFramesThisFile;
      if(d_s!=tr->nFrameSize*fr_n){
	printf("INI Bogus Data Size %li(%lu) in %s\n",d_s,tr->nFrameSize*fr_n,pFS->filename);
	if(d_s<tr->nFrameSize*fr_n){
	  //	  fr_n=(d_s)/tr->nFrameSize;
	  printf("INI Fixing NFrames according DataSize from %li to %lu\n",pSOFTChunk->NFramesThisFile,fr_n);
	}
      }

      // File footer
      rewind(fp);
      if(FindChunkInFileAndOffset(fp,"DATA",TRUE)==~0UL){
	printf("INI Cannot find DATA chunk and jump in file %s\n",cfile);
	return(3);
      }

      if((pFS->nFileFooterSize=f_s-ftell(fp))){
	if(f_s<ftell(fp)){
	  pFS->nFileFooterSize=0;
	  pFS->pFileFooter=NULL;
	  printf("INI File size < ftell(fp) in file %s\n",cfile);
	}
	else{
	  if(!(pFS->pFileFooter=calloc(1,pFS->nFileFooterSize))){
	    printf("INI Cannot allocate for pFileFooter\n");
	    return(1);
	  }
	  
	  if(fread(pFS->pFileFooter,pFS->nFileFooterSize,1,fp) != 1){
	    printf("INI Cannot read filefooter from file %s\n",cfile);
	    return(3);
	  }
	  printf("INI Filled footer(%lu) from %s\n",pFS->nFileFooterSize,rfile);
	}
      }
      else{
	pFS->pFileFooter=NULL;
	printf("INI No footer in %s\n",rfile);
     }
      
    }
    else{
      f_s=pFILEHEADER->sumBytes+4;
      if(f_s!=stat_buf.st_size){
	printf("INI Bogus File Size(%li) field in %s Using filesize(%lu).\n",f_s,pFS->filename,stat_buf.st_size);
	f_s=stat_buf.st_size;
      }
      fr_n=(f_s-tr->size_fileheader)/tr->nFrameSize;

      if(f_s!=fr_n*tr->nFrameSize+tr->size_fileheader && imagetype!=FILE_TYPE_COMPRESSED){
	printf("INI Data loss at file end\nINI Enter 'q' to quit: ");
	fgets(str,4,stdin);
	if(*str == 'q'){
	  return(2);
	}
      }
    }


    pFS->fin_frame=fr_n-1;
    printf("INI %i frames\n",pFS->fin_frame+1);
    
    strcpy(prfile,rfile);

    if(iVersion==VERSION_CHUNK){
      strncpy(rfile,pSOFTChunk->NextFilename,16);
      //      printf("INI P %s\n",pSOFTChunk->PrevFilename);
      //      printf("INI T %s\n",pSOFTChunk->ThisFilename);
      //      printf("INI N %s\n",pSOFTChunk->NextFilename);
    }
    else strcpy(rfile,pFILEHEADER->sumNextFilename);

    for(i=0;i<strlen(rfile);i++){ 
      rfileL[i]=tolower(rfile[i]);
      rfileU[i]=toupper(rfile[i]);
    }    
    rfileL[i]='\0';
    rfileU[i]='\0';
    tr->n_files++;
  }  
  printf("INI Reached tail of series: %s\n",prfile);
  printf("INI Opened %i files\n",tr->n_files);

  if((tr->frame_file_n=(unsigned int*)calloc(tr->n_files,sizeof(unsigned int)))==NULL){
    printf("INI Cannot allocate for frame_file_n\n");
    return(1);
  }

  printf("INI Running sanity check ... ");
  tr->max_n=0;
  tr->iGlobalExperiment=tr->files[0].iExperiment;
  for(i=0;i<tr->n_files;i++){

    if(tr->iGlobalVersion != tr->files[i].iVersion){
      printf("\nINI Version %i(%i) mismatch in %s\nINI Enter 'q' to quit: ",tr->iGlobalVersion,tr->files[i].iVersion,tr->files[i].filename);
      fgets(str,4,stdin);
      if(*str == 'q') return(2);     
    }

    if(tr->iGlobalExperiment != tr->files[i].iExperiment){
      printf("\nINI Experiment type %i(%i) mismatch in %s\nINI Enter 'q' to quit: ",tr->iGlobalExperiment,tr->files[i].iExperiment,tr->files[i].filename);
      fgets(str,4,stdin);
      if(*str == 'q') return(2);     
    }

    if(tr->size_fileheader != tr->files[i].nFileHeaderSize){
      printf("\nINI FileHeader size mismatch in %s\nINI Enter 'q' to quit: ",tr->files[i].filename);
      fgets(str,4,stdin);
      if(*str == 'q') return(2);     
    }

    if(iVersion==VERSION_CHUNK){
      if(!(pSOFTChunk=(SOFT_CHUNK*)FindChunkInBuffer(tr->files[i].pFileHeader,tr->size_fileheader,"SOFT"))){
	printf("\nINI Cannot find SOFT chunk in file %s\n",cfile);
	return(5);
      }
      if(tr->nFrameHeaderSize != pSOFTChunk->FrameHeaderSize){
	printf("\nINI FrameHeader size mismatch in %s\nINI Enter 'q' to quit: ",tr->files[i].filename);
	fgets(str,4,stdin);
	if(*str == 'q') return(2);     
      }
      if(Xdim != pSOFTChunk->XSize || Ydim != pSOFTChunk->YSize){
	printf("\nINI Dimension mismatch in %s\nINI Enter 'q' to quit: ",tr->files[i].filename);
	fgets(str,4,stdin);
	if(*str == 'q') return(2);     
      }
      if(tr->ulDataType != pSOFTChunk->DataType){
	printf("\nINI Data Type mismatch in %s\nINI Enter 'q' to quit: ",tr->files[i].filename);
	fgets(str,4,stdin);
	if(*str == 'q') return(2);     
      }
    }
    else{
      pFILEHEADER=(FILEHEADER*)(tr->files[i].pFileHeader);
      if(Xdim != pFILEHEADER->sumXsize || Ydim != pFILEHEADER->sumYsize){
	printf("\nINI Dimension mismatch in %s\nINI Enter 'q' to quit: ",tr->files[i].filename);
	fgets(str,4,stdin);
	if(*str == 'q') return(2);     
      }
    }

    tr->files[i].offset_frame+=tr->max_n;
    tr->max_n+=tr->files[i].fin_frame+1;
    tr->frame_file_n[i]=tr->max_n;
  }
  printf("OK\n");
  
  if(!(tr->filename=strdup(tr->files->filename))){
    printf("INI Cannot even allocate for file name\n");
    return(1);
  }
  //  tr->filename[8]='\0';
  

  if(!(tr->head=(CAR*)malloc(sizeof(CAR)))){
    printf("INI Cannot allocate for initial car\n");
    return(4);
  } 

  if(!(tr->frameread=(char*)calloc(tr->max_n,1))){
    printf("INI Cannot allocate for frameread\n");
    return(5);
  }   
  tr->n=0;
  tr->replace=tr->tail=tr->head;
  tr->newest=(CAR*)0;
  tr->head->next=(CAR*)0;
  tr->head->prev=(CAR*)0;

  if(!(tr->head->pFrame=malloc(tr->nFrameSize))){
    printf("INI Cannot allocate for initial car frame buffer\n");
    return(6);
  }
  if(!(tr->framemap=(CAR **)calloc(tr->max_n,sizeof(CAR *)))){
    printf("INI Cannot allocate for framemap\n");
    return(7);
  }

  tr->head->iLock=0;
  tr->head->pFrameHeader = tr->head->pFrame;
  tr->head->pFrameImage  = tr->head->pFrame + tr->nFrameHeaderSize;
  tr->head->pimage = (unsigned short *)(tr->head->pFrameImage);

  if(initframe<0){ 
    if(tr->iGlobalExperiment==EXPERIMENT_CONTINUOUS_MODE_EPISODIC){
      initframe=0;
    }
    else{
      initframe=DEFAULT_INITFRAME;
    }
  }

  if(finalframe<0){
    finalframe=tr->max_n-1;
  }

  //Initialize random number generator
  srand(0xFFFFF&time(NULL));
  return(0);
}


int InitializeBuffers(TRAIN *tr){
  int j,k;
  
  if(!(fbuffer=(float*)malloc(XYdim*sizeof(float)))){
    printf("INI Cannot allocate for fbuffer\n");
    return(7);
  }
  if(!(dumbuffer=(unsigned short *)malloc(tr->nFrameImageSize))){
    printf("INI Cannot allocate for dumbuffer\n");
    return(8);
  }
  if(!(mapx=(double *)calloc(XYdim,sizeof(double)))){
    printf("INI Cannot allocate for mapx\n");
    return(9);
  }
  if(!(mapy=(double *)calloc(XYdim,sizeof(double)))){
    printf("INI Cannot allocate for mapy\n");
    return(10);
  }
  if(iDoDivisionByAverage){
    if(!(map0=(double *)calloc(XYdim,sizeof(double)))){
      printf("INI Cannot allocate for map0\n");
      return(11);
    }
  }

  if(do_cocktailaverage){
    if(!(cocktailaverage=(float *)calloc(XYdim,sizeof(float)))){
      printf("INI Cannot allocate for cocktailaverage, can live without it\n");
      do_cocktailaverage=0;
    }
  }

  radius=(int)floor(dradius);
  if(radius!=0){
    if(!(aindex_x=(int *)calloc((1+2*radius)*(1+2*radius),sizeof(int))) || 
       !(aindex_y=(int *)calloc((1+2*radius)*(1+2*radius),sizeof(int)))){
      printf("INI Cannot allocate for aindex\n");
      return(6);
    }
    aindex_n=0;
    for(j=-radius;j<=radius;j++){
      for(k=-radius;k<=radius;k++){
	if(hypot((double)j,(double)k)<=dradius){ 
	  aindex_x[aindex_n]=k;
	  aindex_y[aindex_n++]=j;
	}
      }
    }
    if(!(aindex_x=(int *)realloc(aindex_x,aindex_n*sizeof(int))) ||
       !(aindex_y=(int *)realloc(aindex_y,aindex_n*sizeof(int)))){
      printf("INI Cannot reallocate for aindex\n");
      return(7);
    }
  }
  else{
    if(!(aindex_x=(int *)calloc(1,sizeof(int))) || 
       !(aindex_y=(int *)calloc(1,sizeof(int)))){
      printf("INI Cannot allocate for aindex\n");
      return(6);
    }
    aindex_x[0]=0;
    aindex_y[0]=0;
    aindex_n=1;
  }

  radiusbig=(int)floor(dradiusbig);
  if(radiusbig!=0){
    if(!(aindexbig_x=(int *)calloc((1+2*radiusbig)*(1+2*radiusbig),sizeof(int))) ||
       !(aindexbig_y=(int *)calloc((1+2*radiusbig)*(1+2*radiusbig),sizeof(int)))){
      printf("INI Cannot allocate for aindexbig\n");
      return(8);
    }
    aindexbig_n=0;
    for(j=-radiusbig;j<=radiusbig;j++){
      for(k=-radiusbig;k<=radiusbig;k++){
	if(hypot((double)j,(double)k)<=dradiusbig){
	  aindexbig_x[aindexbig_n]=k;
	  aindexbig_y[aindexbig_n++]=j;
	}
      }
    }
    if(!(aindexbig_x=(int *)realloc(aindexbig_x,aindexbig_n*sizeof(int))) ||
       !(aindexbig_y=(int *)realloc(aindexbig_y,aindexbig_n*sizeof(int)))){
      printf("INI Cannot reallocate for aindexbig\n");
      return(9);
    }
  }
  else{
    if(!(aindexbig_x=(int *)calloc(1,sizeof(int))) ||
       !(aindexbig_y=(int *)calloc(1,sizeof(int)))){
      printf("INI Cannot allocate for aindexbig\n");
      return(8);
    }
    aindexbig_x[0]=0;
    aindexbig_y[0]=0;
    aindexbig_n=0;
  }

  return(0);
}


int InitializeFitting(TRAIN *tr){
  int i,status_OK=1;
  nframes_good=finalframe-initframe+1;
  if(do_remove_curve){
    if(!poly_fitn){ 
      poly_fitn=(int)ceil((double)n_cycles*POLYNOM_COEFF)+POLYNOM_ADD;
      printf("INI Will fit %i smart polies\n",poly_fitn);
    }
    else{
      /*
      if(poly_fitn>2.0*n_cycles){
	printf("INI Too many polies(%i), reducing to 2*n_cycles(%i)\n",poly_fitn,n_cycles);
	poly_fitn=2.0*n_cycles;
      } 
      */
    }

    do_remove_linear=0;
    if( !(remove_cos=(double *)calloc(poly_fitn,sizeof(double))) || 
	!(remove_sin=(double *)calloc(poly_fitn,sizeof(double))) || 
	!(remove_dum=(double *)calloc(poly_fitn,sizeof(double)))){ 
      printf("INI Cannot allocate for remove_cos/sin/dum\n");
      status_OK=0;
    }
    if(status_OK && (remove_curve=(double **)malloc(poly_fitn*sizeof(double*)))){
      for(i=0;i<poly_fitn;i++){
	if(!(remove_curve[i]=(double *)calloc(XYdim,sizeof(double)))){
	  printf("INI Cannot allocate for remove_curve[%i]\n",i);
	  if(i==0) status_OK=0;
	  else{
	    printf("INI Reducing number of polies from %i to %i\n",poly_fitn,i);
	    poly_fitn=i;
	  }
	  break;
	}
      }
    }
    else{
      printf("INI Cannot allocate for remove_curve\n");
      status_OK=0;
    }
    if(!status_OK){
      printf("INI Will try linear fit\n");
      do_remove_curve=0;
      do_remove_linear=1;
    }
    else{
      if(display_removed){ 
	if(!(fbuffer_removed=(float *)calloc(XYdim,sizeof(float)))){
	  printf("INI Cannot allocate for fbuffer_removed, can live without it\n");
	  display_removed=0;
	}
      }
      printf("INI Initialized Poly Fit\n");
    }
  }  

  if(do_remove_linear){
    if(!(remove_sy=(double *)calloc(XYdim,sizeof(double)))){
      printf("INT Cannot allocate for remove_sy\nINT Will keep going\n");
      do_remove_linear=0;
      return(1);
    }
    if(!(remove_sxy=(double *)calloc(XYdim,sizeof(double)))){
      printf("INT Cannot allocate for remove_sxy\nINT Will keep going\n");
      free(remove_sy);
      do_remove_linear=0;
      return(1);
    }
    printf("INI Initialized Linear Fit\n");
  }
  return(0);
}

int InitializeTimeAveraging(TRAIN *tr){
  if(2*timeradius+1 > tr->max_n){
    printf("INI Time average interval(%i) is larger than max frame(%i)\n",2*timeradius+1,tr->max_n);
    do_timeaveraging=0;
    do_remove_curve=1;
    return(1);
  }
  if(!(timeaverage_buffer=(double *)calloc(XYdim,sizeof(double)))){
    printf("INI Cannot allocate for timeaverage_buffer\n");
    do_timeaveraging=0;
    do_remove_curve=1;
    return(1);
  } 
  else{
    do_remove_linear=do_remove_curve=0;
    printf("INI Time averaging over %i+1+%i=%i frames\n",timeradius,timeradius,2*timeradius+1);
  }
  return(0);
}

void ReleaseTrain(TRAIN *tr,int iMode){
  int i;

  printf("INI Releasing train\n");
  for(i=0;i<tr->n_files;i++){
    if(tr->files[i].fp) fclose(tr->files[i].fp);
    if(tr->files[i].pFileFooter) free(tr->files[i].pFileFooter);
    if(tr->files[i].pFileHeader) free(tr->files[i].pFileHeader);
    if(!iMode && compress && !iKeepOriginalFiles){
      if(unlink(tr->files[i].fullfilename)==-1){
	printf("INI Cannot unlink file %s errno=%i\n",tr->files[i].fullfilename,errno);
      }
    }
    if(tr->files[i].fullfilename) free(tr->files[i].fullfilename);
    if(tr->files[i].filename) free(tr->files[i].filename);
  }
  if(tr->pulSynchChannelMax) free(tr->pulSynchChannelMax);
  if(tr->files) free(tr->files);
  if(tr->frame_file_n) free(tr->frame_file_n);
  if(tr->framemap) free(tr->framemap);
}

int AskToQuit(char *pcLineOut){
  char str[8];
  
  printf(pcLineOut);
  if(pcLineOut[strlen(pcLineOut)-1]=='\n') printf("\n");
  printf("Enter 'q' to quit: ");
  fgets(str,4,stdin);
  if(*str == 'q') return(1);
  else return(0);
}
