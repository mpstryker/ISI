#include "iman1.h"

typedef union TwoBytes{
  unsigned short us;
  unsigned char uc[2];
} TwoBytes;

typedef union FourBytes{
  unsigned long ul;
  unsigned char uc[4];
} FourBytes;

typedef FourBytes RECORD;

#define RECORD_N_SIZE (sizeof(unsigned long))
#define EXTRA_RECORD_SIZE (sizeof(RECORD))
#define ORIG_RECORD_SIZE (sizeof(unsigned short))
#define COMP_RECORD_SIZE (sizeof(signed char))

int Compressor(TRAIN *tr,int mode){
  int i,j,k;
  unsigned long ul;
  int diff;
  unsigned int ifr,ffr;
  unsigned long nrecords=0,total_nrecords=0,nrecords_max=0;
  unsigned long compressed_frame_size,original_frame_size;
  unsigned long compressed_file_size,original_file_size,actual_written_size;
  unsigned long ulNCompFrames;
  FILE *fp=NULL;
  char compfile[64],str[64];
  FILESTRUC *pFS=NULL;
  FILEHEADER *pFILEHEADER=NULL;
  void *pCompFileHeader=NULL;
  unsigned long ulCompFileHeaderSize=0;
  unsigned long ulUnCompFileHeaderSize=0;
  SOFT_CHUNK *pSOFTChunk=NULL;
  COMP_CHUNK sCOMPChunk,*pCOMPChunk=NULL;
  CAR *pcar=(CAR*)0;
  RECORD *precord=(RECORD*)0;
  TwoBytes tb;
  unsigned short *pframet=NULL,*pframe0=NULL,*pframe1=NULL;
  signed char *pcframe=NULL;
  unsigned long elements_in_frame;
  struct stat stat_buf;
  int iReturn=0;
  int iContinueFailedCompression=0;

  if(mode!=COMPRESS && mode!=DECOMPRESS){
    printf("COM Unknown compression mode %i\n",mode);
    return(-1);
  }

  original_frame_size=(unsigned long)tr->nFrameSize;
  elements_in_frame=compressed_frame_size=(original_frame_size*COMP_RECORD_SIZE)/ORIG_RECORD_SIZE;

  //  printf("COM OFS=%i C_F_S=%i R_S=%i\n",original_frame_size,compressed_frame_size,EXTRA_RECORD_SIZE);

  if(!(pframe0=(unsigned short *)malloc(elements_in_frame*ORIG_RECORD_SIZE)) ||
     !(pframe1=(unsigned short *)malloc(elements_in_frame*ORIG_RECORD_SIZE)) ||
     !(pcframe=(signed char *)malloc(elements_in_frame*COMP_RECORD_SIZE))){
    printf("COM Cannot allocate for internal buffers\n");
    iReturn=-2;
    goto bailout;
  }

  max_frames_in_cue=2;

  if(mode==COMPRESS){
    // COMPRESS
    for(k=0,pFS=tr->files;k<tr->n_files;k++,pFS=tr->files+k){

      sprintf(compfile,"%sz",pFS->filename);
      if(!stat(compfile,&stat_buf)){
	printf("COM File %s exists. Rewrite (y/n)?",compfile);
	fgets(str,4,stdin);
	if(*str != 'y'){
	  iReturn=-3;
	  goto bailout;
	}
      }
     
      if(!(pCompFileHeader=realloc(pCompFileHeader,pFS->nFileHeaderSize+sizeof(COMP_CHUNK)))){
	printf("COM Cannot reallocate for pFileHeader\n");
	iReturn=-4;
	goto bailout;
      }
      memcpy(pCompFileHeader,pFS->pFileHeader,pFS->nFileHeaderSize);

      ifr=pFS->offset_frame+pFS->ini_frame;
      ffr=pFS->offset_frame+pFS->fin_frame;
      ulNCompFrames=ffr-ifr+1;
      original_file_size=pFS->nFileHeaderSize+original_frame_size*ulNCompFrames+pFS->nFileFooterSize;

      if(pFS->iVersion==VERSION_CHUNK){
	ulCompFileHeaderSize=pFS->nFileHeaderSize+sizeof(COMP_CHUNK);
	if(!(pSOFTChunk = (SOFT_CHUNK*)FindChunkInBuffer(pCompFileHeader,pFS->nFileHeaderSize,"SOFT"))){
	  printf("COM Cannot find SOFT chunk in file header\n");
	  iReturn=5;
	  goto bailout;
	}
	if(pSOFTChunk->Tag[0]=='C'){
	  printf("COM File %s has been compressed already\n",pFS->filename);
	  iReturn=-5;
	  goto bailout;
	}
	pSOFTChunk->Tag[0]='C';
	if(strlen(pSOFTChunk->ThisFilename)) 
	  pSOFTChunk->ThisFilename[strlen(pSOFTChunk->ThisFilename)]='z';
	if(strlen(pSOFTChunk->PrevFilename)) 
	  pSOFTChunk->PrevFilename[strlen(pSOFTChunk->PrevFilename)]='z';
	if(strlen(pSOFTChunk->NextFilename)) 
	  pSOFTChunk->NextFilename[strlen(pSOFTChunk->NextFilename)]='z';
	
	SetChunkNameTag(&sCOMPChunk,sizeof(sCOMPChunk),"COMP","0000");
	sCOMPChunk.CompressedFrameSize=compressed_frame_size;
	sCOMPChunk.CompressedRecordSize=EXTRA_RECORD_SIZE;
	sCOMPChunk.CompressedFrameNumber=ulNCompFrames;
	memset(sCOMPChunk.Free,0,sizeof(sCOMPChunk.Free));

	if(!InsertChunk(pCompFileHeader,pFS->nFileHeaderSize,&sCOMPChunk,sizeof(sCOMPChunk),(void*)pSOFTChunk+sizeof(SOFT_CHUNK))){
	  printf("COM Cannot insert COMP chunk\n");
	  iReturn=-5;
	  goto bailout;
	}
      }
      else{
	ulCompFileHeaderSize=pFS->nFileHeaderSize;
	pFILEHEADER=(FILEHEADER*)pCompFileHeader;
	if(pFILEHEADER->sumTag[0]=='C'){
	  printf("COM File %s has been compressed already\n",pFS->filename);
	  iReturn=-5;
	  goto bailout;
	}
      
	pFILEHEADER->sumTag[0]='C';
	if(strlen(pFILEHEADER->sumOrigFilename)) 
	  pFILEHEADER->sumOrigFilename[strlen(pFILEHEADER->sumOrigFilename)]='z';
	if(strlen(pFILEHEADER->sumPrevFilename)) 
	  pFILEHEADER->sumPrevFilename[strlen(pFILEHEADER->sumPrevFilename)]='z';
	if(strlen(pFILEHEADER->sumNextFilename)) 
	  pFILEHEADER->sumNextFilename[strlen(pFILEHEADER->sumNextFilename)]='z';
	
	pFILEHEADER->sumCompressedFrameSize=compressed_frame_size;
	pFILEHEADER->sumCompressedRecordSize=EXTRA_RECORD_SIZE;      
	pFILEHEADER->sumCompressedFrameNumber=ulNCompFrames;
      }

      printf("COM Compressing file %s(%li frames) to %s\n",pFS->filename,ulNCompFrames,compfile);
      
      if(!(pcar=AddFrame(tr,ifr,AVERAGE_NOT))){
	printf("COM AddFrame return 0 for initframe %i\n",ifr);     
	iReturn=1;
	goto bailout;
      }

      memcpy((void *)pframe0,pcar->pFrame,tr->nFrameSize);
      
      if(!(fp=fopen(compfile,"w"))){
	printf("COM Cannot open file %s\n",compfile);
	iReturn=2;
	goto bailout;
      }
      if(fwrite(pCompFileHeader,ulCompFileHeaderSize,1,fp)!=1){
	printf("COM Cannot write fileheader to %s\n",compfile);
	iReturn=3;
	goto bailout;
      }
      if(fwrite((void*)pframe0,original_frame_size,1,fp)!=1){
	printf("COM Cannot write ini frame to %s\n",compfile);
	iReturn=4;
	goto bailout;
      }

      for(total_nrecords=nrecords=nrecords_max=0,j=ifr+1;j<=ffr;j++){
	if(!(pcar=AddFrame(tr,j,AVERAGE_NOT))){
	  printf("COM AddFrame return 0 for %i\n",j);     
	  iReturn=5;
	  goto bailout;
	}

	memcpy((void *)pframe1,pcar->pFrame,tr->nFrameSize);
	
	for(ul=0;ul<elements_in_frame;ul++){
	  diff=(int)pframe1[ul]-(int)pframe0[ul];
	  if( diff>127 || diff<-128){
	    nrecords++;
	    if(nrecords_max<nrecords){
	      nrecords_max=nrecords;
	      if(!(precord=(RECORD*)realloc((void*)precord,nrecords_max*EXTRA_RECORD_SIZE))){
		printf("COM Cannot realloc for precord\n");
		iReturn=6;
		goto bailout;
	      }
	    }
	    tb.us=pframe1[ul];
	    pcframe[ul]=tb.uc[0];
	    precord[nrecords-1].ul=ul;
	    precord[nrecords-1].uc[3]=tb.uc[1];
	  }
	  else{
	    pcframe[ul]=(signed char)diff;
	  }
	}

	if(!iContinueFailedCompression && nrecords*EXTRA_RECORD_SIZE>original_frame_size-compressed_frame_size){
	  printf("COM WARNING: NO COMPRESSION ACHIEVED, NOISY RECORDS\n");
	  printf("COM WARNING: (extra/orig) elements = %li/%li = %0.2f, waste = %li Bytes\n",nrecords,elements_in_frame,(float)nrecords/elements_in_frame,nrecords*EXTRA_RECORD_SIZE-(original_frame_size-compressed_frame_size));
	  printf("COM WARNING: MAKE SURE HARDWARE GAIN IS SET PROPERLY\n");
	  if(!iKeepOriginalFiles){
	    printf("COM Will keep original files\n");
	    iKeepOriginalFiles=1;
	  }
	  printf("COM Enter <y> to continue ");
	  fgets(str,4,stdin);
	  if(*str != 'y'){
	    iReturn=77;
	    goto bailout;
	  }
       	  iContinueFailedCompression=1;
	}

	//	printf("COM F %i E %i\n",j,nrecords);
	if(fwrite((void*)pcframe,compressed_frame_size,1,fp)!=1){
	  printf("COM Cannot write cframe %i to %s\n",j,compfile);
	  iReturn=7;
	  goto bailout;
	}
	if(fwrite((void*)(&nrecords),RECORD_N_SIZE,1,fp)!=1){
	  printf("COM Cannot write number of records for frame %i to %s\n",j,compfile);
	  iReturn=8;
	  goto bailout;
	}
	if(nrecords){
	  if(fwrite((void*)precord,nrecords*EXTRA_RECORD_SIZE,1,fp)!=1){
	    printf("COM Cannot write records for frame %i to %s\n",j,compfile);
	    iReturn=9;
	    goto bailout;
	  }
	}
	total_nrecords+=nrecords;
	nrecords=0;
	pframet=pframe0;
	pframe0=pframe1;
	pframe1=pframet;    
      }

      if(pFS->nFileFooterSize){
	if(fwrite(pFS->pFileFooter,pFS->nFileFooterSize,1,fp)!=1){
	  printf("COM Cannot write filefooter to %s\n",compfile);
	  iReturn=10;
	  goto bailout;
	}
	else{
	  printf("COM Wrote filefooter to %s, %lu bytes\n",compfile,pFS->nFileFooterSize);	  
	}
      }

      compressed_file_size=ulCompFileHeaderSize+original_frame_size+(RECORD_N_SIZE+compressed_frame_size)*(ulNCompFrames-1)+total_nrecords*EXTRA_RECORD_SIZE+pFS->nFileFooterSize;
      printf("COM Total records %lu\n",total_nrecords);
      printf("COM Original size %lu Compressed size %lu (%.2f%%)\n",original_file_size,compressed_file_size,100.0*((float)original_file_size-(float)compressed_file_size)/(float)original_file_size);

      if(compressed_file_size>original_file_size){
	printf("COM WARNING: NO COMPRESSION ACHIEVED IN FILE %s\n",compfile);
	if(!iKeepOriginalFiles){
	  printf("COM Will keep original files\n");
	  iKeepOriginalFiles=1;
	}
	printf("COM Enter <y> to continue ");
	fgets(str,4,stdin);
	if(*str != 'y'){
	  iReturn=110;
	  goto bailout;
	}
      }
      fclose(fp);    
      fp=NULL;
    }
  }
  else{
    // DECOMPRESS
    printf("COM original_frame_size=%li\n",original_frame_size);
    for(k=0,pFS=tr->files;k<tr->n_files;k++,pFS=tr->files+k){
      actual_written_size=0;

      if(iRemoveTrailingZsFromDecompressedFileName){
	sprintf(compfile,"%s",pFS->filename);
	if(compfile[strlen(compfile)-1]=='z'){
	  compfile[strlen(compfile)-1]='\0';
	}
	else{
	  printf("COM Compressed file name does not have trailing z, adding zz\n");
	  sprintf(compfile,"%szz",pFS->filename);	  
	}
      }
      else{
	sprintf(compfile,"%sz",pFS->filename);
      }
      if(!stat(compfile,&stat_buf)){
	printf("COM File %s exists. Rewrite (y/n)?",compfile);
	fgets(str,4,stdin);
	if(*str != 'y'){
	  iReturn=-3;
	  goto bailout;
	}
      }
      
      if(!(pCompFileHeader=realloc(pCompFileHeader,pFS->nFileHeaderSize))){
	printf("COM Cannot reallocate for pFileHeader\n");
	iReturn=-4;
	goto bailout;
      }
      memcpy(pCompFileHeader,pFS->pFileHeader,pFS->nFileHeaderSize);
      
      ifr=0;

      if(pFS->iVersion==VERSION_CHUNK){
	ulUnCompFileHeaderSize=pFS->nFileHeaderSize-sizeof(COMP_CHUNK);
	if(!(pSOFTChunk = (SOFT_CHUNK*)FindChunkInBuffer(pCompFileHeader,pFS->nFileHeaderSize,"SOFT"))){
	  printf("COM Cannot find SOFT chunk in comp file header\n");
	  iReturn=5;
	  goto bailout;
	}
	if(pSOFTChunk->Tag[0]!='C'){
	  printf("COM File %s is not compressed\n",pFS->filename);
	  iReturn=-6;
	  goto bailout;
	}

	pSOFTChunk->Tag[0]='T';
	if(strlen(pSOFTChunk->ThisFilename)) 
	  pSOFTChunk->ThisFilename[strlen(pSOFTChunk->ThisFilename)-1]='\0';
	if(strlen(pSOFTChunk->PrevFilename)) 
	  pSOFTChunk->PrevFilename[strlen(pSOFTChunk->PrevFilename)-1]='\0';
	if(strlen(pSOFTChunk->NextFilename)) 
	  pSOFTChunk->NextFilename[strlen(pSOFTChunk->NextFilename)-1]='\0';

	if(!(pCOMPChunk = (COMP_CHUNK*)FindChunkInBuffer(pCompFileHeader,pFS->nFileHeaderSize,"COMP"))){
	  printf("COM Cannot find COMP chunk in comp file header\n");
	  iReturn=5;
	  goto bailout;
	}
	ulNCompFrames=pCOMPChunk->CompressedFrameNumber;
	
	if(pCOMPChunk->CompressedRecordSize!=EXTRA_RECORD_SIZE){
	  printf("COM Not supported Comp Record Size %lu(%iu)\n",pCOMPChunk->CompressedRecordSize,EXTRA_RECORD_SIZE);
	  iReturn=6;
	  goto bailout;
	}

	if(pCOMPChunk->CompressedFrameSize!=compressed_frame_size){
	  printf("COM WARNING: Frame Size mismatch %lu (%lu)\n",pCOMPChunk->CompressedFrameSize,compressed_frame_size);
	  getchar();
	}

	if(!RemoveChunk(pCompFileHeader,pFS->nFileHeaderSize,sizeof(sCOMPChunk),pCOMPChunk)){
 	  printf("COM Cannot remove COMP chunk\n");
	  iReturn=-5;
	  goto bailout;
	}
     }
      else{
	ulUnCompFileHeaderSize=pFS->nFileHeaderSize;
	pFILEHEADER=(FILEHEADER*)pCompFileHeader;
	if(pFILEHEADER->sumTag[0]!='C'){
	  printf("COM File %s is not compressed\n",pFS->filename);
	  iReturn=-6;
	  goto bailout;
	}
	ulNCompFrames=pFILEHEADER->sumCompressedFrameNumber;
	
	pFILEHEADER->sumTag[0]='T';

	if(strlen(pFILEHEADER->sumOrigFilename)) 
	  pFILEHEADER->sumOrigFilename[strlen(pFILEHEADER->sumOrigFilename)-1]='\0';
	if(strlen(pFILEHEADER->sumPrevFilename)) 
	  pFILEHEADER->sumPrevFilename[strlen(pFILEHEADER->sumPrevFilename)-1]='\0';
	if(strlen(pFILEHEADER->sumNextFilename)) 
	  pFILEHEADER->sumNextFilename[strlen(pFILEHEADER->sumNextFilename)-1]='\0';

	pFILEHEADER->sumCompressedFrameSize=0;
	pFILEHEADER->sumCompressedRecordSize=0;
	pFILEHEADER->sumCompressedFrameNumber=0;
      }

      ffr=ulNCompFrames-1;
      original_file_size=ulUnCompFileHeaderSize+original_frame_size*ulNCompFrames+pFS->nFileFooterSize;

      printf("COM File %i %s NF=%lu\n",k,compfile,ulNCompFrames);
      printf("COM Decompressing file %s(%lu frames) to %s\n",pFS->filename,ulNCompFrames,compfile);
      
      if(fseek(pFS->fp,pFS->nFileHeaderSize,SEEK_SET)){
	printf("COM Cannot rewind file %s %i (%i %i)\n",pFS->filename,errno,EBADF,EINVAL); 
	iReturn=10;
	goto bailout;
      }

      if(fread((void*)pframe0,original_frame_size,1,pFS->fp)!=1){
	printf("COM Cannot read ini frame from %s\n",pFS->filename);
	iReturn=11;
	goto bailout;
      }

      if(!(fp=fopen(compfile,"w"))){
	printf("COM Cannot open file %s\n",compfile);
	iReturn=12;
	goto bailout;
      }

      if(fwrite(pCompFileHeader,ulUnCompFileHeaderSize,1,fp)!=1){
	printf("COM Cannot write fileheader to %s\n",compfile);
	iReturn=13;
	goto bailout;
      }
      if(fwrite((void*)pframe0,original_frame_size,1,fp)!=1){
	printf("COM Cannot write ini frame to %s\n",compfile);
	iReturn=14;
	goto bailout;
      }

      for(total_nrecords=nrecords=nrecords_max=0,ul=ifr+1;ul<=ffr;ul++){
	if(fread((void*)pcframe,compressed_frame_size,1,pFS->fp)!=1){
	  printf("COM Cannot read frame %lu from %s\n",ul,pFS->filename);
	  iReturn=15;
	  goto bailout;
	}
	
	for(i=0;i<elements_in_frame;i++){
	  pframe1[i]=(unsigned short)((int)pframe0[i]+(int)pcframe[i]);
	}
	
	if(fread((void*)(&nrecords),RECORD_N_SIZE,1,pFS->fp)!=1){
	  printf("COM Cannot read number of records for frame %lu from %s\n",ul,pFS->filename);
	  iReturn=16;
	  goto bailout;
	}
	//	printf("COM F %i R %i\n",ul,nrecords);
	if(nrecords){
	  if(nrecords_max<nrecords){
	    nrecords_max=nrecords;
	    if(!(precord=(RECORD*)realloc((void*)precord,nrecords*EXTRA_RECORD_SIZE))){
	      printf("COM Cannot realloc for precord %lu for frame %lu\n",nrecords,ul);
	      iReturn=17;
	      goto bailout;
	    }
	  }
	  if(fread((void*)precord,nrecords*EXTRA_RECORD_SIZE,1,pFS->fp)!=1){
	    printf("COM Cannot read records for frame %lu from %s\n",ul,pFS->filename);
	    iReturn=18;
	    goto bailout;
	  }
	  for(i=0;i<nrecords;i++){
	    tb.uc[1]=precord[i].uc[3];
	    precord[i].uc[3]=0;
	    tb.uc[0]=pcframe[precord[i].ul];
	    pframe1[precord[i].ul]=tb.us;
	  }
	  total_nrecords+=nrecords;
	}
	
	if(fwrite((void*)pframe1,original_frame_size,1,fp)!=1){
	  printf("COM Cannot write frame %lu to %s\n",ul,compfile);
	  iReturn=19;
	  goto bailout;
	}
	pframet=pframe0;
	pframe0=pframe1;
	pframe1=pframet;    
      }

      if(pFS->nFileFooterSize){
	if(fwrite(pFS->pFileFooter,pFS->nFileFooterSize,1,fp)!=1){
	  printf("COM Cannot write filefooter to %s\n",compfile);
	  iReturn=20;
	  goto bailout;
	}
      }

      compressed_file_size=pFS->nFileHeaderSize+original_frame_size+(RECORD_N_SIZE+compressed_frame_size)*(ffr-ifr)+total_nrecords*EXTRA_RECORD_SIZE+pFS->nFileFooterSize;
      printf("COM Total records %lu\n",total_nrecords);
      printf("COM Original size %lu Compressed file %lu (%.2f%%)\n",original_file_size,compressed_file_size,100.0*((float)original_file_size-(float)compressed_file_size)/(float)original_file_size);
      fclose(fp);   
      fp=NULL;
    }
  }

 bailout:

  if(fp) fclose(fp);
  if(pCompFileHeader) free(pCompFileHeader);
  if(precord) free(precord);
  if(pcframe) free(pcframe);
  if(pframe0) free(pframe0);
  if(pframe1) free(pframe1);
  return(iReturn);
}
