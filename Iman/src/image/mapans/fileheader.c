#include "fileheader.h"

extern int errno; 


const char strListOfKnownChunks[]={"FRAM cost ISOI SOFT DATA COST COMP HARD ROIS SYNC epst EPST GREE "};
const char strListOfKnownExperimentChunks[]={"COST EPST "};
const char strListOfKnownExperimentFrameChunks[]={"cost epst "};

void SetChunkName(void *pChunk,int iSize,char *pcName){
  memcpy(pChunk,pcName,CHUNK_ID_SIZE);
  *((ULONG*)(pChunk+CHUNK_ID_SIZE))=iSize-CHUNK_HEAD_SIZE;
}

void SetChunkNameTag(void *pChunk,int iSize,char *pcName,char *pcTag){
  SetChunkName(pChunk,iSize,pcName);
  if(iSize>=CHUNK_HEAD_SIZE+CHUNK_TAG_SIZE) 
    memcpy(pChunk+CHUNK_HEAD_SIZE,pcTag,CHUNK_TAG_SIZE);
}

int WriteChunk(void *pChunk,int iSize,FILE* fp){
  if(fwrite(pChunk,iSize,1,fp)!=1){
    printf("FIL Cannot write Chunk\n");
    return(FALSE);
  }
  return(TRUE);
}

int ReadChunk(void *pChunk,int iSize,FILE* fp){
  if(fread(pChunk,iSize,1,fp)!=1){
    printf("FIL Cannot read Chunk\n");
    return(FALSE);
  }
  return(TRUE);
}

int ReadChunkHeader(void *pChunk,FILE* fp){
  if(fread(pChunk,CHUNK_HEAD_SIZE,1,fp)!=1){
    printf("FIL Cannot read Chunk header\n");
    return(FALSE);
  }
  if(fseek(fp,-CHUNK_HEAD_SIZE,SEEK_CUR)==-1){
    printf("FIL Cannot fseek file Error=%i\n",errno);
    return(FALSE);
  }
  return(TRUE);
}

int FileOffsetCurrent(long lOffset,FILE* fp){
  if(fseek(fp,lOffset,SEEK_CUR)==-1){
    printf("FIL Cannot fseek file Error=%i\n",errno);
    return(FALSE);
  }
  return(TRUE);
}

int HasName(void *pChunk,char *pcName){
  if(memcmp(pcName,pChunk,CHUNK_ID_SIZE)) return(FALSE);
  else return(TRUE);
}

void PrintChunkName(void *pChunk){
  char str[CHUNK_ID_SIZE+1];
  memcpy(str,pChunk,CHUNK_ID_SIZE);
  str[CHUNK_ID_SIZE]='\0';
  printf("FIL Chunk %s\n",str);
}

int IsChunkKnown(char *pcName){
  char *pc;
  for(pc=(char*)strListOfKnownChunks;*pc;pc+=CHUNK_ID_SIZE+1) 
    if(!strncmp(pc,pcName,CHUNK_ID_SIZE)) return(TRUE);
  return(FALSE);
}

void* FindChunkInBuffer(void *pMultiChunk,ULONG ulMultiChunkSize,char *strName){
  BYTE *pb;
  pb=(BYTE*)pMultiChunk;
  if(!memcmp("ISOI",pb,CHUNK_ID_SIZE)){
    if(!memcmp(strName,"ISOI",CHUNK_ID_SIZE)) return(pMultiChunk);
    pb+=sizeof(ISOI_CHUNK);
  }
  for(;(pb < (BYTE*)pMultiChunk+ulMultiChunkSize) && memcmp(strName,pb,CHUNK_ID_SIZE);
      pb+=CHUNK_HEAD_SIZE+*(ULONG*)(pb+CHUNK_ID_SIZE)){
    //    printf("FIL Found %c%c%c%c %lu\n",*pb,*(pb+1),*(pb+2),*(pb+3),*(ULONG*)(pb+4));
  }
  if(pb < (BYTE*)pMultiChunk+ulMultiChunkSize) return((void*)pb);
  else return(NULL);
}

int FindExperimentChunkInFile(FILE *fp){
  DATA_CHUNK DATAChunk;
  long lOriginalPos;
  void *pv;
  char *pc;
  pv=(void*)&DATAChunk;
  lOriginalPos=ftell(fp);
  rewind(fp);

  while(!feof(fp)){
    if(fread(pv,CHUNK_HEAD_SIZE,1,fp)!=1){
      if(!feof(fp)) printf("FIL Cannot read Chunk header\n");
      fseek(fp,lOriginalPos,SEEK_SET);
      return(-1);
    }

    if(!memcmp("ISOI",pv,CHUNK_ID_SIZE)){
      fseek(fp,CHUNK_TAG_SIZE,SEEK_CUR);
    }
    else{
      for(pc=(char*)strListOfKnownExperimentChunks;*pc;pc+=CHUNK_ID_SIZE+1) 
	if(!strncmp(pc,DATAChunk.ID,CHUNK_ID_SIZE)) 
	  return((pc-strListOfKnownExperimentChunks)/(CHUNK_ID_SIZE+1));
      fseek(fp,DATAChunk.Size,SEEK_CUR);
    }
  }
  fseek(fp,lOriginalPos,SEEK_SET);
  return(-1);
}

unsigned long FindChunkInFileAndOffset(FILE *fp,char *pcName,int iJump){
  DATA_CHUNK DATAChunk;
  long lOriginalPos;
  void *pv;

  if(!IsChunkKnown(pcName)){
    printf("FIL Unknown ini chunk %c%c%c%c\n",pcName[0],pcName[1],pcName[2],pcName[3]);
    return(~0UL);
  }
  else{
    //printf("FIL FindChunkInFileAndOffset %c%c%c%c Jump %i\n",pcName[0],pcName[1],pcName[2],pcName[3],iJump);
  }

  pv=(void*)&DATAChunk;
  lOriginalPos=ftell(fp);
  
  while(!feof(fp)){
    if(fread(pv,CHUNK_HEAD_SIZE,1,fp)!=1){
      printf("FIL Cannot read Chunk header\n");
      fseek(fp,lOriginalPos,SEEK_SET);
      return(~0UL);
    }
    //printf("FIL Chunk size %i\n",DATAChunk.Size);

    if(!IsChunkKnown((char*)pv)){
      printf("FIL Unknown chunk |%c%c%c%c| %i%i%i%i\n",DATAChunk.ID[0],DATAChunk.ID[1],DATAChunk.ID[2],DATAChunk.ID[3],DATAChunk.ID[0],DATAChunk.ID[1],DATAChunk.ID[2],DATAChunk.ID[3]);
      fseek(fp,lOriginalPos,SEEK_SET);
      return(~0UL);
    } 

    if(!memcmp("ISOI",pv,CHUNK_ID_SIZE)){
      if(!memcmp(pcName,"ISOI",CHUNK_ID_SIZE)){
	if(iJump){
	  if(fseek(fp,DATAChunk.Size,SEEK_CUR)==-1){
	    printf("FIL Cannot jump fseek Error=%i\n",errno);
	    fseek(fp,lOriginalPos,SEEK_SET);
	    return(~0UL);
	  }
	}
	return(DATAChunk.Size);
      }
      fseek(fp,CHUNK_TAG_SIZE,SEEK_CUR);
    }
    else{
      if(!memcmp(pcName,pv,CHUNK_ID_SIZE)){
	if(iJump){
	  if(fseek(fp,DATAChunk.Size,SEEK_CUR)==-1){
	    printf("FIL Cannot jump fseek Error=%i\n",errno);
	    fseek(fp,lOriginalPos,SEEK_SET);
	    return(~0UL);
	  }
	}
	return(DATAChunk.Size);
      }
      fseek(fp,DATAChunk.Size,SEEK_CUR);
    }
  }
  fseek(fp,lOriginalPos,SEEK_SET);
  return(FALSE);
}

unsigned long FindChunkInFile(FILE *fp,char *pcName){
  unsigned long ulReturn;
  if((ulReturn=FindChunkInFileAndOffset(fp,pcName,FALSE))!=~0UL) fseek(fp,-CHUNK_HEAD_SIZE,SEEK_CUR);
  return(ulReturn);
}

int InsertChunk(void *pMultiChunk,ULONG ulMultiChunkSize,void *pChunk,ULONG ulChunkSize,void *pMultiChunkEntry){
  long lHeadBytes;
  lHeadBytes=pMultiChunkEntry-pMultiChunk;
  if(lHeadBytes<0){
    printf("FIL Insert entry is higher than chunk head");
    return(FALSE);
  }
  if(lHeadBytes>ulMultiChunkSize){
    printf("FIL Insert entry is lower than chunk tail");
    return(FALSE);
  }
  if(lHeadBytes!=ulMultiChunkSize) 
    memmove(pMultiChunkEntry+ulChunkSize, pMultiChunkEntry, ulMultiChunkSize-lHeadBytes);
  memcpy(pMultiChunkEntry,pChunk,ulChunkSize);
  return(TRUE);
}

int RemoveChunk(void *pMultiChunk,ULONG ulMultiChunkSize,ULONG ulChunkSize,void *pMultiChunkEntry){
  long lHeadBytes;
  lHeadBytes=pMultiChunkEntry-pMultiChunk;
  if(lHeadBytes<0){
    printf("FIL Remove entry is higher than Multichunk head");
    return(FALSE);
  }
  if(lHeadBytes>ulMultiChunkSize){
    printf("FIL Remove entry is lower than Multichunk tail");
    return(FALSE);
  }
  if(lHeadBytes>ulMultiChunkSize-ulChunkSize){ 
    printf("FIL Remove chunk tail is lower than Multichunk tail");
    return(FALSE);
  }
  if(lHeadBytes<ulMultiChunkSize-ulChunkSize)
    memmove(pMultiChunkEntry, pMultiChunkEntry+ulChunkSize, ulMultiChunkSize-lHeadBytes-ulChunkSize);
  return(TRUE);
}

int DataTypeToSizeOf(int iDataType){
  switch(iDataType){
  case DATATYPE_UCHAR:
    return(sizeof(unsigned char));
    break;
  case DATATYPE_USHORT:
    return(sizeof(unsigned short));
    break;
  case DATATYPE_ULONG:
    return(sizeof(ULONG));
    break;
  case DATATYPE_FLOAT:
    return(sizeof(float));
    break;
  default:
    printf("DataTypeToSizeOf: Unknown data type %i",iDataType);
  }
  return(-1);
}
