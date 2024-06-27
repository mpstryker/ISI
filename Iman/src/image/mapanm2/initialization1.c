#include "mapanm21.h"

int InitializeMaps(){
  FILE *fp;
  int i,j,dumi;
  unsigned short dumus;
  double dumd,*pd;
  int map_count;
  char *pcx,*pcy;
  unsigned long ulXYDCrop;
  int iXYShift;
  int iXDCrop,iYDCrop;

  maps=(double**)calloc(MAX_FILE_N,sizeof(double*));

  for(i=0,map_count=0;i<nfiles;i++){
    if((fp=fopen(filenames[i],"r"))==NULL){
      printf("INI Cannot open file %s\n",filenames[i]);
      return(1);
    }
    if(fread((void*)&dumi,sizeof(int),1,fp) != 1){
      printf("INI Cannot read type from file %s\n",filenames[i]);
      return(2);
    }
    if(i==0){
      type=dumi;
    }
    else{
      if((SAVE_TYPE_BIT_DOUBLE_FLOAT&dumi)!=(SAVE_TYPE_BIT_DOUBLE_FLOAT&type)){
	printf("INI Data type mismatch\n");
	return(4);
      }
    }

    if(map_count==0){
      strncpy(savefile,filenames[0],8);
      strncat(savefile,"-",1);
      filestr1=strdup(filenames[0]);
      if((dumi&SAVE_TYPE_BIT_XY_IN_FILES)){
	pcx=strstr(filestr1,"mapraw");
	if(pcx) memmove(pcx,pcx+6,strlen(pcx+6)+1);
      }
      else{
	pcx=strstr(filestr1,"mapxraw");
	pcy=strstr(filestr1,"mapyraw");
	if(pcx) memmove(pcx,pcx+7,strlen(pcx+7)+1);
	if(pcy){ 
	  memmove(pcy,pcy+7,strlen(pcy+7)+1);
	  printf("INI WARNING: Y component is input before X in %s\n",filenames[0]);
	}
      }
    }
    if(map_count==2){
      strncat(savefile,filenames[i]+6,2);
      filestr2=strdup(filenames[i]);
      if((dumi&SAVE_TYPE_BIT_XY_IN_FILES)){
	pcx=strstr(filestr2,"mapraw");
	if(pcx) memmove(pcx,pcx+6,strlen(pcx+6)+1);
      }
      else{
	pcx=strstr(filestr2,"mapxraw");
	pcy=strstr(filestr2,"mapyraw");
	if(pcx) memmove(pcx,pcx+7,strlen(pcx+7)+1);
	if(pcy){
	  memmove(pcy,pcy+7,strlen(pcy+7)+1);
	  printf("INI WARNING: Y component is input before X in %s\n",filenames[2]);
	}
      }
    }
    if(fread((void*)&dumus,sizeof(unsigned short),1,fp) != 1){
      printf("INI Cannot read Xdim from file %s\n",filenames[i]);
      return(2);
    }
    if(i==0){
      Xdim=dumus;
    }
    else{
      if(dumus!=Xdim){
	printf("INI Xdim mismatch %i(%i)\n",Xdim,dumus);
	return(4);
      }
    }
    if(fread((void*)&dumus,sizeof(unsigned short),1,fp) != 1){
      printf("INI Cannot read Ydim from file %s\n",filenames[i]);
      return(2);
    }
    if(i==0){
      Ydim=dumus;
      printf("INI Type=%i Xdim=%i Ydim=%i\n",type,Xdim,Ydim);
      XYdim=(unsigned long)Xdim*(unsigned long)Ydim;
    }
    else{
      if(dumus!=Ydim){
	printf("INI Ydim mismatch %i(%i)\n",Ydim,dumus);
	return(4);
      }
    }
    printf("INI %i\n",dumi&SAVE_TYPE_BIT_XY_IN_FILES);
    if((dumi&SAVE_TYPE_BIT_XY_IN_FILES)==SAVE_XY_IN_ONE_FILE*SAVE_TYPE_BIT_XY_IN_FILES){
      if(!(maps[map_count]=(double *)calloc(XYdim,sizeof(double)))){
	printf("INI Cannot allocate for maps[%i]\n",map_count);
	return(3);
      }
      if(fread((void*)maps[map_count],XYdim*sizeof(double),1,fp) != 1){
	printf("INI Cannot read maps[%i] from file %s\n",map_count,filenames[i]);
	return(2);
      }
      map_count++;
      if(!(maps[map_count]=(double *)calloc(XYdim,sizeof(double)))){
	printf("INI Cannot allocate for maps[%i]\n",map_count);
	return(3);
      }
      if(fread((void*)maps[map_count],XYdim*sizeof(double),1,fp) != 1){
	printf("INI Cannot read maps[%i] from file %s\n",map_count,filenames[i]);
	return(2);
      }
      map_count++;
    }
    else{
      if(!(maps[map_count]=(double *)calloc(XYdim,sizeof(double)))){
	printf("INI Cannot allocate for maps[%i]\n",map_count);
	return(3);
      }
      if(fread((void*)maps[map_count],XYdim*sizeof(double),1,fp) != 1){
	printf("INI Cannot read maps[%i] from file %s\n",map_count,filenames[i]);
	return(2);
      }
      map_count++;
    }

    if((SAVE_TYPE_BIT_ADD_DOUBLE_AT_END&dumi)){
      if(fread((void*)&dumd,sizeof(double),1,fp) != 1){
	printf("INI Cannot read harmonic from file %s\n",filenames[i]);
	type=type&(!((int)(SAVE_TYPE_BIT_ADD_DOUBLE_AT_END)));
	harmonic_d=1.0;
      }
      if(i==0){
	harmonic_d=dumd;
	printf("INI Harmonic %f\n",harmonic_d);
      }
      else{
	if(harmonic_d!=dumd){
	  printf("INI Harmonic mismatch %f %f(%i)\n",harmonic_d,dumd,i);
	  type=type&(!((int)(SAVE_TYPE_BIT_ADD_DOUBLE_AT_END)));
	  harmonic_d=1.0;
	}
      }
    }
    fclose(fp);
  }

  if(iCropROI){
    i=0;
    if(iROIXLeft<0 || iROIXLeft>=Xdim){
      i++;
      printf("GET Incorrect XLeft for ROI %i\n",iROIXLeft);      
    }
    if(iROIXRight<0 || iROIXRight>=Xdim){
      i++;
      printf("GET Incorrect XRight for ROI %i\n",iROIXRight);      
    }
    if(iROIYTop<0 || iROIYTop>=Ydim){
      i++;
      printf("GET Incorrect YTop for ROI %i\n",iROIYTop);
    }
    if(iROIYBottom<0 || iROIYBottom>=Xdim){
      i++;
      printf("GET Incorrect YBottom for ROI %i\n",iROIYBottom);
    }
    if(i){
      printf("GET Incorrect ROI\n");
      iCropROI=0;
    }
    else{
      iXDCrop=iROIXRight-iROIXLeft+1;
      iYDCrop=iROIYBottom-iROIYTop+1;
      ulXYDCrop=(unsigned long)iXDCrop*(unsigned long)iYDCrop;
      iXYShift=iROIXLeft+iROIYTop*Xdim;
      for(j=0;j<map_count;j++){
	pd=maps[j];
	for(i=0;i<ulXYDCrop;i++){
	  dumi=iXYShift+i%iXDCrop+(i/iXDCrop)*Xdim;
	  pd[i]=pd[dumi];
	}
	maps[j]=(double *)realloc(pd,ulXYDCrop*sizeof(double));
      }
      Xdim=iXDCrop;
      Ydim=iYDCrop;
      XYdim=ulXYDCrop;
      printf("INI ROI Xdim=%i Ydim=%i\n",Xdim,Ydim);
    }
  }

  printf("INI Save file - %s\n",savefile);
 
  pminx=0.0;
  pmaxx=Xdim-1.0;
  pmaxy=0.0;
  pminy=Ydim-1.0;

  printf("INI Maps initialized (%i)\n",map_count);
  return(0);
}

int InitializeBuffers(){
  int i,j,k;

  if(!(fbuffer=(float*)malloc(XYdim*sizeof(float)))){
    printf("INI Cannot allocate for fbuffer\n");
    return(1);
  }
  if(!(fbufferr=(float*)malloc(XYdim*sizeof(float)))){
    printf("INI Cannot allocate for fbufferr\n");
    return(1);
  }
  if(!(fbufferi=(float*)malloc(XYdim*sizeof(float)))){
    printf("INI Cannot allocate for fbufferi\n");
    return(1);
  }
  if(!(fbufferx=(float *)calloc(XYdim,sizeof(float)))){
    printf("INI Cannot allocate for fbufferx\n");
    return(2);
  }
  if(!(fbuffery=(float *)calloc(XYdim,sizeof(float)))){
    printf("INI Cannot allocate for fbuffery\n");
    return(3);
  }
  if(!(dbufferx=(double *)calloc(XYdim,sizeof(double)))){
    printf("INI Cannot allocate for dbufferx\n");
    return(4);
  }
  if(!(dbuffery=(double *)calloc(XYdim,sizeof(double)))){
    printf("INI Cannot allocate for dbuffery\n");
    return(5);
  }
  if(!(dbufferphix=(double *)calloc(XYdim,sizeof(double)))){
    printf("INI Cannot allocate for dbufferphix\n");
    return(5);
  }
  if(!(dbufferphiy=(double *)calloc(XYdim,sizeof(double)))){
    printf("INI Cannot allocate for dbufferphiy\n");
    return(5);
  }

  if(!(pdBufferRho1=(double *)calloc(XYdim,sizeof(double)))){
    printf("INI Cannot allocate for pdBufferRho1\n");
    return(5);
  }
  if(!(pdBufferRho2=(double *)calloc(XYdim,sizeof(double)))){
    printf("INI Cannot allocate for pdBufferRho2\n");
    return(5);
  }

  dbuffers=(double**)calloc(sizeof(double*),MAX_FILE_N);
  for(i=0;i<MAX_FILE_N;i++){
    if(!(dbuffers[i]=(double *)calloc(XYdim,sizeof(double)))){
      printf("INI Cannot allocate for dbuffers[%i]\n",i);
      return(4);
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

  printf("INI Buffers initialized\n");
  return(0);
}
