/* Bin analysis                 */
/* Written by V.Kalatsky        */
/* 20Jul2001                    */

#include "binan1.h"

#define SAVE_TYPE_BIT_DOUBLE_FLOAT 1
#define SAVE_TYPE_BIT_ADD_DOUBLE_AT_END 2
#define SAVE_TYPE_BIT_XY_IN_FILES 4

//#define TEST

int InitializeMaps(){
  FILE *fp;
  int i,j,dumi;
  unsigned short dumus;
  double *binbuffer,*pd_s,*pd_d,*pd_map;
  double *pdX=NULL,*pdY=NULL,*pd;
  unsigned short xd,yd;
  unsigned long xy;
  int x,y,xb,yb;
  int iTotalSpatialBinning;
#ifdef TEST
  int k,l;
  double kx,ky,omega;
#endif


  if(nfiles==1){
    for(i=0;i<nfiles;i++){
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
	if(dumi!=type){
	  printf("INI Data type mismatch\n");
	  return(4);
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
      }
      else{
	if(dumus!=Ydim){
	  printf("INI Ydim mismatch %i(%i)\n",Ydim,dumus);
	  return(4);
	}
      }
      if(fread((void*)&nbins,sizeof(int),1,fp) != 1){
	printf("INI Cannot read nbins from file %s\n",filenames[i]);
	return(2);
      }
      printf("INI Nbins=%i\n",nbins);
      if(fread((void*)&harmonic_d,sizeof(double),1,fp) != 1){
	printf("INI Cannot read harmonic from file %s\n",filenames[i]);
      }
      
      maps_orig=(double**)calloc(sizeof(double*),nbins);
      bin_max=(float*)calloc(sizeof(float),nbins);
      bin_min=(float*)calloc(sizeof(float),nbins);
      
      for(j=0;j<nbins;j++){
	if(!(maps_orig[j]=(double *)calloc((unsigned long)Xdim*(unsigned long)Ydim,sizeof(double)))){
	  printf("INI Cannot allocate for maps_orig[%i]\n",j);
	  return(3);
	}
	if(fread((void*)(maps_orig[j]),(unsigned long)Xdim*(unsigned long)Ydim*sizeof(double),1,fp) != 1){
	  printf("INI Cannot read maps_orig[%i] from file %s\n",j,filenames[i]);
	}
      }
      fclose(fp);
    }
  }
  else{
    printf("INI Initializing from %i files\n",nfiles);
    nbins=nfiles;
    maps_orig=(double**)calloc(sizeof(double*),nbins);
    bin_max=(float*)calloc(sizeof(float),nbins);
    bin_min=(float*)calloc(sizeof(float),nbins);

    for(i=0;i<nfiles;i++){
      if((fp=fopen(filenames[i],"r"))==NULL){
	printf("INI Cannot open file %i <%s>\n",i,filenames[i]);
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
	if(dumi!=type){
	  printf("INI Data type mismatch\n");
	  return(4);
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
	xy=(unsigned long)Xdim*(unsigned long)Ydim;
      }
      else{
	if(dumus!=Ydim){
	  printf("INI Ydim mismatch %i(%i)\n",Ydim,dumus);
	  return(4);
	}
      }
            
      if(i==0){
	if(!(pdX=(double *)calloc(xy,sizeof(double)))){
	  printf("INI Cannot allocate for pdX\n");
	  return(3);
	}
	if(!(pdY=(double *)calloc(xy,sizeof(double)))){
	  printf("INI Cannot allocate for pdY\n");
	  return(3);
	}
      }

      if(!(maps_orig[i]=(double *)calloc(xy,sizeof(double)))){
	printf("INI Cannot allocate for maps_orig[%i]\n",i);
	return(3);
      }

      if(fread((void*)(pdX),xy*sizeof(double),1,fp) != 1){
	printf("INI Cannot read pdX from file %s\n",filenames[i]);
      }
      if(fread((void*)(pdY),xy*sizeof(double),1,fp) != 1){
	printf("INI Cannot read pdY from file %s\n",filenames[i]);
      }

      if(i==0){
	if((type&SAVE_TYPE_BIT_ADD_DOUBLE_AT_END)){
	  if(fread((void*)&harmonic_d,sizeof(double),1,fp) != 1){
	    printf("INI Cannot read harmonic from file %s\n",filenames[i]);
	    harmonic_d=1.0;
	  }
	  else harmonic_d=1.0;
	}
	else harmonic_d=1.0;
      }
      
      fclose(fp);

      pd=maps_orig[i];
      for(j=0;j<xy;j++){
	pd[j]=hypot(pdY[j],pdX[j]);
      }
    }
    printf("INI Initialaized from %i maps\n",nfiles);
    SAFE_FREE(pdX);
    SAFE_FREE(pdY);
  }

  if(spatial_binning>1){
    printf("INI Running %ix%i binning\n",spatial_binning,spatial_binning);
    xd=Xdim/spatial_binning;
    yd=Ydim/spatial_binning;
    xy=(unsigned long)xd*(unsigned long)yd;
    iTotalSpatialBinning=spatial_binning*spatial_binning;
    if((binbuffer=(double*)malloc(xy*sizeof(double)))){
      for(j=0;j<nbins;j++){
	memset(binbuffer,0,xy*sizeof(double));
	pd_map=maps_orig[j];
	for(y=0;y<yd;y++){
	  for(x=0;x<xd;x++){
	    pd_s=pd_map+(x+y*Xdim)*spatial_binning;
	    pd_d=binbuffer+(x+y*xd);
	    for(yb=0;yb<spatial_binning;yb++){
	      for(xb=0;xb<spatial_binning;xb++) *pd_d += *(pd_s+xb+yb*Xdim);
	    }
	    if(iNormalizeBinning) *pd_d /= iTotalSpatialBinning;
	  }
	}
	maps_orig[j]=(double*)realloc(maps_orig[j],xy*sizeof(double));
	memcpy(maps_orig[j],binbuffer,xy*sizeof(double));
      }
      Xdim=xd;
      Ydim=yd;
    }
    else{
      printf("INI Cannot allocate for binbuffer. Will go on without binning.\n");
    }
  }

#ifdef TEST

  Xdim=101;
  Ydim=101;
  kx=0.2;
  ky=0.1;
  omega=2.0*M_PI*1.0/(double)nbins;
  for(l=0;l<nbins;l++){
    for(k=0;k<Ydim;k++){
      for(j=0;j<Xdim;j++){
	maps_orig[l][k*Xdim+j]=cos(kx*(double)(j-Xdim/2)+ky*(double)(k-Ydim/2)-omega*(double)l);
      } 
    }
  }

#endif

  XYdim=(unsigned long)Xdim*(unsigned long)Ydim;
  pminx=-0.5;
  pmaxx=(float)Xdim-0.5;
  pmaxy=-0.5;
  pminy=(float)Ydim-0.5;
  //pminy=-0.5;
  //pmaxy=(float)Ydim-0.5;
  
  panelnx=nbins;
  panelny=1;
  for(j=nbins;j<nbins+3;j++){
    for(i=(int)floor(sqrt((double)j));i>0;i--){
      if(!(j%i)){
	if((j/i-i)<(panelnx-panelny)){
	  panelnx=j/i;
	  panelny=i;
	}
      }
    }
  }

  printf("INI Maps initialized\n");
  return(0);
}

int InitializeBuffers(){
  int i;
  unsigned long ul;
  float dumf;

  //  printf("INI InitializeBuffers\n");
  maps=(double**)calloc(sizeof(double*),nbins);
  for(i=0;i<nbins;i++){
    if(!(maps[i]=(double *)calloc(XYdim,sizeof(double)))){
      printf("INI Cannot allocate for maps[%i]\n",i);
      return(4);
    }
    memcpy(maps[i],maps_orig[i],XYdim*sizeof(double));
  }

  fmaps_orig=(float**)calloc(sizeof(float*),nbins);
  for(i=0;i<nbins;i++){
    if(!(fmaps_orig[i]=(float *)calloc(XYdim,sizeof(float)))){
      printf("INI Cannot allocate for fmaps_orig[%i]\n",i);
      return(4);
    }
    Double2FloatShort(XYdim,maps_orig[i],fmaps_orig[i]);
  }

  AverageMaps(Xdim,Ydim,maps,(float**)0,nbins,dradius,dradiusbig,(int)DOUBLE);
  // Val: quick hack, should commented out
  // AverageMaps(Xdim,Ydim,maps_orig,(float**)0,nbins,dradius,dradiusbig,(int)DOUBLE);
  printf("INI Maps averaged\n");

  fmaps=(float**)calloc(sizeof(float*),nbins);
  for(i=0;i<nbins;i++){
    if(!(fmaps[i]=(float *)calloc(XYdim,sizeof(float)))){
      printf("INI Cannot allocate for fmaps[%i]\n",i);
      return(4);
    }
    Double2FloatLong(XYdim,maps[i],fmaps[i],bin_min+i,bin_max+i);
  }

  if(shift_max_to_zero){
    for(ul=0;ul<XYdim;ul++){
      dumf=fmaps[0][ul];
      for(i=1;i<nbins;i++){
	if(dumf<fmaps[i][ul]) dumf=fmaps[i][ul];
      }
      for(i=0;i<nbins;i++){
	fmaps[i][ul]-=dumf;
      }
    } 
    for(i=0;i<nbins;i++){
      FindExtremaEx(XYdim,(unsigned long)Xdim,fmaps[i],bin_min+i,bin_max+i,&ul,&ul);
      if(do_verbose) printf("INI %i Min %f Max %f\n",i,bin_min[i],bin_max[i]);
    }
  }

  if(shift_average_to_zero){
    for(i=0;i<nbins;i++){
      dumf=FindAverage(XYdim,fmaps[i],(float*)0);
      if(do_verbose) printf("INI %i Ave %f\n",i,dumf);
      for(ul=0;ul<XYdim;ul++){
	fmaps[i][ul]-=dumf;
      }      
      FindExtremaEx(XYdim,(unsigned long)Xdim,fmaps[i],bin_min+i,bin_max+i,&ul,&ul);
      if(do_verbose) printf("INI %i Min %f Max %f\n",i,bin_min[i],bin_max[i]);
    }
  }
  
  if(!(fbufferi=(float*)malloc(XYdim*sizeof(float)))){
    printf("INI Cannot allocate for fbufferi\n");
    return(1);
  }

  if(!(fbuffer=(float*)malloc(XYdim*sizeof(float)))){
    printf("INI Cannot allocate for fbuffer\n");
    return(1);
  }

  stenciln=XYdim;
  if(!(stencil=(unsigned long *)malloc(stenciln*sizeof(unsigned long)))){
    printf("INI Cannot allocate for stencil\n");
    return(1);      
  }
  for(ul=0;ul<stenciln;ul++) stencil[ul]=ul;

  smaxx=Xdim-1;
  sminx=0;
  smaxy=Ydim-1;
  sminy=0;

  radius=(int)floor(dradius);
  radiusbig=(int)floor(dradiusbig);
  cradius=(int)floor(dcradius);

  if(iGenerateMaps){
    if(!(pdMapX=(double*)malloc(XYdim*sizeof(double)))){
      printf("INI Cannot allocate for pdMapX\n");
      return(1);
    }
    if(!(pdMapY=(double*)malloc(XYdim*sizeof(double)))){
      printf("INI Cannot allocate for pdMapY\n");
      return(1);
    }
  }

  printf("INI Buffers initialized\n");
  return(0);
}
