#include "iman1.h"

int CocktailAverage(TRAIN *tr){
  int i,j;
  unsigned short us;
  char filen[256];
  FILE *fp;
  struct stat buf;
  double *pca;
  float dumf;
  CAR *pcar;

  sprintf(filen,"%s%s",tr->filename,".ca");
  if(!stat(filen,&buf)){
    printf("Good! Got cocktailaverage saved in file\n");
    if((fp=fopen(filen,"r"))==NULL){
      printf("Cannot open file %s\n",filen);
      return(1);
    }
    if(fread((void*)&us,sizeof(unsigned short),1,fp) != 1){
      printf("Cannot read Xdim from file %s\n",filen);
      return(2);
    }
    if(fread((void*)&us,sizeof(unsigned short),1,fp) != 1){
      printf("Cannot read Ydim from file %s\n",filen);
      return(2);
    }
    if(fread((void*)&camin,sizeof(float),1,fp) != 1){
      printf("Cannot read camin from file %s\n",filen);
      return(2);
    }
    if(fread((void*)&camax,sizeof(float),1,fp) != 1){
      printf("Cannot read camax from file %s\n",filen);
      return(2);
    }
    if(fread((void*)cocktailaverage,XYdim*sizeof(float),1,fp) != 1){
      printf("Cannot read cocktailaverage from file %s\n",filen);
      return(2);
    }
    fclose(fp);
  }
  else{
    if(!(pca=(double *)calloc(XYdim,sizeof(double)))){
      printf("Cannot allocate for pca\n");
      return(3);
    }
    for(i=0;i<tr->max_n;i++){
      pcar=AddFrame(tr,i,AVERAGE);
      for(j=0;j<XYdim;j++){
	pca[j]+=(double)pcar->pimage[j];
      }
      if(!(i%10)) printf("%i\n",i);
    }
    camin=1000000000.0;
    camax=0.0;
    for(j=0;j<XYdim;j++){
      dumf=cocktailaverage[j]=pca[j]/tr->max_n;
      if(camin>dumf) camin=dumf;
      if(camax<dumf) camax=dumf;
    }
    free(pca);
    if((fp=fopen(filen,"w"))==NULL){
      printf("Cannot open file %s\n",filen);
      return(4);
    }
    if(fwrite((void*)&Xdim,sizeof(unsigned short),1,fp) != 1){
      printf("Cannot write Xdim to file %s\n",filen);
      return(5);
    }
    if(fwrite((void*)&Ydim,sizeof(unsigned short),1,fp) != 1){
      printf("Cannot write Ydim to file %s\n",filen);
      return(5);
    }
    if(fwrite((void*)&camin,sizeof(float),1,fp) != 1){
      printf("Cannot write camin to file %s\n",filen);
      return(5);
    }
    if(fwrite((void*)&camax,sizeof(float),1,fp) != 1){
      printf("Cannot write camax to file %s\n",filen);
      return(5);
    }
    if(fwrite((void*)cocktailaverage,XYdim*sizeof(float),1,fp) != 1){
      printf("Cannot write cocktailaverage to file %s\n",filen);
      return(5);
    }
    fclose(fp);
  }
  printf("Completed cocktailaverage\n");
  return(0);
}

void DisplayCocktailAverage(TRAIN *tr){
  cpggray(cocktailaverage,Xdim,Ydim,1,Xdim,1,Ydim,camax,camin,trans);
}

