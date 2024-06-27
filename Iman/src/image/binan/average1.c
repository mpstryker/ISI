/* Map analysis                 */
/* Written by V.Kalatsky        */
/* 20Jul2001                    */

#include "binan1.h"

void AverageMaps(int xd,int yd,double **dm,float **fm,int nn,double dr,double drb,int mode){
  int i,j,k,m,n,x,y,x0,y0;
  int r,rb;
  unsigned long xy;
  double dumd;
  double *localdbuffer;
  float dumf;
  float *localfbuffer;
  int *aindex_x,*aindexbig_x,*aindex_y,*aindexbig_y;
  int aindex_n,aindexbig_n;

  r=(int)dr;
  rb=(int)drb;

  if(r> ((xd>yd) ? xd/2 : yd/2)){
    printf("AVE Low-pass radius is too large (%i)\n",r);
    return;
  }

  xy=(unsigned long)xd*(unsigned long)yd;
  if(r>0){
    if(!(aindex_x=(int *)calloc((1+2*r)*(1+2*r),sizeof(int))) || 
       !(aindex_y=(int *)calloc((1+2*r)*(1+2*r),sizeof(int)))){
      printf("AVE Cannot allocate for aindex\n");
      return;
    }
    aindex_n=0;
    for(j=-r;j<=r;j++){
      for(k=-r;k<=r;k++){
	if(hypot((double)j,(double)k)<=dr){ 
	  aindex_x[aindex_n]=k;
	  aindex_y[aindex_n++]=j;
	}
      }
    }
    printf("AVE Low-pass averaging over %i pixels\n",aindex_n);
    switch(mode){
    case DOUBLE:
      if(!(localdbuffer=(double *)calloc(xy,sizeof(double)))){
	printf("INI Cannot allocate for localdbuffer\n");
	return;
      }
      for(n=0;n<nn;n++){
	for(i=0;i<xy;i++){
	  dumd=0.0;
	  m=0;
	  x0=i%xd;
	  y0=i/xd;
	  for(j=0;j<aindex_n;j++){
	    x=x0+aindex_x[j];
	    y=y0+aindex_y[j];
	    if(x>=0 && x<xd && y>=0 && y<yd){ 
	      dumd+=dm[n][x+y*xd];
	      m++;
	    }
	  }
	  localdbuffer[i]=dumd/(double)m;
	}
	memcpy(dm[n],localdbuffer,xy*sizeof(double));    
      }    
      free(localdbuffer);
      break;
    case FLOAT:
      if(!(localfbuffer=(float *)calloc(xy,sizeof(float)))){
	printf("INI Cannot allocate for localfbuffer\n");
	return;
      }
      for(n=0;n<nn;n++){
	for(i=0;i<xy;i++){
	  dumf=0.0;
	  m=0;
	  x0=i%xd;
	  y0=i/xd;
	  for(j=0;j<aindex_n;j++){
	    x=x0+aindex_x[j];
	    y=y0+aindex_y[j];
	    if(x>=0 && x<xd && y>=0 && y<yd){ 
	      dumf+=fm[n][x+y*xd];
	      m++;
	    }
	  }
	  localfbuffer[i]=dumf/(float)m;
	}
	memcpy(fm[n],localfbuffer,xy*sizeof(float));    
      }    
      free(localfbuffer);      
      break;
    default:
      printf("AVE Unknown mode %i\n",mode);
      return;
    }
    free(aindex_x);
    free(aindex_y);
  }

  if(rb> ((xd>yd) ? xd/2 : yd/2)){
    printf("AVE High-pass radius is too large (%i)\n",rb);
    return;
  }

  if(rb>0){
    if(!(aindexbig_x=(int *)calloc((1+2*rb)*(1+2*rb),sizeof(int))) ||
       !(aindexbig_y=(int *)calloc((1+2*rb)*(1+2*rb),sizeof(int)))){
      printf("AVE Cannot allocate for aindexbig\n");
      return;
    }
    aindexbig_n=0;
    for(j=-rb;j<=rb;j++){
      for(k=-rb;k<=rb;k++){
	if(hypot((double)j,(double)k)<=drb){
	  aindexbig_x[aindexbig_n]=k;
	  aindexbig_y[aindexbig_n++]=j;
	}
      }
    }
    printf("AVE High-pass averaging over %i pixels\n",aindexbig_n);

    switch(mode){
    case DOUBLE:
      if(!(localdbuffer=(double *)calloc(xy,sizeof(double)))){
	printf("INI Cannot allocate for localdbuffer\n");
	return;
      }
      for(n=0;n<nn;n++){
	for(i=0;i<xy;i++){
	  dumd=0;
	  m=0;
	  x0=i%xd;
	  y0=i/xd;
	  for(j=0;j<aindexbig_n;j++){
	    x=x0+aindexbig_x[j];
	    y=y0+aindexbig_y[j];
	    if(x>=0 && x<xd && y>=0 && y<yd){ 
	      dumd+=dm[n][x+y*xd];
	      m++;
	    }
	  }
	  localdbuffer[i]=dumd/(double)m;
	}
	for(i=0;i<xy;i++){
	  dm[n][i]-=localdbuffer[i];
	}
	//	memcpy(dbuffers[n],localdbuffer,xy*sizeof(double));
      }
      free(localdbuffer);
      break;
    case FLOAT:
      if(!(localfbuffer=(float *)calloc(xy,sizeof(float)))){
	printf("INI Cannot allocate for localfbuffer\n");
	return;
      }
      for(n=0;n<nn;n++){
	for(i=0;i<xy;i++){
	  dumf=0;
	  m=0;
	  x0=i%xd;
	  y0=i/xd;
	  for(j=0;j<aindexbig_n;j++){
	    x=x0+aindexbig_x[j];
	    y=y0+aindexbig_y[j];
	    if(x>=0 && x<xd && y>=0 && y<yd){ 
	      dumf+=fm[n][x+y*xd];
	      m++;
	    }
	  }
	  localfbuffer[i]=dumf/(float)m;
	}
	for(i=0;i<xy;i++){
	  fm[n][i]-=localfbuffer[i];
	}
      }
      free(localfbuffer);
      break;
    default:
      printf("AVE Unknown mode %i\n",mode);
      return;
    }
    free(aindexbig_x);
    free(aindexbig_y);
  }
}

void RemoveBias(int xd,int yd,double *db,float *fb,int mode){
  unsigned long xy,i;
  double dbias;
  float fbias;

  xy=(unsigned long)xd*(unsigned long)yd;
  switch(mode){
  case DOUBLE:
    for(dbias=0.0,i=0;i<xy;i++) dbias+=db[i];
    dbias/=(double)xy;
    for(i=0;i<xy;i++) db[i]-=dbias;
    printf("AVE Removed bias: %f\n",dbias);
    break;
  case FLOAT:
    for(fbias=0.0,i=0;i<xy;i++) fbias+=fb[i];
    fbias/=(float)xy;
    for(i=0;i<xy;i++) fb[i]-=fbias;
    printf("AVE Removed bias: %f\n",fbias);
    break;
  default:
    printf("AVE Unknown mode %i\n",mode);
  }
}

float FindAverage(unsigned long n,float *fb,float *faverage){
  unsigned long i;
  float a;
  for(a=0.0,i=0;i<n;i++) a+=fb[i];
  a/=(float)n;
  if(faverage) *faverage=a;
  return(a);
}
