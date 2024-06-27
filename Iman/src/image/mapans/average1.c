/* Map analysis                 */
/* Written by V.Kalatsky        */
/* 06Feb2001                    */
/* Modified 10May2001           */

#include "mapans1.h"

void AverageMaps(int xd,int yd,double **dm,float **fm,int nn,double dr,double drb,int mode){
  int i,j,k,m,n,x,y,x0,y0;
  int r,rb;
  unsigned long xy;
  double dumd;
  double *localdbuffer;
  double *pd;
  float dumf;
  float *localfbuffer;
  float *pf;
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
	pd=dm[n];
	for(i=0;i<xy;i++){
	  dumd=0.0;
	  m=0;
	  x0=i%xd;
	  y0=i/xd;
	  for(j=0;j<aindex_n;j++){
	    x=x0+aindex_x[j];
	    y=y0+aindex_y[j];
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
	  for(j=0;j<aindex_n;j++){
	    x=x0+aindex_x[j];
	    y=y0+aindex_y[j];
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
    free(aindex_x);
    free(aindex_y);
    //    printf("AVE Maps low-pass averaged\n");
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
	pd=dm[n];
	for(i=0;i<xy;i++){
	  dumd=0;
	  m=0;
	  x0=i%xd;
	  y0=i/xd;
	  for(j=0;j<aindexbig_n;j++){
	    x=x0+aindexbig_x[j];
	    y=y0+aindexbig_y[j];
	    if(x>=0 && x<xd && y>=0 && y<yd){ 
	      dumd+=pd[x+y*xd];
	      m++;
	    }
	  }
	  localdbuffer[i]=dumd/(double)m;
	}
	for(i=0;i<xy;i++){
	  pd[i]-=localdbuffer[i];
	}
	memcpy(dbuffers[n],localdbuffer,xy*sizeof(double));
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
	  dumf=0;
	  m=0;
	  x0=i%xd;
	  y0=i/xd;
	  for(j=0;j<aindexbig_n;j++){
	    x=x0+aindexbig_x[j];
	    y=y0+aindexbig_y[j];
	    if(x>=0 && x<xd && y>=0 && y<yd){ 
	      dumf+=pf[x+y*xd];
	      m++;
	    }
	  }
	  localfbuffer[i]=dumf/(float)m;
	}
	for(i=0;i<xy;i++){
	  pf[i]-=localfbuffer[i];
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
    //    printf("AVE Maps high-pass averaged\n");
  }
}

void RemoveBias(int xd,int yd,double *mx,double *my){
  int i;
  unsigned long xy;
  double xbias,ybias;
  xy=(unsigned long)xd*(unsigned long)yd;
  xbias=ybias=0.0;
  for(i=0;i<xy;i++){
    xbias+=mx[i];
    ybias+=my[i];
  }
  xbias/=(double)xy;
  ybias/=(double)xy;
  for(i=0;i<xy;i++){
    mx[i]-=xbias;
    my[i]-=ybias;
  }
  printf("AVE Removed bias: X=%f Y=%f (R=%f P=%f)\n",xbias,ybias,hypot(xbias,ybias),180.0/M_PI*atan2(ybias,xbias));
}
