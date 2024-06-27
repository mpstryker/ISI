#include "mapanm21.h"

void AverageMaps(){
  int i,j,k,m,x,y,x0,y0;
  double dumd;

  if(radius){
    printf("AVE R=%i\n",radius);
    for(k=0;k<MAX_FILE_N;k++){
      for(i=0;i<XYdim;i++){
	dumd=0.0;
	m=0;
	x0=i%Xdim;
	y0=i/Xdim;
	for(j=0;j<aindex_n;j++){
	  x=x0+aindex_x[j];
	  y=y0+aindex_y[j];
	  if(x>=0 && x<Xdim && y>=0 && y<Ydim){ 
	    dumd+=maps[k][x+y*Xdim];
	    m++;
	  }
	}
	dbufferx[i]=dumd/(double)m;
      }
      memcpy(maps[k],dbufferx,XYdim*sizeof(double));    
    }
  }

  if(radiusbig){
    printf("AVE RB=%i\n",radiusbig);
    for(k=0;k<MAX_FILE_N;k++){
      for(i=0;i<XYdim;i++){
	dumd=0;
	m=0;
	x0=i%Xdim;
	y0=i/Xdim;
	for(j=0;j<aindexbig_n;j++){
	  x=x0+aindexbig_x[j];
	  y=y0+aindexbig_y[j];
	  if(x>=0 && x<Xdim && y>=0 && y<Ydim){ 
	    dumd+=maps[k][x+y*Xdim];
	    m++;
	  }
	}
	dbufferx[i]=dumd/(double)m;
      }
      for(i=0;i<XYdim;i++){
	maps[k][i]-=dbufferx[i];
      }
    }
  }
  if(radiusbig || radius) printf("AVE Maps averaged\n");
}
