#include "mapanm21.h"

int AnalyzePW(int xd,int yd,double *mx,double *my){
  int i,j,k,m;
  int ix,iy,xpos,ypos,xrbits,yrbits,xbits,ybits,xangle,yangle;
  int p1,p2;
  double xmapx1,xmapy1,xmapx2,xmapy2;
  double ymapx1,ymapy1,ymapx2,ymapy2;
  double dumd,vp12,vp34,x,y;
  int c;
  int loop[4],loopx[4],loopy[4];

  loop[0]=0;
  loop[1]=1;
  loop[2]=1+xd;
  loop[3]=xd;

  loopx[0]=0;
  loopx[1]=1;
  loopx[2]=1;
  loopx[3]=0;
  loopy[0]=0;
  loopy[1]=0;
  loopy[2]=1;
  loopy[3]=1;

  free(apinwheels);
  apinwheels_n=0;
  for(j=0;j<yd-1;j++){
    for(i=0;i<xd-1;i++){
      k=i+j*xd;

      for(xbits=ybits=ix=iy=m=0;m<4;m++){
	if(mx[k+loop[m]]>0.0){ 
	  ix++;
	  xbits|=1; 
	}
	xbits<<=1;
	if(my[k+loop[m]]>0.0){
	  iy++;
	  ybits|=1; 
	}
	ybits<<=1;
      }
      xbits>>=1;
      ybits>>=1;      

//      printf("%i %i %i %i (%i %i)\n",i,j,xbits,ybits,ix,iy);
      if((ix%4) && (iy%4)){

	if(!Rotator(xbits,&xrbits,&xpos,&xangle)) continue;
	p1=xpos;
	p2=(xpos+3)%4;
	dumd=mx[k+loop[p1]]/(mx[k+loop[p1]]-mx[k+loop[p2]]);
	xmapx1=loopx[p1]+dumd*(loopx[p2]-loopx[p1]);
	xmapy1=loopy[p1]+dumd*(loopy[p2]-loopy[p1]);
	p1=(xpos+(2-xrbits/3))%4;
	p2=(p1+1)%4;
	dumd=mx[k+loop[p1]]/(mx[k+loop[p1]]-mx[k+loop[p2]]);
	xmapx2=loopx[p1]+dumd*(loopx[p2]-loopx[p1]);
	xmapy2=loopy[p1]+dumd*(loopy[p2]-loopy[p1]);

	if(!Rotator(ybits,&yrbits,&ypos,&yangle)) continue;
	p1=ypos;
	p2=(ypos+3)%4;
	dumd=my[k+loop[p1]]/(my[k+loop[p1]]-my[k+loop[p2]]);
	ymapx1=loopx[p1]+dumd*(loopx[p2]-loopx[p1]);
	ymapy1=loopy[p1]+dumd*(loopy[p2]-loopy[p1]);
	p1=(ypos+(2-yrbits/3))%4;
	p2=(p1+1)%4;
	dumd=my[k+loop[p1]]/(my[k+loop[p1]]-my[k+loop[p2]]);
	ymapx2=loopx[p1]+dumd*(loopx[p2]-loopx[p1]);
	ymapy2=loopy[p1]+dumd*(loopy[p2]-loopy[p1]);
	
//	printf("(%f %f)(%f %f)(%f %f)(%f %f)\n",xmapx1,xmapy1,xmapx2,xmapy2,ymapx1,ymapy1,ymapx2,ymapy2);

	if(VECTPROD(xmapx2-xmapx1,xmapy2-xmapy1,ymapx1-xmapx1,ymapy1-xmapy1)*VECTPROD(xmapx2-xmapx1,xmapy2-xmapy1,ymapx2-xmapx1,ymapy2-xmapy1) <= 0.0){

	  if(!(apinwheels=(Pinwheel *)realloc(apinwheels,(++apinwheels_n)*sizeof(Pinwheel)))){
	    printf("Cannot reallocate for apinwheels\n");
	    return(1);
	  }
	  if((dumd=VECTPROD(xmapx2-xmapx1,xmapy2-xmapy1,ymapx2-ymapx1,ymapy2-ymapy1))>0.0) c=1.0;
	  else c=-1.0;

	  //	  c=Charge(xangle,yangle);
	  vp12=VECTPROD(xmapx1,xmapy1,xmapx2,xmapy2);
	  vp34=VECTPROD(ymapx1,ymapy1,ymapx2,ymapy2);
	  x=((ymapx1-ymapx2)*vp12-(xmapx1-xmapx2)*vp34)/dumd;
	  y=((ymapy1-ymapy2)*vp12-(xmapy1-xmapy2)*vp34)/dumd;
	  apinwheels[apinwheels_n-1].x=(x+(double)i);
	  apinwheels[apinwheels_n-1].y=(y+(double)j);

	  apinwheels[apinwheels_n-1].charge=(double)c;
//	  printf("PW: %2i %3i %3i x: %2i %i %i (%4.3f %4.3f)(%4.3f %4.3f) y: %2i %i %i (%4.3f %4.3f)(%4.3f %4.3f)\n",c,i,j,xbits,xrbits,xpos,xmapx1,xmapy1,xmapx2,xmapy2,ybits,yrbits,ypos,ymapx1,ymapy1,ymapx2,ymapy2);
//	  printf("%f %f %f %f %f\n",x,y,dumd,vp12,vp34);
	}
      }
    }
  }
  printf("NA=%i\n",apinwheels_n);
  for(i=0;i<apinwheels_n;i++){
    printf("%i (%f %f)\n",i,apinwheels[i].x,apinwheels[i].y);
  }

  return(0);
}

int Rotator(int bits,int *rbits,int *position,int *angle){
  switch(bits){
  case 0:
    *rbits=0;
    *position=0;
    *angle=-1;
    return(0);
    break;
  case 1:
    *rbits=1;
    *position=0;
    *angle=3;
    return(1);
    break;
  case 2:
    *rbits=1;
    *position=3;
    *angle=1;
    return(1);
    break;
  case 3:
    *rbits=3;
    *position=0;
    *angle=2;
    return(1);
    break;
  case 4:
    *rbits=1;
    *position=2;
    *angle=7;
    return(1);
    break;
  case 5:
    *rbits=5;
    *position=0;
    *angle=-1;
    printf("Bad bits: 5\n");
    return(0);
    break;
  case 6:
    *rbits=3;
    *position=3;
    *angle=0;
    return(1);
    break;
  case 7:
    *rbits=7;
    *position=0;
    *angle=1;
    return(1);
    break;
  case 8:
    *rbits=1;
    *position=1;
    *angle=5;
    return(1);
    break;
  case 9:
    *rbits=3;
    *position=1;
    *angle=4;
    return(1);
    break;
  case 10:
    *rbits=5;
    *position=1;
    *angle=-1;
    printf("Bad bits: 10\n");
    return(0);
    break;
  case 11:
    *rbits=7;
    *position=1;
    *angle=3;
    return(1);
    break;
  case 12:
    *rbits=3;
    *position=2;
    *angle=6;
    return(1);
    break;
  case 13:
    *rbits=7;
    *position=2;
    *angle=5;
    return(1);
    break;
  case 14:
    *rbits=7;
    *position=3;
    *angle=7;
    return(1);
    break;
  case 15:
    *rbits=0;
    *position=0;
    *angle=-1;
    return(0);
    break;
  default:
    printf("Bad bits: %i\n",bits);
    return(0);
    break;    
  }
}

int Charge(int x,int y){
  int i;
  i=y-x;
  if(i<0) i+=8;
  if(!(i%4)) return(0);
  if(i/4) return(-1);
  else return(1);
}

double FindMapAngle(int iX,int iY,double *pdMX1,double *pdMY1,double *pdMX2,double *pdMY2,int iMode){
  unsigned long i,ulXY,ulArea;
  double dX,dY,dXN,dYN,dDum,dDum1,dDum2,dNorm,dNorm1=0.0,dNorm2=0.0;
  double dX1=0.0,dY1=0.0,dX2=0.0,dY2=0.0,dXX,dYY,dDotProd=0;
  //int iIgnoreZeros=1;
  ulXY=(unsigned long)iX*(unsigned long)iY;
  switch(iMode){
  case MODE_INVERSE_NONE:
    for(i=0,dX=dY=dXN=dYN=0.0,ulArea=0;i<ulXY;i++){
      if((dDum1=pdMX1[i]*pdMX1[i]+pdMY1[i]*pdMY1[i])!=0.0 && 
	 (dDum2=pdMX2[i]*pdMX2[i]+pdMY2[i]*pdMY2[i])!=0.0){
	dNorm1+=dDum1;
	dNorm2+=dDum2;
	dNorm=sqrt(dDum1*dDum2);
	dX+= (dDum= pdMX1[i]*pdMX2[i]+pdMY1[i]*pdMY2[i]);
	dXN+=dDum/dNorm;
	dY+= (dDum= pdMX1[i]*pdMY2[i]-pdMX2[i]*pdMY1[i]);
	dYN+=dDum/dNorm;
	ulArea++;

	dX1+=pdMX1[i];
	dY1+=pdMY1[i];
	dX2+=pdMX2[i];
	dY2+=pdMY2[i];
      }
    }
    if(ulArea){
      dXX=(dX2*(dX*ulArea-dY1*dY2)+dX1*(dY2*dY2-dNorm2*ulArea))/((dX1*dX2+dY1*dY2-dX*ulArea)*ulArea);
      dYY=(dY2*(dX*ulArea-dX1*dX2)+dY1*(dX2*dX2-dNorm2*ulArea))/((dX1*dX2+dY1*dY2-dX*ulArea)*ulArea);
    }
    break;
  case MODE_INVERSE_Y:
    for(i=0,dX=dY=dXN=dYN=0.0,ulArea=0;i<ulXY;i++){
      if((dDum1=pdMX1[i]*pdMX1[i]+pdMY1[i]*pdMY1[i])!=0.0 && 
	 (dDum2=pdMX2[i]*pdMX2[i]+pdMY2[i]*pdMY2[i])!=0.0){
	dNorm1+=dDum1;
	dNorm2+=dDum2;
	dNorm=sqrt(dDum1*dDum2);
	dX+= (dDum=  pdMX1[i]*pdMX2[i]-pdMY1[i]*pdMY2[i]);
	dXN+=dDum/dNorm;
	dY+= (dDum=-(pdMX1[i]*pdMY2[i]+pdMX2[i]*pdMY1[i]));
	dYN+=dDum/dNorm;
	ulArea++;
      }
    }
    dXX=0;
    dYY=0;
    break;
  default:
    return(0.0);
  }

  dDum=atan2(dY,dX)*90.0/M_PI;
  dDotProd=dX/sqrt(dNorm1*dNorm2);
  printf("ANA            Angle between maps %f (%f) Norm %f  (DotPr %f)\n",dDum,dDum+180.0,sqrt((dX*dX+dY*dY)/(dNorm1*dNorm2)),dDotProd);
  dDum=atan2(dYN,dXN)*90.0/M_PI;
  if(ulArea){
    dNorm=hypot(dYN,dXN)/(double)ulArea;
    printf("ANA Normalized Angle between maps %f (%f) Norm %f\n",dDum,dDum+180.0,dNorm);
    printf("ANA Area Used %u  All %u\n",ulArea,ulXY);
    if(iMode==MODE_INVERSE_NONE){
      dDotProd=(dX+dXX*dX1+dYY*dY1)/sqrt(dNorm1*(dNorm2+2*dXX*dX2+2*dYY*dY2+(dXX*dXX+dYY*dYY)*ulArea));
      printf("ANA   Best shift for full DotProd (%f,%f) %f\n",dXX,dYY,dDotProd);

      // Best fit
      dXX=-(ulArea*dY*dY*dX2+ulArea*dY1*dY*dNorm2-ulArea*dX1*dX*dNorm2+ulArea*dX*dX*dX2-2.0*dY2*dX1*dY*dX2-dY2*dY2*dY1*dY+dY1*dX2*dX2*dY+dNorm2*dX1*dX1*dX2+dY2*dY2*dX1*dX-dX2*dX2*dX*dX1+dNorm2*dY1*dY1*dX2-2.0*dY2*dY1*dX*dX2)/(dX1*dX1*dY2*dY2+dY1*dY1*dY2*dY2+dY*dY*ulArea*ulArea+dX*dX*ulArea*ulArea-2.0*dX*dY1*ulArea*dY2-2.0*dY*dX1*ulArea*dY2+2.0*dX2*dY*dY1*ulArea-2.0*dX2*dX*dX1*ulArea+dX2*dX2*dY1*dY1+dX2*dX2*dX1*dX1);

      dYY = (dX1*dY*dY2*dY2-dX2*dX2*dY1*dX+dNorm2*dX1*dY*ulArea+dY1*dX*dY2*dY2+2.0*dX2*dX*dX1*dY2-dNorm2*dY1*dY1*dY2-dY*dY*ulArea*dY2-dNorm2*dX1*dX1*dY2+dX*dY1*ulArea*dNorm2-dX*dX*ulArea*dY2-2.0*dX2*dY*dY1*dY2-dX2*dX2*dX1*dY)/(dX1*dX1*dY2*dY2+dY1*dY1*dY2*dY2+dY*dY*ulArea*ulArea+dX*dX*ulArea*ulArea-2.0*dX*dY1*ulArea*dY2-2.0*dY*dX1*ulArea*dY2+2.0*dX2*dY*dY1*ulArea-2.0*dX2*dX*dX1*ulArea+dX2*dX2*dY1*dY1+dX2*dX2*dX1*dX1);

      dDotProd = sqrt( ((dX+dXX*dX1+dYY*dY1)*(dX+dXX*dX1+dYY*dY1)+(dY+dYY*dX1-dXX*dY1)*(dY+dYY*dX1-dXX*dY1))/(dNorm1*(dNorm2+2.0*dXX*dX2+2.0*dYY*dY2+(dXX*dXX+dYY*dYY)*ulArea)) );

      printf("ANA   Best shift  (%f,%f) %f\n",dXX,dYY,dDotProd);
    }
  }
  else{
    printf("ANA Blank image\n");
  }
  return(atan2(dY,dX)*90.0/M_PI);
}

#define N_FOURIER_ITERATIONS 30

// Calculates maximum(ampl. and phase) of a set of Fourier components
int FourierMax(int iNC,double *pdC,int iMode){
  int i;
  int iN;
  double dMin,dMax,dMinX,dMaxX;
  double dIncr,dDum,dX;
  double dDumL,dDumR,dXL,dXR;
  iN=4*iNC;
  dMin=1.0e10;
  dMax=-dMin;
  dMinX=dMaxX=0.0;
  dIncr=2.0*M_PI/iN;
  for(i=0,dX=0.0;i<iN;i++,dX+=dIncr){
    dDum=FourierSum(dX,iNC,pdC,iMode);
    if(dMin>dDum){
      dMin=dDum;
      dMinX=dX;
    }
    if(dMax<dDum){
      dMax=dDum;
      dMaxX=dX;
    }
  }
  dX=dMax>fabs(dMin) ? dMaxX : dMinX;
  dXL=dX-dIncr;
  dXR=dX+dIncr;
  dDumL=FourierSumDir(dXL,iNC,pdC,iMode);
  dDumR=FourierSumDir(dXR,iNC,pdC,iMode);
  for(i=0;i<N_FOURIER_ITERATIONS;i++){
    dDum=FourierSumDir(dX,iNC,pdC,iMode);
    if(dDum==0.0) break;
    if(dDum*dDumL<0.0){
      dXR=dX;
      dDumR=dDum;
      dX=0.5*(dXL+dXR);
    }
    else{
      dXL=dX;
      dDumL=dDum;
      dX=0.5*(dXL+dXR);
    }
  }
  if(iMode){
    pdC[2*iNC]=FourierSum(dX,iNC,pdC,iMode);
    pdC[2*iNC+1]=dX;
  }
  else{
    dDum=FourierSum(dX,iNC,pdC,iMode);
    pdC[2*iNC]=dDum*cos(dX);
    pdC[2*iNC+1]=dDum*sin(dX);
  }
  return(0);
}
// iMode=0: Cartesian, =1: polar
double FourierSum(double x,int iNC,double *pdC,int iMode){
  int i;
  double d=0.0;
  if(iMode)
    for(i=0;i<iNC;i++) d+=pdC[2*i]*cos((i+1)*x-pdC[2*i+1]);
  else
    for(i=0;i<iNC;i++) d+=pdC[2*i]*cos((i+1)*x)+pdC[2*i+1]*sin((i+1)*x);
  return(d);
}

// iMode=0: Cartesian, =1: polar
double FourierSumDir(double x,int iNC,double *pdC,int iMode){
  int i;
  double d=0.0;
  if(iMode)
    for(i=0;i<iNC;i++) d-=(i+1)*(pdC[2*i])*sin((i+1)*x-pdC[2*i+1]);
  else
    for(i=0;i<iNC;i++) 
      d-=(i+1)*(pdC[2*i]*sin((i+1)*x)-pdC[2*i+1]*cos((i+1)*x));
  return(d);
}

