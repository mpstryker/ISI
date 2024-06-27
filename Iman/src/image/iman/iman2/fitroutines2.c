#include "iman1.h"

double PolynomialFitS(int k,int n,double x){
  double pn0_0,pn1_0,a0,a1,p,dk,dumd;
  unsigned long i;
  pn0_0=1.0/sqrt(dk=(double)k);
  if(n==0) return(pn0_0);
  p=2.0*x-(dk-1.0);
  dk*=dk;
  pn1_0=pn0_0*(a0=sqrt(3.0/(dk-1.0)))*p;
  if(n==1) return(pn1_0);
  for(i=2;i<=n;i++){
    dumd=(double)(i*i);
    a1=sqrt((double)(4.0-1.0/dumd)/(dk-dumd));
    dumd=a1*(p*pn1_0-pn0_0/a0);
    pn0_0=pn1_0;
    pn1_0=dumd;
    a0=a1;
  }
  return(pn1_0);
}

// k - number of points
// n - number of polynomials
// x,y - in-out values
// z - polimial coeffs
// mode - 0(DISASSEMBLE), 1(ASSEMBLE)
void PolynomialFitM(int k,int n,double x,double *y,double *z,int mode){
  double pn0_0,pn1_0,p,dumd;
  unsigned long i;
  static double *alpha=(double*)0;
  static int kk=0,nn=0;

  if(kk!=k){
    nn=0;
    kk=k;
  }
  if(nn!=n){
    if(!(alpha=(double *)realloc((void*)alpha,n*sizeof(double)))){
      printf("Cannot allocate for alpha or n=0\n");
      return;
    }
    if(n>nn){
      p=(double)k;
      if(nn==0){ 
	alpha[0]=1.0/sqrt(p);
	nn++;
      }
      p*=p;
      for(i=nn;i<n;i++){ 
	dumd=(double)i;
	dumd*=dumd;
	alpha[i]=sqrt((4.0-1.0/dumd)/(p-dumd));
      }
    }
    nn=n;
  }
  /*
  if(n<1){
    printf("Poly count starts from n=1 - zero order\n");
    return;
  }
  */
  if(mode==DISASSEMBLE){
    *z+=(*y)*(*alpha);
    if(n==1) return;
    z[1]+=(*y)*(pn1_0=(pn0_0=*alpha)*alpha[1]*(p=2.0*x-(double)(k-1)));
    if(n==2) return;
    for(i=2;i<n;i++){
      z[i]+=(*y)*(dumd=alpha[i]*(p*pn1_0-pn0_0/alpha[i-1]));
      pn0_0=pn1_0;
      pn1_0=dumd;
    }
  }
  else{
    *y=(*z)*(*alpha);
    if(n==1) return;
    *y+=z[1]*(pn1_0=(pn0_0=*alpha)*alpha[1]*(p=2.0*x-(double)(k-1)));
    if(n==2) return;
    for(i=2;i<n;i++){
      *y+=z[i]*(dumd=alpha[i]*(p*pn1_0-pn0_0/alpha[i-1]));
      pn0_0=pn1_0;
      pn1_0=dumd;
    }
  }
  return;
}
