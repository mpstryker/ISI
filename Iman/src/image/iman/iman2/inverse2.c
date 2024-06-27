#include "iman1.h"

//#define PLOT_SOL

#define TINY 1.0e-20
#define LARGE 1.0e20
#define EPSILON 1.0e-18
#define N 2
#define ITERATIONM 0.5

void printv(double *a,int n);
double dotprod(double *a,double *b,int n);
double PAP(double *p,int n,double *Ap,double *up,double *vc,double *vs);
int BasisFit(int n,double *co,double *si);

int Inverse(int n,double omega,double alpha,int h){
  int i,j,it=4,k=2*10+1,l,print_p=10,m,o;
  double *vc,*vs;
  double *up,*Ap;
  double *sol[N],*r,*p,a,b,rr;
  double *sol0,*sol1;
  double dumd,co,si,coN,siN,dc,ds,coa,sia;
  int imax[N],imin[N],imax2[N];
  double miny[N],maxy[N],maxy2[N];
#ifdef PLOT_SOL
  int dumid;
  double gminy,gmaxy;
  char str[128];
#endif

  
  if(!(synch_cos=(double *)calloc(n,sizeof(double)))){
    printf("INV Cannot allocate for synch_cos (%luB)\n",(unsigned long)n*sizeof(double));
    return(2);
  }
  if(!(synch_sin=(double *)calloc(n,sizeof(double)))){
    printf("INV Cannot allocate for synch_sin (%luB)\n",(unsigned long)n*sizeof(double));
    return(3);
  }

  if(!do_precise){
    synch_cos[0]=1.0;
    synch_sin[0]=0.0;
    synch_cos[1]=cos(omega*(double)h);
    synch_sin[1]=sin(omega*(double)h);
    for(i=2;i<n;i++){
      synch_cos[i]=synch_cos[i-1]*synch_cos[1]-synch_sin[i-1]*synch_sin[1];
      synch_sin[i]=synch_sin[i-1]*synch_cos[1]+synch_cos[i-1]*synch_sin[1];
    }
    dumd=synch_cos[0]=2.0/(double)n;
    for(i=1;i<n;i++){
      synch_cos[i]*=dumd;
      synch_sin[i]*=dumd;
    }
    BasisFit(n,synch_cos,synch_sin);
    return(0);
  }

  k=2*h-1;
  if(k<1 || k>=n){
    printf("INV Bad harmonic (%i)\n",h);
    free(synch_cos);
    free(synch_sin);
    return(1);
  }

  it=(int)(ITERATIONM*(float)n);
  o=n;

  if(n>10000) print_p=2;

  vc=(double *)calloc(n,sizeof(double));
  vs=(double *)calloc(n,sizeof(double));
  for(i=0;i<N;i++) sol[i]=(double *)calloc(n,sizeof(double));
  up=(double *)calloc(n,sizeof(double));
  r=(double *)calloc(n,sizeof(double));
  p=(double *)calloc(n,sizeof(double));
  Ap=(double *)calloc(n,sizeof(double));

  if(!(vc && vs && up && r && p && Ap && sol[0] && sol[1])){
    printf("INV Cannot allocate\n");
    return(2);
  }

  vc[0]=1.0;
  vs[0]=0.0;
  vc[1]=cos(omega);
  vs[1]=sin(omega);
  for(i=2;i<n;i++){
    vc[i]=vc[i-1]*vc[1]-vs[i-1]*vs[1];
    vs[i]=vs[i-1]*vc[1]+vc[i-1]*vs[1];
  }

  for(m=0;m<N;m++){
    sol0=sol[m];
    for(i=0;i<n;i++) sol0[i]=0.0;
    sol0[k]=2.0/(double)n;
    PAP(sol[m],n,Ap,up,vc,vs);
    for(i=0;i<n;i++) r[i]=p[i]=-Ap[i];
    r[k]+=1.0;
    p[k]+=1.0;

    l=0;
    rr=dotprod(r,r,n);
    while(1){
      a=rr/PAP(p,n,Ap,up,vc,vs);
      for(j=0;j<n;j++) r[j]-=a*Ap[j];
      b=(dumd=dotprod(r,r,n))/rr;
      rr=dumd;
      for(j=0;j<n;j++){ 
	sol0[j]+=a*p[j];
	p[j]=r[j]+b*p[j];
      }
      if(l>it-10) printf("INV %i %f\n",l,rr);
      if(rr<EPSILON || l==it){
	printf("INV It=%i, rr=%e\n",l,rr);
	break;
      }
      l++;
      if(!(l%print_p)) printf("INV %i %e\n",l,rr);
    }
    //  printv(sol[m],n);
    printf("INV First %e\n",sol0[0]);
    
    if(fabs(sol0[0])>fabs(sol0[1])){
      maxy[m]=fabs(sol0[0]);
      imax[m]=0;
      miny[m]=maxy2[m]=fabs(sol0[1]);
      imax2[m]=imin[m]=1;
    }
    else{
      miny[m]=maxy2[m]=fabs(sol0[0]);
      imax2[m]=imin[m]=0;
      maxy[m]=fabs(sol0[1]);    
      imax[m]=1;
    }
    for(i=2;i<n;i++){
      dumd=fabs(sol0[i]);
      if(miny[m]>dumd){ 
	miny[m]=dumd;
	imin[m]=i;
      }
      if(maxy2[m]<dumd){ 
	maxy2[m]=dumd;
	imax2[m]=i;
	if(maxy[m]<maxy2[m]){
	  dumd=maxy[m];
	  maxy[m]=maxy2[m];
	  maxy2[m]=dumd;
	  l=imax[m];
	  imax[m]=imax2[m];
	  imax2[m]=l;
	}
      }
    }

    printf("INV %i MAX %e(%e)(%i) NMAX %e(%i) MIN %e(%i)\n",k,maxy[m],sol[m][imax[m]],imax[m],maxy2[m],imax2[m],miny[m],imin[m]);
#ifdef PLOT_SOL
    if(m==0){ 
      gminy=miny[0];
      gmaxy=maxy[0];
      dumid=cpgopen("/xw");
      cpgpap(14.0,0.5);
      cpgask(0);
    }
    if(gminy>miny[m]) gminy=miny[m];
    if(gmaxy<maxy[m]) gmaxy=maxy[m];
    cpgsci(1);
    cpgenv(0.0,(float)(o),log10(gminy),log10(gmaxy),0,1);
    for(l=0;l<=m;l++){
      cpgsci(2+l);
      cpgmove(0.0,(float)log10(fabs(sol[l][0])));
      for(i=1;i<o;i++){
	cpgdraw((float)i,(float)log10(fabs(sol[l][i])));
      }
    }
#endif
    k++;
  }

#ifdef PLOT_SOL
  sprintf(str,"N=%i %i",n,k);
  cpgsci(1);
  cpglab("","",str);
  cpgask(1);
  cpgclos();
#endif

  free(r);
  free(p);
  free(Ap);
  free(up);


  sol0=sol[0];
  sol1=sol[1];

  coa=cos(alpha);
  sia=sin(alpha);

  for(i=0,dc=sol0[0],ds=sol1[0];i<n/2;i++){ 
    dc+=sol0[2*i+1];
    ds+=sol1[2*i+1];
  }
  synch_cos[0]=dc*coa-ds*sia;
  synch_sin[0]=ds*coa+dc*sia;
  for(j=1;j<n;j++){
    co=vc[j];
    si=vs[j];
    coN=1.0;
    siN=0.0;
    for(i=0,dc=sol0[0],ds=sol1[0];i<(n-1)/2;i++){
      dumd=coN*co-siN*si;
      siN=siN*co+coN*si;
      coN=dumd;
      dc+=sol0[2*i+1]*coN+sol0[2*(i+1)]*siN;
      ds+=sol1[2*i+1]*coN+sol1[2*(i+1)]*siN;
    }
    if(!(n%2)){
      dc+=sol0[n-1]*(coN*co-siN*si);
      ds+=sol1[n-1]*(coN*co-siN*si);
    } 
    synch_cos[j]=dc*coa-ds*sia;
    synch_sin[j]=ds*coa+dc*sia;
  }
  BasisFit(n,synch_cos,synch_sin);

  free(vc);
  free(vs);
  for(i=0;i<N;i++) free(sol[i]);

  return(0);
}

double dotprod(double *a,double *b,int n){
  int i;
  double dumd;
  for(dumd=0.0,i=0;i<n;i++) dumd+=a[i]*b[i];
  return(dumd);
}

double PAP(double *p,int n,double *Ap,double *up,double *vc,double *vs){
  double pap,dumd,co,si,coN,siN,dumd1,dc,ds;
  int i,j;

  for(i=0,dumd=p[0];i<n/2;i++) dumd+=p[2*i+1];
  up[0]=dumd;
  for(j=1;j<n;j++){
    co=vc[j];
    si=vs[j];
    coN=1.0;
    siN=0.0;
    for(i=0,dumd=p[0];i<(n-1)/2;i++){
      dumd1=coN*co-siN*si;
      siN=siN*co+coN*si;
      coN=dumd1;
      dumd+=p[2*i+1]*coN+p[2*(i+1)]*siN;
    }
    if(!(n%2)){
      dumd+=p[n-1]*(coN*co-siN*si);
    } 
    up[j]=dumd;
  }
  pap=dotprod(up,up,n);
  
  for(i=1,dumd=up[0];i<n;i++) dumd+=up[i];
  Ap[0]=dumd;
  Ap[1]=dotprod(up,vc,n);
  Ap[2]=dotprod(up,vs,n);

  for(j=1;j<(n-1)/2;j++){
    dc=up[0];
    ds=0.0;
    co=vc[j+1];
    si=vs[j+1];
    coN=1.0;
    siN=0.0;
    for(i=1;i<n;i++){
      dumd1=coN*co-siN*si;
      siN=siN*co+coN*si;
      coN=dumd1;
      dc+=up[i]*coN;
      ds+=up[i]*siN;
    }
    Ap[2*j+1]=dc;
    Ap[2*j+2]=ds;
  }
  if(!(n%2)){
    dc=up[0];
    co=vc[n/2];
    si=vs[n/2];
    coN=1.0;
    siN=0.0;
    for(i=1;i<n;i++){
      dumd1=coN*co-siN*si;
      siN=siN*co+coN*si;
      coN=dumd1;
      dc+=up[i]*coN;
    }
    Ap[n-1]=dc;
  } 
  
  return(pap);
}

void printv(double *a,int n){
  int i;
  for(i=0;i<n;i++) printf(" INV %i % 8.6f\n",i,a[i]);
  printf("\n");
}

int BasisFit(int n,double *co,double *si){
  int i;
  for(a_cos=a_sin=b_cos=b_sin=0.0,i=0;i<n;i++){
    a_cos+=co[i];
    a_sin+=si[i];
    b_cos+=co[i]*(double)i;
    b_sin+=si[i]*(double)i;
  }
  return(0);
}

#undef PLOT_SOL

