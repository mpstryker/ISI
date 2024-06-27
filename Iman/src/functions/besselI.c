#include <math.h>

/*
#define A0 1.00000000000000000000000000
#define A1 3.51562500000000000000000000
#define A2 3.08990478515625000000000000
#define A3 1.20699405670166015625000000
#define A4 0.26520865503698587417602539
#define A5 0.37294967114576138556003571e-1
#define A6 0.36420866322828260308597237e-2
#define A7 0.26131042482896551560696359e-3
#define A8 0.14354210348270810793839553e-4
#define A9 0.62301260192147616292706391e-6

#define B0 0.39894228040143267793994605
#define B1 0.13298076013381089264664869e-1
#define B2 0.19947114020071633896997303e-2
#define B3 0.55408650055754538602770284e-3
#define B4 0.22625198772766436596131199e-3
#define B5 0.12217607337293875761910848e-3
#define B6 0.82129471545142164843956256e-4
#define B7 0.66094669957757265993469557e-4
#define B8 0.61963753085397436868877710e-4
#define B9 0.66324165339555034278169104e-4
*/

#define SI00  1.0000000
#define SI01  3.5156229
#define SI02  3.0899424
#define SI03  1.2067492
#define SI04  0.2659732
#define SI05  0.0360768
#define SI06  0.0045813

#define BI00  0.39894228
#define BI01  0.01328592
#define BI02  0.00225319
#define BI03 -0.00157565 
#define BI04  0.00916281
#define BI05 -0.02057706
#define BI06  0.02635537
#define BI07 -0.01647633
#define BI08  0.00392377

#define SI10  0.50000000
#define SI11  0.87890594
#define SI12  0.51498869
#define SI13  0.15084934
#define SI14  0.02658733
#define SI15  0.00301532
#define SI16  0.00032411

#define BI10  0.39894228
#define BI11 -0.03988024
#define BI12 -0.00362018
#define BI13  0.00163801
#define BI14 -0.01031555
#define BI15  0.02282967
#define BI16 -0.02895312
#define BI17  0.01787654
#define BI18 -0.00420059
 
double besselI0(double x){
  double ans,y,ax;

  if((ax=fabs(x))<3.75){
    y=x/3.75;
    y*=y;
    //    ans=A0+y*(A1+y*(A2+y*(A3+y*(A4+y*(A5+y*(A6+y*(A7+y*(A8+y*A9))))))));
    ans=SI00+y*(SI01+y*(SI02+y*(SI03+y*(SI04+y*(SI05+y*SI06)))));
  }
  else{
    y=3.75/ax;
    //    ans=(exp(ax)/sqrt(ax))*(B0+y*(B1+y*(B2+y*(B3+y*(B4+y*(B5+y*(B6+y*(B7+y*(B8+y*B9)))))))));
    ans=(exp(ax)/sqrt(ax))*(BI00+y*(BI01+y*(BI02+y*(BI03+y*(BI04+y*(BI05+y*(BI06+y*(BI07+y*BI08))))))));
  }
  return(ans);
}

double besselI1(double x){
  double ans,y,ax;

  if((ax=fabs(x))<3.75){
    y=x/3.75;
    y*=y;
    ans=ax*(SI10+y*(SI11+y*(SI12+y*(SI13+y*(SI14+y*(SI15+y*SI16))))));
  }
  else{
    y=3.75/ax;
    ans=(exp(ax)/sqrt(ax))*(BI10+y*(BI11+y*(BI12+y*(BI13+y*(BI14+y*(BI15+y*(BI16+y*(BI17+y*BI18))))))));
  }
  return(x < 0.0 ? -ans : ans);
}

double besselI0e(double x){
  double ans,y,ax;

  if((ax=fabs(x))<3.75){
    y=x/3.75;
    y*=y;
    ans=exp(-ax)*(SI00+y*(SI01+y*(SI02+y*(SI03+y*(SI04+y*(SI05+y*SI06))))));
  }
  else{
    y=3.75/ax;
    ans=(BI00+y*(BI01+y*(BI02+y*(BI03+y*(BI04+y*(BI05+y*(BI06+y*(BI07+y*BI08))))))))/sqrt(ax);
  }
  return(ans);
}

double besselI1e(double x){
  double ans,y,ax;

  if((ax=fabs(x))<3.75){
    y=x/3.75;
    y*=y;
    ans=exp(-ax)*ax*(SI10+y*(SI11+y*(SI12+y*(SI13+y*(SI14+y*(SI15+y*SI16))))));
  }
  else{
    y=3.75/ax;
    ans=(BI10+y*(BI11+y*(BI12+y*(BI13+y*(BI14+y*(BI15+y*(BI16+y*(BI17+y*BI18))))))))/sqrt(ax);
  }
  return(x < 0.0 ? -ans : ans);
}

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

double besselIn(int n,double x){
  int j;
  double bi,bim,bip,tox,ans;

  if(n==0) return(besselI0(x));
  n=abs(n);
  if(n==1) return(besselI1(x));
  if(x==0.0) return(0.0);
  else{
    tox=2.0/fabs(x);
    bip=ans=0.0;
    bi=1.0;
    for(j=2*(n+(int)sqrt(ACC*(double)n));j>0;j--){
      bim=bip+j*tox*bi;
      bip=bi;
      bi=bim;
      if(fabs(bi)>BIGNO){
	ans*=BIGNI;
	bi*=BIGNI;
	bip*=BIGNI;
      }
      if(j==n) ans=bip;
    }
    ans*=besselI0(x)/bi;
    return(x<0.0 && (n&1) ? -ans : ans);
  }
}

double besselIne(int n,double x){
  int j;
  double bi,bim,bip,tox,ans;

  if(n==0) return(besselI0e(x));
  n=abs(n);
  if(n==1) return(besselI1e(x));
  if(x==0.0) return(0.0);
  else{
    tox=2.0/fabs(x);
    bip=ans=0.0;
    bi=1.0;
    for(j=2*(n+(int)sqrt(ACC*n));j>0;j--){
      bim=bip+j*tox*bi;
      bip=bi;
      bi=bim;
      if(fabs(bi)>BIGNO){
	ans*=BIGNI;
	bi*=BIGNI;
	bip*=BIGNI;
      }
      if(j==n) ans=bip;
    }
    ans*=besselI0e(x)/bi;
    return(x<0.0 && (n&1) ? -ans : ans);
  }
}
/*
double BesselIe(int n,double x,int k){
  int i,j=10;
  double a,x1,ans,st;

  x1=x*8.0;
  a=(2.0*n)*(2.0*n);
  if(k>1) j=k;
  for(ans=st=1.0,i=0;i<=j;i++) ans+=(st*=-(a-(2.0*i+1.0)*(2.0*i+1.0))/(x1*(i+1.0)));

  return(ans/sqrt(2.0*M_PI*x));
}
*/
