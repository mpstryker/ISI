#include "mapanm21.h"

void start_graphics(int paneln){
  int i;
  mainid=cpgopen(device);
  cpgsubp(paneln,2);
  cpgpap(main_win_size,(2.0/(float)paneln)*(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
  cpgsch(0.7);
  //  cpgenv(pminx,pmaxx,pminy,pmaxy,1,1);
  for(i=0;i<2*paneln;i++) cpgenv(pminx+1,pmaxx-1,pminy+1,pmaxy-1,1,-2);
  cpgsch(1.0);
  //  cpglab("","","12345");
  //  cpgpanl(1,1); 
  //  cpgbeg(0,"?",1,1);
  //  cpgenv(pminx,pmaxx,pminy,pmaxy,1,1);
  //  cpglab("","","");
}

void end_graphics(){
  cpgend();
}

void Pallet(int type,float contra,float bright){
  int r=9,h=5,hex=7,dozen=13;
  float rl[]={ 0.0, 0.125,0.25, 0.375,0.5,  0.625,0.75,0.875,1.0};
  float rr[]={ 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
  float rg[]={ 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
  float rb[]={ 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};
  float hl[]={ 0.0, 0.2, 0.4, 0.6, 1.0};
  float hr[]={ 0.0, 0.5, 1.0, 1.0, 1.0};
  float hg[]={ 0.0, 0.0, 0.5, 1.0, 1.0};
  float hb[]={ 0.0, 0.0, 0.0, 0.3, 1.0};
  float hexl[]={ 0.0, 1.0/6.0,      1.0/3.0,  0.5,2.0/3.0, 5.0/6.0,     1.0};
  float hexr[]={ 1.0, 1.0,          0.5,      0.0,0.0,     0.5,         1.0};
  float hexg[]={ 0.0, 0.5,          1.0,      1.0,0.5,     0.0,         0.0};
  float hexb[]={ 0.5, 0.0,          0.0,      0.5,1.0,     1.0,         0.5};
  /*
  float dozenl[]={0.0, 1.0/12.0, 1.0/6.0, 0.25,1.0/3.0, 5.0/12.0, 0.5, 7.0/12.0, 2.0/3.0, 0.75,5.0/6.0, 11.0/12.0, 1.0};
  float dozenr[]={0.0, 0.0,     0.5,      1.0, 1.0,      0.5,     0.0, 0.0,     0.5,       1.0,1.0,     1.0,       0.5};
  float dozeng[]={1.0, 0.5,     0.0,      0.0, 0.5,      1.0,     1.0, 0.5,     0.0,       0.0,0.0,     0.5,       1.0};
  float dozenb[]={0.5, 1.0,     1.0,      0.5, 0.0,      0.0,     0.5, 1.0,     1.0,       0.5,0.5,     0.0,       0.0};
  */
  float dozenl[]={0.0, 1.0/12.0, 1.0/6.0, 0.25,1.0/3.0, 5.0/12.0, 0.5, 7.0/12.0, 2.0/3.0, 0.75,5.0/6.0, 11.0/12.0, 1.0};
  float dozenr[]={0.0, 0.0,     0.5,      1.0, 1.0,      0.5,     0.0, 0.0,     0.5,       1.0,1.0,     0.5,       0.0};
  float dozeng[]={1.0, 0.5,     0.0,      0.0, 0.5,      1.0,     1.0, 0.5,     0.0,       0.0,0.5,     1.0,       1.0};
  float dozenb[]={0.5, 1.0,     1.0,      0.5, 0.0,      0.0,     0.5, 1.0,     1.0,       0.5,0.0,     0.0,       0.5};
  switch(type){
  case 0:
    cpgctab(rl, rr, rg, rb, r, contra, bright);
    break;
  case 1:
    cpgctab(hl, hr, hg, hb, h, contra, bright);
    break;
  case 2:
    cpgctab(hexl, hexr, hexg, hexb, hex, contra, bright);
    break;
  case 3:
    cpgctab(dozenl, dozenr, dozeng, dozenb, dozen, contra, bright);
    break;
  }
}
