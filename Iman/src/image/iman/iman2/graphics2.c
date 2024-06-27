#include "iman1.h"

void Pallet(int type,float contra,float bright){
  int r=9,h=5,hex=7;
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
  }
}


void start_graphics(){
  mainid=cpgopen(device);
  if(n_panels==4){ 
    cpgsubp(2,2);
    cpgpap(12.0,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
  }
  else{
    cpgpap(4.0,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
  }
  cpgsch(0.1);
  cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
  if(n_panels==4){
    cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
    cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
    cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
  }
  //  cpglab("","","12345");
  //  cpgpanl(1,1); 
  //  cpgbeg(0,"?",1,1);
  //  cpgenv(pminx,pmaxx,pminy,pmaxy,1,1);
  //  cpglab("","","");
}

void end_graphics(){
  cpgend();
}
