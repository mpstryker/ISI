/* Bin analysis                 */
/* Written by V.Kalatsky        */
/* 20Jul2001                    */

#include "binan1.h"

void start_graphics(int mode,char *str,int *id){
  int i;

  *id=cpgopen(str);
  cpgsubp(panelnx,panelny);
  cpgpap(window_size,(float)fabs(((double)panelny/(double)panelnx)*(double)((pmaxy-pminy)/(pmaxx-pminx))));
  cpgsch(0.5);
  for(i=0;i<panelnx*panelny;i++) cpgenv(pminx,pmaxx,pminy,pmaxy,1,0);
  cpgsch(3.5);
}

void end_graphics(int id){
  cpgslct(id);
  cpgclos();
}

#define PALLET_SHIFT_SHIFT 1000                                                            

void Pallet(int type,float contra,float bright){
  int r=9,h=5,hex=7,dozen=13,eighteen=19,twentyfour=25;

  float rl[]={ 0.0, 0.125,0.25, 0.375,0.5,  0.625,0.75,0.875,1.0};
  float rr[]={ 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
  float rg[]={ 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
  float rb[]={ 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};

  float hl[]={ 0.0, 0.2, 0.4, 0.6, 1.0};
  float hr[]={ 0.0, 0.5, 1.0, 1.0, 1.0};
  float hg[]={ 0.0, 0.0, 0.5, 1.0, 1.0};
  float hb[]={ 0.0, 0.0, 0.0, 0.3, 1.0};
  
  // RGB Green in the middle perfect hexagon section of the cube
  float hexl1[]={ 0.0, 1.0/6.0, 1.0/3.0,  0.5,2.0/3.0, 5.0/6.0, 1.0};
  float hexr1[]={ 1.0, 1.0,     0.5,      0.0,0.0,     0.5,     1.0};
  float hexg1[]={ 0.0, 0.5,     1.0,      1.0,0.5,     0.0,     0.0};
  float hexb1[]={ 0.5, 0.0,     0.0,      0.5,1.0,     1.0,     0.5};
  
    // RGB Blue in the middle perfect hexagon section of the cube
  float hexl2[]={ 0.0, 1.0/6.0, 1.0/3.0, 0.5, 2.0/3.0, 5.0/6.0, 1.0};
  float hexr2[]={ 1.0, 0.5,     0.0,     0.0, 0.5,     1.0,     1.0};
  float hexg2[]={ 0.5, 1.0,     1.0,     0.5, 0.0,     0.0,     0.5};
  float hexb2[]={ 0.0, 0.0,     0.5,     1.0, 1.0,     0.5,     0.0};
 
  // HSL Blue in the middle S=L=1.0
  float hexl3[]={ 0.0, 1.0/6.0, 1.0/3.0, 0.5, 2.0/3.0, 5.0/6.0, 1.0};
  float hexr3[]={ 1.0, 0.0,     0.0,     0.0, 1.0,     1.0,     1.0};
  float hexg3[]={ 1.0, 1.0,     1.0,     0.0, 0.0,     0.0,     1.0};
  float hexb3[]={ 0.0, 0.0,     1.0,     1.0, 1.0,     0.0,     0.0};

  // HSL Blue @ the left edge S=L=1.0
  float hexr3_1[]={ 0.0, 1.0,     1.0,     1.0, 0.0,     0.0,     0.0};
  float hexg3_1[]={ 0.0, 0.0,     0.0,     1.0, 1.0,     1.0,     0.0};
  float hexb3_1[]={ 1.0, 1.0,     0.0,     0.0, 0.0,     1.0,     1.0};

  float dozenl[]={0.0, 1.0/12.0,1.0/6.0, 0.25,1.0/3.0, 5.0/12.0, 0.5, 7.0/12.0, 2.0/3.0, 0.75,5.0/6.0, 11.0/12.0, 1.0};
  float dozenr[]={0.0, 0.0,     0.5,     1.0, 1.0,     0.5,      0.0, 0.0,      0.5,     1.0,1.0,     0.5,       0.0};
  float dozeng[]={1.0, 0.5,     0.0,     0.0, 0.5,     1.0,      1.0, 0.5,      0.0,     0.0,0.5,     1.0,       1.0};
  float dozenb[]={0.5, 1.0,     1.0,     0.5, 0.0,     0.0,      0.5, 1.0,      1.0,     0.5,0.0,     0.0,       0.5};

  float eighteenl[]={ 0.0, 1.0/18.0,      1.0/9.0,  1.0/6.0,2.0/9.0, 5.0/18.0,     1.0/3.0, 7.0/18.0,      4.0/9.0,  0.5,5.0/9.0, 11.0/18.0,     2.0/3.0, 13.0/18.0,      7.0/9.0,  5.0/6.0,8.0/9.0, 17.0/18.0,     1.0};
  float eighteenr[]={ 1.0, 1.0,          0.5,      0.0,0.0,     0.5,         1.0, 1.0,          0.5,      0.0,0.0,     0.5,         1.0, 1.0,          0.5,      0.0,0.0,     0.5,         1.0};
  float eighteeng[]={ 0.0, 0.5,          1.0,      1.0,0.5,     0.0,         0.0, 0.5,          1.0,      1.0,0.5,     0.0,         0.0, 0.5,          1.0,      1.0,0.5,     0.0,         0.0};
  float eighteenb[]={ 0.5, 0.0,          0.0,      0.5,1.0,     1.0,         0.5, 0.0,          0.0,      0.5,1.0,     1.0,         0.5, 0.0,          0.0,      0.5,1.0,     1.0,         0.5};

  float twentyfourl[]={0.0, 1.0/24.0,1.0/12.0, 0.125,1.0/6.0, 5.0/24.0, 0.25, 7.0/24.0, 1.0/3.0, 0.375,5.0/12.0, 11.0/24.0, 0.5, 13.0/24.0,7.0/12.0, 5.0/8.0,2.0/3.0, 17.0/24.0, 3.0/4.0, 19.0/24.0, 5.0/6.0, 7.0/8.0,11.0/12.0, 23.0/24.0, 1.0};
  float twentyfourr[]={0.0, 0.0,     0.5,     1.0, 1.0,     0.5,      0.0, 0.0,      0.5,     1.0,1.0,     0.5,       0.0, 0.0,     0.5,     1.0, 1.0,     0.5,      0.0, 0.0,      0.5,     1.0,1.0,     0.5,       0.0};
  float twentyfourg[]={1.0, 0.5,     0.0,     0.0, 0.5,     1.0,      1.0, 0.5,      0.0,     0.0,0.5,     1.0,       1.0, 0.5,     0.0,     0.0, 0.5,     1.0,      1.0, 0.5,      0.0,     0.0,0.5,     1.0,       1.0};
  float twentyfourb[]={0.5, 1.0,     1.0,     0.5, 0.0,     0.0,      0.5, 1.0,      1.0,     0.5,0.0,     0.0,       0.5, 1.0,     1.0,     0.5, 0.0,     0.0,      0.5, 1.0,      1.0,     0.5,0.0,     0.0,       0.5};

  switch(type){
  case -2:
    cpgctab(rl, rr, rg, rb, r, contra, bright);
    break;
  case -1:
    cpgctab(hl, hr, hg, hb, h, contra, bright);
    break;
  case 0:
    cpgctab(hexl1, hexr1, hexg1, hexb1, hex, contra, bright);
    break;
  case 1:
    cpgctab(hexl2, hexr2, hexg2, hexb2, hex, contra, bright);
    break;
  case 2:
    cpgctab(hexl3, hexr3, hexg3, hexb3, hex, contra, bright);
    break;
  case 3:
    cpgctab(dozenl, dozenr, dozeng, dozenb, dozen, contra, bright);
    break;
  case 4:
    cpgctab(eighteenl, eighteenr, eighteeng, eighteenb, eighteen, contra, bright);
    break;
  case 5:
    cpgctab(twentyfourl, twentyfourr, twentyfourg, twentyfourb, twentyfour, contra, bright);
    break;
  case 2+PALLET_SHIFT_SHIFT:
    cpgctab(hexl3, hexr3_1, hexg3_1, hexb3_1, hex, contra, bright);
    break;
  default:
    printf("Unknown pallet %i\n",type);
  }
}
