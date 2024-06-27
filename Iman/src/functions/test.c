#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <functions.h>

int main(int argc,char **argv){
  double x,y,z,t;
  x=atof(argv[1]);
  if(argc>2) y=atof(argv[2]);
  if(argc>3) z=atof(argv[3]);
  if(argc>4) t=atof(argv[4]);
  
  //  printf("I%i(%f,%i)=%30.20f\n",(int)y,x,(int)z,BesselIe((int)y,x,(int)z));
  printf("I%i(%f)=%e\n",(int)y,x,besselIne((int)y,x));
  return(0);
}
