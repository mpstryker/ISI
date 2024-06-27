#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void four1(float *data,int nn,int isign);
void realft(float *data,int n,int isign);

int main(int argc,char *argv[]){
  float *pfData,fNorm,fSum,fPower;
  int iNSamples,iNS;
  int i;

  iNSamples=16;
  if(!(pfData=(float *)malloc(iNSamples*sizeof(float)))){
    printf("Cannot malloc for pfData");
    return(1);
  }
  iNS=16;
  fNorm=1.0/sqrt(iNS);
  pfData[0]=-1.0+1;
  pfData[1]=1.0;
  pfData[2]=-1.0-1;
  pfData[3]=1.0;
  pfData[4]=-1.0;
  pfData[5]=1.0;
  pfData[6]=-1.0;
  pfData[7]=1.0;
  for(i=0;i<iNS;i++) pfData[i]=rand()/(double)RAND_MAX;
  fSum=0.0;
  for(i=0;i<iNS;i++){
    fSum+=pfData[i]*pfData[i];
    printf("%.3f ",pfData[i]);
  }
  printf("Sum   %f\n",fSum);
  realft(pfData-1,iNS,1);
  fPower=0.0;
  for(i=0;i<iNS;i++){
    pfData[i]*=fNorm;
    fPower+=pfData[i]*pfData[i];
    printf("%.3f ",pfData[i]);
  }
  fPower*=2.0;
  fPower-=pfData[0]*pfData[0]+pfData[1]*pfData[1];
  printf("Power %f\n",fPower);
  //  realft(pfData-1,iNS,-1);
  //  for(i=0;i<iNS;i++) printf("%.3f ",pfData[i]);
  printf("P/S %f\n",fPower/fSum);
  return(0);
}
