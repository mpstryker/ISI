/* Image analysis               */
/* Written by V.Kalatsky        */
/* 10Aug01                      */

#include "iman1.h"

#define XDIM 100 
#define YDIM 100 

#define TIMER_SHIFT 0
#define TIMER_MULT 50

#define FRAME_PERIOD 100.0
#define FRAME_KX (1.0/(double)XDIM)
#define FRAME_KY 0.0
#define FRAME_WINDOW_LEFT 0.0
#define FRAME_WINDOW_RIGHT 0.1
#define FRAME_SIGMA 10.0
#define FRAME_SIGMAX 30.0
#define FRAME_SIGMAY 4.0
#define FRAME_SHIFT FRAME_SIGMAX
#define FRAME_MULT 4096.0

#define SYNCH_MAX 64
#define SYNCH_SHIFT 0.0
#define SYNCH_MULT ((double)SYNCH_MAX/FRAME_PERIOD)

#define NFRAMES 2000

long synch_simulator(int fr);
unsigned long timer_simulator(int fr);
void frame_simulator(int fr,unsigned short xd,unsigned short yd,unsigned short *buf);

int main(int argc,char *argv[]){
  int i;
  char filename[16]={"fakeimage"};
  FILE *fp;
  FILEHEADER fh;
  FRAMEHEADER frameheader;
  unsigned short *frame_buffer;
  unsigned long XYdim;
  unsigned short Xdim=XDIM,Ydim=YDIM;
  unsigned long Nframes=NFRAMES;

  XYdim=(unsigned long)Xdim*(unsigned long)Ydim;
  if(!(frame_buffer=(unsigned short *)calloc(XYdim,sizeof(unsigned short)))){
    printf("Cannot allocate for frame_buffer");
    exit(-1);
  }

  strncpy(fh.sumTag,"T_01",4);
  fh.sumBytes=sizeof(FILEHEADER)-4+Nframes*(sizeof(FRAMEHEADER)+XYdim*sizeof(unsigned short));	/* byte count of bytes to follow in file */
  fh.sumXsize=Xdim; 
  fh.sumYsize=Ydim; 
  fh.sumNframes=NFRAMES;
  fh.sumNconds=0; 
  fh.sumCondNo=0;
  fh.sumOrder=0;
  fh.sumLoClip=0.0;
  fh.sumHiClip=0.0;
  fh.sumLoPass=0;
  fh.sumHiPass=0;   
  fh.sumLDataType=0;
  fh.sumLFileType=0;
  fh.sumLSizeOf=sizeof(unsigned short);
  //char		fh.sumFree1[2] ;
  fh.sumHeadersize=512;
  fh.sumNreps=10;
  fh.sumRandomize=0; 
  fh.sumXgroup=0;	/* binning done on CCD chip */ // set from xsu.xgroup	// GB: Will be removed !
  fh.sumYgroup=0; 
  fh.sumGoDelay=0;
  fh.sumRatioFirstFrame=0; 
  fh.sumRatioLastFrame=0;
  fh.sumTempBinning=1;   /* set if we collect in full resolution without hardware binning and want to bin */
  fh.sumBit1=0;	/* Offset 64, bit coded, which blocks */ 
  fh.sumBit2=0;   /* bit coded, which frames */
  fh.sumBinning=1;/* binning done on CCD */  // now this is being used for software spatila binning
  fh.sumExposureMsec=50;  // set from xsu.exposure (double)
  fh.sumMagnification=1.0;	/* optical magnification */
  fh.sumGain=1;  
  fh.sumWavelength=610;   
  fh.sumFrameHeaderSize=256;  // used only for stream files, size of the header for each frame
  fh.sumNframeiti=0;   
  fh.sumNframestim=20;   
  //char		fh.sumFree2[20] ;
  //long int	fh.sumSysTime ;		// GB: Will be removed !
  //char		fh.sumSite[4] ;        // GB: Will be removed !
  //short int	fh.sumWinType ;        // GB: Will be removed !
  //short int	fh.sumDay ;            // GB: Will be removed !
  //short int	fh.sumMonth ;          // GB: Will be removed !
  //short int	fh.sumYear ;           // GB: Will be removed !
  //char		fh.sumDateRecorded[6] ;	/* Offset 128, ASCII */ 
  //char		fh.sumUserName[10] ;		/* ASCII */ 
  //char		fh.sumDateCreated[6] ;	/*  ASCII */ 
  strncpy(fh.sumOrigFilename,filename,16) ;	
  fh.sumPrevFilename[0]='\0';
  fh.sumNextFilename[0]='\0';
  fh.sumCompressedRecordSize=2;
  fh.sumCompressedFrameSize=0;
  fh.sumCompressedFrameNumber=0;
  //fh.sumCompressedDumbAss; //To presserve int borders
  //char	        fh.sumFree3[46] ; 
  //char	        fh.sumComments[256] ; /* Offset 256 */



  frameheader.frameheadID=0; // arbitrary for future flexibility
  frameheader.frameheadlength=256;  //8
  frameheader.seqnum=0;  //sequence number 0,1,2,...NOT including paused frames
  frameheader.time=0; //16  from the time call
  frameheader.perfcount.LowPart=0;  // from the performance counter
  frameheader.perfcount.HighPart=0;  // from the performance counter
  frameheader.perfreq.LowPart=1000;
  frameheader.perfreq.HighPart=0;
  frameheader.rep=0; //repetition number of this cycle
  frameheader.trial=0;  // number of stim in this cycle
  frameheader.condition=0;  //stimulus number
  frameheader.frame_of_cond=0; //48  which frame in relation to start of stimulus
  frameheader.frame_paused=0; //  binary variable, marks paused frames(=1[actually !=0])
  frameheader.frame_type=0; //56  frame type, iti and stim for now
  frameheader.seqnum_all=0; //sequence number 0,1,2,...INCLUDING PAUSED FRAMES
  frameheader.synch_in=0; //for synchronization input signals (on port B)

  if(!(fp=fopen(filename,"w"))){
    printf("Cannot open file %s\n",filename);
    exit(1);
  }

  if(fwrite((void*)&fh,sizeof(FILEHEADER),1,fp) != 1){
    printf("Cannot write FILEHEADER to file %s\n",filename);
    exit(2);
  }

  for(i=0;i<NFRAMES;i++){
    frameheader.seqnum=(unsigned long)i;
    frameheader.perfcount.LowPart=timer_simulator(i);
    frameheader.seqnum_all=(long)i;
    frameheader.synch_in=synch_simulator(i);
    frame_simulator(i,Xdim,Ydim,frame_buffer);
    if(fwrite((void*)&frameheader,sizeof(FRAMEHEADER),1,fp) != 1){
      printf("Cannot write FRAMEHEADER to file %s\n",filename);
      exit(3);
    }
    if(fwrite((void*)frame_buffer,XYdim*sizeof(unsigned short),1,fp) != 1){
      printf("Cannot write frame_buffer to file %s\n",filename);
      exit(4);
    }
  }
  fclose(fp);

  return(0);
}


long synch_simulator(int fr){
  return((long)(4*((long)(SYNCH_SHIFT+fr*SYNCH_MULT)%SYNCH_MAX)));
}

unsigned long timer_simulator(int fr){
  return((unsigned long)(TIMER_SHIFT+TIMER_MULT*fr));
}

void frame_simulator(int fr,unsigned short xd,unsigned short yd,unsigned short *buf){
  int i,j;
  double x,y;
  double x1,y1;
  double t;

  t=2.0*M_PI*(double)fr/FRAME_PERIOD;
  for(j=0;j<yd;j++){
    y=(double)j-(double)(yd-1)*0.5;
    for(i=0;i<xd;i++){
      x=(double)i-(double)(xd-1)*0.5;
      //Uniform
      //buf[i+xd*j]=(unsigned short)(4000+40*cos(2.0*M_PI*fr/FRAME_PERIOD)+fr);
      //Fancy
      
      x1=(x*cos(t)-y*sin(t)-FRAME_SHIFT)/FRAME_SIGMAX;
      y1=(x*sin(t)+y*cos(t))/FRAME_SIGMAY;
      //buf[i+xd*j]=(unsigned short)(FRAME_MULT*0.5*(1.0+cos(2.0*M_PI*(x*FRAME_KX+y*FRAME_KY-(double)fr/FRAME_PERIOD))));
      buf[i+xd*j]=(unsigned short)(FRAME_MULT*exp(-(x1*x1+y1*y1)));
      //buf[i+xd*j]=(unsigned short)(0.5*(1.0+cos(2.0*M_PI*fr/FRAME_PERIOD))*FRAME_MULT*((double)rand()/(double)RAND_MAX));
      
    }
  }
  return;
}
