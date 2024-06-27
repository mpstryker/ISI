/* Bin analysis                 */
/* Written by V.Kalatsky        */
/* 20Jul2001                    */

#include "binan1.h"

#define MOVIE_WIN_SIZE 5.0

void DisplayBins(int mode,int xd,int yd,int n,float **bins){
  int callerid;
  int i;
  unsigned long xy;
  float dumfx,dumfy;
  float fmin,fmax;
  unsigned long ulmin,ulmax;
  char ch;

  xy=(unsigned long)xd*(unsigned long)yd;

  cpgqid(&callerid);
  cpgslct(mainid);
  for(i=0;i<n;i++){
    cpgpanl(1+i%panelnx,1+i/panelnx);
    FindExtremaEx(xy,(unsigned long)xd,bins[i],&fmin,&fmax,&ulmin,&ulmax);
    if(do_verbose) printf("DIS B %3i Min=%f Max=%f\n",i,fmin,fmax);
    cpggray(bins[i],xd,yd,1,xd,1,yd,fmax,fmin,trans);
    cpgsci(1);
    cpgwedg("R", -0.0, 0.9, fmax, fmin,"");
    
    cpgsci(2);
    cpgpt1((float)(1+ulmax%Xdim),(float)(1+ulmax/Xdim),2);
    cpgsci(3);
    cpgpt1((float)(1+ulmin%Xdim),(float)(1+ulmin/Xdim),2);
  }

  cpgask(0);
  while(1){
    cpgband(0,0,0.0,0.0,&dumfx,&dumfy,&ch); 
    if(ch=='q' || ch=='D') break;
    printf("DIS X=%f Y=%f\n",dumfx,dumfy);
  }
  if(callerid) cpgslct(callerid);
}

void AnimateBins(int mode,int xd,int yd,int n,float **bins,unsigned long delay){
  int callerid=0;
  int *rv;
  char ch;

  cpgqid(&callerid);

  animatebinsthreadargs.mode=mode;
  animatebinsthreadargs.xd=xd;
  animatebinsthreadargs.yd=yd;
  animatebinsthreadargs.n=n;
  animatebinsthreadargs.bins=bins;
  animatebinsthreadargs.delay=delay;

  pthread_attr_init(&animate_thread_attr);

  pthread_create(&animate_thread,&animate_thread_attr,(void *)&AnimateBinsThread,&animatebinsthreadargs);
  
  ch=(int)getchar();
  
  stop_animation=1;
  rv=&thread_return_value;
  pthread_join (animate_thread,(void**)(&rv));

  if(callerid) cpgslct(callerid);

  printf("DIS Stopped animation\n");
}

void AnimateBinsThread(void *arg){
  int i;
  char str[64];
  unsigned long xy;
  int mode,xd,yd,n;
  float **bins;
  unsigned long delay;
  AnimateBinsThreadArgs *args;
  float fmax,fmin;

  args=(AnimateBinsThreadArgs *)arg;

  mode=args->mode;
  xd=args->xd;
  yd=args->yd;
  n=args->n;
  bins=args->bins;
  delay=args->delay;

  xy=(unsigned long)xd*(unsigned long)yd;
    
  fmin=bin_min[0];
  fmax=bin_max[0];
  for(i=1;i<n;i++){
    if(fmin>bin_min[i]) fmin=bin_min[i];
    if(fmax<bin_max[i]) fmax=bin_max[i];
  }    
  
  cpgopen("/xw");
  cpgpap(MOVIE_WIN_SIZE,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
  cpgsch(1.0);
  cpgenv(pminx,pmaxx,pminy,pmaxy,1,-1);
  cpgask(0);
  while(1){
    for(i=0;i<n;i++){
      if(do_fixed_background){
	cpggray(bins[i],xd,yd,1,xd,1,yd,fmax,fmin,trans);
      }
      else{
	cpggray(bins[i],xd,yd,1,xd,1,yd,bin_max[i],bin_min[i],trans);
      }
      cpgsci(0);
      cpglab("","",str);
      sprintf(str,"%i",i);
      cpgsci(7);
      cpglab("","",str);
      usleep(delay);
      if(stop_animation){
	cpgclos();
	stop_animation=0;
	pthread_exit(0);
      }
    }
  }
}

#define SAVE_TYPE_BIT_ADD_DOUBLE_AT_END 2
#define SAVE_TYPE_BIT_XY_IN_FILES 4

void DisplayMap(int xd,int yd,double **mp){
  int i;
  float dumf;
  float fmaxx,fminx,fmaxy,fminy,fmax,fmin;
  float wedge_offset=0.05,wedge_width=1.4;
  float *fbufferp=NULL,*fbufferr=NULL,*fbufferx=NULL,*fbuffery=NULL;
  int iSaveID;
  char ch;
  int type=6;
  unsigned short dumus;
  unsigned long dumul;
  double dHarmonic=1;
  FILE *fp;
  char str[256];

  if(!(fbufferp=(float*)malloc(XYdim*sizeof(float)))){
    printf("INI Cannot allocate for fbufferp\n");
    return;
  }
  if(!(fbufferr=(float*)malloc(XYdim*sizeof(float)))){
    printf("INI Cannot allocate for fbufferr\n");
    return;
  }
  if(!(fbufferx=(float *)calloc(XYdim,sizeof(float)))){
    printf("INI Cannot allocate for fbufferx\n");
    return;
  }
  if(!(fbuffery=(float *)calloc(XYdim,sizeof(float)))){
    printf("INI Cannot allocate for fbuffery\n");
    return;
  }
  fmax=fmaxy=fmaxx=-(fmin=fminy=fminx=1.0e10);
  for(i=0;i<xd*yd;i++){
    fbufferp[i]=(float)(atan2(mp[1][i],mp[0][i]));
    dumf=fbufferx[i]=(float)(mp[0][i]);
    if(fmaxx<dumf) fmaxx=dumf;
    if(fminx>dumf) fminx=dumf;
    dumf=fbuffery[i]=(float)(mp[1][i]);
    if(fmaxy<dumf) fmaxy=dumf;
    if(fminy>dumf) fminy=dumf;
    dumf=fbufferr[i]=hypot(mp[0][i],mp[1][i]);
    if(fmax<dumf) fmax=dumf;
    if(fmin>dumf) fmin=dumf;
  }
  printf("DIS X(%f %f) Y(%f %f) R(%f %f)\n",fminx,fmaxx,fminy,fmaxy,fmin,fmax);
  cpgqid(&iSaveID);
  cpgopen(device);
  cpgsubp(1,2);
  cpgpap(5.0,(2.0*(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx)))));
  cpgsch(0.7);
  for(i=0;i<2;i++) cpgenv(pminx+1,pmaxx-1,pminy+1,pmaxy-1,1,-2);
  cpgsch(1.0);

  cpgpanl(1,1);
  cpgsch(1.0);
  Pallet(2,1.0,0.5);
  cpgimag(fbufferp,xd,yd,1,xd,1,yd,-M_PI,M_PI,trans);
  cpgsch(2.6);
  cpgwedg("BI",wedge_offset,wedge_width, -180.0, 180.0,"");
  
  cpgpanl(1,2);
  if(fmax==fmin){ fmin=0.0;fmax=2.0;}
  cpggray(fbufferr,xd,yd,1,xd,1,yd,fmax,fmin,trans);
  cpgsci(1);
  cpgsch(2.6);
  cpgwedg("B",wedge_offset,wedge_width, fmax, fmin,"");

  cpgband(7,1,fmax,fmin,&fmax,&fmin,&ch); 

  if(ch=='S'){
    sprintf(str,"%s_map%i_%i",filenames[0],nfiles,iNFComponents);
    printf("DIS Saving curent maps in %s\n",str);
    if((fp=fopen(str,"w"))==NULL){
      printf("DIS Cannot open file %s\n",str);
    }
    i=type&0xE;
    if(fwrite((void*)&i,sizeof(int),1,fp) != 1){
      printf("DIS Cannot write type %i to file %s\n",i,str);
    }
    dumus=(unsigned short)xd;
    if(fwrite((void*)&dumus,sizeof(unsigned short),1,fp) != 1){
      printf("DIS Cannot write Xdim=%i to file %s\n",dumus,str);
    }
    dumus=(unsigned short)yd;
    if(fwrite((void*)&dumus,sizeof(unsigned short),1,fp) != 1){
      printf("DIS Cannot write Ydim=%i to file %s\n",dumus,str);
    }
    dumul=(unsigned long)xd*(unsigned long)yd;
    if(fwrite((void*)mp[0],dumul*sizeof(double),1,fp) != 1){
      printf("DIS Cannot write mp[0] to file %s\n",str);
    }
    if(fwrite((void*)mp[1],dumul*sizeof(double),1,fp) != 1){
      printf("DIS Cannot write mp[1] to file %s\n",str);
    }
    if(fwrite((void*)&dHarmonic,sizeof(double),1,fp) != 1){
      printf("DIS Cannot write harmonic %f to file %s\n",dHarmonic,str);
    }
    fclose(fp);
  }

  cpgclos();
  if(iSaveID) cpgslct(iSaveID);
  free(fbufferp);
  free(fbufferr);
  free(fbufferx);
  free(fbuffery);
}
