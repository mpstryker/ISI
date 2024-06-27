#include "iman1.h"

#define SERIES_WIN_SIZE 14.0

void DisplayImage(TRAIN *tr,int fr,float *imbuffer){
  static char str[10];

  cpgslct(mainid);
  if(n_panels==4){
    cpgpanl(1+panel%2,1+panel/2);
    panel=(1+panel)%4;
  }
  cpgsch(1.4);
  cpggray(imbuffer,Xdim,Ydim,1,Xdim,1,Ydim,fg,bg,trans);
  
  cpgsci(0);
  cpgtext(1.0,pminy*1.01,str);
  sprintf(str,"%i",fr);
  cpgsci(7);
  cpgtext(1.0,pminy*1.01,str);
  
  cpgsci(2);
  cpgpt1((float)(1+ifor%Xdim),(float)(1+ifor/Xdim),2);
  cpgsci(3);
  cpgpt1((float)(1+ibac%Xdim),(float)(1+ibac/Xdim),2);
}

void DisplayImages(TRAIN *tr,int fr1,int fr2,float *imbuffer){
  static char str[10];

  cpgslct(mainid);
  if(n_panels==4){
    cpgpanl(1+panel%2,1+panel/2);
    panel=(1+panel)%4;
  }
  cpgsch(1.5);
  cpggray(imbuffer,Xdim,Ydim,1,Xdim,1,Ydim,fga/image_count,bga/image_count,trans);
  cpgsci(0);
  cpgtext(1.0,pminy*1.015,str);
  sprintf(str,"%i-%i",fr1,fr2);
  cpgsci(7);
  cpgtext(1.0,pminy*1.015,str);
  
  cpgsci(2);
  cpgpt1((float)(1+ifor%Xdim),(float)(1+ifor/Xdim),2);
  cpgsci(3);
  cpgpt1((float)(1+ibac%Xdim),(float)(1+ibac/Xdim),2);
}

void DisplayMap(TRAIN *tr,double *mx,double *my){
  int i,j,m;
  int x,y,x0,y0;
  char str[256],small_str[64],*pc;
  double dumd,*dbuffer;
  float fmaxx,fminx,fmaxy,fminy,*fbufferx,*fbuffery,dumf;
  float fmax,fmin;
  float contour;
  float wedge_offset=0.05,wedge_width=1.4;
  FILE *fp;

  cpgslct(mainid);
  if(n_panels==4){
    cpgpanl(1+panel%2,1+panel/2);
    panel=(1+panel)%4;
  }

  if(!(fbufferx=(float *)calloc(XYdim,sizeof(float)))){
    printf("DIS Cannot allocate for fbufferx\n");
    return;
  }
  if(!(fbuffery=(float *)calloc(XYdim,sizeof(float)))){
    printf("DIS Cannot allocate for fbuffery\n");
    return;
  }

  if(radiusbig!=0){
    if(!(dbuffer=(double *)calloc(XYdim,sizeof(double)))){
      printf("DIS Cannot allocate for dbuffer\n");
      return;
    }
    fminx=10000000000.0;
    fmaxx=-fminx;
    for(i=0;i<XYdim;i++){
      dumd=0;
      m=0;
      x0=i%Xdim;
      y0=i/Xdim;
      for(j=0;j<aindexbig_n;j++){
	x=x0+aindexbig_x[j];
	y=y0+aindexbig_y[j];
	if(x>=0 && x<Xdim && y>=0 && y<Ydim){ 
	  dumd+=mx[x+y*Xdim];
	  m++;
	}
      }
      dumd/=(double)m;
      dbuffer[i]=mx[i]-dumd;
      fbufferx[i]=(float)dumd;
      if(fmaxx<(float)dumd) fmaxx=(float)dumd;
      if(fminx>(float)dumd) fminx=(float)dumd;
    }
    memcpy(mx,dbuffer,XYdim*sizeof(double));

    printf("DIS XBias: min=%f max=%f\n",fminx,fmaxx);
    cpggray(fbufferx,Xdim,Ydim,1,Xdim,1,Ydim,fmaxx,fminx,trans);
    cpgwedg("B", -0.7, 2.3, fmaxx, fminx,"");

    getchar();

    memset(dbuffer,0,XYdim*sizeof(double));
    fminy=10000000000.0;
    fmaxy=-fminy;
    for(i=0;i<XYdim;i++){
      dumd=0;
      m=0;
      x0=i%Xdim;
      y0=i/Xdim;
      for(j=0;j<aindexbig_n;j++){
	x=x0+aindexbig_x[j];
	y=y0+aindexbig_y[j];
	if(x>=0 && x<Xdim && y>=0 && y<Ydim){ 
	  dumd+=my[x+y*Xdim];
	  m++;
	}
      }
      dumd/=(double)m;
      dbuffer[i]=my[i]-dumd;
      fbuffery[i]=(float)dumd;
      if(fmaxy<(float)dumd) fmaxy=(float)dumd;
      if(fminy>(float)dumd) fminy=(float)dumd;
    }
    memcpy(my,dbuffer,XYdim*sizeof(double));

    free(dbuffer);

    printf("DIS YBias: min=%f max=%f\n",fminy,fmaxy);
    cpggray(fbuffery,Xdim,Ydim,1,Xdim,1,Ydim,fmaxy,fminy,trans);

    getchar();

    Pallet(2,1.0,0.5);

    for(i=0;i<XYdim;i++){
      fbuffer[i]=(float)atan2((double)fbuffery[i],(double)fbufferx[i]);
    }
    
    cpgimag(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,-M_PI,M_PI,trans);
    cpgwedg("BI", -0.7, 2.3, -180.0/harmonic, 180.0/harmonic,"");
    
    printf("DIS Save bias images?(y/n): ");
    fgets(str,4,stdin);
    if(*str == 'y'){
      sprintf(str,"%s%s%i%s%i%s%i%s%i%s",tr->filename,".xbias",radius,"_",radiusbig,"_",increment1,"_",initframe,".ps/CPS");
      printf("DIS %s\n",str);
      cpgopen("?");
      cpgpap(4.0,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
      cpgsch(0.1);
      cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
      cpggray(fbufferx,Xdim,Ydim,1,Xdim,1,Ydim,fmaxx,fminx,trans);
      cpgsch(1.0);
      cpgwedg("B", -0.7, 2.3, fmaxx, fminx,"");
      cpgclos();

      sprintf(str,"%s%s%i%s%i%s%i%s%i%s",tr->filename,".ybias",radius,"_",radiusbig,"_",increment1,"_",initframe,".ps/CPS");
      printf("DIS %s\n",str);
      cpgopen("?");
      cpgpap(4.0,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
      cpgsch(0.1);
      cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
      cpggray(fbuffery,Xdim,Ydim,1,Xdim,1,Ydim,fmaxy,fminy,trans);
      cpgsch(1.0);
      cpgwedg("B", -0.7, 2.3, fmaxy, fminy,"");
      cpgclos();
      cpgslct(mainid);
    }
  }
 
  cpgask(0);
  //  cpgpap(12.0,0.36*(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
  cpgpap(10.0,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
  cpgsubp(2,2);
  cpgsch(0.5);
  cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
  cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
  cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
  cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
  cpgsch(1.8);
  Pallet(2,1.0,0.5);

  for(i=0;i<XYdim;i++){
    fbuffer[i]=(float)atan2(my[i],mx[i]);
  }

  cpgpanl(1,1);
  cpgimag(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,-M_PI,M_PI,trans);
  cpgsci(1);
  cpgwedg("BI", wedge_offset,wedge_width, -180.0/harmonic,180.0/harmonic,"");

  fminx=10000000000.0;
  fmaxx=-fminx;
  fminy=10000000000.0;
  fmaxy=-fminy;
  fmin=10000000000.0;
  fmax=-fmin;
  for(i=0;i<XYdim;i++){
    dumf=fbufferx[i]=(float)mx[i];
    if(fmaxx<dumf) fmaxx=dumf;
    if(fminx>dumf) fminx=dumf;
    dumf=fbuffery[i]=(float)my[i];
    if(fmaxy<dumf) fmaxy=dumf;
    if(fminy>dumf) fminy=dumf;
    dumf=fbuffer[i]=(float)hypot(my[i],mx[i]);
    if(fmax<dumf) fmax=dumf;
    if(fmin>dumf) fmin=dumf;
  }
  
  if(do_contours){
    contour=0.0;
    cpgsci(0);
    cpgcont(fbufferx,Xdim,Ydim,1,Xdim,1,Ydim,&contour,-1,trans);
    cpgsci(1);
    cpgcont(fbuffery,Xdim,Ydim,1,Xdim,1,Ydim,&contour,-1,trans);
  }
  
  cpgpanl(1,2);
  cpggray(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,fmax,fmin,trans);
  cpgsci(1);
  cpgwedg("B",wedge_offset,wedge_width, fmax, fmin,"");
  if(do_contours){
    contour=0.0;
    cpgsci(0);
    cpgcont(fbufferx,Xdim,Ydim,1,Xdim,1,Ydim,&contour,-1,trans);
    cpgsci(1);
    cpgcont(fbuffery,Xdim,Ydim,1,Xdim,1,Ydim,&contour,-1,trans);
  }


  cpgpanl(2,1);
  cpggray(fbufferx,Xdim,Ydim,1,Xdim,1,Ydim,fmaxx,fminx,trans);
  cpgsci(1);
  cpgwedg("B",wedge_offset,wedge_width, fmaxx, fminx,"");
  if(do_contours){
    contour=0.0;
    cpgsci(0);
    cpgcont(fbufferx,Xdim,Ydim,1,Xdim,1,Ydim,&contour,-1,trans);
  }
  
  cpgpanl(2,2);
  cpggray(fbuffery,Xdim,Ydim,1,Xdim,1,Ydim,fmaxy,fminy,trans);
  cpgsci(1);
  cpgwedg("B",wedge_offset,wedge_width, fmaxy, fminy,"");
  if(do_contours){
    contour=0.0;
    cpgsci(1);
    cpgcont(fbuffery,Xdim,Ydim,1,Xdim,1,Ydim,&contour,-1,trans);
  }


  printf("DIS Save maps ?(y/n): ");
  fgets(str,4,stdin);
  if(*str == 'y'){
    pc=small_str;
    if(iDoDivisionByAverage){ 
      *pc='d';
      pc++;
    }
    if(do_timeaveraging) sprintf(pc,"t%i",timeaverage_cycles_int);
    else
      if(do_remove_curve) sprintf(pc,"c%i",poly_fitn);
      else
	if(do_remove_linear) sprintf(pc,"l");
	else *pc='\0';

    i=(SAVE_TYPE_BIT_DOUBLE_FLOAT*DOUBLE) | 
      (SAVE_TYPE_BIT_ADD_DOUBLE_AT_END) | 
      (SAVE_TYPE_BIT_XY_IN_FILES*save_in_files);
    if(save_in_files==SAVE_XY_IN_TWO_FILES){
      if(radius!=0 || radiusbig!=0){
	sprintf(str,"%s%s%i%s%i%s%i%s%i%s",tr->filename,".mapx",radius,"_",radiusbig,"_",increment1,"_",initframe,small_str);
      }
      else{
	if(increment1==1) sprintf(str,"%s%s%i%s%i%s%i%s%i%s",tr->filename,".mapxraw",(int)harmonic,"_",n_cycles,"_",initframe,"_",finalframe,small_str);
	else sprintf(str,"%s%s%i%s%i%s%i%s%i%s%i%s",tr->filename,".mapxraw",(int)harmonic,"_",n_cycles,"_",initframe,"_",finalframe,"_",increment1,small_str);
      }
      //      printf("DIS %s\n",str);
      if((fp=fopen(str,"w"))==NULL){
	printf("DIS Cannot open file %s\n",str);
      }
      if(fwrite((void*)&i,sizeof(int),1,fp) != 1){
	printf("DIS Cannot write type to file %s\n",str);
      }
      if(fwrite((void*)&Xdim,sizeof(unsigned short),1,fp) != 1){
	printf("DIS Cannot write Xdim to file %s\n",str);
      }
      if(fwrite((void*)&Ydim,sizeof(unsigned short),1,fp) != 1){
	printf("DIS Cannot write Ydim to file %s\n",str);
      }
      if(fwrite((void*)mx,XYdim*sizeof(double),1,fp) != 1){
	printf("DIS Cannot write mapx to file %s\n",str);
      }
      if(fwrite((void*)&harmonic,sizeof(double),1,fp) != 1){
	printf("DIS Cannot write harmonic to file %s\n",str);
      }
      fclose(fp);
      
      if(radius!=0 || radiusbig!=0){
	sprintf(str,"%s%s%i%s%i%s%i%s%i%s",tr->filename,".mapy",radius,"_",radiusbig,"_",increment1,"_",initframe,small_str);
      }
      else{
	if(increment1==1) sprintf(str,"%s%s%i%s%i%s%i%s%i%s",tr->filename,".mapyraw",(int)harmonic,"_",n_cycles,"_",initframe,"_",finalframe,small_str);
	else sprintf(str,"%s%s%i%s%i%s%i%s%i%s%i%s",tr->filename,".mapyraw",(int)harmonic,"_",n_cycles,"_",initframe,"_",finalframe,"_",increment1,small_str);
      }
      //      printf("DIS %s\n",str);
      if((fp=fopen(str,"w"))==NULL){
	printf("DIS Cannot open file %s\n",str);
      }
      if(fwrite((void*)&i,sizeof(int),1,fp) != 1){
	printf("DIS Cannot write type to file %s\n",str);
      }
      if(fwrite((void*)&Xdim,sizeof(unsigned short),1,fp) != 1){
	printf("DIS Cannot write Xdim to file %s\n",str);
      }
      if(fwrite((void*)&Ydim,sizeof(unsigned short),1,fp) != 1){
	printf("DIS Cannot write Ydim to file %s\n",str);
      }
      if(fwrite((void*)my,XYdim*sizeof(double),1,fp) != 1){
	printf("DIS Cannot write mapy to file %s\n",str);
      }
      if(fwrite((void*)&harmonic,sizeof(double),1,fp) != 1){
	printf("DIS Cannot write harmonic to file %s\n",str);
      }
      fclose(fp);
    }
  
    if(save_in_files==SAVE_XY_IN_ONE_FILE){
      if(radius!=0 || radiusbig!=0){
	sprintf(str,"%s.map%i_%i_%i_%i%s",tr->filename,radius,radiusbig,increment1,initframe,small_str);
      }
      else{
	if(increment1==1) sprintf(str,"%s.mapraw%i_%i_%i_%i%s",tr->filename,(int)harmonic,n_cycles,initframe,finalframe,small_str);
	else sprintf(str,"%s.mapraw%i_%i_%i_%i_%i%s",tr->filename,(int)harmonic,n_cycles,initframe,finalframe,increment1,small_str);
      }
      //      printf("DIS %s\n",str);
      if((fp=fopen(str,"w"))==NULL){
	printf("DIS Cannot open file %s\n",str);
      }
      if(fwrite((void*)&i,sizeof(int),1,fp) != 1){
	printf("DIS Cannot write type to file %s\n",str);
      }
      if(fwrite((void*)&Xdim,sizeof(unsigned short),1,fp) != 1){
	printf("DIS Cannot write Xdim to file %s\n",str);
      }
      if(fwrite((void*)&Ydim,sizeof(unsigned short),1,fp) != 1){
	printf("DIS Cannot write Ydim to file %s\n",str);
      }
      if(fwrite((void*)mx,XYdim*sizeof(double),1,fp) != 1){
	printf("DIS Cannot write mapx to file %s\n",str);
      }
      if(fwrite((void*)my,XYdim*sizeof(double),1,fp) != 1){
	printf("DIS Cannot write mapy to file %s\n",str);
      }
      if(fwrite((void*)&harmonic,sizeof(double),1,fp) != 1){
	printf("DIS Cannot write harmonic to file %s\n",str);
      }
      fclose(fp);      
    }
  }

  if(prompt_to_save_images){
    printf("DIS Save map image?(y/n): ");
    fgets(str,4,stdin);
    if(*str == 'y'){
      sprintf(str,"%s%s%i%s%i%s%i%s%i%s",tr->filename,".map",radius,"_",radiusbig,"_",increment1,"_",initframe,".ps/CPS");
      printf("DIS %s\n",str);
      cpgopen("?");
      cpgpap(4.0,(float)fabs((double)((pmaxy-pminy)/(pmaxx-pminx))));
      cpgsch(0.1);
      cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
      Pallet(2,1.0,0.5);
      cpgimag(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,-M_PI,M_PI,trans);
      cpgsch(1.0);
      cpgwedg("BI", -0.7, 2.3, -180.0/harmonic, 180.0/harmonic,"");
      cpgclos();
      cpgslct(mainid);
    }
  }

  free(fbufferx);
  free(fbuffery);
}

void DisplayVectorField2(int dx,int dy,double *mx,double *my,float *fmx,float *fmy,int mode){
  unsigned long dxy;
  int i,callerid;
  float *fbufx,*fbufy;
  float fmaxx,fminx,fmaxy,fminy,dumf;
  float fmax,fmin;
  char ch;

  cpgqid(&callerid);

  dxy=(unsigned long)dx*(unsigned long)dy;
  if(!(fbufx=(float*)malloc(dxy*sizeof(float)))){
    printf("DIS Cannot allocate for fbufx\n");
    return;
  }
  if(!(fbufy=(float*)malloc(dxy*sizeof(float)))){
    printf("DIS Cannot allocate for fbufy\n");
    free(fbufx);
    return;
  }
  fminx=0.0;
  fmaxx=dx-1.0;
  fmaxy=0.0;
  fminy=dy-1.0+5.0;

  fmin=1.0e10;
  fmax=-fmin;
  if(mode==DOUBLE){
    for(i=0;i<dxy;i++){
      dumf=fbufx[i]=(float)hypot(my[i],mx[i]);
      if(fmax<dumf) fmax=dumf;
      if(fmin>dumf) fmin=dumf;
      fbufy[i]=(float)atan2(my[i],mx[i]);
    }
  }
  if(mode==FLOAT){
    for(i=0;i<dxy;i++){
      dumf=fbufx[i]=(float)hypot((double)fmy[i],(double)fmx[i]);
      if(fmax<dumf) fmax=dumf;
      if(fmin>dumf) fmin=dumf;
      fbufy[i]=(float)atan2((double)fmy[i],(double)fmx[i]);
    }
  }

  cpgopen("/xw");
  cpgask(0);
  cpgpap(10.0,0.5*(float)fabs((fmaxy-fminy)/(fmaxx-fminx)));
  cpgsubp(2,1);
  cpgsch(0.5);
  cpgenv(fminx,fmaxx,fminy,fmaxy,1,-2);
  Pallet(2,1.0,0.5);
  cpgimag(fbufy,dx,dy,1,dx,1,dy,-M_PI,M_PI,trans);
  cpgsch(1.0);
  cpgwedg("BI", -0.7, 2.3, -180.0/harmonic, 180.0/harmonic,"");
  cpgenv(fminx,fmaxx,fminy,fmaxy,1,-2);
  cpggray(fbufx,dx,dy,1,dx,1,dy,fmax,fmin,trans);
  cpgwedg("B", -0.7, 2.3, fmax, fmin,"");
  cpgask(0);
  cpgband(0,0,0.0,0.0,&dumf,&dumf,&ch); 
  cpgclos();
  cpgslct(callerid);
  free(fbufx);
  free(fbufy);
}

void DisplaySeries(TRAIN *tr,int n,double **buf){
  int callerid,seriesid;
  int panelx,panely;
  int i,nn,ns;
  float dumf;
  char ch;

  panelx=n;
  panely=1;
  for(nn=n;nn<n+3;nn++){
    ns=(int)floor(sqrt((double)nn));
    for(i=ns;i>0;i--){
      if(!(nn%i)){
	if((nn/i-i)<(panelx-panely)){
	  panelx=nn/i;
	  panely=i;
	}
      }
    }
  }
  //  printf("DIS Panels X=%i Y=%i\n",panelx,panely);
  cpgqid(&callerid);
  seriesid=cpgopen("/xw");
  cpgsubp(panelx,panely);
  cpgpap(SERIES_WIN_SIZE,(float)fabs(((double)panely/(double)panelx)*(double)((pmaxy-pminy)/(pmaxx-pminx))));
  for(i=0;i<n;i++){
    Double2Float(tr->XYdim,buf[i],fbuffer);
    printf("DIS B %3i Min=%f Max=%f\n",i,bg,fg);
    cpgsch(0.4);
    cpgenv(pminx,pmaxx,pminy,pmaxy,1,-2);
    cpggray(fbuffer,Xdim,Ydim,1,Xdim,1,Ydim,fg,bg,trans);
    cpgsch(4.0);
    cpgsci(1);
    cpgwedg("R", -0.0, 0.9, fg, bg,"");
    
    cpgsci(2);
    cpgpt1((float)(1+ifor%Xdim),(float)(1+ifor/Xdim),2);
    cpgsci(3);
    cpgpt1((float)(1+ibac%Xdim),(float)(1+ibac/Xdim),2);
  }

  cpgask(0);
  cpgband(0,0,0.0,0.0,&dumf,&dumf,&ch); 
  cpgclos();
  cpgslct(callerid);

}
