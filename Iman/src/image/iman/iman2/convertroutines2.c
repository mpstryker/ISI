#include "iman1.h"

int ConvertImage(TRAIN *tr,int fr,int a){
  int i,j;
  double co,si,dumd;
  double dumd1=1.0;
  CAR *pcar;

  if(!(pcar=AddFrame(tr,fr,a))){
    printf("CON AddFrame return 0 fr=%i\n",fr);
    return(1);
  }

  if(iDoDivisionByAverage){
    iDivisionByAverageCount++;
    for(i=0;i<tr->XYdim;i++) map0[i]+=(double)pcar->pimage[i];
  }

  if(synch_sin){
    co=synch_cos[fr-initframe];
    si=synch_sin[fr-initframe];    
  }
  else{
    if(synch_phi) dumd=harmonic*synch_phi[fr-initframe];
    else dumd=phi0+dphi*fr;
    co=cos(dumd);
    si=sin(dumd);
  }

  if(do_timeaveraging){    
    if(AdvanceTimeAverage(tr,fr,a,&timeaverage_buffer_frames_in,&timeaverage_buffer)){
      printf("CON Cannot advance time average for frame %i\n",fr);
      return(5);
    }
    if(!(pcar=AddFrame(tr,fr,a))){
      printf("CON AddFrame return 0 fr=%i\n",fr);
      return(6);
    }

    for(i=0;i<tr->XYdim;i++){
      mapx[i]+=(dumd=(double)pcar->pimage[i]-timeaverage_buffer[i]/timeaverage_buffer_frames_in)*co;
      mapy[i]+=dumd*si;
    }
    return(0);
  }
  else{

    if(do_remove_curve){
      memset(remove_dum,0,poly_fitn*sizeof(double));
      PolynomialFitM(nframes_good,poly_fitn,(double)(fr-initframe),&dumd1,remove_dum,DISASSEMBLE);
      for(j=0;j<poly_fitn;j++){ 
	remove_cos[j]+=remove_dum[j]*co;
	remove_sin[j]+=remove_dum[j]*si;
      }
      for(i=0;i<tr->XYdim;i++){
	mapx[i]+=(dumd=(double)pcar->pimage[i])*co;
	mapy[i]+=dumd*si;
      
	for(j=0;j<poly_fitn;j++){ 
	  remove_curve[j][i]+=remove_dum[j]*dumd;
	}
      }
      return(0);
    }
    else{

      if(do_remove_linear){
	for(i=0;i<tr->XYdim;i++){
	  mapx[i]+=(dumd=(double)pcar->pimage[i])*co;
	  mapy[i]+=dumd*si;

	  remove_sy[i]+=dumd;
	  remove_sxy[i]+=dumd*(double)(fr-initframe);
	}
	return(0);
      }
    }
  }

  for(i=0;i<tr->XYdim;i++){
    mapx[i]+=(dumd=(double)pcar->pimage[i])*co;
    mapy[i]+=dumd*si;   
  }
  
  //  printf("%i ",fr);
  if(do_verbose) printf("CON %4i F %5i(%3i %3i) B %5i(%3i %3i)\n",fr,(int)fg,ifor%Xdim,ifor/Xdim,(int)bg,ibac%Xdim,ibac/Xdim);
  return(0);
}

int ConvertImages(TRAIN *tr,int fr1,int fr2,int a){
  int i,j;
  float dumf;
  double dumd,co,si;
  CAR *pcar1,*pcar2;

  if(fr1<0 || fr1>=tr->max_n){
    printf("CON Frame number(%i) is out of range(%i-%i)\n",fr1,0,tr->max_n-1);
    return(2);
  }
  if(fr2<0 || fr2>=tr->max_n){
    printf("CON Frame number(%i) is out of range(%i-%i)\n",fr2,0,tr->max_n-1);
    return(2);
  }
  if(scheme == SCHEME_DIFF2){
    pcar1=AddFrame(tr,fr1,a);
    pcar2=AddFrame(tr,fr2,a);
    fg=0.0;
    bg=(float)(1<<16);
    dumd=0.0;
    j=fr1>fr2?fr1-fr2:fr1-fr2+tr->max_n;
    for(i=0;i<tr->XYdim;i++){
      if(difference == DIFF_RELATIVE){
	dumf=fbuffer[i]=(float)pcar1->pimage[i]-(float)pcar2->pimage[i];
      }
      else{
	dumf=fbuffer[i]=pcar1->pimage[i]>pcar2->pimage[i]?(float)(pcar1->pimage[i]-pcar2->pimage[i]):(float)(pcar2->pimage[i]-pcar1->pimage[i]);
      }
      if(fg<dumf){ 
	fg=dumf;
	ifor=i;
      }
      if(bg>dumf){ 
	bg=dumf;
	ibac=i;
      }
    }
  }
  else{
    pcar1=AddFrame(tr,fr2,a);
    memset((void*)fbuffer,0,XYdim*sizeof(float));
    for(j=fr2+1;j<=fr1;j++){
      pcar2=AddFrame(tr,j,a);
      for(i=0;i<tr->XYdim;i++){
	fbuffer[i]+=pcar1->pimage[i]>pcar2->pimage[i]?(float)(pcar1->pimage[i]-pcar2->pimage[i]):(float)(pcar2->pimage[i]-pcar1->pimage[i]);
      }
      pcar1=pcar2;
    }
    fg=0.0;
    bg=(float)(1<<16);
    dumd=phi0+dphi*(fr2+fr1)*0.5;
    co=cos(dumd);
    si=sin(dumd);
    j=fr1>fr2?fr1-fr2:fr1-fr2+tr->max_n;
    for(i=0;i<tr->XYdim;i++){
      dumf=(fbuffer[i]/=j);
      if(fg<dumf){ 
	fg=dumf;
	ifor=i;
      }
      if(bg>dumf){ 
	bg=dumf;
	ibac=i;
      }

      mapx[i]+=dumf*co;
      mapy[i]+=dumf*si;

    }
  }
  bga+=bg;
  fga+=fg;
  image_count++;
  printf("CON %.04i-%.04i=%i F %5.2f %5.2f(%.03i %.03i) B %5.2f %5.2f(%.03i %.03i)\n",fr1,fr2,j,fg,fga/image_count,ifor%Xdim,ifor/Xdim,bg,bga/image_count,ibac%Xdim,ibac/Xdim);
  return(0); 
}

int AccumulateTimeAverage(TRAIN *tr,int a,double *t_l,double **t_b){
  int i,j;
  CAR *pcar;
  
  if(!(*t_b)){
    printf("CON AccumulateTimeAverage timeaverage_buffer=NULL\n");
    return(1);
  }

  if(*t_l<0.0){
    *t_l=0.0;
    for(j=initframe_ta;j<initframe+timeradius;j++){
      if(!(pcar=AddFrame(tr,j,a))){
	printf("CON AddFrame return 0 at accumulation for %i\n",j);     
	return(2);
      }
      (*t_l)+=1.0;
      for(i=0;i<tr->XYdim;i++){
	(*t_b)[i]+=(double)pcar->pimage[i];
      }
    }
    printf("CON Accumulated %i frames\n",(int)(*t_l));
  }
  return(0);
}

int AdvanceTimeAverage(TRAIN *tr,int fr,int a,double *t_l,double **t_b){
  CAR *pcar;
  int i;

  if(fr <= finalframe_ta-timeradius){
    if(!(pcar=AddFrame(tr,fr+timeradius,a))){
      printf("CON AddFrame return 0 at time+ for %i\n",fr+timeradius);     
      return(4);
    }
    for(i=0;i<tr->XYdim;i++){
      (*t_b)[i]+=(double)pcar->pimage[i];
    }
    (*t_l)+=1.0;
  }
  if(fr>initframe_ta+timeradius){
    if(!(pcar=AddFrame(tr,fr-timeradius,a))){
      printf("CON AddFrame return 0 at time- for %i\n",fr-timeradius);     
      return(5);
    }
    for(i=0;i<tr->XYdim;i++){
      (*t_b)[i]-=(double)pcar->pimage[i];
    }
    (*t_l)-=1.0;
  }
  return(0);
}
