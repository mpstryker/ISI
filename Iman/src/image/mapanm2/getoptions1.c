#include "mapanm21.h"

void options(int argc,char *argv[]){
  int c;
  int i,j;
  char *pc,*pc1;

  while (1){
    
    if((c=getopt(argc, argv, "D::H:I:L:P:R:S::cd:f:hi:opqsw:"))==-1) break;

    switch (c){

    case 'D':
      iDoDominanceHistogram=1;
      if(optarg) dDominanceThreshold=atof(optarg);
      break;
    case 'H':
      dradiusbig=atof(optarg);
      break;
    case 'I':
      iInverseMode=atoi(optarg);
      break;
    case 'L':
      dradius=atof(optarg);
      break;
    case 'P':
      npanels=atoi(optarg);
      break;

    case 'R':
      for(i=0,j=0;i<strlen(optarg);i++) if(optarg[i]==',') j++;
      if(j==3 || (j==2 && index(optarg,'x'))){
	pc=optarg;
	pc1=index(pc,',');
	*pc1='\0';
	iROIXLeft=atoi(pc);
	pc=pc1+1;
	pc1=index(pc,',');
	*pc1='\0';
	iROIYTop=atoi(pc);
	if(j==3){
	  pc=pc1+1;
	  pc1=index(pc,',');
	  *pc1='\0';
	  iROIXRight=atoi(pc);
	  iROIYBottom=atoi(pc1+1);
	}
	else{
	  pc=pc1+1;
	  pc1=index(pc,'x');
	  *pc1='\0';
	  iROIXRight=atoi(pc);
	  iROIYBottom=atoi(pc1+1);
	  iROIXRight+=iROIXLeft-1;
	  iROIYBottom+=iROIYTop-1;
	}
	if(iROIXRight<iROIXLeft){
	  i=iROIXRight;
	  iROIXRight=iROIXLeft;
	  iROIXLeft=i;
	}
	if(iROIYBottom<iROIYTop){
	  i=iROIYBottom;
	  iROIYBottom=iROIYTop;
	  iROIYTop=i;
	}
	iCropROI=1;
	//	printf("GET ROI %i %i %i %i\n",iROIXLeft,iROIYTop,iROIXRight,iROIYBottom);
      }
      else{
	printf("GET Incorrect ROI format |%s|\n",optarg);
	printf("GET Correct syntax: LT(x,y) RB(x,y)(eg. 10,20,30,40) or LT(x,y) Dim(x,y)(eg. 10,20,20x20)\n");
      }
      break;

    case 'S':
      if(optarg) pcSaveFileSuffix=strdup(optarg);
      break;

    case 'c':
      do_contours=1;
      break;
    case 'f':
      filenames[nfiles]=strdup(optarg);
      nfiles++;
      if(nfiles>MAX_FILE_N){
	printf("Accepts two files only: X map first, Y map second\n");
	exit(1);
      }
      break;
 
    case 'h':
      printf(" -D - plot dominance histogram (default off)\n");
      printf(" -H <float> - high-pass filter domain (default %f)\n",DEFAULT_RADIUSBIG);
      printf(" -I <int> - inversion of  second image (default %i)\n",MODE_INVERSE_DEFAULT);
      printf(" -L <float> - low-pass filter domain (default %f)\n",DEFAULT_RADIUS);
      printf(" -P <int> - number of panels\n");
      printf(" -R <int,int,int,int> or <int,int,intxint> - ROI (XL,YT,XR,YB or XL,YT,XDxYD)\n");
      printf(" -c - draw contour lines\n");
      printf(" -f <filename> - needs %i files (X and Y components)\n",MAX_FILE_N);
      printf(" -h - print this help\n");
      printf(" -i <int> - interactive mode\n");
      printf(" -o - do ocular dominance analysis\n");
      printf(" -p - do pinwheel analysis\n");
      printf(" -q - take sqrt of sum\n");
      printf(" -s - cancel adding Pi to negative phi'sn\n");
      printf(" -w <float> - main window width\n");
      exit(0);
      break;

    case 'i':
      do_interactive=atoi(optarg);
      break;
    case 'o':
      do_od=1;
      break;
    case 'p':
      do_pinwheels=1;
      break;
    case 'q':
      iTakeSQRTOfSum=1;
      break;
    case 's':
      make_phi_positive=0;
      break;
    case 'w':
       main_win_size=atoi(optarg);
      break;
    
    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
      exit(1);
    }
  }
  
  if(optind < argc){
    while (optind < argc){
      filenames[nfiles]=strdup(argv[optind++]);
      nfiles++;
      if(nfiles>MAX_FILE_N){
	printf("Accepts four files at most\n");
	exit(1);
      }
    }
  }  
}
