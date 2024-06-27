/* Map analysis                 */
/* Written by V.Kalatsky        */
/* 06Feb2001                    */
/* Modified 10May2001           */

#include "mapans1.h"

void options(int argc,char *argv[]){
  int c;
  int i,j;
  char *pc,*pc1=NULL,*pc2=NULL;
  
  while (1){
    
    if((c=getopt(argc, argv, "A:BC:D:E:F:G:H:I:J:L:M:O:P:R:S:TVW:Z::ab:c:d:e::f:g::hil::mo:p::rstu:vw:z"))==-1) break;

    switch (c){

    case 'A':
      //      printf("GET ANNOTATION %s\n",optarg);
      if(!index(optarg,',')){
	printf("GET Incorrect annotation string |%s|\n",optarg);
	printf("GET Correct syntax: comma separated values(eg. -A2,32,kHz)\n");
      }
      else{
	iDoWedgeAnnotation=1;
	for(i=0,pc=optarg;*pc;pc++){
	  if(*pc==','){
	    i++;
	    if(i==1){
	      *pc='\0';
	      pc1=pc+1;
	    }
	    if(i==2){ 
	      *pc='\0';
	      pc2=pc+1;
	    }
	    if(i==3){
	      printf("GET Wedge annotation: ignoring trailing characters |%s|\n",pc);
	      *pc='\0';
	      break;
	    }
	  }
	}
	//	printf("GET |%s|%s|%s|\n",optarg,pc1,pc2);
	fWedgeMin=atof(optarg);
	fWedgeMax=atof(pc1);
	if(fWedgeMin==fWedgeMax){
	  printf("GET WARNING Equal limits in the wedge\n");
	  iDoWedgeAnnotation=0;	  
	}
	if(pc2) strWedgeUnits=strdup(pc2);
	else strWedgeUnits=strdup("");
	//	printf("GET |%f|%f|%s|\n",fWedgeMin,fWedgeMax,strWedgeUnits);	
      }
      break;

    case 'B':
      iBlankSource=1;
      iDoNotPlotPhase=1;
      iOverlayContoursColor=1;
      break;

    case 'C':
      color_scheme=atoi(optarg);
      if(color_scheme<COLOR_RGB || color_scheme>COLOR_HLS){
	printf("Bad color scheme (%i), using default (%i)\n",color_scheme,COLOR_DEFAULT);
	color_scheme=COLOR_DEFAULT;
      }
      break;

    case 'D':
      display_mode=atoi(optarg);
      if(display_mode>DISPLAY_MODE6 || display_mode<DISPLAY_MODE2){
	printf("Bad display mode (%i), using default (%i)\n",display_mode,DEFAULT_DISPLAY_MODE);
	display_mode=DEFAULT_DISPLAY_MODE;
      }
      break;

    case 'E':
      dSectionRadius=atof(optarg);
      if(dSectionRadius<=0.0){
	printf("GET Non-positive dSectionRadius");
	dSectionRadius=1.0;
      }
      break;

    case 'F':
      do_fields=1;
      field_mode=atoi(optarg);
      if(field_mode>DISPLAY_MODE6 || field_mode<DISPLAY_MODE2){
	printf("Bad field display mode (%i), using default (%i)\n",field_mode,DEFAULT_FIELD_MODE);
	display_mode=DEFAULT_FIELD_MODE;
      }
      break;
 
    case 'G':
      for(i=0,j=0;i<strlen(optarg);i++) if(optarg[i]==',') j++;
      if(j==1){
	pc=optarg;
	pc1=index(pc,',');
	*pc1='\0';
	pc1++;
	if(*pc!='*'){
	  fFixedGrayScaleLo=atof(pc);
	  iFixedGrayScaleLo=1;
	}
	if(*pc1!='*'){
	  fFixedGrayScaleHi=atof(pc1);
	  iFixedGrayScaleHi=1;
	}
      }
      else{
	printf("GET Incorrect fixed gray-scale wedge format |%s|\n",optarg);
	printf("GET Correct syntax: Low,High(eg. 0.2,2.0)\n");
      }
      break;


    case 'H':
      dradiusbig=atof(optarg);
      break;

    case 'I':
      contours_iniangle=atof(optarg)*M_PI/180.0;
      break;

    case 'J':
      iFlipMaps=atoi(optarg);
      iFlipMaps=iFlipMaps%4;
      if(iFlipMaps<0) iFlipMaps+=4;
      printf("GET INVERSION %i\n",iFlipMaps);
      break;

    case 'L':
      dradius=atof(optarg);
      break;

    case 'M':
      iDoPhaseScale=1;
      dPhaseScale=atof(optarg);
      break;

    case 'O':
      if(optarg){
	if(*optarg=='D'){
	  iDoODAnalysis=1;
	  if(*(optarg+1)) iODAnalysisNBins=atoi(optarg+1);
	  else iODAnalysisNBins=DISTRIBUTION_PHI_OD_BINS_N;
	}
	else dsradius=atof(optarg);
      }
      else dsradius=DEFAULT_SMOOTH_RADIUS;
      break;
      
    case 'P':
      for(i=0,pc=optarg,pc1=pc2=NULL;*pc;pc++){
	if(*pc==','){
	  i++;
	  if(i==1){
	    *pc='\0';
	    pc1=pc+1;
	  }
	  if(i==2){ 
	    *pc='\0';
	    pc2=pc+1;
	  }
	  if(i==3){
	    printf("GET Period: ignoring trailing characters |%s|\n",pc);
	    *pc='\0';
	    break;
	  }
	}
      }
      //printf("GET |%s|%s|%s|\n",optarg,pc1,pc2);
      if(optarg){
	dStimPeriod=atof(optarg);
	if(dStimPeriod<=0.0){
	  printf("GET ERROR dStimPeriod<=0.0\n");
	  break;
	}
	phase_transform=dStimPeriod*harmonic_d/(2.0*M_PI);
	do_phase_transform=1;
      }
      if(pc1){
	iPrePostStimTimeRescale=1;
	dPreStimTime=atof(pc1)/dStimPeriod;
	if(pc2){
	  dPostStimTime=atof(pc2)/dStimPeriod;
	}
	else dPostStimTime=dPreStimTime;
	if(dPreStimTime+dPostStimTime>=1.0){
	  printf("GET WARNING PreStimTime+PostStimTime>=dStimPeriod\n");
	}
      }
      //	printf("GET |%f|%f|%f|\n",dStimPeriod,dPreStimTime,dPostStimTime);
      break;

    case 'R':
      iRotateMaps=atoi(optarg);
      iRotateMaps=iRotateMaps%4;
      if(iRotateMaps<0) iRotateMaps+=4;
      printf("GET ROTATION %i\n",iRotateMaps);
      break;

    case 'S':
      save_mode=atoi(optarg);
      if(save_mode>DISPLAY_MODE6 || save_mode<DISPLAY_MODE2){
	printf("Bad save mode (%i), using default (%i)\n",save_mode,DEFAULT_SAVE_MODE);
	save_mode=DEFAULT_SAVE_MODE;
      }
      break;

    case 'T':
      dradius=5.0;
      dsradius=5.0;
      amplitude_chop_l=0.5;
      depth_lines=3;
      smartcontours_threshold=0.5;
      break;

    case 'V':
      iFullSizeViewPort=1;
      break;

    case 'W':
      for(i=0,j=0;i<strlen(optarg);i++) if(optarg[i]==',') j++;
      if(j==2){
	pc=optarg;
	pc1=index(pc,',');
	*pc1='\0';
	fStim=atof(pc);
	pc=pc1+1;
	pc1=index(pc,',');
	*pc1='\0';
	fPreBlank=atof(pc);
	pc=pc1+1;
	fPostBlank=atof(pc);
	if(fStim<0.0){
	  printf("GET Stim time must be > 0.0\n");
	  iCutBlanksOnWedge=0;
	  break;
	}
	if(fPreBlank<0.0){
	  printf("GET Pre Blank time must be >= 0.0\n");
	  iCutBlanksOnWedge=0;
	  break;
	}
	if(fPostBlank<0.0){
	  printf("GET Post Blank time must be >= 0.0\n");
	  iCutBlanksOnWedge=0;
	  break;
	}	
	iCutBlanksOnWedge=1;
      }
      else{
	printf("GET Incorrect cut wedge format |%s|\n",optarg);
	printf("GET Correct syntax: Stim,Pre Blank,Post Blank(eg. 4,0,2)\n");
      }
      break;

    case 'Z':
      do_pixel_size_transform=1;
      if(optarg) dPixelSizeU=atof(optarg);
      else dPixelSizeU=DEFAULT_PIXEL_SIZE_MICRON;
      break;
      
    case 'a':
      iDoWedgeShift=1;
      break;

    case 'b':
      spatial_binning=atoi(optarg);
      if(spatial_binning>MAX_SPATIAL_BINNING) spatial_binning=MAX_SPATIAL_BINNING;
      if(spatial_binning<1) spatial_binning=1;
      break;

    case 'c':
      contour_scheme=atoi(optarg);
      if(contour_scheme>10 || contour_scheme<0 || contour_scheme==3 || contour_scheme==7){
	printf("Bad contour scheme %i, using default %i\n",contour_scheme,CONTOUR_DEFAULT);
	contour_scheme=CONTOUR_DEFAULT;
      }
      contour_scheme_phi=contour_scheme & 3;
      contour_scheme_rho=contour_scheme >> 2;
      break;
     
    case 'd':
      depth_lines=atoi(optarg);
      if(depth_lines<0){
	depth_lines=0;
      }
      
      for(pc1=NULL,pc=optarg;*pc;pc++){
	if(*pc==','){
	  *pc='\0';
	  pc1=pc+1;
	}
      }
      depth_lines=atof(optarg);
      if(depth_lines<0){
	depth_lines=0;
      }
      if(pc1) fContourSpacingAngle=atof(pc1);
      break;

    case 'e':
      iBlendBorder=1;
      if(optarg){
	iBlendBorderPixel=atoi(optarg);
	if(iBlendBorderPixel<=0) iBlendBorderPixel=DEFAULT_BLEND_BORDER;
      }
      break;

    case 'f':
      filenames[nfiles]=strdup(optarg);
      nfiles++;
      if(nfiles>MAX_FILE_N){
	printf("Accepts two(X map first, Y map second) or one(XY in one file) files only\n");
	exit(1);
      }
      break;
 
    case 'g':
      iDoAmplitudeNorm=1;
      if(optarg){
	dAmplitudeNormCoeff=atof(optarg);
      }
      else dAmplitudeNormCoeff=-1.0;
      break;

    case 'h':
      printf(" -A <float(min),float(max),[char*(nits)]> - do wedge annotation (default off)\n");
      printf(" -B - load blank image\n");
      printf(" -C [%i,%i] -  color scheme (default %i)\n",COLOR_RGB,COLOR_HLS,COLOR_DEFAULT);
      printf(" -D [%i,%i,%i] - display mode (default %i)\n",DISPLAY_MODE2,DISPLAY_MODE4,DISPLAY_MODE6,DEFAULT_DISPLAY_MODE);
      printf(" -E <float> - sectioning smooth radius (default %f)\n",DEFAULT_SECTION_RADIUS);
      printf(" -F [%i,%i,%i] - field display mode (default off)\n",DISPLAY_MODE2,DISPLAY_MODE4,DISPLAY_MODE6);
      printf(" -G <float,float> - fixed gray-scale wedge (Low,High)\n");
      printf(" -H <float> - high-pass filter domain (default %f)\n",DEFAULT_RADIUSBIG);
      printf(" -I <float> - contours initial angle in degrees (default %f)\n",DEFAULT_CONTOURS_INI_ANGLE_RAD/M_PI*180.0);
      printf(" -J <int> - inverse maps X %i Y %i (default off)\n",INVERSION_X,INVERSION_Y);
      printf(" -L <float> - low-pass filter domain (default %f)\n",DEFAULT_RADIUS);
      printf(" -M <float> - phase scale factor (default none=1)\n");
      printf(" -O <float> - smooth contours radius (default %f)\n",DEFAULT_SMOOTH_RADIUS);
      printf(" -OD <int> - number of bins for OD analysis (default %i)\n",DISTRIBUTION_PHI_OD_BINS_N);
      printf(" -P <float,[float,float]> - period (deg) for phase transformation (default off)\n");
      printf(" -R <int> - rotate maps by <int>*90deg (default off)\n");
      printf(" -S [%i,%i,%i] - save mode (default %i)\n",DISPLAY_MODE2,DISPLAY_MODE4,DISPLAY_MODE6,DEFAULT_SAVE_MODE);
      printf(" -W <float,float,float> - cut wedge (Stim,Pre Blank,Post Blank)\n");
      printf(" -Z [float] pixel size in microns (default transform off, default value %f)\n",DEFAULT_PIXEL_SIZE_MICRON);
      printf(" -a do wedge shift (used for auditory) (default off)\n");
      printf(" -b <int> spatial binning (default 1)\n");
      printf(" -c <int> countour scheme (default %i)\n",CONTOUR_DEFAULT);
      printf(" -d <int[,float]> number of countour lines and contour spacing (default 0)\n");
      printf(" -e blend border (default off, default border %i)\n",DEFAULT_BLEND_BORDER);
      printf(" -f <filename> - needs two (X and Y components) or one (XY in one file)files \n");
      printf(" -g [float] - normalize map by [float] (default off, default0 1)\n");
      printf(" -h - print this help\n");
      printf(" -i interactive image manipulation (default off)\n");
      printf(" -l[radius] - lens correction (default off)\n");
      printf(" -m - mark min/max sites (default off)\n");
      printf(" -o <float> - chop radius (default %f)\n",DEFAULT_CHOP_RADIUS);
      printf(" -p - do pinwheel analysis (default off)\n");
      printf(" -r - remove bias (default off)\n");
      printf(" -s - save response dump (default off)\n");
      printf(" -t - do discontinuous phase (default off)\n");
      printf(" -u <string> - attach suffix <string> to saved file (default _a)\n");
      printf(" -v - verbose mode (default off)\n");
      printf(" -w <float[,float]> - window size [, interactive window size] (default %f %f)\n",WINDOW_SIZE,INTERACTIVE_WINDOW_SIZEX);
      printf(" -z - draw zero lines (default off)\n");
      exit(0);
      break;

    case 'i':
      int_mode=INT_MODE_DEFAULT;
      break;
     
    case 'l':
      iLensCorrection=1;
      if(optarg) dLensRadius=atof(optarg);
      break;
     
    case 'm':
      iMarkMinMax=1;
      break;
     
    case 'o':
      dcradius=atof(optarg);
      break;
      
    case 'p':
      do_pinwheels=1;
      if(optarg) iSaverPinwheels=1;
      break;
     
    case 'r':
      removebias=1;
      break;
     
    case 's':
      do_dump=1;
      break;
     
    case 't':
      do_discontinuous_phase=1;
      break;

    case 'u':
      SAFE_FREE(pcSaveFileSuffix);
      if(optarg) pcSaveFileSuffix=strdup(optarg);
      break;

    case 'v':
      do_verbose=1;
      break;
     
    case 'w':
      if(optarg){
	pc1=optarg;
	pc2=NULL;
	for(pc=optarg;*pc;pc++){
	  if(*pc==','){
	    *pc='\0';
	    pc2=pc+1;
	    break;
	  }
	}
	window_size=atof(pc1);
	if(pc2) fInteractiveWindowWidth=atof(pc2);
      }
      break;
     
    case 'z':
      do_zerolines=1;
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
	printf("Accepts two(X map first, Y map second) or one(XY in one file) files only\n");
	exit(1);
      }
    }
  }  
}
