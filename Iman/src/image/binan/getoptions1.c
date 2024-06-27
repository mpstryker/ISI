/* Bin analysis                 */
/* Written by V.Kalatsky        */
/* 20Jul2001                    */

#include "binan1.h"

void options(int argc,char *argv[]){
  int c;
  
  while (1){
    
    if((c=getopt(argc, argv, "AFH:L:M::a:b:himo:rt:vw:z"))==-1) break;

    switch (c){

    case 'A':
      shift_average_to_zero=1;
      break;

    case 'F':
      do_fixed_background=1;
      break;

    case 'H':
      dradiusbig=atof(optarg);
      break;

    case 'L':
      dradius=atof(optarg);
      break;

    case 'M':
      iMakeMaps=1;
      if(optarg) iNFComponents=atoi(optarg);
      break;

    case 'a':
      animation_delay=(unsigned long)atoi(optarg);
      do_animation=1;
      break;

    case 'b':
      spatial_binning=atoi(optarg);
      if(spatial_binning>MAX_SPATIAL_BINNING) spatial_binning=MAX_SPATIAL_BINNING;
      if(spatial_binning<1) spatial_binning=1;
      break;

    case 'h':
      printf(" -A - shift image average to zero (default off)\n");
      printf(" -F - fixed background (default floating)\n");
      printf(" -H <float> - high-pass filter domain (default %f)\n",DEFAULT_RADIUSBIG);
      printf(" -L <float> - low-pass filter domain (default %f)\n",DEFAULT_RADIUS);
      printf(" -M <int> - make maps using <int> Fourier components (default off)\n");
      printf(" -a <int> animation delay (default %lu)\n",DEFAULT_ANIMATION_DELAY);
      printf(" -b <int> spatial binning (default 1, max %i)\n",MAX_SPATIAL_BINNING);
      printf(" -h - print this help\n");
      printf(" -i interactive image manipulation (default off)\n");
      printf(" -m - mark max site (default off)\n");
      printf(" -o <float> - chop radius (default %f)\n",DEFAULT_CHOP_RADIUS);
      printf(" -r - remove bias (default off)\n");
      printf(" -t - bin size in seconds (default %0.1f)\n",DEFAULT_BIN_SIZE_SEC);
      printf(" -v - verbose mode (default off)\n");
      printf(" -w <float> - window size (default %f)\n",WINDOW_SIZE);
      printf(" -z - shift maxima to zero (default off)\n");
      exit(0);
      break;

    case 'i':
      int_mode=INT_MODE_DEFAULT;
      break;
     
    case 'm':
      do_marking=1;
      break;
     
    case 'o':
      dcradius=atof(optarg);
      break;
      
    case 'r':
      removebias=1;
      break;
     
    case 't':
      fBinSizeSec=atof(optarg);
      break;
      
    case 'v':
      do_verbose=1;
      break;
     
    case 'w':
      window_size=atof(optarg);
      break;
     
    case 'z':
      shift_max_to_zero=1;
      break;

    default:
      printf ("?? getopt returned character code 0%o ??\n", c);
      exit(1);
    }
  }
  
  if(optind < argc){
    while (optind < argc){
      filenames=(char**)realloc(filenames,(nfiles+1)*sizeof(char*));
      filenames[nfiles]=strdup(argv[optind++]);
      nfiles++;
      /*
      if(nfiles>MAX_FILE_N){
	printf("Accepts %i file(s) only\n",MAX_FILE_N);
	exit(1);
      }
      */
    }
  }  
}
