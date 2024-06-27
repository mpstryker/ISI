#include "iman1.h"

void options(int argc,char *argv[]){
  int c;
  char *pc;
  
  while (1){
    
    if((c=getopt(argc, argv, "1A:B:C:DE:F:H:I:L:N::PR::QS:T:Uabc::f:hil:kmp:qr::s::tu:vw:"))==-1) break;

    switch (c){

    case '1':
      do_remove_linear=1;
      break;
    case 'A':
      additionalframes=atoi(optarg);
      iDrawAxesAndLabelOnGreen=0;
      break;
    case 'B':
      initframe=atoi(optarg);
      break;
    case 'C':
      n_cycles=atoi(optarg);
      break;
    case 'D':
      iDoDivisionByAverage=1;
      break;
    case 'E':
      finalframe=atoi(optarg);
      break;
    case 'F':
      harmonic=atof(optarg);
      break;
    case 'H':
      dradiusbig=atof(optarg);
      break;
    case 'I':
      stimulus_period=atof(optarg);
      do_intervals|=1;
      ulSynchType|=SYNCHRONIZATION_QUASI_TIME;
      break;
    case 'L':
      dradius=atof(optarg);
      break;
    case 'N':
      if(optarg) n_bins=atoi(optarg);
      else n_bins=1;
      break;
    case 'P':
      do_precise=1;
      break;
    case 'Q':
      compress=DECOMPRESS;
      break;
    case 'R':
      if(optarg){
	frame_rule=atoi(optarg);
	if(frame_rule>=FRAME_RULE_N){
	  printf("GET Unknown Frame Rule %i. Reverting to default %i\n",frame_rule,FRAME_RULE_DEFAULT);
	  frame_rule=FRAME_RULE_DEFAULT;
	}
      }
      else frame_rule=FRAME_RULE_DEFAULT;
      break;
    case 'S':
      do_statistics=1;
      statistics_num_bins=atoi(optarg);
      break;
    case 'T':
      //      timeradius=atoi(optarg);
      do_timeaveraging=1;
      timeaverage_blocks=atof(optarg);
      //      timeaverage_cycles=ceil(timeaverage_blocks);
      timeaverage_cycles=timeaverage_blocks;
      timeaverage_cycles_int=(int)(timeaverage_cycles);
      break;
    case 'U':
      sacrifice_boundaries=0;
      break;
    case 'a':
      average=AVERAGE;
      break;
    case 'b':
      find_best_cycle=1;
      break;
    case 'c':
      do_contours=1;
      if(optarg) fContrast=atof(optarg);
      break;
    case 'f':
      nfiles++;
      if(nfiles>MAX_FILE_N){
	printf("GET Accepts %i file(s) only\n",MAX_FILE_N);
	exit(1);
      }
      filenames[nfiles-1]=strdup(optarg);
      break;
 
    case 'h':
      printf(" -1 - do linear fit, (default off, -r and -R supersede this option)\n");
      printf(" -A <int> - additional frames (default 0)\n");
      printf(" -B <int> - initial frame (default EPST=0 COST=%i)\n",DEFAULT_INITFRAME);
      printf(" -C <int> - number of cycles (default %i)\n",NCYCLES);
      printf(" -D - divide by average image(evens illumination non-uniformities) (default off)\n");
      printf(" -E <int> - final frame (default %i+%i*NFRAMES/CYCLE)\n",DEFAULT_INITFRAME,NCYCLES);
      printf(" -F <float> - harmonic (default %.2f)\n",HARMONIC);
      printf(" -H <float> - high-pass filter domain (default %.1f) (size of green)\n",DRADIUSBIG);
      printf(" -I <float> - period of stimulus in seconds, uses interframe times (default off)\n");
      printf(" -L <float> - low-pass filter domain (default %.1f)\n",DRADIUS);
      printf(" -N [int(default 1)] - number of bins for temporal frame binning, requires -T option (default off)\n");
      printf(" -P - do precise analysis (default off)\n");
      printf(" -Q - decompress file\n");
      printf(" -R [int(default %i)] - rule number for binning\n",FRAME_RULE_DEFAULT);
      printf(" -S <int> - number of bins, do statistics (default off)\n");
      printf(" -T <float> - 1/2 of time averaging window (cycles/condition blocks) (default off, trace default %i frames)\n",TIMERADIUS);
      printf(" -U - sacrifice boundaries off, works with -T (default on) (format green output)\n");
      printf(" -a - do averaging (default off)\n");
      printf(" -b - find best cycle (default off)\n");
      printf(" -c - draw contour lines (default off)\n");
      printf(" -f <filename> - needs %i files\n",MAX_FILE_N);
      printf(" -h - print this help\n");
      printf(" -i - show interframe timing (default off)\n");
      printf(" -k - keep original (compressed/decompressed) files (default %i)\n",DEFAULT_KEEP_ORIGINAL_FILES);
      printf(" -l <int> - Fourie series lower cut-off (default %i)\n",CUT_OFF_L);
      printf(" -m - multiple frame mode (default single)\n");
      printf(" -p <float> - period of stimulus (frames) (default frames/cycle, can be bogus)\n");
      printf(" -q - compress file\n");
      printf(" -r [int(default Ncycles*%.2f+%i)] - degree of fitting polynomial,  (default off)\n",POLYNOM_COEFF,POLYNOM_ADD);
      printf(" -s [int(default 0)] - use synchronization channel (default off)\n");
      printf(" -t - do trace analysis (default off)\n");
      printf(" -u <int> - Fourie series upper cut-off (default %i)\n",CUT_OFF_U);
      printf(" -v - verbose mode, prints out a lot of crap (default off)\n");
      printf(" -w<float>/-w<float,float> - size of green save image/interactive,plot trace windows (default %f/%f,%f)\n",GREEN_IMAGE_SIZE_SAVE,TRACE_INTERACTIVE_WINDOW_SIZEX,TRACE_PLOT_WINDOW_SIZEX);
      exit(0);
      break;

    case 'i':
      do_intervals|=2;
      break;

    case 'k':
      iKeepOriginalFiles=!DEFAULT_KEEP_ORIGINAL_FILES;
      break;

    case 'l':
      cut_offl=atoi(optarg);
      break;

    case 'm':
      single_multiple=DO_MULTIPLE;
      break;

    case 'p':
      period=atof(optarg);
      ulSynchType|=SYNCHRONIZATION_QUASI_USER;
      break;

    case 'q':
      compress=COMPRESS;
      break;

    case 'r':
      do_remove_curve=1;
      if(optarg) poly_fitn=atoi(optarg);
      break;

    case 's':
      do_synch=1;
      if(optarg) iSynchChannel=atoi(optarg);
      else iSynchChannel=0;
      ulSynchType|=SYNCHRONIZATION_TRUE_EXTERNAL;
      break;

    case 't':
      do_trace=1;
      break;

    case 'u':
      cut_offu=atoi(optarg);
      break;

    case 'v':
      do_verbose=1;
      break;

    case 'w':
      if(optarg){
	if((pc=strtok(optarg,","))){
	  if(*pc){
	    fSizeOfGreenSave=fSizeOfTraceInteractiveWindow=atof(pc);    
	    if((pc=strtok(NULL,","))){
	      fSizeOfTracePlotWindow=atof(pc);
	    }
	  }
	}
      }
      break;

    default:
      printf ("GET getopt returned character code 0%o\n", c);
      exit(1);
    }
  }
  
  if(optind < argc){
    while (optind < argc){
      nfiles++;
      if(nfiles>MAX_FILE_N){
	printf("GET Accepts %i file(s) only\n",MAX_FILE_N);
	exit(1);
      }
      filenames[nfiles-1]=strdup(argv[optind++]);
    }
  }  
}
