// fileheader.h: collection of headers and chunks
/* Image analysis               */
/* Written by V.Kalatsky        */
/* 16Mar02                      */

#include <stdio.h>
#include <string.h> 
#include <errno.h>

#ifndef _BYTE
#define _BYTE
typedef unsigned char BYTE;
#endif
#ifndef _BOOL
#define _BOOL
typedef int BOOL;
#endif
#ifndef _USHORT
#define _USHORT
typedef unsigned short USHORT;
#endif
#ifndef _ULONG
#define _ULONG
typedef unsigned long ULONG;
#endif

typedef struct FILEHEADER {
  char		sumTag[4] ;	/* S?01, ? = F,S, M, R for expt location for sum_file Q for sumsq, R for ratio*/
  long int	sumBytes ;	/* byte count of bytes to follow in file */
  short int	sumXsize ;   // set from xsu.xstart,xstop
  short int	sumYsize ;   // set from xsu.ystart,ystop
  short int	sumNframes ; // set from xsu.sequentil
  short int	sumNconds ;  // <-----must be set for expt to run
  short int	sumCondNo	;	// six entries set for compatibility with IDL version
  short int	sumOrder ;
  float 	sumLoClip ;
  float		sumHiClip ;
  short int	sumLoPass ;
  short int	sumHiPass ;   
  long		sumLDataType ;
  long		sumLFileType ;
  long		sumLSizeOf ;
  char		sumFree1[2] ;
  short int	sumHeadersize ;
  short int	sumNreps ;    // <-----must be set for expt to run
  short int	sumRandomize ; 
  short int	sumXgroup ;	/* binning done on CCD chip */ // set from xsu.xgroup	// GB: Will be removed !
  short int	sumYgroup ;   // set from xsu.xgroup								// GB: Will be removed !
  short int	sumGoDelay ;
  short int	sumRatioFirstFrame ; 
  short int	sumRatioLastFrame ;
  short int	sumTempBinning ;   /* set if we collect in full resolution without hardware binning and want to bin */
  long int	sumBit1 ;	/* Offset 64, bit coded, which blocks */ 
  long int	sumBit2 ;   /* bit coded, which frames */
  short int	sumBinning ;/* binning done on CCD */  // now this is being used for software spatila binning
  short int	sumExposureMsec ;  // set from xsu.exposure (double)
  float		sumMagnification ;	/* optical magnification */
  short int	sumGain ;          // set from xsu.gain
  short int	sumWavelength ;   
  long		sumFrameHeaderSize ;  // used only for stream files, size of the header for each frame
  short int	sumNframeiti ;   
  short int	sumNframestim ;   
  char		sumFree2[20] ;
  long int	sumSysTime ;		// GB: Will be removed !
  char		sumSite[4] ;        // GB: Will be removed !
  short int	sumWinType ;        // GB: Will be removed !
  short int	sumDay ;            // GB: Will be removed !
  short int	sumMonth ;          // GB: Will be removed !
  short int	sumYear ;           // GB: Will be removed !
  char		sumDateRecorded[6] ;	/* Offset 128, ASCII */ 
  char		sumUserName[10] ;		/* ASCII */ 
  char		sumDateCreated[6] ;	/*  ASCII */ 
  char		sumOrigFilename[16] ;	
/* encoded as Sddmy.ea, Gddmy.ea, or Rddmy.ea, where S=sum, R=red, G=green
date of experiment start: dd = day, m=month in hex, y = last digit of year, 
e = experiment number (0, 1, 2, ...), a = run number (a, b, c, d, ... ) */
  char		sumPrevFilename[16];
  char		sumNextFilename[16];
  unsigned short sumCompressedRecordSize;
  unsigned long sumCompressedFrameSize;
  unsigned long sumCompressedFrameNumber;
  unsigned short sumCompressedDumbAss; //To presserve int borders
  char	        sumFree3[46] ; // Val: used 44 bytes out of 90 in sumFree3 for PrevFilename, NextFilename, sumCompressedFrameSize, sumCompressedRecordSize and sumCompressedFrameNumber
  char	        sumComments[256] ;		/* Offset 256 */
  /* DATA (unlimited) starts at offset 512, or 508 bytes after tag & bytecount */

} FILEHEADER ;

typedef struct I64{
  unsigned long LowPart; 
  unsigned long HighPart;
} I64;

typedef struct FRAMEHEADER {
  unsigned long frameheadID; // arbitrary for future flexibility
  unsigned long frameheadlength;  //8
  unsigned long seqnum;  //sequence number 0,1,2,...NOT including paused frames
  unsigned long time; //16  from the time call
  I64 perfcount;  // from the performance counter
  I64 perfreq; //32
  long rep; //repetition number of this cycle
  long trial;  // number of stim in this cycle
  long condition;  //stimulus number
  long frame_of_cond; //48  which frame in relation to start of stimulus
  long frame_paused; //  binary variable, marks paused frames(=1[actually !=0])
  long frame_type; //56  frame type, iti and stim for now
  long seqnum_all; //sequence number 0,1,2,...INCLUDING PAUSED FRAMES
  long synch_in; //for synchronization input signals (on port B)
  long stim_params[48];  // extra parameters
} FRAMEHEADER; //256 bytes total


//Val: all records are multiples of 4

#define CHUNK_ID_SIZE 4
#define CHUNK_SIZE_SIZE sizeof(ULONG)
#define CHUNK_HEAD_SIZE (CHUNK_ID_SIZE+CHUNK_SIZE_SIZE)
#define CHUNK_TAG_SIZE 4

typedef struct FRAM_CHUNK{		// Experiment independent chunk sizeof(FRAM_CHUNK)=64
  char		ID[4];			// "FRAM"=FRAMe header
  ULONG		Size;			// size of the frameheader to follow sizeof(FRAM_CHUNK)-8 = 56
  char		Tag[4];			//01 01
  ULONG		FrameSeqNumber;	        //02 02 Frame sequential number as reported by frame grabber
  ULONG		FrameRingNumber;	//03 03 Frame ring number as reported by frame grabber
  ULONG		TimeArrivalUsecLo;	//04 04 Arrival time(microsec) as reported by frame grabber, low part
  ULONG		TimeArrivalUsecHi;	//05 05 Arrival time(microsec) as reported by frame grabber, high part
  ULONG		TimeDelayUsec;		//06 06 Difference between arrival and WaitEvent time
  ULONG		PotentiallyBad;		//07 07 As returned by GrabWaitFrameEx
  ULONG		Locked;			//08 08 As returned by GrabWaitFrameEx
  ULONG		FrameSeqCount;		//09 09 Frame sequential number (should be: FrameSeqCount+1=FrameSeqNumber)
  ULONG		CallbackResult;		//10 10 Return value of the experiment callback
  ULONG		Free[4];		//11 14
} FRAM_CHUNK;

typedef struct FRAM_COST_CHUNK{
  char		ID[4];			// "cost" COST experiment frame header. sizeof(FRAM_COST_CHUNK)=64
  ULONG		Size;			// size of the frameheader and frame data to follow(sizeof(FRAM_COST_CHUNK)-8+frameSize)
  char		Tag[4];			//01 01
  ULONG		HeaderSize;		//02 02 size of the frameheader to follow(sizeof(FRAM_COST_CHUNK)-16)
  ULONG		SynchChannel[4];	//02 06 Up to 4 channels. Increase the arrays size if needed or use SYNC_CHUNK
  ULONG		SynchChannelDelay[4];	//06 10 Relative to arrival time (microseconds)
  ULONG		Free[4];		//10 14
} FRAM_COST_CHUNK;

typedef struct FRAM_COST_HEADER{	// 128Bytes
  FRAM_CHUNK		framChunk;
  FRAM_COST_CHUNK	framCOSTChunk;
} FRAM_COST_HEADER;

// EPST frame types
#define EPST_FRAME_ITI          0
#define EPST_FRAME_STIM         1
#define EPST_FRAME_PRE_STIM     2
#define EPST_FRAME_POST_STIM    3
#define EPST_FRAME_PAUSE        4

typedef struct FRAM_EPST_CHUNK{
  char		ID[4];			// "epst" EPST experiment frame header. sizeof(EPST_FRAM_CHUNK)=64
  ULONG		Size;			// size of the frameheader and frame data to follow(sizeof(FRAM_EPST_CHUNK)-8+frameSize)
  char		Tag[4];			//01 01
  ULONG		HeaderSize;		//02 02 size of the frameheader to follow(sizeof(FRAM_EPST_CHUNK)-16)
  ULONG		SeqNumber;		//03 03 sequence number 0,1,2,...NOT including paused frames
  ULONG		Repetition;		//04 04 repetition number of this cycle
  ULONG		Trial;			//05 05 number of stim in this cycle
  ULONG		Condition;		//06 06 stimulus number
  ULONG		FrameOfCondition;	//07 07 which frame in relation to start of stimulus
  ULONG		FramePaused;		//08 08 BOOL variable, marks paused frames(=1[actually !=0])
  ULONG		FrameType;		//09 09 frame type, iti and stim for now
  ULONG		Free[5];		//10 14
} FRAM_EPST_CHUNK;

typedef struct FRAM_EPST_HEADER{	// 128Bytes
  FRAM_CHUNK		framChunk;
  FRAM_EPST_CHUNK	framEPSTChunk;
} FRAM_EPST_HEADER;


typedef struct ISOI_CHUNK{
  char		ID[4];			// "ISOI"=Intrinsic Signal Optical Imaging
  ULONG		Size;			// size of the file to follow, 0 if unknown
  char		Tag[4];			// "COIM"=COntinuous IMaging 
} ISOI_CHUNK;

typedef struct HEAD_CHUNK{
  char		ID[4];			// "HEAD"
  ULONG		Size;			// Size of the chunk (HEAD) to follow. Fixed to 512 bytes
  char		Tag[4];			//01 01
  char		DateTimeRecorded[24];	//02 07 ASCII UNIX time
  char		UserName[16];		//03 11 ASCII 
  char		SubjectID[16];		//04 15 ASCII 
  char		ThisFilename[16];	//05 19 Encoded in MakeFileName
  char		PrevFilename[16];	//06 23
  char		NextFilename[16];	//07 27
  ULONG		DataType;		//08 28
  ULONG		FileType;		//09 29
  ULONG		SizeOfDataType;		//10 30
  ULONG		XSize;			//11 31 X size of stored frames 
  ULONG		YSize;			//12 32 Y size of stored frames
  ULONG		XCCDSize;		//13 33 X size of HW binned and not SW binned frame or CCD chip pixel size 
  ULONG		YCCDSize;		//14 34 Y size of HW binned and not SW binned frame or CCD chip pixel size 
  ULONG		ROIXPosition;		//15 35 X coordinate of upper-left conner of ROI (before binning)
  ULONG		ROIYPosition;		//16 36 Y coordinate of upper-left conner of ROI (before binning)
  ULONG		ROIXSize;		//17 37 X size of ROI (before binning) 
  ULONG		ROIYSize;		//18 38 Y size of ROI (before binning) 
  ULONG		ROINumber;		//19 39 Sequential ROI number. Set to 0 for main ROI
  ULONG		TemporalBinning;	//20 40
  ULONG		SpatialBinning;		//21 41 union with two USHORTs, first = X, second = Y
  ULONG		NConditions;		//22 42 Number of conditions (for episodic stimuli)
  ULONG		NRepetitions;		//23 43 Number of repetitions (for episodic stimuli)
  ULONG		Randomize;		//24 44 (for episodic stimuli)
  ULONG		IntegrationTimeUsec;	//25 45
  ULONG		InterFrameTimeUsec;	//26 46
  ULONG		HardwareBinning;	//27 47 union with two USHORTs, first = X(V), second = Y(H)
  ULONG		HardwareGain;		//28 48
  long		HardwareOffset;		//29 49
  ULONG		OpticalMagnification;	//30 50 union with two USHORTs, first=top lens focal lenght, second = bottom
  ULONG		WaveLengthAndFilterWidth;//31 51 union with two USHORTs, first = wave lenght, second = filter width
  ULONG		FrameHeaderSize;	//32 52 Size of the header for each frame
  ULONG		NFramesTotal;		//33 53 Expected number of frames 
  ULONG		NFrames;		//34 54 Number of frames in this file
  ULONG		NFramesITI;		//35 55 (for episodic stimuli)
  ULONG		NFRramesStim;		//36 56 (for episodic stimuli)
  ULONG		NSynchChannels;		//37 57 Up to 4 channels. Increase size of SynchChannelMax if needed or use SYNC_CHUNK
  ULONG		SynchChannelMax[4];	//38 61
  char		Free[12];		//42 64
  char		Comments[256];		//43 128
} HEAD_CHUNK;

typedef struct HARD_CHUNK{
  char		ID[4];			// "HARD" - hardware chunk
  ULONG		Size;			// size of the chunk (HARD) to follow = 256
  char		Tag[4];			//01 01
  char		CameraName[16];		//02 05 
  ULONG		CameraType;		//03 06 
  ULONG		ResolutionX;		//04 07 Physical X resolution of CCD chip  
  ULONG		ResolutionY;		//05 08 Physical Y resolution of CCD chip 
  ULONG		PixelSizeX;		//06 09 micrometers
  ULONG		PixelSizeY;		//07 10 micrometers
  ULONG		CCDApertureX;		//08 11 micrometers
  ULONG		CCDApertureY;		//09 12 micrometers
  ULONG		IntegrationTime;	//10 13 microseconds
  ULONG		InterFrameTime;		//11 14 microseconds
  ULONG		HardwareBinningV;	//12 15 Vertical   or X binning
  ULONG		HardwareBinningH;	//13 16 Horizontal or Y binning
  ULONG		HardwareGain;		//14 17
  long		HardwareOffset;		//15 18
  ULONG		CCDSizeX;		//16 19 Hardware binned CCD X size 
  ULONG		CCDSizeY;		//17 20 Hardware binned CCD Y size 
  ULONG		DynamicRange;		//18 21
  ULONG		OpticsFocalLengthTop;	//19 22 Top lens focal lenght (the one closer to the camera), millimeters
  ULONG		OpticsFocalLengthBottom;//20 23 Bottom lens focal lenght, millimeters
  char		Comments[164];		//21 64
} HARD_CHUNK;

typedef struct SOFT_CHUNK{
  char		ID[4];			// "SOFT" - software chunk
  ULONG		Size;			// size of the chunk (SOFT) to follow = 256 
  char		Tag[4];			//01 01
  char		DateTimeRecorded[24];	//02 07 ASCII UNIX time
  char		UserName[16];		//03 11 ASCII 
  char		SubjectID[16];		//04 15 ASCII 
  char		ThisFilename[16];	//05 19 Encoded in MakeFileName
  char		PrevFilename[16];	//06 23
  char		NextFilename[16];	//07 27
  ULONG		DataType;		//08 28
  ULONG		FileType;		//09 29
  ULONG		SizeOfDataType;		//10 30
  ULONG		XSize;			//11 31 X size of stored frames 
  ULONG		YSize;			//12 32 Y size of stored frames
  ULONG		ROIXPosition;		//13 33 X coordinate of upper-left conner of ROI (before binning)
  ULONG		ROIYPosition;		//14 34 Y coordinate of upper-left conner of ROI (before binning)
  ULONG		ROIXSize;		//15 35 X size of ROI (before binning) 
  ULONG		ROIYSize;		//16 36 Y size of ROI (before binning) 
  ULONG		ROIXPositionAdjusted;	//17 37 X coordinate of upper-left conner of adjusted ROI (before binning)
  ULONG		ROIYPositionAdjusted;	//18 38 Y coordinate of upper-left conner of adjusted ROI (before binning)
  ULONG		ROINumber;		//19 39 Sequential ROI number. Set to 0 for main ROI
  ULONG		TemporalBinning;	//20 40
  ULONG		SpatialBinningX;	//21 41 X
  ULONG		SpatialBinningY;	//22 42 Y
  ULONG		FrameHeaderSize;	//23 43 Size of the header for each frame
  ULONG		NFramesTotal;		//24 44 Expected number of frames 
  ULONG		NFramesThisFile;	//25 45 Number of frames in this file
  ULONG		WaveLength;		//26 46 Wave lenght, nanometers
  ULONG		FilterWidth;		//27 47 Filter width, nanometers
  char		Comments[68];		//28 64
} SOFT_CHUNK;

typedef struct COST_CHUNK{
  char		ID[4];			// "COST" - Continuous STimulation chunk
  ULONG		Size;			// size of the chunk (COST) to follow = 64
  char		Tag[4];			//01 01
  ULONG		NSynchChannels;		//02 02 Up to 4 channels. Add new values bellow if needed or use SYNC_CHUNK
  ULONG		SynchChannelMax[4];	//03 06
  ULONG		NStimulusChanels;	//07 07 Used if no Sync Channels present
  ULONG		StimulusPeriod[4];	//08 11 milliseconds
  char		Comments[20];		//12 16
} COST_CHUNK;

typedef struct EPST_CHUNK{
  char		ID[4];			// "EPST" - EPisodic STimulation chunk
  ULONG		Size;			// size of the chunk (EPST) to follow = 64
  char		Tag[4];			//01 01
  ULONG		NConditions;		//02 02 Number of conditions (for episodic stimuli)
  ULONG		NRepetitions;		//03 03 Number of repetitions (for episodic stimuli)
  ULONG		Randomize;		//04 04 (for episodic stimuli)
  ULONG		NFramesITI;		//05 05 (for episodic stimuli)
  ULONG		NFramesStim;		//06 06 (for episodic stimuli)
  ULONG         NFramesBlankPre;        //07 07 Pre-stimulus frames (for episodic stimuli)
  ULONG         NFramesBlankPost;       //08 08 Post-stimulus frames (for episodic stimuli)
  char          Comments[32];           //09 16 
} EPST_CHUNK;

typedef struct GREE_CHUNK{
  char		ID[4];			// "GREE" - GREEn chunk
  ULONG		Size;			// size of the chunk (GREE) to follow = 32
  char		Tag[4];			//01 01
  float 	LoClip;			//02 02
  float		HiClip;			//03 03
  long		LoPass;			//04 04
  long		HiPass;   		//05 05
  char		Comments[12];		//06 08
} GREE_CHUNK;

typedef struct DATA_CHUNK{
  char		ID[4];		// "DATA"
  ULONG		Size;		// size of the chunk (DATA) to follow
} DATA_CHUNK;

typedef struct SYNC_CHUNK{
  char		ID[4];		// "SYNC"
  ULONG		Size;		// size of the chunk (SYNC) to follow (sizeof(SYNC_CHUNK)-8)
  char		Tag[4];		//01 01
  ULONG		RecordSize;	//02 02
  ULONG		RecordType;	//03 03 BYTE, USHORT,...
  ULONG		RecordCount;	//04 04
  ULONG		ChannelNumber;	//05 05
  ULONG		ChannelMaxValue;//06 06
  ULONG		StimulusPeriod; //07 07 milliseconds
  ULONG		Free[3];        //08 10
} SYNC_CHUNK;

typedef struct ROIS_CHUNK{
  char		ID[4];		// "ROIS"
  ULONG		Size;		// size of the chunk (ROIS) to follow (sizeof(ROIS_CHUNK)-8)
  char		Tag[4];		//
  ULONG		RecordSize;	//
  ULONG		RecordCount;	//
} ROIS_CHUNK;

typedef struct COMP_CHUNK{
  char		ID[4];			// "COMP" - COMPressed chunk
  ULONG		Size;			// size of the chunk (COMP) to follow = 28
  char		Tag[4];			//01 01
  ULONG		CompressedRecordSize;	//02 02 
  ULONG		CompressedFrameSize;	//03 03
  ULONG		CompressedFrameNumber;	//04 04 
  ULONG		Free[3];		//05 07
} COMP_CHUNK;

typedef struct ROISRecord{
  int		iNumber;
  ULONG	left;
  ULONG	top;
  ULONG	right;
  ULONG	bottom;
}ROISRecord;


typedef struct FILE_HEADER{		// ?Bytes
  ISOI_CHUNK		ISOIChunk;
  SOFT_CHUNK		SOFTChunk;
  HARD_CHUNK		HARDChunk;
  DATA_CHUNK		DATAChunk;
} FILE_HEADER;

typedef struct FILE_COST_HEADER{	// ?Bytes
  ISOI_CHUNK		ISOIChunk;
  SOFT_CHUNK		SOFTChunk;
  HARD_CHUNK		HARDChunk;
  COST_CHUNK		COSTChunk;
  DATA_CHUNK		DATAChunk;
} FILE_COST_HEADER;

typedef struct FILE_EPST_HEADER{	// ?Bytes
  ISOI_CHUNK		ISOIChunk;
  SOFT_CHUNK		SOFTChunk;
  HARD_CHUNK		HARDChunk;
  EPST_CHUNK		EPSTChunk;
  DATA_CHUNK		DATAChunk;
} FILE_EPST_HEADER;

#define TRUE  1
#define FALSE 0

#define DATATYPE_UCHAR		11L
#define DATATYPE_USHORT		12L
#define DATATYPE_ULONG		13L
#define DATATYPE_FLOAT		14L
#define DATATYPE_DOUBLE		15L
#define DEFAULT_DATATYPE	DATATYPE_USHORT

#define VERSION_UNKNOWN 0
#define VERSION_HEADER  1
#define VERSION_CHUNK   2

#define FILE_TYPE_UNKNOWN    0
#define FILE_TYPE_GREEN      1
#define FILE_TYPE_STREAM     2
#define FILE_TYPE_COMPRESSED 3
#define FILE_TYPE_ANALYSIS   4


void SetChunkName(void *pChunk,int iSize,char *pcName);
void SetChunkNameTag(void *pChunk,int iSize,char *pcName,char *pcTag);
int WriteChunk(void *pChunk,int iSize,FILE* fp);
int ReadChunk(void *pChunk,int iSize,FILE* fp);
int ReadChunkHeader(void *pChunk,FILE* fp);
int FindExperimentChunkInFile(FILE *fp);
int FileOffsetCurrent(long lOffset,FILE* fp);
int HasName(void *pChunk,char *pcName);
void PrintChunkName(void *pChunk);

void* FindChunkInBuffer(void *pMultiChunk,ULONG ulMultiChunkSize,char *strName);
unsigned long FindChunkInFileAndOffset(FILE *fp,char *pcName,int iJump);
unsigned long FindChunkInFile(FILE *fp,char *strName);
int InsertChunk(void *pMultiChunk,ULONG ulMultiChunkSize,void *pChunk,ULONG ulChunkSize,void *pMultiChunkEntry);
int RemoveChunk(void *pMultiChunk,ULONG ulMultiChunkSize,ULONG ulChunkSize,void *pMultiChunkEntry);


int DataTypeToSizeOf(int iDataType);

