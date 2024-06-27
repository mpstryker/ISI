#include <stdio.h>
#include <stdlib.h>

#define BUFSIZE 1024
#define MULT 1 

int main(int argc,char **argv){
  char str[BUFSIZE*MULT];
  char filen1[256]={"/home/dump"},filen[256]={"dump"};
  char ch;
  FILE *fp;
  unsigned long i,n;

  if(argc!=3){
    printf("Usage: mdtest int char\n");
    exit(1);
  }
  n=(unsigned long)atoi(argv[1])/MULT;
  ch=argv[2][0];

  if(ch=='w'){
    if(!(fp=fopen(filen,"w"))){
      printf("Cannot open file %s for writing\n",filen);
      exit(1);
    }
    for(i=0;i<n;i++) fwrite((void *)str,BUFSIZE*MULT,1,fp);
  }
  else{
    if(!(fp=fopen(filen,"r"))){
      printf("Cannot open file %s\n",filen);
      exit(1);
    }
    for(i=0;i<n;i++) fread((void *)str,BUFSIZE*MULT,1,fp);
  }

  fclose(fp);
}
