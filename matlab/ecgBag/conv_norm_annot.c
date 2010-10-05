//-------------------------------------------------------------------------
// filename:    conv_norm_annot.c
// version:             1.0
// authors:     G. Clifford, with much help from P. Hayton, N. Townsend 2001
// (C) authors, under the GNU pulic license
//-------------------------------------------------------------------------
// description: Program to convert the MIT Polysomnographic normal DB 
// annotation files produced by rdann into floating point 
// numbers that represent the beat type and (maybe) patient ailments.

// Expects 2 arguments - the /data/smp/nsrdb file you want to convert 
// and the annotation file you want to record the converted annotations in. 
//
// note that statistics are produced (to stdout) every WINDOW minutes on 
// how many normals, ectopics atefactual and other beats are recorded.
// A the output file is two columns : <time> <beat type (int)>  
// The beat type is an integer such that 
// 1 == normal
// 2 -> 19 == abnormal/ectopics 
// 30 and above == artefact
// 20 -- unspecified
// 25 -- paced 
// 0 -- error - classification not recognised
//
// To compile : gcc -Wall conv_norm_annot.c -o conv_norm_annot -lm
// e.g.
//  rdann -r 16265 -a atr > 16265.ascii  
//  conv_norm_annot 16265.ascii file.out 
//
//
// Reads in MIT polysomographic DB annotations and writes them out
// in a new format: (float)time float(sleep_ann) float(Label)


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

// user defined variables
#define sample_freq 128
#define WINDOW 5

  // globals
  int normals_in_5_min=0;
  int ectopics_in_5_min=0;
  int unspecified_in_5_min=0;
  int paced_in_5_min=0;
  int artefact_in_5_min=0;
  float fivemin_HR=0.0;
  float start_of_wind=0.0;

// file pointers declared
FILE *fptr_in;
FILE *fptr_out;

// define limits and defaults
#define MAX_LINE_LEN 460000 
#define TRUE 1
#define FALSE 0
#define YES 1
#define NO 0

char* scanline (char **ptr, char *buffer, int max)
{
  int n=0;
  char *start = buffer;

  if (**ptr == '\0') return NULL;

  while (**ptr != '\n' && n < max) 
    {
    *buffer = **ptr;
    if (**ptr == '\0') return start;

    (*ptr)++;
    buffer++;
    n++;
    }
  *buffer = **ptr;
  (*ptr)++;
  buffer++;  
  *buffer='\0';  

  return start;
}

int scanstring (char **ptr, char *buffer)
{
  int n=0;

  if (**ptr == '\0') return FALSE;

  // Skip whitespace
  while (**ptr == ' ' || **ptr == '\n') (*ptr)++;
  // Copy until the next space is found
  while (**ptr != '\n' && **ptr != ' ' && **ptr != '\0') 
    {
    *buffer = **ptr;
    if (**ptr == '\0') return TRUE;

    (*ptr)++;
    buffer++;
    n++;
    }
  (*ptr)++;
  *buffer='\0';  

  return TRUE;
}

char *parse_line(FILE *fp, char *line, int max_len, char *parts, char **part, int *n_parts, char** text_stream)
{
char *cp;
char *result;
int    in_white;
int    in_quote;
int    k;

// If fp == NULL, we scan from the text stream instead 
if (fp)                       result = fgets(line, max_len, fp);
else if (text_stream != NULL) result = scanline(text_stream, line, max_len);
else result = 0; 

if (result) {
    cp       = &line[0];
    in_white = YES;
    in_quote = NO;
    *n_parts = 0;
    k        = 0;
    while (*cp != '\0') {
	if (in_quote) {
	    if (*cp == '"') {  
		in_quote = NO;
		in_white = YES;
		parts[k++] = '\0';
		}
	    else parts[k++] = *cp;
	    }
	/* check for white space */
/* 	else if ((*cp == '"') && in_white) {  */
/* 	    in_quote = YES; */
/* 	    part[(*n_parts)++] = &parts[k]; */
/* 	    } */
	else if (*cp == '%') { /* ignore comments */
	    *cp = '\0';
	    continue;
	    }
	else if ((*cp == ' ') || (*cp == '\t') || (*cp == ':') || (*cp == '\n')) {
	    if (!in_white) {
		in_white = YES;
		if (*n_parts) parts[k++] = '\0';
		}
	    }
	else {
	    if (in_white) part[(*n_parts)++] = &parts[k];
	    parts[k++] = *cp;
	    in_white = NO;
	    }
	cp++;
	}

    if (!in_white && *n_parts) parts[k] = '\0';
    }
return result;
}

void reset_globals()
{
 normals_in_5_min=0;
 ectopics_in_5_min=0;
 unspecified_in_5_min=0;
 paced_in_5_min=0;
 artefact_in_5_min=0;
}

int main(int argc, char *argv[])
{
  char  *part[MAX_LINE_LEN/2];
  char   line[MAX_LINE_LEN];
  char   parts[3*MAX_LINE_LEN/2+1];
  int    n_parts;
  float time_s;
  int time_int;

  // if we have less than 2 arguments (executable+file to read)
  if (argc < 3)
  {
    printf("\n Usage: %s <record in> <record_out> \n", argv[0]);
    exit(0);
  }
  // open record (1st argument of executable) for input and check
  // if it can be opened

  if ( !(fptr_out = fopen(argv[2],"w")) )
    {
      fprintf(stderr, "\n Error: cannot open record %s\n",argv[2]);
      exit(0);
    }
  if ( !(fptr_in = fopen(argv[1],"r")) )
    {
      fprintf(stderr, "\n Error: cannot open record %s\n",argv[1]);
      exit(0);
    }

  // loop over the data until the end of the file
  while (parse_line(fptr_in, line, MAX_LINE_LEN, parts, part, &n_parts, NULL))
    {
/* fprintf(stderr,"%d parts : last part is %s\n",n_parts,part[n_parts-1]); */
    time_int = atoi(part[3]);
    // time_float = atof(part[5]);
    time_s = (float)(time_int)/sample_freq;

    //    printf("%i %f\n",time_int,time_s);
    //        printf("%i %i %i %i %i %i %i %i\n",atoi(part[0]),atoi(part[1]),atoi(part[2]),atoi(part[3]),atoi(part[4]),atoi(part[5]),atoi(part[6]),atoi(part[7]));

	//    printf("%i %f\n",n_parts,time_s);
    if (n_parts < 8) fprintf(stderr," -20 "); // error
    else if (n_parts-1 == 7) 
      { // label state - if there is no other label third column is 0
	if (!strcmp(part[4],"N")) {
	  fprintf(fptr_out,"%f 1\n",time_s);
	  normals_in_5_min +=1;
	}
	else if (!strcmp(part[4],"L")){
	  fprintf(fptr_out,"%f 2\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"R")){
	  fprintf(fptr_out,"%f 3\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"B")){
	  fprintf(fptr_out,"%f 4\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"A")){
	  fprintf(fptr_out,"%f 5\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"a")){
	  fprintf(fptr_out,"%f 6\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"J")){
	  fprintf(fptr_out,"%f 7\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"S")){
	  fprintf(fptr_out,"%f 8\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"V")){
	  fprintf(fptr_out,"%f 9\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"r")){
	  fprintf(fptr_out,"%f 10\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"F")){
	  fprintf(fptr_out,"%f 11\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"e")){
	  fprintf(fptr_out,"%f 12\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"j")){
	  fprintf(fptr_out,"%f 13\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"n")){
	  fprintf(fptr_out,"%f 14\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"E")){
	  fprintf(fptr_out,"%f 15\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"f")){
	  fprintf(fptr_out,"%f 16\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"Q")){
	  fprintf(fptr_out,"%f 17\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"\?")){
	  fprintf(fptr_out,"%f 20\n",time_s);
	  unspecified_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"P")){
	  fprintf(fptr_out,"%f 25\n",time_s);
	  paced_in_5_min +=1;
	} 	
	else if (!strcmp(part[4],"|")){
	  fprintf(fptr_out,"%f 30\n",time_s);
	  artefact_in_5_min +=1;
	} 
	else if (!strcmp(part[4],"~")){
	  fprintf(fptr_out,"%f 31\n",time_s);
	  artefact_in_5_min +=1; 
	} 
	else if (!strcmp(part[4],"{")){
	  fprintf(fptr_out,"%f 32\n",time_s);
	  artefact_in_5_min +=1; 
	} 
	else if (!strcmp(part[4],"}")){
	  fprintf(fptr_out,"%f 33\n",time_s);
	  artefact_in_5_min +=1; 
	} 
	else if (!strcmp(part[4],"\\")){
	  fprintf(fptr_out,"%f 33\n",time_s);
	  artefact_in_5_min +=1; 
	} 
	else fprintf(fptr_out,"%f 0\n",time_s); // error
      } 
    else if (n_parts-1 == 8) // we've clicked over onto the next day
      {
	time_int = atoi(part[4]);
	// time_float = atof(part[5]);
	time_s = (float)(time_int)/sample_freq;

	if (!strcmp(part[5],"N")) {
	  fprintf(fptr_out,"%f 1\n",time_s);
	  normals_in_5_min +=1;
	}
	else if (!strcmp(part[5],"L")){
	  fprintf(fptr_out,"%f 2\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"R")){
	  fprintf(fptr_out,"%f 3\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"B")){
	  fprintf(fptr_out,"%f 4\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"A")){
	  fprintf(fptr_out,"%f 5\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"a")){
	  fprintf(fptr_out,"%f 6\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"J")){
	  fprintf(fptr_out,"%f 7\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"S")){
	  fprintf(fptr_out,"%f 8\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"V")){
	  fprintf(fptr_out,"%f 9\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"r")){
	  fprintf(fptr_out,"%f 10\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"F")){
	  fprintf(fptr_out,"%f 11\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"e")){
	  fprintf(fptr_out,"%f 12\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"j")){
	  fprintf(fptr_out,"%f 13\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"n")){
	  fprintf(fptr_out,"%f 14\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"E")){
	  fprintf(fptr_out,"%f 15\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"f")){
	  fprintf(fptr_out,"%f 16\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"Q")){
	  fprintf(fptr_out,"%f 17\n",time_s);
	  ectopics_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"\?")){
	  fprintf(fptr_out,"%f 20\n",time_s);
	  unspecified_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"P")){
	  fprintf(fptr_out,"%f 25\n",time_s);
	  paced_in_5_min +=1;
	} 	
	else if (!strcmp(part[5],"|")){
	  fprintf(fptr_out,"%f 30\n",time_s);
	  artefact_in_5_min +=1;
	} 
	else if (!strcmp(part[5],"~")){
	  fprintf(fptr_out,"%f 31\n",time_s);
	  artefact_in_5_min +=1; 
	} 
	else if (!strcmp(part[5],"{")){
	  fprintf(fptr_out,"%f 32\n",time_s);
	  artefact_in_5_min +=1; 
	} 
	else if (!strcmp(part[5],"}")){
	  fprintf(fptr_out,"%f 33\n",time_s);
	  artefact_in_5_min +=1; 
	} 
	else if (!strcmp(part[5],"\\")){
	  fprintf(fptr_out,"%f 33\n",time_s);
	  artefact_in_5_min +=1; 
	} 
	else fprintf(fptr_out,"%f 0\n",time_s); // error
      }
    else fprintf(fptr_out,"%f -100\n",time_s); // big error
    //fprintf(stdout,"%i %s %i\n",time_int,part[n_parts-1], n_parts-1);
    
    //printf("%i ",((int)(time_s*1000))%(WINDOW*60*1000) );
    //    if( ((int)(time_s*1000))%(WINDOW*60*1000)<50){ 
    if( (time_s-start_of_wind)> (60*WINDOW) ){ 
      start_of_wind=time_s;
      fivemin_HR = normals_in_5_min/WINDOW;
      printf("%f %f %i %i %i %i %i\n",time_s,fivemin_HR,normals_in_5_min,ectopics_in_5_min,artefact_in_5_min,unspecified_in_5_min,paced_in_5_min);
      reset_globals();
    }
    
  }
  return(1);
} 








