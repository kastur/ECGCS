//-------------------------------------------------------------------------
// filename:    conv_sleep_annot.c
// version:             1.0
// authors:     G. Clifford, with much help from P. Hayton, N. Townsend 2001
// (C) authors, under the GNU pulic license
//-------------------------------------------------------------------------
// description: Program to convert the MIT Polysomnographic slp annotations 
// files produced by rdann into floating point numbers that represent the 
// sleep continuum, and patient ailments.
//
// Expects 2 arguments - the slp file you want to convert and the annotation
// file you want to record the converted              
//
// To compile : gcc -Wall conv_sleep_annot.c -o conv_sleep_annot -lm
// e.g.
//  rdann -r slp01a -a st > slp01a.ascii  
//  conv_sleep_annot slp01a.ascii file.out 
//
//
// Reads in MIT polysomographic DB annotations and writes them out
// in a new format: (float)time float(sleep_ann) float(arousal/condition)
// where sleep_ann and arousal/condition = 
/*  -1  W   subject is awake
     0  R   REM sleep
     1  1   sleep stage 1
     2  2   sleep stage 2
     2  2   sleep stage 3
     4  4   sleep stage 4

  -0.1     No label in 3rd column
    -2  H   Hypopnea
  -0.2  HA  Hypopnea with arousal
    -3  OA  Obstructive apnea
  -0.3  X   Obstructive apnea with arousal
    -4  CA  Central apnea
  -0.4  CAA Central apnea with arousal
    -5  L   Leg movements
  -0.5  LA  Leg movements with arousal
  -0.6  A   Unspecified arousal
   -10  MT  Movement time

  In theory, the arousals will only manifest in the 3rd column
  - this is to show that the sleep stage may give misleading
  information about the arousal (plot the third column and 
  one sees that any values between -0.1 and -0.6 mean that the 
  patient experienced an arousal of between 5-30s in this 30s 
  window.

  Classification: Based on Mayela's feedback: 

  Classifying outside stages 1-4, awake and 
  REM sleep is more difficult:
  REM sleep is indistinguishable from Awake 
  when using EEG only and has therefore been 
  given a classification between stage 1 and W. 
  Movement in MT usually causes a strong EEG 
  artefact so this is placed at the end of the list. 
  Assuming that the analysis is based on the
  EEG, only the categories 'with arousal' are 
  those that can be associated with a specific 
  sleep stage: Awake. Therefore these are clustered 
  together and placed near 0.
  An hypopnea or an apnoea can occur at any sleep stage, 
  including REM sleep. 
  Leg movements are only noticeable in EMG signals from the leg. 
  An arousal can occur with no apparent reason. 
  Normal subjects experience about 15 arousals during the night. 
  Arousals caused by OSA are short in duration (3-20 seconds).

  The labels are time stamped in seconds (the period is 30 seconds
  _FROM_ the time)
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

// user defined values
#define sample_freq 250

// file pointers declared
FILE *fptr_in;
FILE *fptr_out;

#define MAX_LINE_LEN 460000
#define TRUE 1
#define FALSE 0
#define YES 1
#define NO 0
#define samp_freq 250

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
    time_int = atoi(part[4]);
    // time_float = atof(part[5]);
    time_s = time_int/samp_freq;
    // printf("%i\n",time_int);

    if (n_parts < 10) fprintf(stderr,"-20"); // error
    else if (n_parts-1 == 9) 
      { // label state - if there is no other label third column is 0
	if (!strcmp(part[9],"W")) fprintf(fptr_out,"%f -1 -0.1\n",time_s);
	else if (!strcmp(part[9],"1")) fprintf(fptr_out,"%f 1 -0.1\n",time_s);
	else if (!strcmp(part[9],"2")) fprintf(fptr_out,"%f 2 -0.1\n",time_s);
	else if (!strcmp(part[9],"3")) fprintf(fptr_out,"%f 3 -0.1\n",time_s);
	else if (!strcmp(part[9],"4")) fprintf(fptr_out,"%f 4 -0.1\n",time_s);
	else if (!strcmp(part[9],"R")) fprintf(fptr_out,"%f 0 -0.1\n",time_s);
	else if (!strcmp(part[9],"H")) fprintf(fptr_out,"%f -2 -0.1\n",time_s);
	else if (!strcmp(part[9],"HA")) fprintf(fptr_out,"%f -0.2 -0.1\n",time_s);
	else if (!strcmp(part[9],"OA")) fprintf(fptr_out,"%f -3 -0.1\n",time_s);
	else if (!strcmp(part[9],"X")) fprintf(fptr_out,"%f -0.3 -0.1\n",time_s);
	else if (!strcmp(part[9],"CA")) fprintf(fptr_out,"%f -4 -0.1\n",time_s);
	else if (!strcmp(part[9],"CAA")) fprintf(fptr_out,"%f -0.4 -0.1\n",time_s);
	else if (!strcmp(part[9],"L")) fprintf(fptr_out,"%f -5 -0.1\n",time_s);
	else if (!strcmp(part[9],"LA")) fprintf(fptr_out,"%f -0.5 -0.1\n",time_s);
	else if (!strcmp(part[9],"A")) fprintf(fptr_out,"%f -0.6 -0.1\n",time_s);
	else if (!strcmp(part[9],"MT")) fprintf(fptr_out,"%f -10 -0.1\n",time_s);
	else fprintf(fptr_out,"%f -21\n",time_s); // error
      }
    else if (n_parts-1 > 9) 
      { // OR label sleep state and arousal level  //* == arousal -- only relevant annot
	if (!strcmp(part[9],"W")) fprintf(fptr_out,"%f 0 ",time_s);
	else if (!strcmp(part[9],"1")) fprintf(fptr_out,"%f 1 ",time_s);
	else if (!strcmp(part[9],"2")) fprintf(fptr_out,"%f 2 ",time_s);
	else if (!strcmp(part[9],"3")) fprintf(fptr_out,"%f 3 ",time_s);
	else if (!strcmp(part[9],"4")) fprintf(fptr_out,"%f 4 ",time_s);
	else if (!strcmp(part[9],"R")) fprintf(fptr_out,"%f 0 ",time_s);
	else if (!strcmp(part[9],"H")) fprintf(fptr_out,"%f -2 ",time_s);
	else if (!strcmp(part[9],"HA")) fprintf(fptr_out,"%f -0.1 ",time_s);//*
	else if (!strcmp(part[9],"OA")) fprintf(fptr_out,"%f -3 ",time_s);
	else if (!strcmp(part[9],"X")) fprintf(fptr_out,"%f -0.3 ",time_s); //*
	else if (!strcmp(part[9],"CA")) fprintf(fptr_out,"%f -4 ",time_s);  
	else if (!strcmp(part[9],"CAA")) fprintf(fptr_out,"%f -0.4 ",time_s);//*
	else if (!strcmp(part[9],"L")) fprintf(fptr_out,"%f -5 ",time_s);
	else if (!strcmp(part[9],"LA")) fprintf(fptr_out,"%f -0.5 ",time_s);//*
	else if (!strcmp(part[9],"A")) fprintf(fptr_out,"%f -0.6 ",time_s); //*
	else if (!strcmp(part[9],"MT")) fprintf(fptr_out,"%f -10 ",time_s); 
	else fprintf(fptr_out,"%f -21 ",time_s); // error
	// label specific awareness level 
	if (!strcmp(part[10],"W")) fprintf(fptr_out,"-1\n");
	else if (!strcmp(part[10],"1")) fprintf(fptr_out,"1\n");
	else if (!strcmp(part[10],"2")) fprintf(fptr_out,"2\n");
	else if (!strcmp(part[10],"3")) fprintf(fptr_out,"3\n");
	else if (!strcmp(part[10],"4")) fprintf(fptr_out,"4\n");
	else if (!strcmp(part[10],"R")) fprintf(fptr_out,"0\n");
	else if (!strcmp(part[10],"H")) fprintf(fptr_out,"-2\n");
	else if (!strcmp(part[10],"HA")) fprintf(fptr_out,"-0.2\n"); //*
	else if (!strcmp(part[10],"OA")) fprintf(fptr_out,"-3\n");
	else if (!strcmp(part[10],"X")) fprintf(fptr_out,"-0.3\n");  //*
	else if (!strcmp(part[10],"CA")) fprintf(fptr_out,"-4\n");
	else if (!strcmp(part[10],"CAA")) fprintf(fptr_out,"0.4\n");//*
	else if (!strcmp(part[10],"L")) fprintf(fptr_out,"-5\n");
	else if (!strcmp(part[10],"LA")) fprintf(fptr_out,"-0.5\n"); //*
	else if (!strcmp(part[10],"A")) fprintf(fptr_out,"-0.6\n");  //*
	else if (!strcmp(part[10],"MT")) fprintf(fptr_out,"-10\n");
	else fprintf(fptr_out,"%f -5\n",time_s); // error
	//	else fprintf(fptr_out,"%f -22\n",time_s); // error
      }
    //    else if (n_parts-1 > 10) 
    //  {
    //  }
    else fprintf(fptr_out,"%f -100\n",time_s); // big error
    
    //fprintf(stdout,"%i %s %i\n",time_int,part[n_parts-1], n_parts-1);
    }
  return(1);
}








