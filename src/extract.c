/*
 * extract.c -- extract phase from input data
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "extract.h"
#include "util.h"
#include "pi.h"
/*
 *  Get the phase data from the input data.  Allocate the temporary
 *  memory required and read the input file.
 *
 *  in_format = 0 for 8-byte complex data,
 *              1 for 4-byte complex data,
 *              2 for 1-byte phase data,
 *              3 for 4-byte float phase data
 *
 *  Output: normalized (0 - 1) phase values in array "phase" 
 */
void GetPhase(int in_format, FILE *ifp, const char *infile, float *phase,
              int xsize, int ysize)
{
  void  *in_data;
    
    
    printf("in_format = %d\n\n", in_format);
               /* floating-pt phase values */
    AllocateFloat((float **)(&in_data), xsize*ysize, "input data");

  printf("Reading input data...\n");
    
  
  if (in_format==3) /* floating-pt phase */
    ReadFloat(ifp, (float *)in_data, xsize*ysize, infile);
  

  ExtractPhase(in_format, in_data, phase, xsize, ysize, 1);
  free(in_data);
}

/* Extract the phase from the input data.  If status_flag is 0 */
/* then do not print any status.                               */
void ExtractPhase(int in_format, void *in_data, float *phase,
                  int xsize, int ysize, int status_flag)
{
    int            j;
    static double  one_over_twopi = 1.0/TWOPI;
    static double  scale;
    float          *float_phase = (float *)in_data;
    
    
    
    /* 4-byte float phase */
    scale = one_over_twopi;
    for (j=0; j<xsize*ysize; j++) {
        /* re-scale phase to interval (0,1) */
        phase[j] = float_phase[j];
    }
    
}
