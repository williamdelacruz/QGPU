/*
 * getqual.c -- generate or process quality map and set quality mode
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "getqual.h"
#include "file.h"
#include "util.h"
#include "pi.h"
#include "qualvar.h"
#include "qualpseu.h"
#include "qualgrad.h"
#include "definitions.h"

/* Generate a new quality map from the phase data, or process */
/* the quality map that was input.                            */
void GetQualityMap(int mode, float *qual_map, float *phase, 
                   unsigned char *bitflags,
                   int tsize, int xsize, int ysize)
{

    float  *temp;
    double rmin, rmax, rscale;
    int    i, j, k;
    
    /* process phase gradients */
    if (mode==variance)
    {
        AllocateFloat(&temp, xsize*ysize, "temp data");
        PhaseVariance(phase, qual_map, bitflags, temp,
                      tsize, xsize, ysize);
        free(temp);
    
        /* convert from cost to quality, and scale to interval (0,1) */
        for (rmin = rmax = qual_map[0], k=0; k<xsize*ysize; k++)
        {
            if (rmin > qual_map[k]) rmin = qual_map[k];
            if (rmax < qual_map[k]) rmax = qual_map[k];
        }
        
        printf("Min & max phase derivative variance = %lf, %lf\n",
               rmin, rmax);
        rscale = (rmin != rmax) ? 1.0/(rmax - rmin) : 0.0;
        
        for (k=0; k<xsize*ysize; k++)
        {
            qual_map[k] = (rmax - qual_map[k])*rscale;
            if (bitflags && bitflags[k]==BORDER)
                qual_map[k] = 0.0;
        }
    }
    else if (mode==gradient)
    {
        AllocateFloat(&temp, xsize*ysize, "temp data");
        MaxPhaseGradients(phase, qual_map, bitflags, temp,
                          tsize, xsize, ysize);
        free(temp);
    
        /* convert from cost to quality, and scale to interval (0,1) */
        for (rmin = rmax = qual_map[0], k=0; k<xsize*ysize; k++)
        {
            if (rmin > qual_map[k]) rmin = qual_map[k];
            if (rmax < qual_map[k]) rmax = qual_map[k];
        }
        
        printf("Min&max 'max phase gradient' = %lf %lf\n", rmin,rmax);
        rscale = (rmin != rmax) ? 1.0/(rmax - rmin) : 0.0;
        
        for (k=0; k<xsize*ysize; k++)
        {
            qual_map[k] = (rmax - qual_map[k])*rscale;
            if (bitflags && bitflags[k]==BORDER)
                qual_map[k] = 0.0;
        }
    }
    else if (mode==pseudocorr)
    {
        AllocateFloat(&temp, xsize*ysize, "temp data");
        PseudoCorrelation(phase, qual_map, bitflags, temp,
                          tsize, xsize, ysize);
    
        for (k=0; k<xsize*ysize; k++)
        {
            if (bitflags && bitflags[k]==BORDER)
                qual_map[k] = 0.0;
        }
        free(temp);
    }
    else if (mode==none)
    {
        for (k=0; k<xsize*ysize; k++)
        {
            qual_map[k] = 1.0;
            if (bitflags && bitflags[k]==BORDER)
                qual_map[k] = 0.0;
        }
    }
    else /* quality map was input */
    {
        /* scale to interval (0,1) */
        for (rmin = rmax = qual_map[0], k=0; k<xsize*ysize; k++)
        {
            if (rmin > qual_map[k]) rmin = qual_map[k];
            if (rmax < qual_map[k]) rmax = qual_map[k];
        }
    
        printf("Min & max corr. coeff. = %lf, %lf\n", rmin, rmax);
        rscale = (rmin != rmax) ? 1.0/(rmax - rmin) : 0.0;
    
        for (k=0; k<xsize*ysize; k++)
        {
            qual_map[k] = (qual_map[k] - rmin)*rscale;
            if (bitflags && bitflags[k]==BORDER)
                qual_map[k] = 0.0;
        }
    }
}

