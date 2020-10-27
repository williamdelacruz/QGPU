/*
 * getmask.c -- read and process mask input
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "getmask.h"
#include "util.h"
#include "file.h"
#include "maskfat.h"
#include "definitions.h"


/* Read mask file and process matrix */
void GetMask(FILE *mfp, char *bmaskfile, unsigned char **bitflags, int xsize, int ysize)
{
    unsigned char *in_data;
    int k;
    
    /*  temporary data array */
    AllocateByte(&in_data, xsize*ysize, "bitflag data");

    /* border mask data */
    printf("Processing border mask data...\n");
    
    if (mfp!=NULL)
        ReadByte(mfp, in_data, xsize*ysize, bmaskfile);
    else
    {
        for (k=0; k<xsize*ysize; k++)
            in_data[k] = 255;
    }
    
    for (k=0; k<xsize*ysize; k++)
        if (in_data[k]==0)
            in_data[k] = BORDER;
        else in_data[k] = 0;
    
    FattenMask(in_data, 1, xsize, ysize);
    
    // copy in_data to bitflags
    for (k=0; k<xsize*ysize; k++)
        (*bitflags)[k] = in_data[k];
    free(in_data);
}

