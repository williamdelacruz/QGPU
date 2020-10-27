#ifndef __GETQUAL
#define __GETQUAL

typedef enum {none, gradient, variance, pseudocorr,
              corr_coeffs, dxdygrad
} UnwrapMode;

void GetQualityMap(int mode, float *qual_map, float *phase,
                   unsigned char *bitflags,
                   int tsize, int xsize, int ysize);
#endif