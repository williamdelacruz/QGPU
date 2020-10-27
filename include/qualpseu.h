#ifndef __QUALPSEU
#define __QUALPSEU
void PseudoCorrelation(float *phase, float *result,
                     unsigned char *bitflags,
                     float *temp1, int tsize, int xsize, int ysize);
void SqrAvgFilter(float *in, float *out, int xsize, int ysize,
                  int size, unsigned char *bitflags,
                  int add_flag);
#endif