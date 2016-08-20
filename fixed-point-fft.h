#ifndef _AUDIO_GOLENC_FIXED_POINT_FFT_H_
#define _AUDIO_GOLENC_FIXED_POINT_FFT_H_
#include <iostream>
#include "fixedpoint.h"

//x:实部  y:虚部  n：fft长度
int fixed_point_fft(int* x, int* y, int n);

#endif
