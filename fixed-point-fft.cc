#include "fixed-point-fft.h"

//Fixed-point Fast Fourier Transform
//=======================
#define M_PI 3.141592653589793238462643383279502
static void make_fixedpoint_sintbl(int n, int* sintbl)
{
    int i, n2, n4, n8;
    float c, s, dc, ds, t;

    n2 = n / 2;  n4 = n / 4;  n8 = n / 8;
    t = sin(M_PI / n);
    dc = 2 * t * t;  ds = sqrt(dc * (2 - dc));
    t = 2 * dc;  c = 1; s = 0;
    sintbl[n4] = Quantification(1, FFT_SAMP_MAX);
    sintbl[0] = Quantification(0, FFT_SAMP_MAX);
    for (i = 1; i < n8; i++) {
        c -= dc;  dc += t * c;
        s += ds;  ds -= t * s;
        sintbl[i] = Quantification(s, FFT_SAMP_MAX);
        sintbl[n4 - i] = Quantification(c, FFT_SAMP_MAX);
    }
    if (n8 != 0) sintbl[n8] = Quantification(sqrt(0.5), FFT_SAMP_MAX);
    for (i = 0; i < n4; i++)
        sintbl[n2 - i] = sintbl[i];
    for (i = 0; i < n2 + n4; i++)
        sintbl[i + n2] = - sintbl[i];
}

static void make_fixedpoint_bitrev(int n, int* bitrev)
{
    int i, j, k, n2;

    n2 = n / 2;  i = j = 0;
    for ( ; ; ) {
        bitrev[i] = j;
        if (++i >= n) break;
        k = n2;
        while (k <= j) { j -= k; k /= 2; }
        j += k;
    }
}

int fixed_point_fft(int* x, int* y, int n)
{
    static int    fixed_last_n = 0;     /* previous n */
    static int    *fixed_bitrev = NULL; /* bit reversal table */
    static int *fixed_sintbl = NULL; /* trigonometric function table */
    int i, j, k, ik, h, d, k2, n4;
    int t, s, c, dx, dy;

    /* preparation */
    n4 = n / 4;
    if (n != fixed_last_n || n == 0) {
        fixed_last_n = n;
        if (fixed_sintbl != NULL) free(fixed_sintbl);
        if (fixed_bitrev != NULL) free(fixed_bitrev);
        if (n == 0) return 0;             /* free the memory */
        fixed_sintbl = (int*)malloc((n + n4) * sizeof(int));
        fixed_bitrev = (int*)malloc(n * sizeof(int));
        if (fixed_sintbl == NULL || fixed_bitrev == NULL) {
            std::cerr << "out of memory\n";
            return 1;
        }
        make_fixedpoint_sintbl(n, fixed_sintbl);
        make_fixedpoint_bitrev(n, fixed_bitrev);
    }

    /* bit reversal */
    for (i = 0; i < n; i++) {
        j = fixed_bitrev[i];
        if (i < j) {
            t = x[i];  x[i] = x[j];  x[j] = t;
            t = y[i];  y[i] = y[j];  y[j] = t;
        }
    }

    /* transformation */
    for (k = 1; k < n; k = k2) {
        h = 0;  k2 = k + k;  d = n / k2;
        for (j = 0; j < k; j++) {
            c = fixed_sintbl[h + n4];
            s = fixed_sintbl[h];
            for (i = j; i < n; i += k2) {
                ik = i + k;
                dx = FixedRound(FixedMul(s, y[ik]) + FixedMul(c, x[ik]), FFT_FRAC_BITS);
                dy = FixedRound(FixedMul(c, y[ik]) - FixedMul(s, x[ik]), FFT_FRAC_BITS);
                x[ik] = x[i] - dx;  x[i] += dx;
                y[ik] = y[i] - dy;  y[i] += dy;
                /*
                dx = s * y[ik] + c * x[ik];
                dy = c * y[ik] - s * x[ik];
                x[ik] = x[i] - dx;  x[i] += dx;
                y[ik] = y[i] - dy;  y[i] += dy;
                */
            }
            h += d;
        }
    }
    return 0;  /* finished successfully */
}
//=======================

