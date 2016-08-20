#ifndef _HR_FIXED_POINT_H_
#define _HR_FIXED_POINT_H_
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include "common.h"
#define M_2PI 6.283185307179586476925286766559005

typedef int64_t Int64;
typedef unsigned long Uint64;
typedef int Int32;
typedef unsigned int Uint32;
typedef short Int16;
typedef unsigned short Uint16;

const Int32 COMMON_FRAC_BITS = 16;
const Int32 COMMON_SAMP_MAX = 65535;

const Int32 WAVE_FRAC_BITS = 0;
const Int32 WAVE_SCALE_FACTOR = 1;

const Int32 FFT_FRAC_BITS = 16;
const Int32 FFT_SAMP_MAX = 65535;

const Int32 WINDOW_FRAC_BITS = 16;
const Int32 WINDOW_SAMP_MAX = 65535;

const Int32 LOG_FRAC_BITS = 16;
const Int32 LOG_SAMP_MAX = 65536;

const Int32 DELTA_FRAC_BITS = 16;
const Int32 DELTA_SAMP_MAX = 65535;

const Int32 NORMALIZE_FRAC_BITS = LOG_FRAC_BITS;
const Int32 NORMALIZE_SAMP_MAX = LOG_SAMP_MAX;

inline Int32 Quantification(float val, Int32 SAMP_MAX) {
    if (val + 1E-6 > 1.0f)
        return SAMP_MAX;
    else
        return floor(val*SAMP_MAX+0.5);
}

inline Uint32 FixedSqrt(Uint64 x) {
    Uint64 res = 0;
    Uint64 op  = x;
    Uint64 one = 1uL << 48; // The second-to-top bit is set:
                              //use 1u << 14 for uint16_t type;
                              //use 1uL<<30 for uInt32_t type

    // "one" starts at the highest power of four <= than the argument.
    while (one > op)
    {
        one >>= 2;
    }

    while (one != 0)
    {
        if (op >= res + one)
        {
            op = op - (res + one);
            res = res +  2 * one;
        }
        res >>= 1;
        one >>= 2;
    }
    return res;
}

inline Uint16 FixedSqrt32(Uint32 x) {
    Uint64 res = 0;
    Uint64 op  = x;
    Uint64 one = 1uL << 30; // The second-to-top bit is set:
                              //use 1u << 14 for uint16_t type;
                              //use 1uL<<30 for uInt32_t type

    // "one" starts at the highest power of four <= than the argument.
    while (one > op)
    {
        one >>= 2;
    }

    while (one != 0)
    {
        if (op >= res + one)
        {
            op = op - (res + one);
            res = res +  2 * one;
        }
        res >>= 1;
        one >>= 2;
    }
    return res;
}

inline Int32 FindHighestOrderBit(Uint64 n) {
    n |= (n >>  1);
    n |= (n >>  2);
    n |= (n >>  4);
    n |= (n >>  8);
    n |= (n >> 16);
    n |= (n >> 32);
    n = n - (n >> 1);
    Int32 index = 0;
    while (n >>= 1) { index++; }
    return index;
}

inline Int16 FindHighestOrderBit(Uint32 n) {
    n |= (n >>  1);
    n |= (n >>  2);
    n |= (n >>  4);
    n |= (n >>  8);
    n |= (n >> 16);
    n = n - (n >> 1);
    Int32 index = 0;
    while (n >>= 1) { index++; }
    return index;
}

inline Int16 FindNormBitNum(Uint32 a) {
    int16_t zeros;
              
    if (a == 0) {
        return 0;
    } 

    if (!(0xFFFF8000 & a)) {
        zeros = 16;
    } else {
        zeros = 0;
    }     
    if (!(0xFF800000 & (a << zeros))) zeros += 8;
    if (!(0xF8000000 & (a << zeros))) zeros += 4;
    if (!(0xE0000000 & (a << zeros))) zeros += 2;
    if (!(0xC0000000 & (a << zeros))) zeros += 1;

    return zeros;
}

inline Int16 GetEffectiveBitsNum(Uint32 n) {
    return FindHighestOrderBit(n) + 1;
}

inline Int16 GetEffectiveBitsNum(Uint64 n) {
    return FindHighestOrderBit(n) + 1;
}

inline Int32 FixedRound(Int64 val, Int32 FRAC_BITS) {
    return (val + (1L<<(FRAC_BITS-1))) >> FRAC_BITS;
}

inline Int32 FixedRoundSqrt(Int64 val, Int32 FRAC_BITS) {
    return FixedSqrt((val + (1L<<(FRAC_BITS-1))) >> FRAC_BITS);
}

inline Int64 FixedMul(Int32 a, Int32 b) {
    Int64 prod = (Int64)(a)*b;
    return prod;
}

inline Uint64 FixedMul(Uint32 a, Uint32 b) {
    Uint64 prod = (Int64)(a)*b;
    return prod;
}

inline Int32 FixedRoundMul(Int32 a, Int32 b, Int32 FRAC_BITS) {
    return (FixedMul(a, b) + (1L<<(FRAC_BITS-1))) >> FRAC_BITS;
}

inline Uint32 FixedRoundMul(Uint32 a, Uint32 b, Int32 FRAC_BITS) {
    return (FixedMul(a, b) + (1L<<(FRAC_BITS-1))) >> FRAC_BITS;
}

inline Int32 FixedDiv(Int32 a, Int32 b, Int32 SAMP_MAX, Int32 FRAC_BITS) {
    Int64 prod = (Int64)(a)*(SAMP_MAX/b);
    return (prod + (1L<<(FRAC_BITS-1))) >> FRAC_BITS;
}

inline Int32 FixedHalfOf(Int32 a) {
    return (a>>1);
}

inline Int32 FixedLog2(Int32 frac_part_plus_one, Int32 FRAC_BITS) {
    static Int32 *log2_lookup_table;
    Int32 table_size = 1L<<FRAC_BITS;
    if (log2_lookup_table == NULL) {
        Malloc(&log2_lookup_table, table_size);
        for (Int32 i = 0; i < table_size; ++i) {
            log2_lookup_table[i] = log2((1024 + i + 0.5)/1024.0) * (LOG_SAMP_MAX - 1) + 0.5;
        }
    }
    Int32 index = (frac_part_plus_one & 0x3ff);
    return log2_lookup_table[index];

}

inline Int32 FixedLog(Uint64 x) {
    Int32 highest_order = FindHighestOrderBit(x);
    const Int32 FRAC_BITS=10;
    const Int32 inv_log2e = Quantification(0.6931471805599453, LOG_SAMP_MAX);
    Int32 frac_part_plus_one;
    if (FRAC_BITS <= highest_order) {
        frac_part_plus_one = (x >> (highest_order - FRAC_BITS));
    } else {
        frac_part_plus_one = (x << (FRAC_BITS - highest_order));
    }
    Int32 log_value = FixedLog2(frac_part_plus_one, FRAC_BITS);
    log_value += (highest_order*LOG_SAMP_MAX);
    return FixedRoundMul(log_value, inv_log2e, LOG_FRAC_BITS);
}

inline Int32 FixedLog2(Uint64 x) {
    Int32 highest_order = FindHighestOrderBit(x);
    const Int32 FRAC_BITS=10;
    Int32 frac_part_plus_one;
    if (FRAC_BITS <= highest_order) {
        frac_part_plus_one = (x >> (highest_order - FRAC_BITS));
    } else {
        frac_part_plus_one = (x << (FRAC_BITS - highest_order));
    }
    Int32 log_value = FixedLog2(frac_part_plus_one, FRAC_BITS);
    log_value += (highest_order*LOG_SAMP_MAX);
    return log_value;
}

inline Int32 FixedLogAndAddLn2(Uint64 x) {
    Int32 highest_order = FindHighestOrderBit(x);
    const Int32 FRAC_BITS=10;
    const Int32 inv_log2e = Quantification(0.6931471805599453, LOG_SAMP_MAX);
    Int32 frac_part_plus_one;
    if (FRAC_BITS <= highest_order) {
        frac_part_plus_one = (x >> (highest_order - FRAC_BITS));
    } else {
        frac_part_plus_one = (x << (FRAC_BITS - highest_order));
    }
    Int32 log_value = FixedLog2(frac_part_plus_one, FRAC_BITS);
    log_value += ((highest_order+1)*LOG_SAMP_MAX);
    return FixedRoundMul(log_value, inv_log2e, LOG_FRAC_BITS);
}

inline Int32 FixedLog32(Uint32 x) {
    Int32 highest_order = FindHighestOrderBit(x);
    const Int32 FRAC_BITS=10;
    const Int32 inv_log2e = Quantification(0.6931471805599453, LOG_SAMP_MAX);
    Int32 frac_part_plus_one;
    if (FRAC_BITS <= highest_order) {
        frac_part_plus_one = (x >> (highest_order - FRAC_BITS));
    } else {
        frac_part_plus_one = (x << (FRAC_BITS - highest_order));
    }
    Int32 log_value = FixedLog2(frac_part_plus_one, FRAC_BITS);
    log_value += (highest_order*LOG_SAMP_MAX);
    return FixedRoundMul(log_value, inv_log2e, LOG_FRAC_BITS);
}

#endif
