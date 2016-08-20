#include "audio-process.h"

static const Int16 kAllPassCoefsQ15[2] = { 20972, 5571 };
static const Int16 kAllPassCoefsQ13[2] = { 5243, 1392 };

static Int32 CalScalingSquare(const Int16* wave,
                              Int32 wave_len) {
    Int16 nbits = GetEffectiveBitsNum((Uint32)wave_len);
    size_t i;
    Int16 smax = -1;
    Int16 sabs;
    const Int16 *sptr = wave;
    Int16 t;
    size_t looptimes = wave_len;

    for (i = looptimes; i > 0; i--) {
        sabs = (*sptr > 0 ? *sptr++ : -*sptr++);
        smax = (sabs > smax ? sabs : smax);
    }
    t = FindNormBitNum((Int32)(smax)*smax);

    if (smax == 0) {
        return 0; // Since norm(0) returns 0
    } else {
        return (t > nbits) ? 0 : nbits - t;
    }
}

Uint32 CalSquareSum(const Int16* vec,
                    Int32 vec_len,
                    Int16* scale_factor) {
    Uint32 square_sum = 0;           
    Int32 scaling = CalScalingSquare(vec, vec_len);
    const Int16 *ptr = vec;   

    for (Int32 i = 0; i < vec_len; i++) {  
        square_sum += ((*ptr * *ptr) >> scaling);
        ptr++;            
    }
    *scale_factor = scaling;  

    return square_sum;
}

Uint32 CalNormSquareSum(const Int16* vec,
                        Int32 vec_len) {
    Uint32 normalized_square_sum = 0;           
    Int32 scaling = CalScalingSquare(vec, vec_len);
    const Int16 *ptr = vec;   

    for (Int32 i = 0; i < vec_len; i++) {  
        normalized_square_sum += ((*ptr * *ptr) >> scaling);
        ptr++;            
    }
    normalized_square_sum = ((normalized_square_sum/vec_len) << scaling) 
                            + ((normalized_square_sum%vec_len)<<scaling)/vec_len;

    return normalized_square_sum;
}

Int32 CalNormInnerProd(const Int16* vec1,
                        const Int16* vec2,
                        Int32 vec_len) {
    Int32 inner_prod = 0;           
    Int16 nbits = GetEffectiveBitsNum((Uint32)vec_len);

    for (Int32 i = 0; i < vec_len; i++) {  
        inner_prod += ((vec1[i]*vec2[i]) >> nbits);
    }
    return ((inner_prod/vec_len)<<nbits) 
                + ((inner_prod%vec_len)<<nbits)/vec_len;
}

Uint32 GetWaveEnergy(const Int16* wave,
                    Int32 wave_len,
                    Int16* scale_factor) {
    return CalSquareSum(wave, wave_len, scale_factor);
}

Uint32 GetWaveNormEnergy(const Int16* wave,
                        Int32 wave_len) {
    return CalNormSquareSum(wave, wave_len);
}

void ApplyWindow(Int16* wave,
                 Int32 wave_len,
                 const Uint16* window,
                 Int32 window_len) {
    for (Int32 i = 0; i < wave_len; ++i) {
        wave[i] = FixedRoundMul(wave[i], window[i], WINDOW_FRAC_BITS);
    }
}

Int16 CalCrossZeroRate(const Int16* wave,
                       Int32 wave_len) {
    Int16 cross_zero_num = 0;
    for (Int32 i = 0; i < wave_len; ++i) {
        Uint16 a1 = wave[i], a2 = wave[i+1];
        cross_zero_num += ((a1^a2)>>15);
    }
    return cross_zero_num; 
}

void HighPassFilter(const Int16* signal_in,
                    Int32 signal_length,
                    Int16* signal_out,
                    Int32* filter_state) {
    const float b0 = 0.97138, b1 = -0.97138;
    const float a1 = -0.94276;
    
    //signal_out[0] = b0*signal[0] - a1*filter_state;
    signal_out[0] = b0*signal_in[0];
    float last_y = signal_out[0];

    for (Int32 i = 1; i < signal_length; ++i) {
        float tmp = b0*signal_in[i] + b1*signal_in[i-1];
        signal_out[i] = tmp - a1*last_y;
        last_y = signal_out[i];
    }
}

void HighPassFilterFixedPoint(const Int16* signal_in,
                            Int32 signal_length,
                            Int16* signal_out,
                            Int32* filter_state) {
    const Int16 b0_Q14 = 15915, b1_Q14 = -15915;
    const Int16 a1_Q14 = -15446;
    
    //signal_out[0] = b0*signal[0] - a1*filter_state;
    filter_state[1] += b0_Q14*signal_in[0];
    filter_state[0] = filter_state[1] - a1_Q14*(filter_state[0]>>14);
    signal_out[0] = (filter_state[0]>>14);

    for (Int32 i = 1; i < signal_length; ++i) {
        Int32 tmp = b0_Q14*signal_in[i] + b1_Q14*signal_in[i-1];
        filter_state[0] = tmp - a1_Q14 * (filter_state[0]>>14);
        signal_out[i] = (filter_state[0]>>14);
    }
    filter_state[1] = b1_Q14*signal_in[signal_length-1];
}

void HighPassFilterFixedPointMoreSafe(const Int16* signal_in,
                            Int32 signal_length,
                            Int16* signal_out,
                            Int32* filter_state) {
    const Int16 b0_Q14 = 15915, b1_Q14 = -15915;
    const Int16 a1_Q13 = -7723;
    
    //signal_out[0] = b0*signal[0] - a1*filter_state;
    *filter_state += ((b0_Q14*signal_in[0])>>2);
    signal_out[0] = (*filter_state)>>12;

    for (Int32 i = 1; i < signal_length; ++i) {
        Int32 tmp = b0_Q14*signal_in[i] + b1_Q14*signal_in[i-1];
        *filter_state = ((tmp>>3) - (a1_Q13 * (*filter_state>>14)))<<1;
        signal_out[i] = (*filter_state>>12);
    }
    *filter_state += ((b1_Q14*signal_in[signal_length-1])>>2);
}

void LowPassFilterCutOff500(const Int16* signal_in,
                            Int32 signal_length,
                            Int16* signal_out,
                            Int32* filter_states) {
    assert(signal_length >= 2);
    //filter coefficient
    const float b0 = 0.00225158, b1 = 0.00450317, b2 = 0.00225158;
    const float a1 = -1.86136115, a2 = 0.87036748;

    //butterworth 2 second lowpass filter, cutoff freq 500Hz
    Int32 y_1 = filter_states[0], y_2 = filter_states[1];
    signal_out[0] = b0*signal_in[0] - a1*y_1 - a2*y_2;
    y_2 = y_1;
    y_1 = signal_out[0];
    signal_out[1] = b0*signal_in[1] + b1*signal_in[0] - a1*y_1 - a2*y_2;
    y_2 = y_1;
    y_1 = signal_out[1];

    for (Int32 i = 2; i < signal_length; ++i) {
        float tmp = b0*signal_in[i] + b1*signal_in[i-1] + b2*signal_in[i-2];
        signal_out[i] = tmp - a1*y_1 - a2*y_2;
        y_2 = y_1;
        y_1 = signal_out[i];
    }
    filter_states[0] = y_1;
    filter_states[1] = y_2;
}

void LowPassFilterCutOff500FixedPoint(const Int16* signal_in,
                                      Int32 signal_length,
                                      Int16* signal_out,
                                      Int32* filter_states) {
    assert(signal_length >= 2);
    //filter coefficient
    const Int32 b0_Q14 = 37, b1_Q14 = 74, b2_Q14 = 37;
    const Int32 a1_Q12 = -7624, a2_Q12 = 3565; //half of a1 and a2, beacause |a1| > 1.0
    
    //butterworth 2 second lowpass filter, cutoff freq 500Hz
    Int32 y_1 = filter_states[0], y_2 = filter_states[1];
    Int32 tmp = ((b0_Q14*signal_in[0])>>2) - a1_Q12*(y_1>>12) - a2_Q12*(y_2>>12);
    y_2 = y_1;
    y_1 = tmp;
    //if (tmp >= (1<<27)) signal_out[0] = 32767;
    //else if (tmp <= -(1<<27)) signal_out[0] = -32768;
    //else signal_out[0] = tmp>>12;
    signal_out[0] = tmp>>12;

    tmp = ((b0_Q14*signal_in[1])>>2) + ((b1_Q14*signal_in[0])>>2)
            - a1_Q12*(y_1>>12) - a2_Q12*(y_2>>12);

    //if (tmp >= (1<<27)) signal_out[1] = 32767;
    //else if (tmp <= -(1<<27)) signal_out[1] = -32768;
    //else signal_out[1] = tmp>>12;
    signal_out[1] = tmp>>12;
    y_2 = y_1;
    y_1 = tmp;

    for (Int32 i = 2; i < signal_length; ++i) {
        Int32 t = (b0_Q14*signal_in[i]) + (b1_Q14*signal_in[i-1])
                       + (b2_Q14*signal_in[i-2]);
        t = (t>>2) - a1_Q12*(y_1>>12) - a2_Q12*(y_2>>12);
        signal_out[i] = t>>12;
        y_2 = y_1;
        y_1 = t;

        //if (t >= (1<<31)) y_1 = (1<<31)-1;
        //else if (t <= -(1<<31)) y_1 = -(1<<31);
        //else y_1 = t;
    }
    filter_states[0] = y_1;
    filter_states[1] = y_2;
}

static void AllPassFilter(const Int16* data_in,
                          Int32 data_length,
                          Int16 filter_coefficient,
                          Int16* filter_state,
                          Int16* data_out) {
    // The filter can only cause overflow (in the w16 output variable)
    // if more than 4 consecutive input numbers are of maximum value and
    // has the the same sign as the impulse responses first taps.
    // First 6 taps of the impulse response:
    // 0.6399 0.5905 -0.3779 0.2418 -0.1547 0.0990

    Int16 tmp16 = 0;
    Int32 tmp32 = 0;
    Int32 state32 = ((Int32) (*filter_state) << 16);  // Q15

    for (Int32 i = 0; i < data_length; i++) {
        tmp32 = state32 + filter_coefficient * *data_in;
        tmp16 = (Int16) (tmp32 >> 16);  // Q(-1)
        *data_out++ = tmp16;
        state32 = (*data_in << 14) - filter_coefficient * tmp16;  // Q14
        state32 <<= 1;  // Q15.
        data_in += 2;
    }

    *filter_state = (Int16) (state32 >> 16);  // Q(-1)
}

void Downsampling(const Int16* signal_in,
                  Int16* signal_out,
                  Int32* filter_state,
                  Int32 in_length) {
    Int16 tmp16_1 = 0, tmp16_2 = 0;
    Int32 tmp32_1 = filter_state[0];
    Int32 tmp32_2 = filter_state[1];
    // Downsampling by 2 gives half length.
    Int32 half_length = (in_length >> 1); 

    // Filter coefficients in Q13, filter state in Q0.
    for (Int32 n = 0; n < half_length; n++) {
        // All-pass filtering upper branch.
        tmp16_1 = (Int16) ((tmp32_1 >> 1) +
                ((kAllPassCoefsQ13[0] * *signal_in) >> 14));
        *signal_out = tmp16_1;
        tmp32_1 = (Int32)(*signal_in++) - ((kAllPassCoefsQ13[0] * tmp16_1) >> 12);

        // All-pass filtering lower branch.
        tmp16_2 = (Int16) ((tmp32_2 >> 1) +
                ((kAllPassCoefsQ13[1] * *signal_in) >> 14));
        *signal_out++ += tmp16_2;
        tmp32_2 = (Int32)(*signal_in++) - ((kAllPassCoefsQ13[1] * tmp16_2) >> 12);
    }
    // Store the filter states.
    filter_state[0] = tmp32_1;
    filter_state[1] = tmp32_2;
}

static void SplitFilter(const Int16* data_in,
                        Int32 data_length,
                        Int16* upper_state,
                        Int16* lower_state,
                        Int16* hp_data_out,
                        Int16* lp_data_out) {
    Int32 half_length = data_length >> 1;  // Downsampling by 2.
    Int16 tmp_out;

    // All-pass filtering upper branch.
    AllPassFilter(&data_in[0], half_length, kAllPassCoefsQ15[0], upper_state,
            hp_data_out);

    // All-pass filtering lower branch.
    AllPassFilter(&data_in[1], half_length, kAllPassCoefsQ15[1], lower_state,
            lp_data_out);

    // Make LP and HP signals.
    for (Int32 i = 0; i < half_length; i++) {
        tmp_out = *hp_data_out;
        *hp_data_out++ -= *lp_data_out;
        *lp_data_out++ += tmp_out;
    }
}

void ComputeLinearEnergy(Int64* spectrum_power,
                         Int32 spec_dim,
                         Int32 num_linear_bins,
                         Int64* linear_energies) {
    const Int32 freq_width = (spec_dim + num_linear_bins - 1)/num_linear_bins;
    Int32 idx = 0;
    for (Int32 i = 0; i < spec_dim; i += freq_width) {
        Int64 sum = 0;
        for (Int32 j = 0; j < freq_width; ++j) {
            sum += spectrum_power[i + j]/freq_width;
        }
        linear_energies[idx++] = sum;
    }
}

Int32 GetLargestComponentIndex(Int64* spectrum_power,
                               Int32 spec_dim) {
    Int64 largest = -FLT_MAX;
    Int32 largest_idx = -1;
    for (Int32 i = 0; i < spec_dim; ++i) {
        if (largest < spectrum_power[i]) {
            largest = spectrum_power[i];
            largest_idx = i;
        }
    }
    return largest_idx;
}

float CalInterSpectrumFlatnessFloat(float* spectrum_power,
                                    Int32 spec_dim) {
    assert(((spec_dim - 1)&spec_dim) == 0);//spec_dim == 2^n
    
    float mean = 0, var = 0;
    for (Int32 i = 0; i < spec_dim; ++i) {
        float t = spectrum_power[i];
        mean += t;
        var += t*t;
    }
    mean /= spec_dim;
    var = var/spec_dim - mean*mean;
    return log2(var + 1E-6);
}

Int32 CalInterSpectrumFlatness(Int64* spectrum_power,
                               Int32 spec_dim) {
    assert(((spec_dim - 1)&spec_dim) == 0);//spec_dim == 2^n

    Int16 nbits = GetEffectiveBitsNum((Uint32)spec_dim) - 1;//8
    Int64 mean = 0;
    for (Int32 i = 0; i < spec_dim; ++i) {
        mean += (spectrum_power[i]>>nbits);
    }
    Int64 max = 0, abs = 0;
    for (Int32 i = 0; i < spec_dim; ++i) {
        abs = spectrum_power[i] - mean; 
        if (abs < 0) abs = -abs;
        if (max < abs) max = abs;
    }

    Int16 max_nbits = GetEffectiveBitsNum((Uint64)max);
    Int16 right_shift = (max_nbits>23) ? max_nbits - 23 : 0;
    Int64 var = 0;
    for (Int32 i = 0; i < spec_dim; ++i) {
        Int64 t = (spectrum_power[i] - mean)>>right_shift;
        var += ((t*t)>>nbits);
    }
    return FixedLog2(var + 1) + (2*right_shift<<16);
}

Int32 CalLongTermSpectrumFlatness(Int64* spectrum_powers, //2d array [num_frame * spec_dim]
                                  Int32 num_frame,
                                  Int32 spec_dim) {
    Int32 tot_flatness = 0;
    if (num_frame == 1) return 0;
    for (Int32 s = 0; s < spec_dim; ++s) {
        //cal mean
        Int64 mean = 0;
        for (Int32 f = 0; f < num_frame; ++f) {
            mean += (spectrum_powers[f*spec_dim + s]/num_frame);
        }

        //get max value in s-th col of spectrum_powers
        Int64 max = 0, abs = 0;
        for (Int32 f = 0; f < num_frame; ++f) {
            abs = spectrum_powers[f*spec_dim + s] - mean; 
            if (abs < 0) abs = -abs;
            if (max < abs) max = abs;
        }
        Int16 max_nbits = GetEffectiveBitsNum((Uint64)max);
        Int16 right_shift = (max_nbits>23) ? max_nbits - 23 : 0; //48bit
        Int64 var = 0;
        for (Int32 f = 0; f < num_frame; ++f) {
            Int64 t = (spectrum_powers[f*spec_dim + s] - mean)>>right_shift;
            var += ((t*t)/num_frame);
        }

        Int32 cur_flatness = FixedLog2(var + 1) + (2*right_shift<<16);
        tot_flatness += cur_flatness;
    }
    return tot_flatness/spec_dim;
}

float CalLongTermSpectrumFlatnessFloat(float* spectrum_powers, //2d array [num_frame * spec_dim]
                                  Int32 num_frame,
                                  Int32 spec_dim) {
    float tot_flatness = 0;
    if (num_frame == 1) return 0;
    for (Int32 s = 0; s < spec_dim; ++s) {
        float mean = 0, var = 0;
        for (Int32 f = 0; f < num_frame; ++f) {
            float t = spectrum_powers[f*spec_dim + s];
            mean += t;
            var += t*t;
        }
        mean /= num_frame;
        var = var/num_frame - mean*mean;
        if (var < 0) var = 0;
        float cur_flatness = log2(var + 1E-6);
        tot_flatness += cur_flatness;
    }
    return tot_flatness/spec_dim;
}

Int32 CalPeriodicity(Int64* spectrum_power,
                     Int32 spec_dim) {
    assert(spec_dim == 256);
    const Int32 low_freq = 2, high_freq = 16;
    const Int32 num_freq_copies = spec_dim/16;
    Int32 max = 0;
    for (Int32 i = low_freq; i < high_freq; ++i) {
        Int32 sum = 0;
        for (Int32 j = 0; j < num_freq_copies; ++j) {
            sum += FixedLog2(spectrum_power[j*num_freq_copies + i] + 1);
        }
        if (max < sum) max = sum;
    }
    return max/num_freq_copies;
}

void CalHarmonicityAndClarity(Int16* frame_wave, // apply hamming windows, 1kHz sampling rate
                               Int16* hamming_window,
                               Int32 frame_len,
                               Int32* harmonicity,
                               Int32* clarity) {
    Int32 corr_max = -1E7, corr_min = 1E7;
    Int32 sum_numerator = 0, sum_denominator = 0;

    Int32 corr_0 = 0;
    sum_numerator = CalNormSquareSum(frame_wave, frame_len);
    sum_denominator = CalNormSquareSum(hamming_window, frame_len);
    //Q.15
    sum_denominator >>= 15;
    //sum_denominator += (1<<15);
    corr_0 = ((sum_numerator/sum_denominator)<<15)
                + ((sum_numerator%sum_denominator)<<15)/sum_denominator;

    for (Int32 offset = 2; offset < frame_len/2; ++offset) {
        sum_numerator = CalNormInnerProd(frame_wave, &frame_wave[offset], frame_len - offset); 
        sum_denominator = CalNormInnerProd(hamming_window, &hamming_window[offset], frame_len - offset);
        //Q.15
        sum_denominator >>= 15;
        //sum_denominator += (1<<15);
        Int32 corr = ((sum_numerator/sum_denominator)<<15)
                        + ((sum_numerator%sum_denominator)<<15)/sum_denominator;
        if (corr_max < corr) {
            corr_max = corr;
        } 
        if (corr_min > corr) {
            corr_min = corr;
        }
    }


    assert(corr_max >= corr_min);
    //clarity = 1.0 - sqrt((corr_0 - corr_max)/(corr_0 - corr_min))
    //Q.15
    Int32 numer = (corr_0 - corr_max);
    if (numer <= 0) numer = -numer;

    Int32 denom = corr_0 - corr_min;
    if (denom < 0) denom = -denom;
    if (denom < numer || denom <= 1) {
        *clarity = 0;
    } else {
        Int32 val = FixedDiv(numer, 1 + denom, 1<<30, 16);
        *clarity = (1<<7) - FixedSqrt32(val);
    }
    //printf("fixed corr_0 %.2f max %.2f min %.2f\n", corr_0/1.0, corr_max/1.0, corr_min/1.0);

    *harmonicity = FixedLog2((Uint32)(corr_max + 1));
}

void CalHarmonicityAndClarity(Int16* frame_wave, // apply hamming windows, 1kHz sampling rate
                               Int16* hamming_window,
                               Int32 frame_len,
                               Int32* hamming_corr_inv,
                               Int32* harmonicity,
                               Int32* clarity) {
    Int32 corr_max = -1E7, corr_min = 1E7;
    Int32 sum_numerator = 0;

    Int32 corr_0 = 0;
    sum_numerator = CalNormSquareSum(frame_wave, frame_len);
    corr_0 = FixedRoundMul(sum_numerator, hamming_corr_inv[0], 15);

    for (Int32 offset = 2; offset < frame_len/2; ++offset) {
        sum_numerator = CalNormInnerProd(frame_wave, &frame_wave[offset], frame_len - offset); 
        Int32 corr = FixedRoundMul(sum_numerator, hamming_corr_inv[offset], 15);

        if (corr_max < corr) {
            corr_max = corr;
        } 
        if (corr_min > corr) {
            corr_min = corr;
        }
    }


    assert(corr_max >= corr_min);
    //clarity = 1.0 - sqrt((corr_0 - corr_max)/(corr_0 - corr_min))
    //Q.15
    Int32 numer = (corr_0 - corr_max);
    if (numer <= 0) numer = -numer;

    Int32 denom = corr_0 - corr_min;
    if (denom < 0) denom = -denom;
    if (denom < numer || denom <= 1) {
        *clarity = 0;
    } else {
        Int32 val = FixedDiv(numer, 1 + denom, 1<<30, 16);
        *clarity = (1<<7) - FixedSqrt32(val);
    }
    //printf("fixed corr_0 %.2f max %.2f min %.2f\n", corr_0/1.0, corr_max/1.0, corr_min/1.0);

    *harmonicity = FixedLog2((Uint32)(corr_max + 1));
}

void CalHarmonicityAndClarityFloat(float* frame_wave, // apply hamming windows, 1kHz sampling rate
                               float* hamming_window,
                               Int32 frame_len,
                               float* harmonicity,
                               float* clarity) {
    float corr_max = -1E7, corr_min = 1E7;
    float sum_numerator = 0, sum_denominator = 0;
    for (Int32 i = 0; i < frame_len; ++i) {
        sum_numerator += frame_wave[i]*frame_wave[i];
        sum_denominator += 1 + hamming_window[i]*hamming_window[i];
    }
    float corr_0 = sum_numerator/sum_denominator;


    for (Int32 offset = 2; offset < frame_len/2; ++offset) {
        sum_numerator = sum_denominator = 0;
        for (Int32 idx = 0; idx < frame_len - offset; ++idx) {
            sum_numerator += frame_wave[idx]*frame_wave[offset + idx];
            sum_denominator += 1 + hamming_window[idx]*hamming_window[offset + idx];
        }
        float corr = sum_numerator/sum_denominator;
        if (corr_max < corr) {
            corr_max = corr;
        }
        if (corr_min > corr) {
            corr_min = corr;
        }
    }

    assert(corr_max >= corr_min);
    float tmp = fabs((corr_0 - corr_max)/(corr_0 - corr_min + 1));
    //if (fabs(tmp) <= 1E-6) tmp = 1;
    *clarity = 1.0 - sqrt(tmp);

    *harmonicity = log2(corr_max + 1);
}
