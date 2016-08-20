
#ifndef _AUDIO_GOLDENC_AUDIO_PROCESS_H_
#define _AUDIO_GOLDENC_AUDIO_PROCESS_H_
#include <iostream>
#include <cfloat>
#include <cassert>
#include <cstdio>
#include "fixedpoint.h"

Uint32 GetWaveEnergy(Int16* wave, Int32 wave_len);
Uint32 GetWaveNormEnergy(const Int16* wave, Int32 wave_len);

void HighPassFilter(const Int16* signal_in,
                    Int32 signal_length,
                    Int16* signal_out,
                    Int32* filter_state);

void HighPassFilterFixedPoint(const Int16* signal_in,
                            Int32 signal_length,
                            Int16* signal_out,
                            Int32* filter_state);

void HighPassFilterFixedPointMoreSafe(const Int16* signal_in,
                            Int32 signal_length,
                            Int16* signal_out,
                            Int32* filter_state);

void LowPassFilterCutOff500(const Int16* signal_in,
                            Int32 signal_length,
                            Int16* signal_out,
                            Int32* filter_states);
void LowPassFilterCutOff500FixedPoint(const Int16* signal_in,
                                      Int32 signal_length,
                                      Int16* signal_out,
                                      Int32* filter_states);

inline void GenHammingWindow(float* hamming_window,
                             Int32 window_len) {
    for (Int32 i = 0; i < window_len; ++i) {
        hamming_window[i] = 0.54 - 0.46*cos(M_2PI * i / (window_len-1));
    }
}

void ApplyWindow(Int16* wave,
                 Int32 wave_len,
                 const Uint16* window,
                 Int32 window_len);

inline void ApplyHammingWindow(Int16* wave,
                               Int32 wave_len,
                               const Uint16* hamming_window,
                               Int32 window_len) {
    ApplyWindow(wave, wave_len, hamming_window, window_len);
}

inline void ApplySineWindow(Int16* wave,
                     Int32 wave_len,
                     const Uint16* sine_window,
                     Int32 window_len) {
    ApplyWindow(wave, wave_len, sine_window, window_len);
}

Int16 CalCrossZeroRate(const Int16* wave,
                       Int32 wave_len);

void ComputeLinearEnergy(Int64* spectrum_power,
                         Int32 spec_dim,
                         Int32 num_linear_bins,
                         Int64* linear_energies);
Int32 GetLargestComponentIndex(Int64* spectrum_power, Int32 spec_dim);

float CalInterSpectrumFlatnessFloat(float* spectrum_power, Int32 spec_dim);
Int32 CalInterSpectrumFlatness(Int64* spectrum_power, Int32 spec_dim);
Int32 CalLongTermSpectrumFlatness(Int64* spectrum_powers, //2d array [num_frame * spec_dim]
                                  Int32 num_frame,
                                  Int32 spec_dim);
float CalLongTermSpectrumFlatnessFloat(float* spectrum_powers, //2d array [num_frame * spec_dim]
                                  Int32 num_frame,
                                  Int32 spec_dim);

Int32 CalPeriodicity(Int64* spectrum_power, Int32 spec_dim);

void CalHarmonicityAndClarityFloat(float* frame_wave, // apply hamming windows, 1kHz sampling rate
                               float* hamming_window,
                               Int32 frame_len,
                               float* harmonicity,
                               float* clarity);
void CalHarmonicityAndClarity(Int16* frame_wave, // apply hamming windows, 1kHz sampling rate
                               Int16* hamming_window,
                               Int32 frame_len,
                               Int32* harmonicity,
                               Int32* clarity);
void CalHarmonicityAndClarity(Int16* frame_wave, // apply hamming windows, 1kHz sampling rate
							Int16* hamming_window,
							Int32 frame_len,
							Int32* hamming_corr_inv,
							Int32* harmonicity,
							Int32* clarity);

#endif
