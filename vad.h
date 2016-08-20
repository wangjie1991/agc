#ifndef _AUDIO_GOLDENC_VAD_H_
#define _AUDIO_GOLDENC_VAD_H_
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cassert>
#include "fixedpoint.h"
#include "fixed-point-fft.h"
#include "audio-process.h"

typedef enum { kInSilence, kInSpeech, kUnknown } VadState;

const Int32 sil_to_speech_hangover = 7;
const Int32 speech_to_sil_hangover = 20;

const Int32 energy_vad_frame_length = 240; // 16kHz, 15ms
const Int32 energy_vad_frame_shift = 160; // 16kHz, 10ms
const Int32 energy_vad_num_frames_in_buf = 20;
const Int32 energy_vad_lookback_frames = 4;

const Int32 energy_vad_leading_sil_frames = 20;

//const Uint32 energy_smooth_factor = 32113;// 32768*0.98
const Uint32 energy_smooth_factor = 29491;// 32768*0.9

const Int32 energy_vad_buf_size = energy_vad_frame_length+(energy_vad_num_frames_in_buf-1)*energy_vad_frame_shift;

typedef struct {
    VadState state;

    Uint32 buf_tot_energy_hi, buf_tot_energy_lo;//hi save sum of upper 16bit energy
    Uint32 energy_buf[energy_vad_lookback_frames];
    //vad energy threshold
    Uint32 energy_threshold;

    //hangover
    Int16 consecutive_voice_num; //consecutive num of voice frame
    Int16 consecutive_silence_num;

    //butterworth highpass filter state
    Int32 highpass_state[2];

    //statistic variable
    Int32 voice_frame_num;

    Int32 buf_tot_crosszero_rate;
    Int16 crosszero_rate_buf[energy_vad_lookback_frames];
    //cross zero rate range
    Int32 low_crosszero_rate;
    Int32 high_crosszero_rate;

    Int32 start_pos;
    Int32 end_pos;
    Int32 elems_in_buf;
    Int16 buf[energy_vad_buf_size];
}EnergyVad;

const Int32 spectrum_vad_frame_length = 240; // 16kHz, 15ms
const Int32 spectrum_vad_frame_shift = 160; // 16kHz, 10ms
const Int32 spectrum_vad_fft_size = 512;
const Int32 spectrum_vad_num_linear_bins = 40;
const Int32 spectrum_vad_num_frames_in_buf = energy_vad_num_frames_in_buf - 1;
const Int32 spectrum_vad_lookback_frames = 4;

const Int32 spectrum_vad_leading_sil_frames = 20;
const Int32 spectrum_vad_feature_num = 5;

const Int32 spectrum_vad_buf_size = spectrum_vad_frame_length+(spectrum_vad_num_frames_in_buf-1)*spectrum_vad_frame_shift;

typedef struct {
    VadState state;

    Int32 spectrum_buf_elems, spectrum_buf_idx;
    Int64 spectrum_buf[spectrum_vad_lookback_frames*spectrum_vad_fft_size/2];
    //vad energy threshold
    Int32 thresholds[spectrum_vad_feature_num];//inter-flatness, longterm-flatness, periodicity, hamonicity, clarity
    Int16 vote_weights[spectrum_vad_feature_num];//[2, 2, 2, 3, 3]

    //gaussian pdf
    Int32 sil_mean[spectrum_vad_feature_num];
    Int32 sil_var[spectrum_vad_feature_num]; // std variance
    Int32 voice_mean[spectrum_vad_feature_num];
    Int32 voice_var[spectrum_vad_feature_num]; // std variance

    //hangover
    Int16 consecutive_voice_num; //consecutive num of voice frame
    Int16 consecutive_silence_num;

    //butterworth highpass filter(cutoff 300Hz) state
    Int32 highpass_state[2];
    //butterworth lowpass filter(cutoff 500Hz) state
    Int32 lowpass_state[2];
    Int16 filter_buf[spectrum_vad_frame_length]; // lowpass filter out buf
    Int16 rate1k_wave[spectrum_vad_frame_length/16]; // 1KHz sampling

    Int16 hamming_window[spectrum_vad_frame_length];
    Int16 lowpass_hamming[spectrum_vad_frame_length/16];
    Int32 hamming_corr_invert[spectrum_vad_frame_length/32];

    Int32 fft_real_part[spectrum_vad_fft_size];
    Int32 fft_img_part[spectrum_vad_fft_size];
    Int64 power_spectrum[spectrum_vad_fft_size/2];
    Int64 linear_energies[spectrum_vad_num_linear_bins];

    //statistic variable
    Int32 voice_frame_num;

    Int32 start_pos;
    Int32 end_pos;
    Int32 elems_in_buf;
    Int16 buf[spectrum_vad_buf_size];
}SpectrumEnergyVad;

void InitVad(EnergyVad* vad_inst);

bool SimpleEnergyVoiceDetect(const Int16* wave_frame, Int32 frame_len, EnergyVad* vad, Uint32 mean_energy);
bool MultiChannelEnergyVoiceDetect();
bool SpectrumFlatnessVoiceDetect();

bool IsVoice(const Int16* wave_frame, Int32 frame_len, Int32 cur_frame, EnergyVad* vad_inst);

void CalculateVadEnergyThreshold(const Int16* wave,
	Int32 wave_len,
	Int32 frame_len,
	Int32 frame_shift,
	EnergyVad* vad_inst);


Int32 DetectWaveVAD(const Int16* wave, Int32 wave_len, Int32 frame_len,
        Int32 frame_shift, Int32 leading_silence_frames, std::vector<float>& voice_seg);

Int32 DetectWaveVADTwoStage(const Int16* wave,
                    Int32 wave_len,
                    Int32 tempo_frame_len,
                    Int32 tempo_frame_shift,
                    Int32 spec_frame_len,
                    Int32 spec_frame_shift,
                    Int32 leading_sil_frames,
                    std::vector<float>& voice_seg);

Int32 DetectWaveVADFinal(const Int16* wave,
                        Int32 wave_len,
                        Int32 tempo_frame_len,
                        Int32 tempo_frame_shift,
                        Int32 spec_frame_len,
                        Int32 spec_frame_shift,
                        Int32 leading_sil_frames);

#endif
