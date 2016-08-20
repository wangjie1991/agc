#include "vad.h"

inline static void UpdateEnergyThreshold(EnergyVad* vad_inst,
                                    Uint32 mean_energy) {
    Uint32 threshold = vad_inst->energy_threshold;
    threshold = FixedRoundMul(threshold, energy_smooth_factor, 15) \
                                 + FixedRoundMul(mean_energy, 32768 - energy_smooth_factor, 15);
    vad_inst->energy_threshold = threshold > 10000 ? threshold : 10000;
}

bool SimpleEnergyVoiceDetect(const Int16* wave_frame,
                             Int32 frame_len,
                             EnergyVad* vad_inst,
                             Uint32* cur_energy) { //output cur mean energy
    const Int16* wave = wave_frame;
    Uint32 energy = GetWaveNormEnergy(wave, frame_len);
    Int32 num_frames = energy_vad_lookback_frames;

    Uint32 energy_remove = vad_inst->energy_buf[0];
    vad_inst->buf_tot_energy_hi -= (energy_remove>>16);
    vad_inst->buf_tot_energy_lo -= (energy_remove & 0xffff);
    vad_inst->buf_tot_energy_hi += (energy>>16);
    vad_inst->buf_tot_energy_lo += (energy & 0xffff);

    for (Int32 i = 1; i < num_frames; ++i) {
        vad_inst->energy_buf[i-1] = vad_inst->energy_buf[i];
    }
    vad_inst->energy_buf[num_frames-1] = energy;

    Uint32 hi = vad_inst->buf_tot_energy_hi, lo = vad_inst->buf_tot_energy_lo;
    Uint32 mean_energy = (hi/num_frames)<<16;
    mean_energy += (((hi%num_frames)<<16) + lo)/num_frames;
    
    bool is_voice = false;
    Uint32 threshold = vad_inst->energy_threshold;
    if (threshold*2 < mean_energy) {
        is_voice = true;
    } else {
        UpdateEnergyThreshold(vad_inst, mean_energy);
    }
    *cur_energy = mean_energy;

    return is_voice;
}

bool CrossZeroRateVoiceDetect(const Int16* wave_frame,
                              Int32 frame_len,
                              EnergyVad* vad_inst,
                              Int32* cur_rate) {
    Int16 rate = CalCrossZeroRate(wave_frame, frame_len);
    Int16 rate_remove = vad_inst->crosszero_rate_buf[0];
    vad_inst->buf_tot_crosszero_rate += (rate - rate_remove);

    Int32 num_frames = energy_vad_lookback_frames;
    for (Int32 i = 1; i < num_frames; ++i) {
        vad_inst->crosszero_rate_buf[i-1] = vad_inst->crosszero_rate_buf[i];
    }
    vad_inst->crosszero_rate_buf[num_frames-1] = rate;
    
    Int16 mean_rate = vad_inst->buf_tot_crosszero_rate/num_frames;
    bool is_voice = false;
    if (vad_inst->low_crosszero_rate < mean_rate
            && vad_inst->high_crosszero_rate > mean_rate) {
        is_voice = true;
    }
    *cur_rate = mean_rate;

    return is_voice;
}

void ComputePowerSpectrumFixedPoint(const Int32 *real_part, //input, use first 256 elems
	const Int32 *img_part, //input, use first 256 elems
	const Int32 num_elems, //256   
	Int64 *power_spectrum //output
	) {

	const Int32 half_dim = num_elems;
	Int64 first_energy = FixedRoundMul(real_part[0], real_part[0], 1);
	Int64 last_energy = FixedRoundMul(img_part[1], img_part[1], 1);
	for (Int32 i = 1; i < half_dim; i++) {
		Int32 real = real_part[i], img = img_part[i];
		power_spectrum[i] = FixedRoundMul(real, real, 1) + FixedRoundMul(img, img, 1);
	}
	power_spectrum[0] = first_energy;
	power_spectrum[half_dim] = last_energy;
}


bool SpectrumFlatnessVoiceDetect(const Int16* wave_frame,
                                 Int32 frame_len,
                                 SpectrumEnergyVad* vad_inst,
                                 Int32* cur_flatness) {
    assert(frame_len <= spectrum_vad_fft_size);

    Int32* real_part = vad_inst->fft_real_part;
    Int32* img_part = vad_inst->fft_img_part;
    Int16* hamming_window = vad_inst->hamming_window;
    const Int32 fft_size = spectrum_vad_fft_size;
    for (Int32 i = 0; i < frame_len; ++i) {
        real_part[i] = FixedRoundMul(wave_frame[i], hamming_window[i], 15);
    }
    memset(real_part + frame_len, 0, (fft_size - frame_len)*sizeof(Int32));
    memset(img_part, 0, fft_size*sizeof(Int32));
    fixed_point_fft(real_part, img_part, fft_size);

    const Int32 half_dim = fft_size/2;
    ComputePowerSpectrumFixedPoint(real_part, img_part,
                                   half_dim, vad_inst->power_spectrum);

    GetLargestComponentIndex(vad_inst->power_spectrum, half_dim);
    ComputeLinearEnergy(vad_inst->power_spectrum, half_dim,
                        spectrum_vad_num_linear_bins, vad_inst->linear_energies);

    return false;
}

void ComputeVADMeasure(const Int16* wave_frame,
                       Int32 frame_len,
                       SpectrumEnergyVad* vad_inst,
                       Int32* cur_measures) {
    //windowing
    assert(frame_len <= spectrum_vad_fft_size);
    assert(frame_len == spectrum_vad_frame_length);

    Int32* real_part = vad_inst->fft_real_part;
    Int32* img_part = vad_inst->fft_img_part;
    Int16* hamming_window = vad_inst->hamming_window;
    const Int32 fft_size = spectrum_vad_fft_size;
    for (Int32 i = 0; i < frame_len; ++i) {
        real_part[i] = FixedRoundMul(wave_frame[i], hamming_window[i], 15);
    }
    memset(real_part + frame_len, 0, (fft_size - frame_len)*sizeof(Int32));
    memset(img_part, 0, fft_size*sizeof(Int32));
    fixed_point_fft(real_part, img_part, fft_size);
    const Int32 half_dim = fft_size/2;
    ComputePowerSpectrumFixedPoint(real_part, img_part,
                                   half_dim, vad_inst->power_spectrum);

    Int32 inter_frame_flat = CalInterSpectrumFlatness(vad_inst->power_spectrum,
                                                      half_dim);
    {
        Int32 idx = vad_inst->spectrum_buf_idx,
              num_elems = vad_inst->spectrum_buf_elems;
        memcpy(vad_inst->spectrum_buf + idx*half_dim,
                vad_inst->power_spectrum, half_dim*sizeof(float));
        vad_inst->spectrum_buf_idx = (idx + 1)%spectrum_vad_lookback_frames;
        if (num_elems < spectrum_vad_lookback_frames) {
            vad_inst->spectrum_buf_elems = num_elems + 1;
        }
    }
    Int32 longterm_flat = CalLongTermSpectrumFlatness(vad_inst->spectrum_buf,
                                                      vad_inst->spectrum_buf_elems,
                                                      half_dim);
    cur_measures[0] = inter_frame_flat;
    cur_measures[1] = longterm_flat;

    Int32 periodicity = CalPeriodicity(vad_inst->power_spectrum, half_dim);
    cur_measures[2] = periodicity;

    //calculate pitch feature
    Int16* filter_buf = vad_inst->filter_buf;
    LowPassFilterCutOff500FixedPoint(wave_frame, frame_len,
                                     filter_buf, vad_inst->lowpass_state);
    Int32 lowpass_frame_len = frame_len/16;
    for (Int32 i = 0; i < lowpass_frame_len; ++i) {
        vad_inst->rate1k_wave[i] = FixedRoundMul(filter_buf[i*16], vad_inst->lowpass_hamming[i], 15);
    }
    Int32 harmonicity = 0, clarity = 0;
    CalHarmonicityAndClarity(vad_inst->rate1k_wave, vad_inst->lowpass_hamming,
                             lowpass_frame_len, vad_inst->hamming_corr_invert,
                             &harmonicity, &clarity);
    cur_measures[3] = harmonicity;
    cur_measures[4] = clarity;
}

bool SpectrumAndHarmonDetect(const Int16* wave_frame,
                             Int32 frame_len,
                             SpectrumEnergyVad* vad_inst,
                             Int32* cur_measures) {
    ComputeVADMeasure(wave_frame, frame_len, vad_inst, cur_measures);

    Int32 voice_vote = 0, sum_weights = 0;
    for (Int32 i = 0; i < spectrum_vad_feature_num; ++i) {
        if (cur_measures[i] > vad_inst->thresholds[i]*1.1 + 1) {
            voice_vote += vad_inst->vote_weights[i];
        }
        else if  (cur_measures[i] < vad_inst->thresholds[i]*1.1 + 1) {
            const Int32 smooth_factor = 29491; // 32768*0.9
            Int32 thresh = vad_inst->thresholds[i];
            vad_inst->thresholds[i] = FixedRoundMul(thresh, smooth_factor, 15)
                                + FixedRoundMul(cur_measures[i], 32768 - smooth_factor, 15);
        }
        sum_weights += vad_inst->vote_weights[i];
    }
    if (voice_vote > sum_weights/2) {
        return true;
    } else {
        return false;
    }
}

bool SpectrumAndHarmonDetectBaseGaussian(const Int16* wave_frame,
                                         Int32 frame_len,
                                         SpectrumEnergyVad* vad_inst,
                                         Int32* cur_measures) {
    ComputeVADMeasure(wave_frame, frame_len, vad_inst, cur_measures);

    Int32 voice_prob = 0, sil_prob = 0;
    for (Int32 i = 0; i < spectrum_vad_feature_num; ++i) {
        Int32 tmp = cur_measures[i] - vad_inst->voice_mean[i];
        if (tmp < 0) tmp = -tmp;
        Int32 z_voice = FixedDiv(tmp, vad_inst->voice_var[i], 1<<30, 26);
        voice_prob -= z_voice;

        //silence probability
        tmp = cur_measures[i] - vad_inst->sil_mean[i];
        if (tmp < 0) tmp = -tmp;
        Int32 z_sil = FixedDiv(tmp, vad_inst->sil_var[i], 1<<30, 26);
        sil_prob -= z_sil;
    }
    //printf("voice %d sil %d\n", voice_prob, sil_prob);
    const Int32 kfactor = 16;
    if (voice_prob > sil_prob + kfactor) {
        return true;
    } else {
        return false;
    }
}

bool IsVoice(const Int16* wave_frame,
             Int32 frame_len,
             Int32 cur_frame,
             EnergyVad* vad_inst) {
    bool is_voice_frame = false;
    Uint32 mean_energy = 0;
    is_voice_frame = SimpleEnergyVoiceDetect(wave_frame, frame_len, vad_inst, &mean_energy);

    //std::cout << "mean_enrgy " << mean_energy << " threshold " << vad_inst->energy_threshold << std::endl;
    if (is_voice_frame) {
        if (vad_inst->state == kInSilence) {
            vad_inst->consecutive_voice_num++;
            if (vad_inst->consecutive_voice_num >= sil_to_speech_hangover) {
                vad_inst->state = kInSpeech;
                vad_inst->consecutive_voice_num = 0;
                vad_inst->voice_frame_num = 0;
            }
        } else if (vad_inst->state == kInSpeech) {
            vad_inst->consecutive_silence_num = 0;

            vad_inst->voice_frame_num++;
            const Int32 Too_Long_Voice_segment = 1000;
            if (vad_inst->voice_frame_num >=  Too_Long_Voice_segment) {
                fprintf(stderr, "too long voice segment, begin at %d frame\n", cur_frame - Too_Long_Voice_segment);
            }
        } else {
            fprintf(stderr, "unkown vad state\n");
        }
    } else { // is silence
        if (vad_inst->state == kInSilence) {
            vad_inst->consecutive_voice_num = 0;
        } else if (vad_inst->state == kInSpeech) {
            vad_inst->consecutive_silence_num++;
            if (vad_inst->consecutive_silence_num >= speech_to_sil_hangover) {
                vad_inst->state = kInSilence;
                vad_inst->consecutive_silence_num = 0;
            }
        } else {
            fprintf(stderr, "unkown vad state\n");
        }
    }

    if (vad_inst->state == kInSpeech) return true;
    else return false;
}

bool IsVoice(const Int16* wave_frame,
             Int32 frame_len,
             Int32 cur_frame,
             SpectrumEnergyVad* vad_inst) {
    Int32 spectrum_measures[spectrum_vad_feature_num];
    bool is_voice_frame = SpectrumAndHarmonDetect(wave_frame, frame_len,
                                                  vad_inst, spectrum_measures);
    /*
    for (Int32 i = 0; i < spectrum_vad_feature_num - 1; ++i) {
        printf("[T %.2f M %.2f] ", vad_inst->thresholds[i]/65536.0, spectrum_measures[i]/65536.0);
    }
    printf("[T %.2f M %.2f]\n", vad_inst->thresholds[4]/128.0, spectrum_measures[4]/128.0);
    */

    if (is_voice_frame) {
        if (vad_inst->state == kInSilence) {
            vad_inst->consecutive_voice_num++;
            if (vad_inst->consecutive_voice_num >= sil_to_speech_hangover) {
                vad_inst->state = kInSpeech;
                vad_inst->consecutive_voice_num = 0;
                vad_inst->voice_frame_num = 0;
            }
        } else if (vad_inst->state == kInSpeech) {
            vad_inst->consecutive_silence_num = 0;

            vad_inst->voice_frame_num++;
            const Int32 Too_Long_Voice_segment = 1000;
            if (vad_inst->voice_frame_num >=  Too_Long_Voice_segment) {
                fprintf(stderr, "too long voice segment, begin at %d frame\n", cur_frame - Too_Long_Voice_segment);
            }

        } else {
            fprintf(stderr, "unkown vad state\n");
        }
    } else { // is silence
        if (vad_inst->state == kInSilence) {
            vad_inst->consecutive_voice_num = 0;
        } else if (vad_inst->state == kInSpeech) {
            vad_inst->consecutive_silence_num++;
            if (vad_inst->consecutive_silence_num >= speech_to_sil_hangover) {
                vad_inst->state = kInSilence;
                vad_inst->consecutive_silence_num = 0;
            }
        } else {
            fprintf(stderr, "unkown vad state\n");
        }
    }

    if (vad_inst->state == kInSpeech) return true;
    else return false;
}

void InitVad(SpectrumEnergyVad* vad_inst) {
    vad_inst->state = kInSilence;

    vad_inst->spectrum_buf_elems = vad_inst->spectrum_buf_idx = 0;

    memset(vad_inst->thresholds, 0, sizeof(vad_inst->thresholds));
    //init gaussian mean and variance
    {
        Int32 sil_mean[spectrum_vad_feature_num] = { 2493710, 971882, 943086, 794977, 59};
        Int32 sil_var[spectrum_vad_feature_num] = { 340940, 200610, 204087, 53633, 30};
        Int32 voice_mean[spectrum_vad_feature_num] = { 3324520, 1309320, 1267790, 1008580, 77};
        Int32 voice_var[spectrum_vad_feature_num] = { 376120, 178552, 214402, 161976, 24};
        memcpy(vad_inst->sil_mean, sil_mean, sizeof(sil_mean));
        memcpy(vad_inst->sil_var, sil_var, sizeof(sil_var));
        memcpy(vad_inst->voice_mean, voice_mean, sizeof(voice_mean));
        memcpy(vad_inst->voice_var, voice_var, sizeof(voice_var));
    }


    memset(vad_inst->spectrum_buf, 0, sizeof(vad_inst->spectrum_buf));

    vad_inst->consecutive_voice_num = vad_inst->consecutive_silence_num = 0;
    Int16 weights[] = {2, 2, 2, 3, 3};
    memcpy(vad_inst->vote_weights, weights, sizeof(weights));

    vad_inst->highpass_state[0] = 0;
    vad_inst->highpass_state[1] = 0;
    memset(vad_inst->lowpass_state, 0, sizeof(vad_inst->lowpass_state));
    for (Int32 i = 0; i < spectrum_vad_frame_length; ++i) {
        float val = 0.54 - 0.46*cos(M_2PI * i / (spectrum_vad_frame_length-1));
        vad_inst->hamming_window[i] = Quantification(val, 32767);
    }

    {
        const int frame_len = spectrum_vad_frame_length/16;
        for (Int32 i = 0; i < frame_len; ++i) {
            float val = 0.54 - 0.46*cos(M_2PI * i / (frame_len-1));
            vad_inst->lowpass_hamming[i] = Quantification(val, 32767);
        }
        const Int32 hamming_corr_inv[] = { 3433, 3467, 3566, 3738,
                                           3993, 4349, 4830, 5475,
                                           6335, 7488, 9050, 11194 };
        memcpy(vad_inst->hamming_corr_invert, hamming_corr_inv, sizeof(hamming_corr_inv));
    }

    vad_inst->voice_frame_num = 0;

    vad_inst->start_pos = vad_inst->end_pos = 0;
    vad_inst->elems_in_buf = 0;
    memset(vad_inst->buf, 0, sizeof(vad_inst->buf));
}

void InitVad(EnergyVad* vad_inst) {
    vad_inst->state = kInSilence;

    vad_inst->buf_tot_energy_hi = vad_inst->buf_tot_energy_lo = 0;
    memset(vad_inst->energy_buf, 0, sizeof(vad_inst->energy_buf));

    vad_inst->consecutive_voice_num = vad_inst->consecutive_silence_num = 0;

    vad_inst->voice_frame_num = 0;

    vad_inst->highpass_state[0] = 0;
    vad_inst->highpass_state[1] = 0;

    vad_inst->energy_threshold = 0;
    vad_inst->low_crosszero_rate = 5;
    vad_inst->high_crosszero_rate = 20;

    vad_inst->start_pos = vad_inst->end_pos = 0;
    vad_inst->elems_in_buf = 0;
    memset(vad_inst->buf, 0, sizeof(vad_inst->buf));
}

void ResetVadState(SpectrumEnergyVad* vad_inst) {
    vad_inst->state = kInSilence;

    vad_inst->spectrum_buf_elems = vad_inst->spectrum_buf_idx = 0;

    memset(vad_inst->spectrum_buf, 0, sizeof(vad_inst->spectrum_buf));

    vad_inst->consecutive_voice_num = vad_inst->consecutive_silence_num = 0;

    vad_inst->highpass_state[0] = 0;
    vad_inst->highpass_state[1] = 0;
    memset(vad_inst->lowpass_state, 0, sizeof(vad_inst->lowpass_state));

    vad_inst->voice_frame_num = 0;

    vad_inst->start_pos = vad_inst->end_pos = 0;
    vad_inst->elems_in_buf = 0;
    memset(vad_inst->buf, 0, sizeof(vad_inst->buf));
}

void CalculateVadEnergyThreshold(const Int16* wave,
                                 Int32 wave_len,
                                 Int32 frame_len,
                                 Int32 frame_shift,
                                 EnergyVad* vad_inst) {
  
	Int16 *wave_frame = (Int16*)malloc(frame_len*sizeof(Int16));
    Int32 overlap = frame_len - frame_shift;

    Uint32 total_energy_high = 0, total_energy_low = 0;
    Int32 num_frames = (wave_len - frame_len)/frame_shift + 1;

    memcpy(wave_frame + frame_shift, wave, sizeof(Int16)*overlap);

    for (Int32 i = 0; i < num_frames; ++i) {
        Int32 begin = i*frame_shift + overlap;
        memmove(wave_frame, wave_frame + frame_shift, overlap*sizeof(Int16));
        memcpy(wave_frame + overlap, wave + begin, sizeof(Int16)*frame_shift);

        Uint32 energy = GetWaveNormEnergy(wave_frame, frame_len);
        total_energy_high += (energy>>16);
        total_energy_low += (energy & 0xffff);

        if (i + energy_vad_num_frames_in_buf >= num_frames) {
            if (vad_inst->elems_in_buf == 0) {
                memcpy(vad_inst->buf, wave_frame, frame_len*sizeof(Int16));
            } else if (vad_inst->elems_in_buf < energy_vad_num_frames_in_buf) {
                Int32 offset = overlap + frame_shift*vad_inst->elems_in_buf;
                memcpy(vad_inst->buf + offset, wave_frame + overlap, frame_shift*sizeof(Int16));
            }
            vad_inst->elems_in_buf++;
        }
    }
    assert(vad_inst->elems_in_buf == energy_vad_num_frames_in_buf);

    Uint32 thresh = (total_energy_high/num_frames)<<16;
    thresh += (((total_energy_high%num_frames)<<16) + total_energy_low)/num_frames;
    if (thresh < 10000) thresh = 10000;
    vad_inst->energy_threshold = thresh;
	free(wave_frame);
}

void CalculateVadEnergyThreshold(const Int16* wave,
                                 Int32 wave_len,
                                 Int32 frame_len,
                                 Int32 frame_shift,
                                 SpectrumEnergyVad* vad_inst) {
	Int16 *wave_frame = (Int16*)malloc(frame_len*sizeof(Int16));
    Int32 overlap = frame_len - frame_shift;

    Int32 num_frames = (wave_len - frame_len)/frame_shift + 1;

    memcpy(wave_frame + frame_shift, wave, sizeof(Int16)*overlap);

    for (Int32 i = 0; i < num_frames; ++i) {
        Int32 begin = i*frame_shift + overlap;
        memmove(wave_frame, wave_frame + frame_shift, overlap*sizeof(Int16));
        memcpy(wave_frame + overlap, wave + begin, sizeof(Int16)*frame_shift);

        Int32 measures[spectrum_vad_feature_num];
        ComputeVADMeasure(wave_frame, frame_len, vad_inst, measures);

        for (Int32 f = 0; f < spectrum_vad_feature_num; ++f) {
            vad_inst->thresholds[f] += measures[f];
        }

        //save wave to buffer
        if (i + spectrum_vad_num_frames_in_buf >= num_frames) {
            if (vad_inst->elems_in_buf == 0) {
                memcpy(vad_inst->buf, wave_frame, frame_len*sizeof(Int16));
            } else if (vad_inst->elems_in_buf < spectrum_vad_num_frames_in_buf) {
                Int32 offset = overlap + frame_shift*vad_inst->elems_in_buf;
                memcpy(vad_inst->buf + offset, wave_frame + overlap, frame_shift*sizeof(Int16));
            }
            vad_inst->elems_in_buf++;
        }
    }
    assert(vad_inst->elems_in_buf == spectrum_vad_num_frames_in_buf);
    
    for (Int32 f = 0; f < spectrum_vad_feature_num; ++f) {
        vad_inst->thresholds[f] /= num_frames;
    }
    for (Int32 f = 0; f < spectrum_vad_feature_num; ++f) {
        if (vad_inst->thresholds[f] > vad_inst->voice_mean[f]) {
            continue;
        }
        vad_inst->sil_mean[f] = vad_inst->thresholds[f];
    }
	free(wave_frame);
}

void DetectWaveVADLearnPhase(const Int16* wave,
                              Int32 wave_len,
                              EnergyVad &tempo_vad_inst,
                              SpectrumEnergyVad &spec_vad_inst,
                              Int32 tempo_frame_len,
                              Int32 tempo_frame_shift,
                              Int32 spec_frame_len,
                              Int32 spec_frame_shift,
                              Int32 learn_frames_num) {
    Int32 num_frames = (wave_len - tempo_frame_len)/tempo_frame_shift + 1;

    assert(num_frames >= learn_frames_num);
    assert(tempo_frame_len <= spec_frame_len);

    //tempo vad learning phase
    Int32 tempo_learn_samples = tempo_frame_len + (learn_frames_num-1)*tempo_frame_shift;
    CalculateVadEnergyThreshold(wave, tempo_learn_samples,
            tempo_frame_len, tempo_frame_shift, &tempo_vad_inst);

    //learning phase
    Int32 spec_learn_samples = spec_frame_len + (learn_frames_num-1)*spec_frame_shift;
    CalculateVadEnergyThreshold(wave, spec_learn_samples,
                spec_frame_len, spec_frame_shift, &spec_vad_inst);
}

void DetectWaveVADPredictPhase(EnergyVad &tempo_vad_inst,
                                SpectrumEnergyVad &spec_vad_inst,
                                const Int16* tempo_wave_frame,
                                Int32 tempo_frame_len,
                                Int32 tempo_frame_shift,
                                Int16* spec_wave_frame,
                                Int32 spec_frame_len,
                                Int32 spec_frame_shift,
                                const Int32 cur_frame,
                                Int16* output_voice_buf,
                                Int32* voice_length) {

    const Int32 tempo_overlap = tempo_frame_len - tempo_frame_shift,
                spec_overlap = spec_frame_len - spec_frame_shift;
    { // full buf already
        Int32 offset = energy_vad_buf_size - tempo_frame_shift;
        memmove(tempo_vad_inst.buf, tempo_vad_inst.buf + tempo_frame_shift,
                offset*sizeof(Int16));
        memcpy(tempo_vad_inst.buf + offset, tempo_wave_frame + tempo_overlap,
                tempo_frame_shift*sizeof(Int16));
    }

    *voice_length = 0;
    VadState tempo_last_state = tempo_vad_inst.state;
    if (IsVoice(tempo_wave_frame, tempo_frame_len, cur_frame, &tempo_vad_inst)) {
        //invoke spectrum Energy VAD
        //send wave buffer to spectrum energy VAD
        assert(energy_vad_buf_size == spectrum_vad_buf_size);
        memcpy(spec_vad_inst.buf, tempo_vad_inst.buf,
                energy_vad_buf_size*sizeof(Int16));

        if (tempo_last_state != kInSpeech && tempo_vad_inst.state == kInSpeech) { //first time to wakeup spectrum vad
            const Int32 frames_in_buf = spectrum_vad_num_frames_in_buf;
            bool in_voice = false;
            for (Int32 i = 0; i < frames_in_buf - 1; ++i) {
                Int16 *history_wave_frame = (Int16*)malloc(spec_frame_len*sizeof(Int16));
                const Int32 begin = i*spec_frame_shift,
                            history_frame = cur_frame - (frames_in_buf - i);
                memcpy(history_wave_frame, spec_vad_inst.buf + begin,
                        spec_frame_len*sizeof(Int16));

                //check buffered wave and update VAD threhold
                if (IsVoice(history_wave_frame, spec_frame_len, history_frame, &spec_vad_inst)) {
                    in_voice = true;
                }
				free(history_wave_frame);
            }
            //send wave buffer to asr wake-up system
            if (in_voice) {
                memcpy(output_voice_buf, spec_vad_inst.buf, spectrum_vad_buf_size*sizeof(Int16));
                *voice_length += energy_vad_buf_size;
            }
        }

        { // fill spectrum vad frame
            Int32 offset = spectrum_vad_buf_size - (spec_frame_len - tempo_frame_len),
                  diff = spec_frame_len - tempo_frame_len;
            memcpy(spec_wave_frame, spec_vad_inst.buf + offset, diff*sizeof(Int16));
            memcpy(spec_wave_frame + diff, tempo_wave_frame, tempo_frame_len*sizeof(Int16));
        }

        { // full buf already
            Int32 offset = spectrum_vad_buf_size - spec_frame_shift;
            memmove(spec_vad_inst.buf, spec_vad_inst.buf + spec_frame_shift,
                    offset*sizeof(Int16));
            memcpy(spec_vad_inst.buf + offset, spec_wave_frame + spec_overlap,
                    spec_frame_shift*sizeof(Int16));
        }

        //
        VadState spec_last_state = spec_vad_inst.state;
        if (IsVoice(spec_wave_frame, spec_frame_len, cur_frame, &spec_vad_inst)) {
            //send wave buffer to asr wake-up system
            if (spec_last_state != kInSpeech && spec_vad_inst.state == kInSpeech) { //begin voice segment
                memcpy(output_voice_buf + *voice_length, spec_vad_inst.buf,
                        spectrum_vad_buf_size*sizeof(Int16));
                *voice_length += energy_vad_buf_size;
            } else {
                memcpy(output_voice_buf + *voice_length, spec_wave_frame + spec_overlap,
                        spec_frame_shift*sizeof(Int16));
                *voice_length += spec_frame_shift;
            }
        }
    } else { // if tempo_vad.IsVoice
        ResetVadState(&spec_vad_inst);
    }
}

Int32 DetectWaveVADFinal(const Int16* wave,
                        Int32 wave_len,
                        Int32 tempo_frame_len,
                        Int32 tempo_frame_shift,
                        Int32 spec_frame_len,
                        Int32 spec_frame_shift,
                        Int32 leading_sil_frames) {
    if (wave_len < tempo_frame_len) { // too short
        return 0;
    }
    Int32 num_frames = (wave_len - tempo_frame_len)/tempo_frame_shift + 1;

    assert(tempo_frame_len <= spec_frame_len);

    EnergyVad tempo_vad_inst;
    InitVad(&tempo_vad_inst);
    SpectrumEnergyVad spec_vad_inst;
    InitVad(&spec_vad_inst);

    Int16* highpass_wave = NULL;
    Malloc(&highpass_wave, wave_len);

    HighPassFilterFixedPoint(wave, wave_len, highpass_wave, (tempo_vad_inst.highpass_state));

    DetectWaveVADLearnPhase(wave, wave_len, tempo_vad_inst, spec_vad_inst,
                tempo_frame_len, tempo_frame_shift, spec_frame_len, spec_frame_shift,
                leading_sil_frames);

    //predict phase wave frame buffer
    const Int32 tempo_overlap = tempo_frame_len - tempo_frame_shift;
	Int16 *tempo_wave_frame = (Int16*)malloc(tempo_frame_len*sizeof(Int16));
	Int16 *spec_wave_frame = (Int16*)malloc(spec_frame_len*sizeof(Int16));


    Int32 tempo_learn_samples = tempo_frame_len + (leading_sil_frames-1)*tempo_frame_shift;
    memcpy(tempo_wave_frame + tempo_frame_shift, highpass_wave + tempo_learn_samples,
                    sizeof(Int16)*tempo_overlap);

	Int16 *voice_buf = (Int16*)malloc((energy_vad_buf_size + spec_frame_len)*sizeof(Int16));
    Int32 voice_length = 0;
    for (Int32 f = leading_sil_frames; f < num_frames; ++f) {

        Int32 begin = (f-leading_sil_frames)*tempo_frame_shift + tempo_learn_samples + tempo_overlap;
        memmove(tempo_wave_frame, tempo_wave_frame + tempo_frame_shift,
                                    tempo_overlap*sizeof(Int16));
        memcpy(tempo_wave_frame + tempo_overlap, highpass_wave + begin,
                sizeof(Int16)*tempo_frame_shift);

        voice_length = 0;
        DetectWaveVADPredictPhase(tempo_vad_inst, spec_vad_inst,
                    tempo_wave_frame, tempo_frame_len, tempo_frame_shift,
                    spec_wave_frame, spec_frame_len, spec_frame_shift,
                    f, voice_buf, &voice_length);
        if (voice_length > 0) {
            //invoke asr wakeup system
            std::cout << f << "(" << voice_length << ") ";//<< std::endl;
        }
    }
	free(tempo_wave_frame);
	free(spec_wave_frame);
	free(voice_buf);

    return 0;
}
