#ifndef _DIGITAL_AGC_H_
#define _DIGITAL_AGC_H_
#include "vad.h"

const static int num_learn_frame = 21;
const static int frame_len = energy_vad_frame_length;          //fs=16000, 15ms
const static int frame_shift = energy_vad_frame_shift;         //10ms
const static int frame_overlap = frame_len - frame_shift;
const static int num_sub_frames = 10;      //每帧分成10小帧，即1ms
const static int samples_sub_frames = 16;  //每小帧有16个采样点
const static int lookback_len = sil_to_speech_hangover*frame_shift + frame_overlap;

typedef struct DigitalAgc_
{
	int        gain;
	int        gains[num_sub_frames + 1];
	const int  *gain_table;  //Q16
	EnergyVad  vad_inst;
	int        count_load_frames;
	short      frame_buff[frame_len];
	short      *learn_buff;    //缓存前20帧用作噪声估计
	int        is_first_speech_frame;//判断是否是语音帧的第一帧，1是，0否
	int        pre_isVoice;                //保存前1帧是否为语音帧, 1是语音帧，0不是
	int        capacitorSlow;
	int        capacitorFast;
} DigitalAgc;

//audio_data:帧移数据   frame_buff：一整帧数据
int DigitalAgc_Init(DigitalAgc *agc_inst);
int DigitalAgc_Process(DigitalAgc *agc_inst, short *audio_data, unsigned int data_len);
int DigitalAgc_Process_Use_VAD(DigitalAgc *agc_inst, short *audio_data, unsigned int data_len, int &vad_state);

#endif