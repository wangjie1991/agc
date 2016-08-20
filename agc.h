#ifndef _DIGITAL_AGC_H_
#define _DIGITAL_AGC_H_
#include "vad.h"

const static int num_learn_frame = 21;
const static int frame_len = energy_vad_frame_length;          //fs=16000, 15ms
const static int frame_shift = energy_vad_frame_shift;         //10ms
const static int frame_overlap = frame_len - frame_shift;
const static int num_sub_frames = 10;      //ÿ֡�ֳ�10С֡����1ms
const static int samples_sub_frames = 16;  //ÿС֡��16��������
const static int lookback_len = sil_to_speech_hangover*frame_shift + frame_overlap;

typedef struct DigitalAgc_
{
	int        gain;
	int        gains[num_sub_frames + 1];
	const int  *gain_table;  //Q16
	EnergyVad  vad_inst;
	int        count_load_frames;
	short      frame_buff[frame_len];
	short      *learn_buff;    //����ǰ20֡������������
	int        is_first_speech_frame;//�ж��Ƿ�������֡�ĵ�һ֡��1�ǣ�0��
	int        pre_isVoice;                //����ǰ1֡�Ƿ�Ϊ����֡, 1������֡��0����
	int        capacitorSlow;
	int        capacitorFast;
} DigitalAgc;

//audio_data:֡������   frame_buff��һ��֡����
int DigitalAgc_Init(DigitalAgc *agc_inst);
int DigitalAgc_Process(DigitalAgc *agc_inst, short *audio_data, unsigned int data_len);
int DigitalAgc_Process_Use_VAD(DigitalAgc *agc_inst, short *audio_data, unsigned int data_len, int &vad_state);

#endif