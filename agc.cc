#include <stdlib.h>
#include "agc.h"
//#include "vad.h"
//generator gain table
//; MinGain = 0; CompRatio = 3; Knee = 1;
const static int genGainTable[32] = 
{
	//MaxGain = 15
	52016, 65536, 82570, 104032,
	131070, 165123, 207897, 260622,
	318099, 356074, 366600, 368270,
	368500, 368531, 368535, 368536,
	368536, 368536, 368536, 368536,
	368536, 368536, 368536, 368536,
	368536, 368536, 368536, 368536,
	368536, 368536, 368536, 368536
};

//2^n~(2^(n+1)-1)-->zeros=31-n
static int DigitalAgc_NormU32(unsigned int value)
{
	int zeros = 0;

	if (0 == value)
		return 0;

	if (!(0xFFFF0000 & value))
		zeros = 16;
	if (!(0xFF000000 & (value << zeros)))
		zeros += 8;
	if (!(0xF0000000 & (value << zeros)))
		zeros += 4;
	if (!(0xC0000000 & (value << zeros)))
		zeros += 2;
	if (!(0x80000000 & (value << zeros)))
		zeros += 1;

	return zeros;
}

int DigitalAgc_Init(DigitalAgc *agc_inst)
{
	InitVad(&(agc_inst->vad_inst));
	agc_inst->gain_table          = genGainTable;
	agc_inst->gain                = 65536;//Q16
	agc_inst->count_load_frames = 0;
	agc_inst->learn_buff = (short*)malloc(sizeof(short)*(num_learn_frame*frame_shift));
	agc_inst->is_first_speech_frame = 1;
	agc_inst->pre_isVoice = 0;
	agc_inst->capacitorFast = 0;
	agc_inst->capacitorSlow = 0;
	memset(agc_inst->learn_buff, 0, sizeof(short)*(num_learn_frame*frame_shift));
	memset(agc_inst->gains, 0, sizeof(int)*(num_sub_frames + 1));
	memset(agc_inst->frame_buff, 0, sizeof(short)*frame_len);
	
	return 0;
}

/*
every frames(10ms) will separate into 10 sub_frames(1ms).
only suit for fs=16000Hz, bit per samples is 16
audio_data is the input audio data, and the output data which after agc will store in audio_data
*/
int DigitalAgc_Process(DigitalAgc *agc_inst, short *audio_data, unsigned int data_len)
{

	if (data_len > num_sub_frames*samples_sub_frames)
	{
		//数据不能超过160个采样点
		return -1;
	}

	// find max amplitude per sub frame
	// iterate over sub frames
	int max_amplitude_per_sub_frames[num_sub_frames] = { 0 };
	int max_amplitude;
	int temp;
	int decay = -65;
	for (int i = 0; i < num_sub_frames; i++)
	{
		max_amplitude = 0;
		// iterate over samples
		for (int j = 0; j < samples_sub_frames; j++)
		{
			temp = (int)(audio_data[i*samples_sub_frames + j] * audio_data[i*samples_sub_frames + j]);
			if (temp > max_amplitude)
			{
				max_amplitude = temp;
			}
		}
		max_amplitude_per_sub_frames[i] = max_amplitude;
	}

	// calculate gain per sub frame
	// Translate signal level into gain, using a piecewise linear approximation
	// find number of leading zeros
	agc_inst->gains[0] = agc_inst->gain;
	for (int i = 0; i < 10; i++)
	{
		// Fast envelope follower
		//  decay time = -131000 / -1000 = 131 (ms)
		agc_inst->capacitorFast = agc_inst->capacitorFast + ((agc_inst->capacitorFast) >> 16)*(-1000) +
			(((0x0000FFFF & (agc_inst->capacitorFast))*(-1000)) >> 16);
		if (max_amplitude_per_sub_frames[i] > agc_inst->capacitorFast)
		{
			agc_inst->capacitorFast = max_amplitude_per_sub_frames[i];
		}
		// Slow envelope follower
		if (max_amplitude_per_sub_frames[i] > agc_inst->capacitorSlow)
		{
			// increase capacitorSlow
			agc_inst->capacitorSlow = agc_inst->capacitorSlow + ((max_amplitude_per_sub_frames[i] - agc_inst->capacitorSlow) >> 16)*(500) +
				(((0x0000FFFF & (max_amplitude_per_sub_frames[i] - agc_inst->capacitorSlow))*(500)) >> 16);
		}
		else
		{
			// decrease capacitorSlow
			agc_inst->capacitorSlow = agc_inst->capacitorSlow + ((agc_inst->capacitorSlow) >> 16)*(decay)+
				(((0x0000FFFF & (agc_inst->capacitorSlow))*(decay)) >> 16);
		}

		// use maximum of both capacitors as current level
		int cur_level;
		if (agc_inst->capacitorFast > agc_inst->capacitorSlow)
		{
			cur_level = agc_inst->capacitorFast;
		}
		else
		{
			cur_level = agc_inst->capacitorSlow;
		}
		// Translate signal level into gain, using a piecewise linear approximation
		// find number of leading zeros
		int zeros = DigitalAgc_NormU32((unsigned int)cur_level);
		if (cur_level == 0)
		{
			zeros = 31;
		}
		//利用线性插值获取当前的增益
		int temp32 = (cur_level << zeros) & 0x7FFFFFFF;
		short frac = (short)(temp32 >> 19); //Q12
		temp32 = (agc_inst->gain_table[zeros - 1] - agc_inst->gain_table[zeros]) * frac;
		agc_inst->gains[i + 1] = agc_inst->gain_table[zeros] + (temp32 >> 12);
	}

	// Limit gain to avoid overload distortion
	int gain32;
	for (int i = 0; i < 10; i++)
	{
		// To prevent wrap around
		int zeros = 10;
		if (agc_inst->gains[i + 1] > 47453132)
		{
			zeros = 16 - DigitalAgc_NormU32((unsigned int)agc_inst->gains[i + 1]);
		}
		gain32 = (agc_inst->gains[i + 1] >> zeros) + 1;
		gain32 = gain32 * gain32;
		// check for overflow
		while ((((gain32) >> 13)*(((max_amplitude_per_sub_frames[i] >> 12) + 1)) + (((0x00001FFF & (gain32))*(((max_amplitude_per_sub_frames[i] >> 12) + 1))) >> 13))
				> ((2 * (1 - zeros + 10) >= 0) ? ((int)(32767) << (2 * (1 - zeros + 10))) : ((int)(32767) >> (-(2 * (1 - zeros + 10))))))
		{
			// multiply by 253/256 ==> -0.1 dB
			if (agc_inst->gains[i + 1] > 8388607)
			{
				// Prevent wrap around
				agc_inst->gains[i + 1] = (agc_inst->gains[i + 1] >> 8) * 253;
			}
			else
			{
				agc_inst->gains[i + 1] = (agc_inst->gains[i + 1] * 253) >> 8;
			}
			gain32 = (agc_inst->gains[i + 1] >> zeros) + 1;
			gain32 = gain32 * gain32;
		}
	}
	
	// gain reductions should be done 1 ms earlier than gain increases. ????
	for (int i = 1; i < num_sub_frames; i++)
	{
		if (agc_inst->gains[i] > agc_inst->gains[i + 1])
		{
			agc_inst->gains[i] = agc_inst->gains[i + 1];
		}
	}

	// save start gain for next frame
	agc_inst->gain = agc_inst->gains[num_sub_frames];

	// apply gain
	// handle first sub frame separately
	int delta = agc_inst->gains[1] - agc_inst->gains[0];
	gain32 = (agc_inst->gains[0] << 4);;
	// iterate over samples
	for (int i = 0; i < samples_sub_frames; i++)
	{
		//????
		int temp32 = (int)audio_data[i] * ((gain32 + 127) >> 7);
		int out_temp = (temp32 >> 16);
		if (out_temp > 4095)
		{
			audio_data[i] = 32767;
		}
		else if (out_temp < -4096)
		{
			audio_data[i] = -32767;
		}
		else
		{
			temp32 = (int)audio_data[i] * (gain32 >> 4);
			audio_data[i] = (short)(temp32 >> 16);
		}
		gain32 += delta;
	}

	// iterate over subframes
	for (int i = 1; i < num_sub_frames; i++)
	{
		delta = agc_inst->gains[i + 1] - agc_inst->gains[i];
		gain32 = (agc_inst->gains[i] << 4);
		// iterate over samples
		for (int j = 0; j < samples_sub_frames; j++)
		{
			int temp32 = (int)audio_data[i * samples_sub_frames + j] * (gain32 >> 4);
			audio_data[i * samples_sub_frames + j] = (short)(temp32 >> 16);
			gain32 += delta;
		}
	}

	return 0;

}

int DigitalAgc_Process_Use_VAD(DigitalAgc *agc_inst, short *audio_data, unsigned int data_len, int &vad_state)
{
	//前20帧用来估计噪声
	if (agc_inst->count_load_frames < num_learn_frame)
	{
		memcpy(agc_inst->learn_buff, audio_data, sizeof(short)*frame_shift);
		agc_inst->learn_buff += frame_shift;
		agc_inst->count_load_frames++;
	}
	else
	{
		if (agc_inst->count_load_frames == num_learn_frame)
		{
			HighPassFilterFixedPoint(agc_inst->learn_buff-frame_overlap, frame_overlap, agc_inst->frame_buff + frame_shift, (agc_inst->vad_inst.highpass_state));
			agc_inst->learn_buff = agc_inst->learn_buff - (num_learn_frame*frame_shift);
			CalculateVadEnergyThreshold(agc_inst->learn_buff, num_learn_frame*frame_shift, frame_len, frame_shift, &agc_inst->vad_inst);

			//HighPassFilterFixedPoint(audio_data, frame_shift, agc_inst->frame_buff + frame_shift, (agc_inst->vad_inst.highpass_state));
			agc_inst->count_load_frames++;
			free(agc_inst->learn_buff);
		}

		memmove(agc_inst->frame_buff, agc_inst->frame_buff + frame_shift, frame_overlap*sizeof(short));

		HighPassFilterFixedPoint(audio_data, frame_shift, agc_inst->frame_buff + frame_overlap, (agc_inst->vad_inst.highpass_state));

		{ // has been full buf
			Int32 offset = energy_vad_buf_size - frame_shift;
			memmove(agc_inst->vad_inst.buf, agc_inst->vad_inst.buf + frame_shift, offset*sizeof(short));
			memcpy(agc_inst->vad_inst.buf + offset, agc_inst->frame_buff + frame_overlap, frame_shift*sizeof(short));
		}

		if (IsVoice(agc_inst->frame_buff, frame_len, 0, &agc_inst->vad_inst)) {
			if (0 == agc_inst->pre_isVoice)
			{
				agc_inst->is_first_speech_frame = 1;
			}
			else
			{
				agc_inst->is_first_speech_frame = 0;
			}
			vad_state = agc_inst->is_first_speech_frame;
			agc_inst->pre_isVoice = 1;
			
			if (agc_inst->is_first_speech_frame)
			{
				for (int i = 0; i < energy_vad_num_frames_in_buf; i++)
				{
					//DigitalAgc_Process(agc_inst, agc_inst->vad_inst.buf + i*frame_shift, frame_shift);//xxx,
					DigitalAgc_Process(agc_inst, agc_inst->vad_inst.buf + frame_overlap + i*frame_shift, frame_shift);//xxx,
				}
			}
			else
			{
				//memcpy(audio_data, agc_inst->frame_buff, frame_shift*sizeof(short));
				memcpy(audio_data, agc_inst->frame_buff + frame_overlap, frame_shift*sizeof(short));
				DigitalAgc_Process(agc_inst, audio_data, data_len);//xxx,
			}
			
		}
		else
		{
			agc_inst->pre_isVoice = 0;
		}
	}
	//
	return 0;
}

