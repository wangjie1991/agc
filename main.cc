#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "agc.h"

//���ַ���str1�����str2����ӽ��������res��
void addStr(char* res, char *str1, char *str2)
{
	int res_index;
	int str1_index;
	for (str1_index = 0, res_index = 0; str1[str1_index] != 0; str1_index++, res_index++)
	{
		res[res_index] = str1[str1_index];
	}

	int str2_index;
	for (str2_index = 0; str2[str2_index] != 0; str2_index++, res_index++)
	{
		res[res_index] = str2[str2_index];
	}
	res[res_index] = 0;
}

//���ַ���str1+str2+str3����ӽ��������res��
void addStr(char* res, char *str1, char *str2, char *str3)
{
	int res_index;
	int str1_index;
	for (str1_index = 0, res_index = 0; str1[str1_index] != 0; str1_index++, res_index++)
	{
		res[res_index] = str1[str1_index];
	}

	int str2_index;
	for (str2_index = 0; str2[str2_index] != 0; str2_index++, res_index++)
	{
		res[res_index] = str2[str2_index];
	}

	int str3_index;
	for (str3_index = 0; str3[str3_index] != 0; str3_index++, res_index++)
	{
		res[res_index] = str3[str3_index];
	}

	res[res_index] = 0;
}

const int samples_per_frame = 160; //16000Hz, 10ms

int main(int argc, char* argv[])
{
	FILE *fp_src, *fp_dst;
	int actual_len;
	//�����ļ�·����
	DigitalAgc agc_inst;
	DigitalAgc_Init(&agc_inst);

  if (argc != 3)
  {
    printf("Usage:%s input.wav output.wav\n", argv[0]);
    return -1;
  }

	fp_src = fopen(argv[1], "rb");
	if (!fp_src)
	{
		printf("failed to open src wav file, exit...\n");
		exit(1);
	}

	fp_dst = fopen(argv[2], "wb");
	if (!fp_dst)
	{
		printf("failed to open dst wav file, exit...\n");
		exit(1);
	}

	//�����ڴ�ռ�
	//short *src_all;
	int readLen = sizeof(short)*samples_per_frame;
	short *src_all = (short*)malloc(readLen);
	//short src_all[160];
	short *temp = (short*)malloc(sizeof(short)*(energy_vad_buf_size - frame_overlap));
	memset(temp, 0, sizeof(short)*(energy_vad_buf_size - frame_overlap));
	//���ζ�ȡpcm���ݡ��������
	//while (feof(fp_src)==0) //feof����ֵ0�����ļ�û�������ļ��������ط�0
	short data[3360] = { 0 };
	while (1)
	{
		//���ļ���ȡ����
		actual_len = fread(src_all, 1, readLen, fp_src);
		if (actual_len == 0)
		{
			//�����Ѿ���ȡ���
			break;
		}
		int vad_state;
		DigitalAgc_Process_Use_VAD(&agc_inst, src_all, samples_per_frame, vad_state);
		if (1 == vad_state)
		{
			fseek(fp_dst, -sizeof(short)*(energy_vad_buf_size - frame_shift), 1);
			fwrite(agc_inst.vad_inst.buf + frame_overlap, 1, sizeof(short)*(energy_vad_buf_size - frame_overlap), fp_dst);
		}
		else
		{
			fwrite(src_all, 1, actual_len, fp_dst);
		}
	}

	//�ͷſռ�
	free(src_all);
	//�ر��ļ�
	fclose(fp_src);
	fclose(fp_dst);

	return 0;
}

