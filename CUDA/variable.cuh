#include "common.cuh"
#include"grid.cuh"
#include "particle.cuh"
int step = 20000;//2016_11_23���Ӳ���
int lizi_count = step * 30;
//DEBUG�޴� RELEASE�д�
//�����ڴ� DEBUG�޴� �ȶ�
//��������� DEBUG�ȶ�
//�ش�BUG���£���ⳡ��tidȫ��д����tid=gridIdx.x*blockDim.x+threadIdx.x�޸�Ϊ
//tid=blockIdx.x*blockDim.x+threadIdx.x�޸���release�ȶ�
//***********device**********
Grid *device_G = NULL;
Grid *device_Gn = NULL;
Grid *dev_Gtemp = NULL;
float *d_rds = NULL, *d_rds1 = NULL;
float *d_stac_Bx = NULL, *d_stac_Bz = NULL, *d_stac_ex = NULL, *d_stac_ez = NULL;
Paticle *d_pat_ion = NULL;
Pre_Paticle *d_pre_ion = NULL;
Paticle *d_pat_elc = NULL;
Pre_Paticle *d_pre_elc = NULL;
//float *d_ex1 = NULL, *d_iex1 = NULL, *d_ey1 = NULL, *d_iey1 = NULL, *d_ez1 = NULL, *d_iez1 = NULL,
//*d_hx1 = NULL, *d_ihx1 = NULL, *d_hy1 = NULL, *d_ihy1 = NULL, *d_hz1 = NULL, *d_ihz1 = NULL,
float *d_sigmaz1 = NULL, *d_sigmaz = NULL;
float *DS_number = NULL;

//************Host***********
Paticle *pat_elc = NULL, *pat_ion = NULL;
Paticle wuchafenxi[30];//������������ɾ��
int paticle_count = 0;
int host_temptail = 0;
float *HS_number = NULL;
Grid *G_gpu;