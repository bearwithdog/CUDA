#include "lty_paralle_common.cuh"
//#include "common.cuh"
#include <stdio.h>
#include<iostream>
#include<time.h>
#include"device_functions.h"
#include "device_launch_parameters.h"
#include <math.h>
#include "sm_20_atomic_functions.hpp"
using namespace std;
extern int Ne,int nxx,int nzz,int kongxue,int tail,int enlarge;
extern  float dr,float dz,float ca_posz,float KTe;
//extern float *d_ex1,float *d_iex1,float *d_ey1,float *d_iey1,float *d_ez1,float *d_iez1,
//	   float *d_hx1,float *d_ihx1,float *d_hy1,float *d_ihy1,float *d_hz1,float *d_ihz1,
//	   float *d_sigmaz1,float *d_sigmaz,float *sigmaz1,float *sigmaz;
extern Paticle *pat_elc;
__device__ float d_ex1[51 * 101], d_iex1[51 * 101], d_ey1[51 * 101], d_iey1[51 * 101], d_ez1[51 * 101], d_iez1[51 * 101],
d_hx1[51 * 101], d_ihx1[51 * 101], d_hy1[51 * 101], d_ihy1[51 * 101], d_hz1[51 * 101], d_ihz1[51 * 101];
//extern float *d_sigmaz1 = NULL, *d_sigmaz = NULL;
__constant__ float D_parameter[12]={1.0e-003,1.0e-003,5.0e-002,1.0e-001,1.66782e-12,-5.5594e-15,5.5594e-17,0.005,0.03,
									4.095e-16,2.4033e-20,9.1e-31};
//0dr 1dz 2R 3L 4dt 5qe 6qi 7ca_posr 8ca_posz 9KTe 10KTi 11Me /*3.2044e-19*/
__global__ void device_initialchang(Grid* device_g,Grid* device_gn)//x是行 y是列
{
	int x=51,y=101;
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	while(tid<x*y)
	{
		(device_g+tid)->ex=0;
		(device_g+tid)->ey=0;
		(device_g+tid)->ez=0;
		(device_g+tid)->hx=0;
		(device_g+tid)->hy=0;
		(device_g+tid)->hz=0;

		(device_g+tid)->ave_ex=0;
		(device_g+tid)->ave_ey=0;
		(device_g+tid)->ave_ez=0;
		(device_g+tid)->ave_hx=0;
		(device_g+tid)->ave_hy=0;
		(device_g+tid)->ave_hz=0;

		(device_g+tid)->ne[0]=0;
		(device_g+tid)->ne[1]=0;
		(device_g+tid)->jr=0.0;
		(device_g+tid)->jz=0.0;
		(device_g+tid)->jy=0.0;
		(device_g+tid)->jr_ion=0.0;
		(device_g+tid)->jz_ion=0.0;
		(device_g+tid)->jy_ion=0.0;
		(device_g+tid)->Pmax=0.0;

		(device_gn+tid)->ex=0;
		(device_gn+tid)->ey=0;
		(device_gn+tid)->ez=0;
		(device_gn+tid)->hx=0;
		(device_gn+tid)->hy=0;
		(device_gn+tid)->hz=0;
		(device_gn+tid)->ne[2]=0;
		tid+=gridDim.x*blockDim.x;
	}
	//printf("success,device_initial_chang");
	
}

//void device_initial_always(float *device_rds,float *device_rds1,Paticle *p_pat_elc,Grid *device_G,float ca_posz,int Ne,int nx,int nz,float dr,float dz,float KTe)
__global__ void initial_always(Paticle *pat_elc,int n,float *S_number,int tmp)//n是step*30
{
	
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	while(tid<n)
	{
				int i=tid%tmp;
				float afa=0;//rds5*(pi/12);
				float cita=0;//rds6*pi/2;
				float vv = 1e6;
				pat_elc[tid].pz=D_parameter[8];
				pat_elc[tid].pr =0.04*S_number[i]+0.001;//发射位置垫高
				pat_elc[tid].py=2*pi*S_number[i];
				pat_elc[tid].vr = vv*sin(afa)*cos(cita);
				pat_elc[tid].vy = vv*sin(afa)*sin(cita);
				pat_elc[tid].vz=vv;
				//pat_elc[tid].blei = pat_elc[tid].pr/D_parameter[0];
//				pat_elc[tid].blek = pat_elc[tid].pz/D_parameter[1];
				tid+=gridDim.x*blockDim.x;
     }	
	
}

__global__ void device_ave_field(int x,int y,Grid* device_g,Grid* device_gn) //x是行 y是列
{
	int u=0,v=0;
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	while(tid<x*(y))
	{

		//bool bv=tid%y;//v;
		//bool bu=tid/(y);//u;
	/*G[u][v].ave_ex=(bv*G[u][v].ex+bu*G[u-1+1-bu][v].ex+(1-bv)*(1-bu)*G[u][v].ex+(1-bv)*bu*G[u][v].ex)/(bv+bu+(1-bv)*(1-bu)+(1-bv)*bu);
	G[u][v].ave_ez=(G[u][v].ez+bv*G[u][v-1+1-bv].ez)/(2-(1-bv));
	G[u][v].ave_hy=(G[u][v].hy+bu*G[u-1+1-bu][v].hy+bv*G[u][v-1+1-bv].hy+bu*bv*G[u-1+1-bu][v-1+1-bv].hy
					+Gn[u][v].hy+bu*Gn[u-1+1-bu][v].hy+bv*Gn[u][v-1+1-bv].hy+bu*bv*Gn[u-1+1-bu][v-1+1-bv].hy)/(8/(3-bu-bv+(1-bu)*(1-bv)));

	G[u][v].ave_ey=G[u][v].ey;
	G[u][v].ave_hx=(G[u][v].hx+Gn[u][v].hx+bv*G[u][v-1+1-bv].hx+bv*Gn[u][v-1+1-bv].hx)/(2+2*bv);
	G[u][v].ave_hz=(G[u][v].hz+Gn[u][v].hz+bu*G[u-1+1-bu][v].hz+bu*Gn[u-1+1-bu][v].hz)/(2+2*bu);*/
	/*device_g[tid].ave_ex=(bv*device_g[tid].ex+bu*device_g[tid-bu*y].ex+(1-bv)*(1-bu)*device_g[tid].ex+(1-bv)*bu*device_g[tid].ex)/(bv+bu+(1-bv)*(1-bu)+(1-bv)*bu);
	device_g[tid].ave_ez=(device_g[tid].ez+bv*device_g[tid-bv].ez)/(2-(1-bv));
	device_g[tid].ave_hy=(device_g[tid].hy+bu*device_g[tid-bu*y].hy+bv*device_g[tid-bv].hy+bu*bv*device_g[tid-bu*y-bv].hy
					+device_gn[tid].hy+bu*device_gn[tid-bu*y].hy+bv*device_gn[tid-bv].hy+bu*bv*device_gn[tid-bu*y-bv].hy)/(8/(3-bu-bv+(1-bu)*(1-bv)));

	device_g[tid].ave_ey=device_g[tid].ey;
	device_g[tid].ave_hx=(device_g[tid].hx+device_gn[tid].hx+bv*device_g[tid-bv].hx+bv*device_gn[tid-bv].hx)/(2+2*bv);
	device_g[tid].ave_hz=(device_g[tid].hz+device_gn[tid].hz+bu*device_g[tid-bu*y].hz+bu*device_gn[tid-bu*y].hz)/(2+2*bu);*/
		u=tid/(y);
		v=tid%y;
		int tid_temp=u*(y)+v;
		if(u==0&&v!=0)
				{
				(*(device_g+tid_temp)).ave_ex=((*(device_g+tid_temp)).ex);
				(*(device_g+tid_temp)).ave_ez=((*(device_g+tid_temp)).ez+(*(device_g+tid_temp-1)).ez)/2;
				(*(device_g+tid_temp)).ave_hy=((*(device_g+tid_temp)).hy+(*(device_g+tid_temp-1)).hy+(*(device_gn+tid_temp)).hy+(*(device_gn+tid_temp-1)).hy)/4;

				(*(device_g+tid_temp)).ave_ey=(*(device_g+tid_temp)).ey;
				(*(device_g+tid_temp)).ave_hx=((*(device_g+tid_temp)).hx+(*(device_gn+tid_temp)).hx+(*(device_g+tid_temp-1)).hx+(*(device_gn+tid_temp-1)).hx)/4;
				(*(device_g+tid_temp)).ave_hz=((*(device_g+tid_temp)).hz+(*(device_gn+tid_temp)).hz)/2;
				}
		else if(v==0&&u!=0)
			 {
				(*(device_g+tid_temp)).ave_ex=((*(device_g+tid_temp)).ex+(*(device_g+tid_temp-y)).ex)/2;
				(*(device_g+tid_temp)).ave_ez=((*(device_g+tid_temp)).ez);
				(*(device_g+tid_temp)).ave_hy=((*(device_g+tid_temp)).hy+(*(device_g+tid_temp-y)).hy+(*(device_gn+tid_temp)).hy+(*(device_gn+tid_temp-y)).hy)/4;
				(*(device_g+tid_temp)).ave_ey=(*(device_g+tid_temp)).ey;
				(*(device_g+tid_temp)).ave_hx=((*(device_g+tid_temp)).hx+(*(device_gn+tid_temp)).hx)/2;
				(*(device_g+tid_temp)).ave_hz=((*(device_g+tid_temp)).hz+(*(device_gn+tid_temp)).hz+(*(device_g+tid_temp-y)).hz+(*(device_gn+tid_temp-y)).hz)/4;
			 }
		else if(v==0&&u==0)
			 {
				(*(device_g+tid_temp)).ave_ex=(*(device_g+tid_temp)).ex;
				(*(device_g+tid_temp)).ave_ez=(*(device_g+tid_temp)).ez;
				(*(device_g+tid_temp)).ave_hy=((*(device_g+tid_temp)).hy+(*(device_gn+tid_temp)).hy)/2;
				(*(device_g+tid_temp)).ave_ey=(*(device_g+tid_temp)).ey;
				(*(device_g+tid_temp)).ave_hx=((*(device_g+tid_temp)).hx+(*(device_gn+tid_temp)).hx)/2;
				(*(device_g+tid_temp)).ave_hz=((*(device_g+tid_temp)).hz+(*(device_gn+tid_temp)).hz)/2;
			 }
		 else{
				(*(device_g+tid_temp)).ave_ex=((*(device_g+tid_temp)).ex+(*(device_g+tid_temp-y)).ex)/2;
				(*(device_g+tid_temp)).ave_ez=((*(device_g+tid_temp)).ez+(*(device_g+tid_temp-1)).ez)/2;
				(*(device_g+tid_temp)).ave_hy=((*(device_g+tid_temp)).hy+(*(device_g+tid_temp-y)).hy+(*(device_g+tid_temp-1)).hy+(*(device_g+tid_temp-y-1)).hy
					+(*(device_gn+tid_temp)).hy+(*(device_gn+tid_temp-y)).hy+(*(device_gn+tid_temp-1)).hy+(*(device_gn+tid_temp-y-1)).hy)/8;
				(*(device_g+tid_temp)).ave_ey=(*(device_g+tid_temp)).ey;
				(*(device_g+tid_temp)).ave_hx=((*(device_g+tid_temp)).hx+(*(device_gn+tid_temp)).hx+(*(device_g+tid_temp-1)).hx+(*(device_gn+tid_temp-1)).hx)/4;
				(*(device_g+tid_temp)).ave_hz=((*(device_g+tid_temp)).hz+(*(device_gn+tid_temp)).hz+(*(device_g+tid_temp-y)).hz+(*(device_gn+tid_temp-y)).hz)/4;
				 }
		   tid+=gridDim.x*blockDim.x;
		  // printf("success,device_initial_chang");
	}		
	
}
__device__ float d_min_1(float x,float y)
{
	return(x<y?x:y);
}
__device__ float d_max_1(float x,float y)
{
	return(x>y?x:y); 
}

__global__ void device_update_last(float *device_static_magneticX, float *device_static_magneticZ, float *device_static_electricX, float *device_static_electricZ, Paticle *p_pat_elc, Pre_Paticle *d_pre_elc, Grid *device_G, float *S_number, int tail)
{
	//printf("%d",tail);
	float  Qm_ion = 7.33945e+5;
	float  Qm = -1.7588e+11;
	float mur = 4.0*pi*1.0e-7;
	float E[3] = { 0, 0, 0 }, B[3] = { 0, 0, 0 };
	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	int threadtid = threadIdx.x;
	//Paticle p_prepat_elc;
	__shared__ float elcVr[thread]; __shared__ float elcVy[thread]; __shared__ float elcVz[thread];
	__shared__ float elcPr[thread]; __shared__ float elcPy[thread]; __shared__ float elcPz[thread];
	__shared__ float PrelcPr[thread];  __shared__ float  PrelcPz[thread];
	while (tid<tail)
	{
		elcVr[threadIdx.x] = p_pat_elc[tid].vr;
		elcVy[threadIdx.x] = p_pat_elc[tid].vy;
		elcVz[threadIdx.x] = p_pat_elc[tid].vz;
		elcPr[threadIdx.x] = p_pat_elc[tid].pr;
		elcPy[threadIdx.x] = p_pat_elc[tid].py;
		elcPz[threadIdx.x] = p_pat_elc[tid].pz;
		
		d_pre_elc[tid].pr = elcPr[threadIdx.x];//
		d_pre_elc[tid].pz = elcPz[threadIdx.x];//
		int ii = (int)(elcPr[threadIdx.x] / D_parameter[0]);               //云方程
		int kk = (int)(elcPz[threadIdx.x] / D_parameter[1]);
		float wr = (elcPr[threadIdx.x] / D_parameter[0]) - ii;
		float wz = (elcPz[threadIdx.x] / D_parameter[1]) - kk;
		float s1 = (1 - wr)*(1 - wz);
		float s2 = (1 - wr)*wz;
		float s3 = wr*(1 - wz);
		float s4 = wr*wz;
		int grid_temp1 = ii*(nz + 1) + kk;//ii kk
		int grid_temp2 = ii*(nz + 1) + kk + 1;//ii kk+1
		int grid_temp3 = (ii + 1)*(nz + 1) + kk;//ii+1 kk
		int grid_temp4 = (ii + 1)*(nz + 1) + kk + 1;//ii+1 kk+1
		//printf("u:%d\t", tid);
		E[0] = ((*(device_G + grid_temp1)).ave_ex*s1 + (*(device_G + grid_temp3)).ave_ex*s3 + (*(device_G + grid_temp2)).ave_ex*s2 + (*(device_G + grid_temp4)).ave_ex*s4)
			+((*(device_static_electricX+grid_temp1))*s1+(*(device_static_electricX+grid_temp3))*s3+(*(device_static_electricX+grid_temp2))*s2+(*(device_static_electricX+grid_temp4))*s4);//+stac_ex[ii][kk];
		//printf("u:%d\tE0=%f\n", tid,E[0]);
		E[2] = ((*(device_G + grid_temp1)).ave_ez*s1 + (*(device_G + grid_temp3)).ave_ez*s3 + (*(device_G + grid_temp2)).ave_ez*s2 + (*(device_G + grid_temp4)).ave_ez*s4)
			+((*(device_static_electricZ+grid_temp1))*s1+(*(device_static_electricZ+grid_temp3))*s3+(*(device_static_electricZ+grid_temp2))*s2+(*(device_static_electricZ+grid_temp4))*s4);//+(*(device_static_electricZ+grid_temp1));
		//printf("u:%d\tE2=%f\n", tid,E[2]);
		B[1] = (((*(device_G + grid_temp1)).ave_hy*s1 + (*(device_G + grid_temp3)).ave_hy*s3 + (*(device_G + grid_temp2)).ave_hy*s2 + (*(device_G + grid_temp4)).ave_hy*s4)*mur);
		//printf("u:%d\tB1=%f\n", tid, B[1]);
		E[1] = ((*(device_G + grid_temp1)).ave_ey*s1 + (*(device_G + grid_temp3)).ave_ey*s3 + (*(device_G + grid_temp2)).ave_ey*s2 + (*(device_G + grid_temp4)).ave_ey*s4);
		//printf("u:%d\tE1=%f\n", tid,E[1]);
		B[0] = (((*(device_G + grid_temp1)).ave_hx*s1 + (*(device_G + grid_temp3)).ave_hx*s3 + (*(device_G + grid_temp2)).ave_hx*s2 + (*(device_G + grid_temp4)).ave_hx*s4)*mur)
		+(*(device_static_magneticX+grid_temp1))*s1+(*(device_static_magneticX+grid_temp3))*s3+((*(device_static_magneticX+grid_temp2))*s2+(*(device_static_magneticX+grid_temp4))*s4);//+(*(device_static_magneticX+grid_temp3));
		//B[0]=((G[ii][kk].ave_hx*s1+G[ii+1][kk].ave_hx*s3+G[ii][kk+1].ave_hx*s2+G[ii+1][kk+1].ave_hx*s4)*mur)
		//+(stac_Bx[ii][kk] * s1 + stac_Bx[ii + 1][kk] * s3 + stac_Bx[ii][kk + 1] * s2 + stac_Bx[ii + 1][kk + 1] * s4);
		//printf("u:%d\tB0=%f\n", tid, B[0]);
		B[2] = (((*(device_G + grid_temp1)).ave_hz*s1 + (*(device_G + grid_temp3)).ave_hz*s3 + (*(device_G + grid_temp2)).ave_hz*s2 + (*(device_G + grid_temp4)).ave_hz*s4)*mur)
		+((*(device_static_magneticZ+grid_temp1))*s1+(*(device_static_magneticZ+grid_temp3))*s3+(*(device_static_magneticZ+grid_temp2))*s2+(*(device_static_magneticZ+grid_temp4))*s4);//+(*(device_static_magneticZ+grid_temp3));
		//printf("B2=%f\n", B[2]);

		float u1[3] = { 0 }, u2[3] = { 0 }, u3[3] = { 0 };//u_n-1/2,u-,,u+,u_n+1/2    //    electron 分步求解
		float t[3] = { 0 }, s[3] = { 0 };
		float pp[3][3] = { 0 };

		u1[0] = elcVr[threadIdx.x] + (D_parameter[4] / 2)*Qm*E[0];
		u1[1] = elcVy[threadIdx.x] + (D_parameter[4] / 2)*Qm*E[1];
		u1[2] = elcVz[threadIdx.x] + (D_parameter[4] / 2)*Qm*E[2];

		for (int m = 0; m<3; m++)
		{
			t[m] = (B[m] * Qm*D_parameter[4]) / 2;      // 鲍尔斯旋转，求t
			s[m] = (2 * t[m]) / (1 + t[m] * t[m]);        //求s
		}
		pp[0][0] = 1 - s[2] * t[2] - s[1] * t[1];      //3*3矩阵元素
		pp[0][1] = s[1] * t[0] + s[2];
		pp[0][2] = s[2] * t[0] - s[1];
		pp[1][0] = s[0] * t[1] - s[2];
		pp[1][1] = 1 - s[2] * t[2] - s[0] * t[0];
		pp[1][2] = s[0] + s[2] * t[1];
		pp[2][0] = s[0] * t[2] + s[1];
		pp[2][1] = s[1] * t[2] - s[0];
		pp[2][2] = 1 - s[1] * t[1] - s[0] * t[0];

		u2[0] = pp[0][0] * u1[0] + pp[0][1] * u1[1] + pp[0][2] * u1[2];
		u2[1] = pp[1][0] * u1[0] + pp[1][1] * u1[1] + pp[1][2] * u1[2];
		u2[2] = pp[2][0] * u1[0] + pp[2][1] * u1[1] + pp[2][2] * u1[2];

		for (int m = 0; m<3; m++)
			u3[m] = u2[m] + (D_parameter[4] / 2)*Qm*E[m];

		elcVr[threadIdx.x] = u3[0];
		elcVy[threadIdx.x] = u3[1];
		elcVz[threadIdx.x] = u3[2];

		float cit1a = 0;
		if ((elcPr[threadIdx.x]) == 0)
			cit1a = 0;
		else  cit1a = atan((elcVy[threadIdx.x] * D_parameter[4]) / (elcPr[threadIdx.x] + elcVr[threadIdx.x] * D_parameter[4]));
		float temp_x = elcPr[threadIdx.x] + elcVr[threadIdx.x] * D_parameter[4];
		if (temp_x >= 0)
			elcPr[threadIdx.x] = sqrt((elcPr[threadIdx.x] + elcVr[threadIdx.x] * D_parameter[4])*(elcPr[threadIdx.x] + elcVr[threadIdx.x] * D_parameter[4]) + elcVy[threadIdx.x] * D_parameter[4] * elcVy[threadIdx.x] * D_parameter[4]);
		else
		{
			elcPr[threadIdx.x] = -temp_x;
			elcVr[threadIdx.x] = -elcVr[threadIdx.x];
		}
		elcPz[threadIdx.x] = elcPz[threadIdx.x] + elcVz[threadIdx.x] * D_parameter[4];
		elcPy[threadIdx.x] = elcPy[threadIdx.x] + cit1a;
		if ((elcPr[threadIdx.x]) == 0)
		{
			elcVr[threadIdx.x] = elcVr[threadIdx.x];
			elcVy[threadIdx.x] = elcVy[threadIdx.x];
		}
		elcVr[threadIdx.x] = cos(cit1a)*elcVr[threadIdx.x] + sin(cit1a)*elcVy[threadIdx.x];
		elcVy[threadIdx.x] = -sin(cit1a)*elcVr[threadIdx.x] + cos(cit1a)*elcVy[threadIdx.x];
		p_pat_elc[tid].vr = elcVr[threadIdx.x];
		p_pat_elc[tid].vy = elcVy[threadIdx.x];
		p_pat_elc[tid].vz = elcVz[threadIdx.x];
		p_pat_elc[tid].pr = elcPr[threadIdx.x];
		p_pat_elc[tid].py = elcPy[threadIdx.x];
		p_pat_elc[tid].pz = elcPz[threadIdx.x];
		tid += gridDim.x*blockDim.x;
	}
}
__global__ void device_update_ion(float *device_static_magneticX, float *device_static_magneticZ, float *device_static_electricX, float *device_static_electricZ, Paticle *p_pat_elc, Pre_Paticle *d_pre_elc, Grid *device_G, float *S_number, int tail)
{
	//printf("%d",tail);
	float  Qm_ion = 7.33945e+5;
	float  Qm = -1.7588e+11;
	float mur = 4.0*pi*1.0e-7;
	float E[3] = { 0, 0, 0 }, B[3] = { 0, 0, 0 };
	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	int threadtid = threadIdx.x;
	//Paticle p_prepat_elc;
	__shared__ float elcVr[thread]; __shared__ float elcVy[thread]; __shared__ float elcVz[thread];
	__shared__ float elcPr[thread]; __shared__ float elcPy[thread]; __shared__ float elcPz[thread];
	__shared__ float PrelcPr[thread];  __shared__ float  PrelcPz[thread];
	while (tid<tail)
	{
		elcVr[threadIdx.x] = p_pat_elc[tid].vr;
		elcVy[threadIdx.x] = p_pat_elc[tid].vy;
		elcVz[threadIdx.x] = p_pat_elc[tid].vz;
		elcPr[threadIdx.x] = p_pat_elc[tid].pr;
		elcPy[threadIdx.x] = p_pat_elc[tid].py;
		elcPz[threadIdx.x] = p_pat_elc[tid].pz;

		d_pre_elc[tid].pr = elcPr[threadIdx.x];//
		d_pre_elc[tid].pz = elcPz[threadIdx.x];//
		int ii = (int)(elcPr[threadIdx.x] / D_parameter[0]);               //云方程
		int kk = (int)(elcPz[threadIdx.x] / D_parameter[1]);
		float wr = (elcPr[threadIdx.x] / D_parameter[0]) - ii;
		float wz = (elcPz[threadIdx.x] / D_parameter[1]) - kk;
		float s1 = (1 - wr)*(1 - wz);
		float s2 = (1 - wr)*wz;
		float s3 = wr*(1 - wz);
		float s4 = wr*wz;
		int grid_temp1 = ii*(nz + 1) + kk;//ii kk
		int grid_temp2 = ii*(nz + 1) + kk + 1;//ii kk+1
		int grid_temp3 = (ii + 1)*(nz + 1) + kk;//ii+1 kk
		int grid_temp4 = (ii + 1)*(nz + 1) + kk + 1;//ii+1 kk+1
		//printf("u:%d\t", tid);
		E[0] = ((*(device_G + grid_temp1)).ave_ex*s1 + (*(device_G + grid_temp3)).ave_ex*s3 + (*(device_G + grid_temp2)).ave_ex*s2 + (*(device_G + grid_temp4)).ave_ex*s4)
			+((*(device_static_electricX+grid_temp1))*s1+(*(device_static_electricX+grid_temp3))*s3+(*(device_static_electricX+grid_temp2))*s2+(*(device_static_electricX+grid_temp4))*s4);//+stac_ex[ii][kk];
		//printf("u:%d\tE0=%f\n", tid,E[0]);
		E[2] = ((*(device_G + grid_temp1)).ave_ez*s1 + (*(device_G + grid_temp3)).ave_ez*s3 + (*(device_G + grid_temp2)).ave_ez*s2 + (*(device_G + grid_temp4)).ave_ez*s4)
			+((*(device_static_electricZ+grid_temp1))*s1+(*(device_static_electricZ+grid_temp3))*s3+(*(device_static_electricZ+grid_temp2))*s2+(*(device_static_electricZ+grid_temp4))*s4);//+(*(device_static_electricZ+grid_temp1));
		//printf("u:%d\tE2=%f\n", tid,E[2]);
		B[1] = (((*(device_G + grid_temp1)).ave_hy*s1 + (*(device_G + grid_temp3)).ave_hy*s3 + (*(device_G + grid_temp2)).ave_hy*s2 + (*(device_G + grid_temp4)).ave_hy*s4)*mur);
		//printf("u:%d\tB1=%f\n", tid, B[1]);
		E[1] = ((*(device_G + grid_temp1)).ave_ey*s1 + (*(device_G + grid_temp3)).ave_ey*s3 + (*(device_G + grid_temp2)).ave_ey*s2 + (*(device_G + grid_temp4)).ave_ey*s4);
		//printf("u:%d\tE1=%f\n", tid,E[1]);
		B[0] = (((*(device_G + grid_temp1)).ave_hx*s1 + (*(device_G + grid_temp3)).ave_hx*s3 + (*(device_G + grid_temp2)).ave_hx*s2 + (*(device_G + grid_temp4)).ave_hx*s4)*mur)
		+(*(device_static_magneticX+grid_temp1))*s1+(*(device_static_magneticX+grid_temp3))*s3+((*(device_static_magneticX+grid_temp2))*s2+(*(device_static_magneticX+grid_temp4))*s4);//+(*(device_static_magneticX+grid_temp3));
		//B[0]=((G[ii][kk].ave_hx*s1+G[ii+1][kk].ave_hx*s3+G[ii][kk+1].ave_hx*s2+G[ii+1][kk+1].ave_hx*s4)*mur)
		//+(stac_Bx[ii][kk] * s1 + stac_Bx[ii + 1][kk] * s3 + stac_Bx[ii][kk + 1] * s2 + stac_Bx[ii + 1][kk + 1] * s4);
		//printf("u:%d\tB0=%f\n", tid, B[0]);
		B[2] = (((*(device_G + grid_temp1)).ave_hz*s1 + (*(device_G + grid_temp3)).ave_hz*s3 + (*(device_G + grid_temp2)).ave_hz*s2 + (*(device_G + grid_temp4)).ave_hz*s4)*mur)
		+((*(device_static_magneticZ+grid_temp1))*s1+(*(device_static_magneticZ+grid_temp3))*s3+(*(device_static_magneticZ+grid_temp2))*s2+(*(device_static_magneticZ+grid_temp4))*s4);//+(*(device_static_magneticZ+grid_temp3));
		//printf("B2=%f\n", B[2]);

		float u1[3] = { 0 }, u2[3] = { 0 }, u3[3] = { 0 };//u_n-1/2,u-,,u+,u_n+1/2    //    electron 分步求解
		float t[3] = { 0 }, s[3] = { 0 };
		float pp[3][3] = { 0 };

		u1[0] = elcVr[threadIdx.x] + (D_parameter[4] / 2)*Qm_ion*E[0];
		u1[1] = elcVy[threadIdx.x] + (D_parameter[4] / 2)*Qm_ion*E[1];
		u1[2] = elcVz[threadIdx.x] + (D_parameter[4] / 2)*Qm_ion*E[2];

		for (int m = 0; m<3; m++)
		{
			t[m] = (B[m] * Qm_ion*D_parameter[4]) / 2;      // 鲍尔斯旋转，求t
			s[m] = (2 * t[m]) / (1 + t[m] * t[m]);        //求s
		}
		pp[0][0] = 1 - s[2] * t[2] - s[1] * t[1];      //3*3矩阵元素
		pp[0][1] = s[1] * t[0] + s[2];
		pp[0][2] = s[2] * t[0] - s[1];
		pp[1][0] = s[0] * t[1] - s[2];
		pp[1][1] = 1 - s[2] * t[2] - s[0] * t[0];
		pp[1][2] = s[0] + s[2] * t[1];
		pp[2][0] = s[0] * t[2] + s[1];
		pp[2][1] = s[1] * t[2] - s[0];
		pp[2][2] = 1 - s[1] * t[1] - s[0] * t[0];

		u2[0] = pp[0][0] * u1[0] + pp[0][1] * u1[1] + pp[0][2] * u1[2];
		u2[1] = pp[1][0] * u1[0] + pp[1][1] * u1[1] + pp[1][2] * u1[2];
		u2[2] = pp[2][0] * u1[0] + pp[2][1] * u1[1] + pp[2][2] * u1[2];

		for (int m = 0; m<3; m++)
			u3[m] = u2[m] + (D_parameter[4] / 2)*Qm_ion*E[m];

		elcVr[threadIdx.x] = u3[0];
		elcVy[threadIdx.x] = u3[1];
		elcVz[threadIdx.x] = u3[2];

		float cit1a = 0;
		if ((elcPr[threadIdx.x]) == 0)
			cit1a = 0;
		else  cit1a = atan((elcVy[threadIdx.x] * D_parameter[4]) / (elcPr[threadIdx.x] + elcVr[threadIdx.x] * D_parameter[4]));
		float temp_x = elcPr[threadIdx.x] + elcVr[threadIdx.x] * D_parameter[4];
		if (temp_x >= 0)
			elcPr[threadIdx.x] = sqrt((elcPr[threadIdx.x] + elcVr[threadIdx.x] * D_parameter[4])*(elcPr[threadIdx.x] + elcVr[threadIdx.x] * D_parameter[4]) + elcVy[threadIdx.x] * D_parameter[4] * elcVy[threadIdx.x] * D_parameter[4]);
		else
		{
			elcPr[threadIdx.x] = -temp_x;
			elcVr[threadIdx.x] = -elcVr[threadIdx.x];
		}
		elcPz[threadIdx.x] = elcPz[threadIdx.x] + elcVz[threadIdx.x] * D_parameter[4];
		elcPy[threadIdx.x] = elcPy[threadIdx.x] + cit1a;
		if ((elcPr[threadIdx.x]) == 0)
		{
			elcVr[threadIdx.x] = elcVr[threadIdx.x];
			elcVy[threadIdx.x] = elcVy[threadIdx.x];
		}
		elcVr[threadIdx.x] = cos(cit1a)*elcVr[threadIdx.x] + sin(cit1a)*elcVy[threadIdx.x];
		elcVy[threadIdx.x] = -sin(cit1a)*elcVr[threadIdx.x] + cos(cit1a)*elcVy[threadIdx.x];
		p_pat_elc[tid].vr = elcVr[threadIdx.x];
		p_pat_elc[tid].vy = elcVy[threadIdx.x];
		p_pat_elc[tid].vz = elcVz[threadIdx.x];
		p_pat_elc[tid].pr = elcPr[threadIdx.x];
		p_pat_elc[tid].py = elcPy[threadIdx.x];
		p_pat_elc[tid].pz = elcPz[threadIdx.x];
		tid += gridDim.x*blockDim.x;
	}
}
__global__ void current_ion(Paticle *p_pat_elc, Pre_Paticle *d_pre_elc, Grid *device_G, float *S_number, int tail, int t)
{

	__shared__ float elcVr[thread]; __shared__ float elcVy[thread]; __shared__ float elcVz[thread];
	__shared__ float elcPr[thread]; __shared__ float elcPy[thread]; __shared__ float elcPz[thread];
	__shared__ float PrelcPr[thread]; __shared__ float  PrelcPz[thread];

	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	while (tid < tail)
	{

		elcVy[threadIdx.x] = p_pat_elc[tid].vy;

		elcPr[threadIdx.x] = p_pat_elc[tid].pr;

		elcPz[threadIdx.x] = p_pat_elc[tid].pz;

		PrelcPr[threadIdx.x] = d_pre_elc[tid].pr;//
		PrelcPz[threadIdx.x] = d_pre_elc[tid].pz;//
		if ((elcPz[threadIdx.x] >= D_parameter[3]) || (elcPz[threadIdx.x] < 0) || (elcPr[threadIdx.x] >= D_parameter[2]))
		{
			//int i = tid %30;
			//float afa = 0;//rds5*(pi/12);
			//float cita = 0;//rds6*pi/2;
			//float vv = sqrt(D_parameter[9] * 2 / D_parameter[11]);
			//p_pat_elc[tid].pz = D_parameter[8];
			//p_pat_elc[tid].pr = 0.005*S_number[i];
			//p_pat_elc[tid].py = 2 * pi*S_number[i];
			//p_pat_elc[tid].vr = vv*sin(afa)*cos(cita);
			//p_pat_elc[tid].vy = vv*sin(afa)*sin(cita);
			//p_pat_elc[tid].vz = 3e7;
			//p_pat_elc[tid].blei = p_pat_elc[tid].pr / D_parameter[0];
			//p_pat_elc[tid].blek = p_pat_elc[tid].pz / D_parameter[1];
		}
		else
		{
			int i = (int)(elcPr[threadIdx.x] / D_parameter[0]);
			int j = (int)(elcPz[threadIdx.x] / D_parameter[1]);
			int ii = (int)(PrelcPr[threadIdx.x] / D_parameter[0]);               //云方程
			int kk = (int)(PrelcPz[threadIdx.x] / D_parameter[1]);
			float wrr = PrelcPr[threadIdx.x] / D_parameter[0] - ii;
			float wzz = PrelcPz[threadIdx.x] / D_parameter[1] - kk;
			float newwrr = elcPr[threadIdx.x] / D_parameter[0] - i;
			float newwzz = elcPz[threadIdx.x] / D_parameter[1] - j;
			float V = abs(pi*((ii + 1)*D_parameter[0] * (ii + 1)*D_parameter[0] - ii*D_parameter[0] * ii*D_parameter[0])*D_parameter[1]);
			float V1 = abs(pi*((i + 1)*D_parameter[0] * (i + 1)*D_parameter[0] - i*D_parameter[0] * i*D_parameter[0])*D_parameter[1]);


			int grid_1 = i*(nz + 1) + j;
			int grid_2 = i*(nz + 1) + j + 1;
			int grid_3 = (i + 1)*(nz + 1) + j;
			int grid_4 = (i + 1)*(nz + 1) + j + 1;
			int grid_temp1 = ii*(nz + 1) + kk;//ii kk
			int grid_temp2 = ii*(nz + 1) + kk + 1;//ii kk+1
			int grid_temp3 = (ii + 1)*(nz + 1) + kk;//ii+1 kk
			int grid_temp4 = (ii + 1)*(nz + 1) + kk + 1;//ii+1 kk+1
			/*(*(device_G + grid_1)).Q += D_parameter[6] * (1 - wrr)*(1 - wzz);
			(*(device_G + grid_2)).Q += D_parameter[6] * (1 - wrr)*wzz;
			(*(device_G + grid_3)).Q += D_parameter[6] * wrr*(1 - wzz);
			(*(device_G + grid_4)).Q += D_parameter[6] * wrr*wzz;*/
			//float area = pi*((ii*D_parameter[0] + D_parameter[0])*(ii*D_parameter[0] + D_parameter[0]) - ii*D_parameter[0] * ii*D_parameter[0])*D_parameter[1];
			/*(*(device_G + grid_1)).den = (*(device_G + grid_1)).Q / area;
			(*(device_G + grid_2)).den = (*(device_G + grid_2)).Q / area;
			(*(device_G + grid_3)).den = (*(device_G + grid_3)).Q / area;
			(*(device_G + grid_4)).den = (*(device_G + grid_4)).Q / area;*/

			//粒子穿越网格，四种情况

			float xp = d_min_1(d_min_1(ii*D_parameter[0], i*D_parameter[0]) + D_parameter[0], d_max_1(d_max_1(ii*D_parameter[0], i*D_parameter[0]), (elcPr[threadIdx.x] + PrelcPr[threadIdx.x]) / 2));
			float zp = d_min_1(d_min_1(kk*D_parameter[1], j*D_parameter[1]) + D_parameter[1], d_max_1(d_max_1(kk*D_parameter[1], j*D_parameter[1]), (elcPz[threadIdx.x] + PrelcPz[threadIdx.x]) / 2));
			float fr1 = D_parameter[6] * (xp - PrelcPr[threadIdx.x]) / D_parameter[4];
			float fz1 = D_parameter[6] * (zp - PrelcPz[threadIdx.x]) / D_parameter[4];

			float fr2 = D_parameter[6] * (elcPr[threadIdx.x] - xp) / D_parameter[4];
			float fz2 = D_parameter[6] * (elcPz[threadIdx.x] - zp) / D_parameter[4];

			float wr1 = (xp + PrelcPr[threadIdx.x]) / 2 / D_parameter[1] - ii;
			float wz1 = (zp + PrelcPz[threadIdx.x]) / 2 / D_parameter[1] - kk;
			float wr2 = (xp + elcPr[threadIdx.x]) / 2 / D_parameter[0] - i;
			float wz2 = (zp + elcPz[threadIdx.x]) / 2 / D_parameter[1] - j;
			//printf("%d,%e,%e\n",tail, wr2, wz2);
			/*float da = (*(device_G + grid_1)).den*p_pat_elc[tid].vz;
			float da1 = (*(device_G + grid_3)).den*p_pat_elc[tid].vz;*/
			///////////////////////////////////////////////////电子电流密度	1
			//int logic = ii*i / ((ii - 0.001)*(i - 0.001));//逻辑体积 分支优化
			//float logicV = logic*V1 + (1 - logic)*V;
			if (ii == 0 || i == 0)
			{
				atomicAdd(&((*(device_G + grid_temp1)).jr_ion), (fr1*(1 - wz1) / V));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_temp2)).jr_ion), (fr1*(wz1) / V));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_temp1)).jz_ion), (fz1*(1 - wr1) / V));//注意是否除以i
				__syncthreads();
				atomicAdd(&((*(device_G + grid_temp3)).jz_ion), (fz1*(wr1) / V));
				__syncthreads();

				atomicAdd(&((*(device_G + grid_1)).jr_ion), (fr2*(1 - wz2) / V));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_2)).jr_ion), ((fr2*wz2) / V));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_1)).jz_ion), (fz2*(1 - wr2) / V));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_3)).jz_ion), (fz2*wr2 / V));
				__syncthreads();

			}
			else{
				atomicAdd(&((*(device_G + grid_temp1)).jr_ion), (fr1*(1 - wz1) / V1));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_temp2)).jr_ion), (fr1*(wz1) / V1));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_temp1)).jz_ion), (fz1*(1 - wr1) / V1));//注意是否除以i
				__syncthreads();
				atomicAdd(&((*(device_G + grid_temp3)).jz_ion), (fz1*(wr1) / V1));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_1)).jr_ion), (fr2*(1 - wz2) / V1));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_2)).jr_ion), ((fr2*wz2) / V1));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_1)).jz_ion), (fz2*(1 - wr2) / V1));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_3)).jz_ion), (fz2*wr2 / V1));
				__syncthreads();

			}

			float qedens = D_parameter[6] / V;   //电流密度Jy方向
			float Jc = qedens*elcVy[threadIdx.x];
			float mid1, mid2;
			mid1 = xp / D_parameter[0] - i;
			mid2 = zp / D_parameter[1] - j;
			float A1 = (1 - wrr)*(1 - wzz);
			float A3 = wrr*(1 - wzz);
			float A2 = (1 - wrr)*wzz;
			float A4 = wrr*wzz;

			float M1 = (1 - mid1)*(1 - mid2);
			float M3 = mid1*(1 - mid2);
			float M2 = (1 - mid1)*mid2;
			float M4 = mid1*mid2;

			float B1 = (1 - newwrr)*(1 - newwzz);
			float B3 = newwrr*(1 - newwzz);
			float B2 = (1 - newwrr)*newwzz;
			float B4 = newwrr*newwzz;

			atomicAdd(&(*(device_G + grid_temp1)).jy_ion, (Jc*(A1 + M1)) / 4 - 0.0);
			__syncthreads();
			atomicAdd(&(*(device_G + grid_temp2)).jy_ion, (Jc*(A2 + M2)) / 4 - 0.0);
			__syncthreads();
			atomicAdd(&(*(device_G + grid_temp3)).jy_ion, (Jc*(A3 + M3)) / 4 - 0.0);
			__syncthreads();
			atomicAdd(&(*(device_G + grid_temp4)).jy_ion, (Jc*(A4 + M4)) / 4 - 0.0);
			__syncthreads();

			atomicAdd(&(*(device_G + grid_1)).jy_ion, (Jc*(B1 + M1)) / 4 - 0.0);
			__syncthreads();
			atomicAdd(&(*(device_G + grid_2)).jy_ion, (Jc*(B2 + M2)) / 4 - 0.0);
			__syncthreads();
			atomicAdd(&(*(device_G + grid_3)).jy_ion, (Jc*(B3 + M3)) / 4 - 0.0);
			__syncthreads();
			atomicAdd(&(*(device_G + grid_4)).jy_ion, (Jc*(B4 + M4)) / 4 - 0.0);
			__syncthreads();
		}
		tid += gridDim.x*blockDim.x;
	}
}
__global__ void current(Paticle *p_pat_elc, Pre_Paticle *d_pre_elc, Grid *device_G, float *S_number, int tail,int t)
{

	__shared__ float elcVr[thread]; __shared__ float elcVy[thread]; __shared__ float elcVz[thread];
	__shared__ float elcPr[thread]; __shared__ float elcPy[thread]; __shared__ float elcPz[thread];
	__shared__ float PrelcPr[thread];__shared__ float  PrelcPz[thread];

	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	while (tid < tail)
	{
		
		elcVy[threadIdx.x] = p_pat_elc[tid].vy;
		
		elcPr[threadIdx.x] = p_pat_elc[tid].pr;

		elcPz[threadIdx.x] = p_pat_elc[tid].pz;

		PrelcPr[threadIdx.x] = d_pre_elc[tid].pr;//
		PrelcPz[threadIdx.x] = d_pre_elc[tid].pz;//
		if ((elcPz[threadIdx.x] >= D_parameter[3]) || (elcPz[threadIdx.x] < 0) || (elcPr[threadIdx.x] >= D_parameter[2]))
		{
			//int i = tid %30;
			//float afa = 0;//rds5*(pi/12);
			//float cita = 0;//rds6*pi/2;
			//float vv = sqrt(D_parameter[9] * 2 / D_parameter[11]);
			//p_pat_elc[tid].pz = D_parameter[8];
			//p_pat_elc[tid].pr = 0.005*S_number[i];
			//p_pat_elc[tid].py = 2 * pi*S_number[i];
			//p_pat_elc[tid].vr = vv*sin(afa)*cos(cita);
			//p_pat_elc[tid].vy = vv*sin(afa)*sin(cita);
			//p_pat_elc[tid].vz = 3e7;
			//p_pat_elc[tid].blei = p_pat_elc[tid].pr / D_parameter[0];
			//p_pat_elc[tid].blek = p_pat_elc[tid].pz / D_parameter[1];
		}
		else
		{
			int i = (int)(elcPr[threadIdx.x] / D_parameter[0]);
			int j = (int)(elcPz[threadIdx.x] / D_parameter[1]);
			int ii = (int)(PrelcPr[threadIdx.x] / D_parameter[0]);               //云方程
			int kk = (int)(PrelcPz[threadIdx.x] / D_parameter[1]);
			float wrr = PrelcPr[threadIdx.x] / D_parameter[0] - ii;
			float wzz = PrelcPz[threadIdx.x] / D_parameter[1] - kk;
			float newwrr = elcPr[threadIdx.x] / D_parameter[0] - i;
			float newwzz = elcPz[threadIdx.x] / D_parameter[1] - j;
			float V = abs(pi*((ii + 1)*D_parameter[0] * (ii + 1)*D_parameter[0] - ii*D_parameter[0] * ii*D_parameter[0])*D_parameter[1]);
			float V1 = abs(pi*((i + 1)*D_parameter[0] * (i + 1)*D_parameter[0] - i*D_parameter[0] * i*D_parameter[0])*D_parameter[1]);
			
			
			int grid_1 = i*(nz + 1) + j;
			int grid_2 = i*(nz + 1) + j + 1;
			int grid_3 = (i + 1)*(nz + 1) + j;
			int grid_4 = (i + 1)*(nz + 1) + j + 1;
			int grid_temp1 = ii*(nz + 1) + kk;//ii kk
			int grid_temp2 = ii*(nz + 1) + kk + 1;//ii kk+1
			int grid_temp3 = (ii + 1)*(nz + 1) + kk;//ii+1 kk
			int grid_temp4 = (ii + 1)*(nz + 1) + kk + 1;//ii+1 kk+1
			/*(*(device_G + grid_1)).Q += D_parameter[5] * (1 - wrr)*(1 - wzz);
			(*(device_G + grid_2)).Q += D_parameter[5] * (1 - wrr)*wzz;
			(*(device_G + grid_3)).Q += D_parameter[5] * wrr*(1 - wzz);
			(*(device_G + grid_4)).Q += D_parameter[5] * wrr*wzz;*/
			//float area = pi*((ii*D_parameter[0] + D_parameter[0])*(ii*D_parameter[0] + D_parameter[0]) - ii*D_parameter[0] * ii*D_parameter[0])*D_parameter[1];
			/*(*(device_G + grid_1)).den = (*(device_G + grid_1)).Q / area;
			(*(device_G + grid_2)).den = (*(device_G + grid_2)).Q / area;
			(*(device_G + grid_3)).den = (*(device_G + grid_3)).Q / area;
			(*(device_G + grid_4)).den = (*(device_G + grid_4)).Q / area;*/

			//粒子穿越网格，四种情况

			float xp = d_min_1(d_min_1(ii*D_parameter[0], i*D_parameter[0]) + D_parameter[0], d_max_1(d_max_1(ii*D_parameter[0], i*D_parameter[0]), (elcPr[threadIdx.x] + PrelcPr[threadIdx.x]) / 2));
			float zp = d_min_1(d_min_1(kk*D_parameter[1], j*D_parameter[1]) + D_parameter[1], d_max_1(d_max_1(kk*D_parameter[1], j*D_parameter[1]), (elcPz[threadIdx.x] + PrelcPz[threadIdx.x]) / 2));
			float fr1 = D_parameter[5] * (xp - PrelcPr[threadIdx.x]) / D_parameter[4];
			float fz1 = D_parameter[5] * (zp - PrelcPz[threadIdx.x]) / D_parameter[4];
			
			float fr2 = D_parameter[5] * (elcPr[threadIdx.x] - xp) / D_parameter[4];
			float fz2 = D_parameter[5] * (elcPz[threadIdx.x] - zp) / D_parameter[4];
			
			float wr1 = (xp + PrelcPr[threadIdx.x]) / 2 / D_parameter[1] - ii;
			float wz1 = (zp + PrelcPz[threadIdx.x]) / 2 / D_parameter[1] - kk;
			float wr2 = (xp + elcPr[threadIdx.x]) / 2 / D_parameter[0] - i;
			float wz2 = (zp + elcPz[threadIdx.x]) / 2 / D_parameter[1] - j;
			//printf("%d,%e,%e\n",tail, wr2, wz2);
			/*float da = (*(device_G + grid_1)).den*p_pat_elc[tid].vz;
			float da1 = (*(device_G + grid_3)).den*p_pat_elc[tid].vz;*/
			///////////////////////////////////////////////////电子电流密度	1
			//int logic = ii*i / ((ii - 0.001)*(i - 0.001));//逻辑体积 分支优化
			//float logicV = logic*V1 + (1 - logic)*V;
			if (ii == 0||i==0 )
			{
				atomicAdd(&((*(device_G + grid_temp1)).jr), (fr1*(1 - wz1) / V));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_temp2)).jr), (fr1*(wz1) / V));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_temp1)).jz), (fz1*(1 - wr1) / V));//注意是否除以i
				__syncthreads();
				atomicAdd(&((*(device_G + grid_temp3)).jz), (fz1*(wr1) / V));
				__syncthreads();
				
				atomicAdd(&((*(device_G + grid_1)).jr), (fr2*(1 - wz2) / V));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_2)).jr), ((fr2*wz2) / V));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_1)).jz), (fz2*(1 - wr2) / V));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_3)).jz), (fz2*wr2 / V));
				__syncthreads();

			}
			else{
				atomicAdd(&((*(device_G + grid_temp1)).jr), (fr1*(1 - wz1) / V1));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_temp2)).jr), (fr1*(wz1) / V1));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_temp1)).jz), (fz1*(1 - wr1) / V1));//注意是否除以i
				__syncthreads();
				atomicAdd(&((*(device_G + grid_temp3)).jz), (fz1*(wr1) / V1));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_1)).jr), (fr2*(1 - wz2) / V1));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_2)).jr), ((fr2*wz2) / V1));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_1)).jz), (fz2*(1 - wr2) / V1));
				__syncthreads();
				atomicAdd(&((*(device_G + grid_3)).jz), (fz2*wr2 / V1));
				__syncthreads();
				
			}
	
			float qedens = D_parameter[5] / V;   //电流密度Jy方向
			float Jc = qedens*elcVy[threadIdx.x];
			float mid1, mid2;
			mid1 = xp / D_parameter[0] - i;
			mid2 = zp / D_parameter[1] - j;
			float A1 = (1 - wrr)*(1 - wzz);
			float A3 = wrr*(1 - wzz);
			float A2 = (1 - wrr)*wzz;
			float A4 = wrr*wzz;

			float M1 = (1 - mid1)*(1 - mid2);
			float M3 = mid1*(1 - mid2);
			float M2 = (1 - mid1)*mid2;
			float M4 = mid1*mid2;

			float B1 = (1 - newwrr)*(1 - newwzz);
			float B3 = newwrr*(1 - newwzz);
			float B2 = (1 - newwrr)*newwzz;
			float B4 = newwrr*newwzz;
			
			atomicAdd(&(*(device_G + grid_temp1)).jy, (Jc*(A1 + M1)) / 4-0.0);
			__syncthreads();
			atomicAdd(&(*(device_G + grid_temp2)).jy, (Jc*(A2 + M2)) / 4-0.0);
			__syncthreads();
			atomicAdd(&(*(device_G + grid_temp3)).jy, (Jc*(A3 + M3)) / 4-0.0);
			__syncthreads();
			atomicAdd(&(*(device_G + grid_temp4)).jy, (Jc*(A4 + M4)) / 4-0.0);
			__syncthreads();
		
			atomicAdd(&(*(device_G + grid_1)).jy, (Jc*(B1 + M1)) / 4-0.0);
			__syncthreads();
			atomicAdd(&(*(device_G + grid_2)).jy, (Jc*(B2 + M2)) / 4-0.0);
			__syncthreads();
			atomicAdd(&(*(device_G + grid_3)).jy, (Jc*(B3 + M3)) / 4-0.0);
			__syncthreads();
			atomicAdd(&(*(device_G + grid_4)).jy, (Jc*(B4 + M4)) / 4-0.0);
			__syncthreads();
		}
		tid += gridDim.x*blockDim.x;
	}
}

__global__ void device_define_G(int n, int m, Grid *device_gn, Grid *device_g)
{
	
		int tid=blockIdx.x*blockDim.x+threadIdx.x;
		while(tid<n*m)
		{
			/*int i=tid/m;
			int j=tid%m;
			int tid_temp=i*(m+1)+j;*/
			(device_gn+tid)->ey=(device_g+tid)->ey;
			(device_gn+tid)->hx=(device_g+tid)->hx;
			(device_gn+tid)->hz=(device_g+tid)->hz;
			(device_gn+tid)->ex=(device_g+tid)->ex;
			(device_gn+tid)->hy=(device_g+tid)->hy;
			(device_gn+tid)->ez=(device_g+tid)->ez;
			(device_g+tid)->ne[2]=0;

			(device_g+tid)->jr=0.0;
			(device_g+tid)->jz=0.0;
			(device_g+tid)->jy=0.0;
			(device_g+tid)->Q=0.0;
			(device_g+tid)->den=0.0;

			(device_g+tid)->jr_ion=0.0;
			(device_g+tid)->jz_ion=0.0;
			(device_g+tid)->jy_ion=0.0;
			(device_g+tid)->Q_ion=0.0;
			(device_g+tid)->den_ion=0.0;
			tid+=gridDim.x*blockDim.x;
		}
		
		/* printf("success,device_define_G");*/
}

__device__ void L_InitialPML(int tid)
{
	
	    d_ex1[tid] = 0;
	    d_iex1[tid] = 0;
		d_ey1[tid] = 0;
		d_iey1[tid] = 0;
		d_ez1[tid] = 0;
		d_iez1[tid] = 0;
		d_hx1[tid] = 0;
		d_ihx1[tid] = 0; 
		d_hy1[tid]=0; 
		d_ihy1[tid]=0; 
		d_hz1[tid]=0; 
		d_ihz1[tid]=0;
		
}

__global__ void kernel_L_InitialPML(int nxx, int nzz)
{
	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	while (tid<nxx*nzz)
	{
		L_InitialPML(tid);
		tid += gridDim.x*blockDim.x;
	}
	
}
__global__ void cacuchang_hx(Grid *device_G,Grid *device_Gn,float *d_sigmaz1,float *d_sigmaz,float dt,float dr,float dz,int nxx,int nzz)              
{
	
	float epsl=(8.854e-12);
	float mur=4.0*pi*1.0e-7;
	//float t0=120*dt,T=100*dt;
	float /*ca=0.0,*//*cb=0.0,*/ca1=0.0,cb1=0.0/*,ca2=0,cb2=0*/;
	float /*da=0.0,db=0.0,*/da1=0.0,db1=0.0,da2=0,db2=0;
	int tid=blockDim.x*blockIdx.x+threadIdx.x;
	/////********************   TE01    ******************/////
	
	while (tid<nxx*nzz)
	{
		
				int i=tid/nzz;
				int k=tid%nzz;
				if (i != nx&&k != nz)
				{
					da2 = 1;//(1+dt*sigmaz1[k]/epsl/2);
					db2 = 1;//(1-dt*sigmaz1[k]/epsl/2);
					da1 = (2 * epsl - dt*d_sigmaz[k]) / (2 * epsl + dt*d_sigmaz[k]);
					db1 = (2 * epsl) / (2 * epsl + d_sigmaz[k] * dt);
					d_hx1[tid] = da1*d_hx1[tid] + db1*dt*((*(device_Gn + tid+1)).ey - (*(device_Gn + tid)).ey) / dz;
					(*(device_G + tid)).hx = (*(device_Gn + tid)).hx + (d_hx1[tid] - d_ihx1[tid]) / mur;
					d_ihx1[tid] = d_hx1[tid];
				}
			tid+=gridDim.x*blockDim.x;
	}
	
}
__global__ void cacuchang_hy(Grid *device_G, Grid *device_Gn, float *d_sigmaz1, float *d_sigmaz, float dt, float dr, float dz, int nxx, int nzz)
{

	float epsl = (8.854e-12);
	float mur = 4.0*pi*1.0e-7;
	//float t0=120*dt,T=100*dt;
	float /*ca=0.0,*//*cb=0.0,*/ca1 = 0.0, cb1 = 0.0/*,ca2=0,cb2=0*/;
	float /*da=0.0,db=0.0,*/da1 = 0.0, db1 = 0.0, da2 = 0, db2 = 0;
	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	/////********************   TE01    ******************/////

	while (tid<nxx*nzz)
	{

		int i = tid / nzz;
		int k = tid%nzz;
		if (i != nx&&k != nz)
		{
			da2 = 1;//(1+dt*sigmaz1[k]/epsl/2);
			db2 = 1;//(1-dt*sigmaz1[k]/epsl/2);
			ca1 = (2 * epsl - dt*d_sigmaz[k]) / (2 * epsl + dt*d_sigmaz[k]);
			cb1 = 2 * epsl / (2 * epsl + dt*d_sigmaz[k]);
			d_hy1[tid] = d_hy1[tid] + dt*(((*(device_Gn + tid+nzz)).ez - (*(device_Gn + tid)).ez) / dr - ((*(device_Gn + tid+1)).ex - (*(device_Gn + tid)).ex) / dz);
			(*(device_G + tid)).hy = ca1*(*(device_Gn + tid)).hy + cb1*(d_hy1[tid] - d_ihy1[tid]) / mur;
			d_ihy1[tid] = d_hy1[tid];
		}
		tid += gridDim.x*blockDim.x;
	}

}
__global__ void cacuchang_hz(Grid *device_G, Grid *device_Gn, float *d_sigmaz1, float *d_sigmaz, float dt, float dr, float dz, int nxx, int nzz)
{

	float epsl = (8.854e-12);
	float mur = 4.0*pi*1.0e-7;
	//float t0=120*dt,T=100*dt;
	float /*ca=0.0,*//*cb=0.0,*/ca1 = 0.0, cb1 = 0.0/*,ca2=0,cb2=0*/;
	float /*da=0.0,db=0.0,*/da1 = 0.0, db1 = 0.0, da2 = 0, db2 = 0;
	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	/////********************   TE01    ******************/////
	while (tid<nxx*nzz)
	{

		int i = tid / nzz;
		int k = tid%nzz;
		if (i != nx&&k != nz)
		{
			da2 = 1;//(1+dt*sigmaz1[k]/epsl/2);
			db2 = 1;//(1-dt*sigmaz1[k]/epsl/2);
			d_hz1[tid] = d_hz1[tid] - dt*((*(device_Gn + tid+nzz)).ey - (*(device_Gn + tid)).ey) / dr - dt*((*(device_Gn + tid+nzz)).ey + (*(device_Gn + tid)).ey) / (2 * (i + 0.5)*dr);
			(*(device_G + tid)).hz = (*(device_Gn + tid)).hz + (da2*d_hz1[tid] - db2*d_ihz1[tid]) / mur;
			d_ihz1[tid] = d_hz1[tid];
		}
		tid += gridDim.x*blockDim.x;
	}

}

__global__ void cacuchang_ex(Grid *device_G,Grid *device_Gn,float *d_sigmaz1,float *d_sigmaz,float dt,float dr,float dz,int nxx,int nzz)   
{
	
	float epsl=(8.854e-12);
	const float cgm=0.0;
	float mur=4.0*pi*1.0e-7;
	float A,B,C,D,v,a6=0.0;
	//float bate=0;
	float omega=0,freq=0;
	float t0=120*dt,T=100*dt;
	float ca=0.0,cb=0.0,ca1=0.0,cb1=0.0,ca2=0,cb2=0;
	float da=0.0,db=0.0,da1=0.0,db1=0.0,da2=0,db2=0;

    v=1/sqrt(mur*epsl);
	A=(2*epsl-cgm*dt)/(2*epsl-cgm*dt);
	B=(2*dt)/(2*epsl+cgm*dt);
	C=1;
	D=dt/mur;	
	freq=1e8;
	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	while (tid < nxx*nzz)
	{
		int i = tid / nzz;
		int k = tid%nzz;
		if (i != nx && k != nz)
		{
			if (k == 0)
			{
				float v = 1 / sqrt(mur*epsl);
				float coff = (v*dt - dz) / (v*dt + dz);		
				device_G[tid].ex = device_Gn[tid + 1].ex + coff*(device_Gn[tid + 1].ex - device_Gn[tid].ex);
			}
			else
			{
				ca = (2 * epsl - dt*d_sigmaz1[k]) / (2 * epsl + dt*d_sigmaz1[k]);
				cb = 2 * dt / (2 * epsl + dt*d_sigmaz1[k]);
				d_ex1[tid] = ca*d_ex1[tid] - cb*(device_G[tid].hy - device_G[tid - 1].hy) / dz;
				device_G[tid].ex = device_Gn[tid].ex + (d_ex1[tid] - d_iex1[tid]) - B*(device_G[tid].jr + device_G[tid].jr_ion);
				d_iex1[tid] = d_ex1[tid];
			}
		}
		tid += gridDim.x*blockDim.x;
	}
		
	
}

__global__ void cacuchang_ey(Grid *device_G,Grid *device_Gn,float *d_sigmaz1,float *d_sigmaz,float dt,float dr,float dz,int nxx,int nzz)   
{
	
	float epsl=(8.854e-12);
	const float cgm=0.0;
	float mur=4.0*pi*1.0e-7;
	float A,B,C,D,v,a6=0.0;
	float omega=0,freq=0;
	float t0=120*dt,T=100*dt;
	float ca=0.0,cb=0.0,ca1=0.0,cb1=0.0,ca2=0,cb2=0;
	float da=0.0,db=0.0,da1=0.0,db1=0.0,da2=0,db2=0;

    v=1/sqrt(mur*epsl);
	A=(2*epsl-cgm*dt)/(2*epsl-cgm*dt);
	B=(2*dt)/(2*epsl+cgm*dt);
	C=1;
	D=dt/mur;	
	freq=1e8;
	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	while(tid<nxx*nzz)
	{
		int i=tid/nz;
		int k=tid%nz;
		if (i != nx&&k != nz)
		{
			int tid1 = i*nzz + k;
			int tid2 = i*nzz + k + 1;
			int tid3 = (i + 1)*nzz + k;
			int tid4 = (i + 1)*nzz + k + 1;
			if (i == 0)
			{
				device_G[tid].ey = 0;//te01
			}
			else  if (k == 0)
			{
				float v = 1 / sqrt(mur*epsl);
				float coff = (v*dt - dz) / (v*dt + dz);
				device_G[tid].ey = device_Gn[tid + 1].ey + coff*(device_Gn[tid + 1].ey - device_Gn[tid].ey);//TE01

			}
			else
			{
				da = (2 * epsl - dt*d_sigmaz1[k]) / (2 * epsl + dt*d_sigmaz1[k]); //TE01
				db = 2 / (2 * epsl + dt*d_sigmaz1[k]);
				d_ey1[tid] = d_ey1[tid] + dt*((device_G[tid].hx - device_G[tid - 1].hx) / dz - (device_G[tid].hz - device_G[tid-nzz].hz) / dr);//此处严重错误
				device_G[tid].ey = da*device_Gn[tid].ey + db*(d_ey1[tid] - d_iey1[tid]) - B*(device_G[tid].jy + device_G[tid].jy_ion) * 10;
				d_iey1[tid] = d_ey1[tid];
			}
		}
			 tid+=gridDim.x*blockDim.x;
	}
	
}
__global__ void cacuchang_ez(Grid *device_G,Grid *device_Gn,float *d_sigmaz1,float *d_sigmaz,float dt,float dr,float dz,int nxx,int nzz)
{
	
	float epsl=(8.854e-12);
	const float cgm=0.0;
	float mur=4.0*pi*1.0e-7;
	float A,B,C,D,v,a6=0.0;
	//float bate=0;
	float omega=0,freq=0;
	float t0=120*dt,T=100*dt;
	float ca=0.0,cb=0.0,ca1=0.0,cb1=0.0,ca2=0,cb2=0;
	float da=0.0,db=0.0,da1=0.0,db1=0.0,da2=0,db2=0;

    v=1/sqrt(mur*epsl);
	A=(2*epsl-cgm*dt)/(2*epsl-cgm*dt);
	B=(2*dt)/(2*epsl+cgm*dt);
	//printf("%e\n",B);
	C=1;
	D=dt/mur;	
	freq=1e8;
	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	while(tid<nxx*nzz)
	{
		int i=tid/nzz;
		int k=tid%nzz;
		if (i != nx&&k != nz)
		{
			if (i == 0)
			{
				ca2 = 1;//(2*epsl+dt*sigmaz[k])/(2*epsl);
				cb2 = 1;//(2*epsl-dt*sigmaz[k])/(2*epsl);
				d_ez1[tid] = d_ez1[tid] + dt * 4 * (*(device_G + tid)).hy / dr / epsl;
				(*(device_G + tid)).ez = (*(device_Gn + tid)).ez + (ca2*d_ez1[tid] - cb2*d_iez1[tid]) - B*((*(device_G + tid)).jz + (*(device_G + tid)).jz_ion);
				d_iez1[tid] = d_ez1[tid];
			}
			else
			{
				ca2 = 1;//(2*epsl+dt*sigmaz[k])/(2*epsl);//TM01
				cb2 = 1;//(2*epsl-dt*sigmaz[k])/(2*epsl);
				d_ez1[tid] = d_ez1[tid] + dt*((1 / (2 * i*dr) + 1 / dr)*(*(device_G + tid)).hy + ((1 / (2 * i*dr)) - 1 / dr)*(*(device_G + tid - nzz)).hy);
				(*(device_G + tid)).ez = (*(device_Gn + tid)).ez + (ca2*d_ez1[tid] - cb2*d_iez1[tid]) / epsl - B*((*(device_G + tid)).jz + (*(device_G + tid)).jz_ion);
				d_iez1[tid] = d_ez1[tid];
			}
		}
	 tid+=gridDim.x*blockDim.x;
	}
}
