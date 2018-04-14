#ifndef LTY_PARALLE_COMMON_CUH
#define LTY_PARALLE_COMMON_CUH
#include"grid.cuh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda.h"
#include "particle.cuh"
#include<curand_kernel.h>
#define block 18
#define thread 128
#define nx 50
#define nz 100
#define pi 3.1415926
//	还没拆分
 void curandnumber(int n);//随机数生成函数,n是随机数的个数
__device__ float d_min_1(float x,float y);
__device__ float d_max_1(float x,float y);
__global__ void kernel_L_InitialPML(int nxx,int nzz);
__device__ void L_InitialPML(int tid);
__global__ void device_update_last(float *device_static_magneticX, float *device_static_magneticZ, float *device_static_electricX, float *device_static_electricZ, Paticle *p_pat_elc, Pre_Paticle *d_pre_elc, Grid *device_G, float *S_number, int tail);
__global__ void device_update_ion(float *device_static_magneticX, float *device_static_magneticZ, float *device_static_electricX, float *device_static_electricZ, Paticle *p_pat_elc, Pre_Paticle *d_pre_elc, Grid *device_G, float *S_number, int tail);
__global__ void initial_always(Paticle *p_pat_elc,int n,float *S_number,int tmp);//n是step*30
__global__ void device_define_G(int n, int m, Grid *Gn, Grid *G);
__global__ void device_initialchang(Grid* device_g,Grid* device_gn);
__global__ void device_ave_field(int x,int y,Grid* device_g,Grid* device_gn);
__global__ void current(Paticle *p_pat_elc, Pre_Paticle *d_pre_elc, Grid *device_G, float *S_number, int tail,int t);
__global__ void current_ion(Paticle *p_pat_elc, Pre_Paticle *d_pre_elc, Grid *device_G, float *S_number, int tail, int t);
__global__ void cacuchang_hx(Grid *device_G,Grid *device_Gn,float *d_sigmaz1,float *d_sigmaz,float dt,float dr,float dz,int nxx,int nzz);
__global__ void cacuchang_hy(Grid *device_G,Grid *device_Gn,float *d_sigmaz1,float *d_sigmaz,float dt,float dr,float dz,int nxx,int nzz);
__global__ void cacuchang_hz(Grid *device_G,Grid *device_Gn,float *d_sigmaz1,float *d_sigmaz,float dt,float dr,float dz,int nxx,int nzz);
__global__ void cacuchang_ex(Grid *device_G,Grid *device_Gn,float *d_sigmaz1,float *d_sigmaz,float dt,float dr,float dz,int nxx,int nzz);
__global__ void cacuchang_ey(Grid *device_G,Grid *device_Gn,float *d_sigmaz1,float *d_sigmaz,float dt,float dr,float dz,int nxx,int nzz);
__global__ void cacuchang_ez(Grid *device_G,Grid *device_Gn,float *d_sigmaz1,float *d_sigmaz,float dt,float dr,float dz,int nxx,int nzz);
#endif

