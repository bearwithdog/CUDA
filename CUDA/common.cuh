#pragma once
#ifndef _COMMON_CUH
#define _COMMON_CUH
#include "parameter.cuh"
//#include "grid.h"


const float pi=3.1415926;        

const float c=2.99792458e8;      //光速   

const float q=-1.60e-19;         // 电子带电荷量 

const float boci=1.6022e-19;     //波尔茨曼常数，转换后的1ev=boci

const float Me=9.1e-31;          //电子质量

const float Mi=2.18e-25;         //离子质量

const float  Qm=-1.7588e+11;     //电子荷质比

const float  Qm_ion=7.33945e+5;     //离子荷质比

const int filename_size=1024;

float epsl=(8.854e-12);

const float sbepsl=8.854e-8;

const float cgm=0.0;

const int pml=10;

float mur=4.0*pi*1.0e-7;

const float cof=1.0/100;

int nz;
int nx;
int nzz;
int nxx;

float dz;  //轴向空间步长
float dr; //径向空间步长

//float dt=p.dr/2/c;  //时间步长
float dt;  //时间步长

float L;
float R;

float ca_posr;     //阴极半径
float ca_posz;    //阴极长度

float wall_p;    //阳极电势
float scr_p;    //屏栅电势
float ca_p;    //阴极电势

float ni;             //原子密度

float KTe;
float KTi;

int Ne;
float qe;
float qi;
int total_time;
int time_inter;
int times=0;
float wegith;
//const int N=10000;
float dti;	
	float *  stac_ex;
	float *  stac_ez;
	float *  stac_Bx;//改成一维的指针 一定要注意
	float *  stac_Bz;
	float * arr; 
	float * arr1;
	float ** ex1;
	float ** iex1;
	float ** ey1;
	float ** iey1;
	float ** ez1;
	float ** iez1;

	float ** hz1;
	float ** ihz1;
	float ** hy1;
	float ** ihy1;
	float ** hx1;
	float ** ihx1;

void initial_para(parameter p)
{

 nz=p.L/p.dz;
 nx=p.R/p.dr;
 nzz=nz+1;
 nxx=nx+1;

 dz=p.dz;  //轴向空间步长
 dr=p.dr; //径向空间步长
 

 dt=p.dr/2/c;  //时间步长

 //dt=5e-11;  //时间步长,变态eps用这个
 //dt=1/sqrt(1/p.dr/p.dr+1/p.dz/p.dz)/c*0.646;
 dti=dt;
 L=p.L;
 R=p.R;

 ca_posr=p.ca_rd;     //阴极半径
 ca_posz=p.ca_len;    //阴极长度

 wall_p=p.anode_p;    //阳极电势
 scr_p=p.screen_p;    //屏栅电势
 ca_p=p.cathode_p;    //阴极电势

 ni=p.ni;             //原子密度
 
 KTe=p.energy_e*boci;
 KTi=p.energy_ion*boci;

 times=1;
 Ne=30;
 //Ne=int(-p.Ib*dt*10/q/p.wegith);
 qe=-p.Ib*times*dt/Ne;
 //qe=p.wegith*q;
 qi=-qe;

 wegith=qe/q;

 total_time=p.total_time;
 time_inter=p.time_inter;
}
float *rds;
float *rds1;


#endif