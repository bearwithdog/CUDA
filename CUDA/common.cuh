#pragma once
#ifndef _COMMON_CUH
#define _COMMON_CUH
#include "parameter.cuh"
//#include "grid.h"


const float pi=3.1415926;        

const float c=2.99792458e8;      //����   

const float q=-1.60e-19;         // ���Ӵ������ 

const float boci=1.6022e-19;     //��������������ת�����1ev=boci

const float Me=9.1e-31;          //��������

const float Mi=2.18e-25;         //��������

const float  Qm=-1.7588e+11;     //���Ӻ��ʱ�

const float  Qm_ion=7.33945e+5;     //���Ӻ��ʱ�

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

float dz;  //����ռ䲽��
float dr; //����ռ䲽��

//float dt=p.dr/2/c;  //ʱ�䲽��
float dt;  //ʱ�䲽��

float L;
float R;

float ca_posr;     //�����뾶
float ca_posz;    //��������

float wall_p;    //��������
float scr_p;    //��դ����
float ca_p;    //��������

float ni;             //ԭ���ܶ�

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
	float *  stac_Bx;//�ĳ�һά��ָ�� һ��Ҫע��
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

 dz=p.dz;  //����ռ䲽��
 dr=p.dr; //����ռ䲽��
 

 dt=p.dr/2/c;  //ʱ�䲽��

 //dt=5e-11;  //ʱ�䲽��,��̬eps�����
 //dt=1/sqrt(1/p.dr/p.dr+1/p.dz/p.dz)/c*0.646;
 dti=dt;
 L=p.L;
 R=p.R;

 ca_posr=p.ca_rd;     //�����뾶
 ca_posz=p.ca_len;    //��������

 wall_p=p.anode_p;    //��������
 scr_p=p.screen_p;    //��դ����
 ca_p=p.cathode_p;    //��������

 ni=p.ni;             //ԭ���ܶ�
 
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