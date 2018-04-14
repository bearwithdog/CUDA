#ifndef PARAMETER_CUH
#define PARAMETER_CUH

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "readfile.cuh"

/////////////////////////////////////
class parameter  {
	
private:

  readfile  rf;

public:
	float  dr;
	float  dz;
	float  R;            //�ŵ��Ұ뾶
	float  L;            //�ŵ��ҳ���
	float  ca_rd;        //�����뾶
	float  ca_len;       //��������

	float  anode_p;         //��ʼ��������
	float  cathode_p;       //��ʼ��������
	float  screen_p;        //��ʼ��դ����

	float  ni,Ib;        //ԭ���ܶȣ��������
	float  energy_e,energy_ion;  //��ʼ������������ʼ��������
	float  wegith;                   //����Ȩ��
	float  time_inter;                //
    float  total_time ;               

  char      *input_file_name;        // command line input or default value
  char      *output_file_name;

  parameter(char* inputdir);
  void      initial_particle( void );
  void      initial_grid( void );
  void      initial_static_potential(void);
  void      inital_time_set(void);
  

};

#endif