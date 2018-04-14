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
	float  R;            //放电室半径
	float  L;            //放电室长度
	float  ca_rd;        //阴极半径
	float  ca_len;       //阴极长度

	float  anode_p;         //初始阳极电势
	float  cathode_p;       //初始阴极电势
	float  screen_p;        //初始屏栅电势

	float  ni,Ib;        //原子密度，发射电流
	float  energy_e,energy_ion;  //初始电子能量，初始离子能量
	float  wegith;                   //电子权重
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