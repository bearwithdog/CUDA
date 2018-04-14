#include "iostream"
#include "math.h"
#include "fstream"
#include "ostream"
#include "sstream"
#include <iomanip>
#include <vector>
#include <algorithm>
#include <time.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda.h"
#include"device_functions.h"
#include "common.cuh"
#include"grid.cuh"
#include "particle.cuh"
#include "UPML.cuh"
#include "inteper.cuh"
#include "static_electric.cuh"
#include "static_magnetic.cuh"
#include "lty_paralle_common.cuh"
#include "curand_kernel.h"
#include <stdlib.h>
#include <stdio.h>
#include "curand.h"
#include"variable.cuh"
using namespace std;

int main()
{
	int tmp = 5000;
	cudaDeviceReset();
	curandnumber(tmp);
	void data_save(Paticle *pat_elc,Paticle *pat_ion);
	void data_save_ion(Paticle *pat_elc, Paticle *pat_ion);//最后时刻粒子信息
	//void data_save(Paticle *pat_elc);
	void data_save(Grid* G_GPU);
	void current_save(Grid* G);
	char* indir="..\\data\\input.txt";
	parameter p(indir);	
	initial_para(p);//初始化参数
	//cout<<qe<<" "<<Me<<endl;
	 sigmaz1=one_array_malloc(nzz);
	 sigmaz =one_array_malloc(nzz);
	 G_gpu = new Grid[nxx*nzz];//仅仅用于验证场的正确性 数据回拷
	 pat_elc=new Paticle[lizi_count];//主机上的粒子信息，用于接受GPU算完之后的结果
	 pat_ion=new Paticle[lizi_count];//主机上的粒子信息，用于接受GPU算完之后的结果
	 cudaError_t cudaStatus;
	cudaMalloc((void**)&device_G,nxx*nzz*sizeof(Grid));
	cudaMalloc((void**)&dev_Gtemp, nxx*nzz*sizeof(Grid));
	cudaMalloc((void**)&device_Gn,nxx*nzz*sizeof(Grid));
	//cudaMalloc((void**)&device_tail,sizeof(int));//用于后移线程索引
	cudaStatus=cudaMalloc((void**)&d_pat_elc,lizi_count*sizeof(Paticle));
	cudaStatus = cudaMalloc((void**)&d_pre_elc, lizi_count*sizeof(Pre_Paticle));//总的模拟粒子数目 GPU上
	cudaStatus=cudaMalloc((void**)&d_pat_ion,lizi_count*sizeof(Paticle));
	cudaStatus = cudaMalloc((void**)&d_pre_ion, lizi_count*sizeof(Pre_Paticle));//总的模拟粒子数目 GPU上
		if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMemcpy failed!");
       }
	cudaStatus=cudaMalloc((void**)&d_stac_Bx,nxx*nzz*sizeof(float));
	cudaStatus=cudaMalloc((void**)&d_stac_Bz,nxx*nzz*sizeof(float));
	cudaStatus=cudaMalloc((void**)&d_stac_ex,nxx*nzz*sizeof(float));
	cudaStatus=cudaMalloc((void**)&d_stac_ez,nxx*nzz*sizeof(float));
	read_cross("E:\\cross_section_net.txt");
	 static_magnetic();
	 static_electric();
	 initialsigma();
	cudaStatus=cudaMemcpy(d_stac_ex,stac_ex,nxx*nzz*sizeof(float),cudaMemcpyHostToDevice);
	cudaStatus=cudaMemcpy(d_stac_ez,stac_ez,nxx*nzz*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(d_stac_Bx,stac_Bx,nxx*nzz*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(d_stac_Bz,stac_Bz,nxx*nzz*sizeof(float),cudaMemcpyHostToDevice);
	device_initialchang<<<block,thread>>>(device_G,device_Gn);
	kernel_L_InitialPML<<<block,thread>>>(nxx, nzz);
	cudaStatus = cudaMalloc((void**)&d_sigmaz1, nzz*sizeof(float));
	cudaStatus = cudaMalloc((void**)&d_sigmaz, nzz*sizeof(float));
	cudaStatus = cudaMemcpy(d_sigmaz1, sigmaz1, nzz*sizeof(float), cudaMemcpyHostToDevice);
	cudaStatus = cudaMemcpy(d_sigmaz, sigmaz, nzz*sizeof(float), cudaMemcpyHostToDevice);
	clock_t start=0, end=0;
	
	 ofstream time("D:\\PIC\\time.txt",ios::app);
	 
	initial_always<<<block,thread>>>(d_pat_elc,lizi_count,DS_number,tmp);//电子初始化
	initial_always << <block, thread >> >(d_pat_ion, lizi_count, DS_number,tmp);//离子初始化
	device_define_G << <block,thread >> >(nxx, nzz, device_Gn, device_G);
	start = clock();
	 for(int t=0;t<1;t++)
	 {
		 host_temptail+=30;
		//device_ave_field<<<block,thread>>>(nxx,nzz,device_G,device_Gn);
		device_update_last << <block, thread >> >(d_stac_Bx, d_stac_Bz, d_stac_ex, d_stac_ez, d_pat_elc,d_pre_elc,device_G,DS_number,host_temptail);
		//device_update_ion << <block, thread >> >(d_stac_Bx, d_stac_Bz, d_stac_ex, d_stac_ez, d_pat_ion,d_pre_ion,device_G,DS_number,host_temptail);
		//cudaMemcpy(wuchafenxi, (d_pat_elc), 30 * sizeof(Paticle), cudaMemcpyDeviceToHost);
		//data_save(wuchafenxi);//每个步长考出30个粒子信息作为比对，实际可以不要
		current << <block, thread >> >(d_pat_elc, d_pre_elc, device_G, DS_number, host_temptail,t);
		//current_ion << <block, thread >> >(d_pat_ion, d_pre_ion, device_G, DS_number, host_temptail,t);
		/*cudaMemcpy(G_gpu, device_G, nxx*nzz*sizeof(Grid), cudaMemcpyDeviceToHost);
		current_save(G_gpu);*/
		cacuchang_hx<<<block,thread>>>(device_G,device_Gn,d_sigmaz1,d_sigmaz,dt,dr,dz,nxx,nzz);
		cacuchang_hy << <block,thread >> >(device_G, device_Gn, d_sigmaz1, d_sigmaz, dt, dr, dz, nxx, nzz);
		cacuchang_hz<< <block, thread >> >(device_G, device_Gn, d_sigmaz1, d_sigmaz, dt, dr, dz, nxx, nzz);
		cacuchang_ex << <block, thread >> >(device_G, device_Gn, d_sigmaz1, d_sigmaz, dt, dr, dz, nxx, nzz);
		cacuchang_ey << <block, thread >> >(device_G, device_Gn, d_sigmaz1, d_sigmaz, dt, dr, dz, nxx, nzz);
		cacuchang_ez << <block, thread >> >(device_G, device_Gn, d_sigmaz1, d_sigmaz, dt, dr, dz, nxx, nzz);
	    device_define_G<<<block,thread>>>(nxx,nzz,device_Gn,device_G);

		/*end = clock();
		cout<<t<<"\t"<<"\t"<<(float)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;
        time<<t<<"  Szabo1 Run time:  "<<(float)(end - start) / CLOCKS_PER_SEC<<"S"<<endl;*/
	 }
	 end = clock();
	 cout<< (float)(end - start) / CLOCKS_PER_SEC << "S" << endl;
	 cudaMemcpy(G_gpu,device_G, nxx*nzz*sizeof(Grid), cudaMemcpyDeviceToHost);
	 data_save(G_gpu);//最后时刻场信息
	 cudaStatus=cudaMemcpy(pat_elc,d_pat_elc,lizi_count*sizeof(Paticle),cudaMemcpyDeviceToHost);
	 cudaStatus=cudaMemcpy(pat_ion,d_pat_ion,lizi_count*sizeof(Paticle),cudaMemcpyDeviceToHost);
	 data_save(pat_elc,pat_ion);//最后时刻粒子信息
	 data_save_ion(pat_ion,pat_elc);//最后时刻粒子信息
	cudaFree(device_G);
	cudaFree(device_Gn);
	cudaFree(d_pat_elc);
	cudaFree(d_pat_ion);
	cudaFree(d_rds);
	cudaFree(d_rds1);
	cudaFree(d_stac_Bx);
	cudaFree(d_stac_Bz);
	cudaFree(d_stac_ex);
	cudaFree(d_stac_ez);
	system("pause");
	return 0;
}

void data_save(Paticle *pat_elc,Paticle *pat_ion)
{
	
	  /*************************************************/
		time_inter=1000;
		//if(t%time_inter==0)//time_inter
		
		char s_ele[1000];  
		sprintf_s(s_ele,"D:\\PIC\\输出e\\bingxing1.txt");
		ofstream on_ele(s_ele);	
	    for(int i=0;i<host_temptail;i++)
 		 {		
	        on_ele<<pat_elc[i].pr<<"\t"<<pat_elc[i].py<<"\t"<<pat_elc[i].pz<<"\t"<<pat_elc[i].vr<<"\t"<<pat_elc[i].vy<<"\t"<<pat_elc[i].vz<<endl;
		 }  	
}
void data_save_ion(Paticle *pat_elc, Paticle *pat_ion)
{

	/*************************************************/
	time_inter = 1000;
	//if(t%time_inter==0)//time_inter
	//{
	char s_ele[1000];
	sprintf_s(s_ele, "D:\\PIC\\输出e\\bingxing2.txt");
	ofstream on_ele(s_ele);
	for (int i = 0; i<host_temptail; i++)
	{
		on_ele << pat_elc[i].pr << "\t" << pat_elc[i].py << "\t" << pat_elc[i].pz << "\t" << pat_elc[i].vr << "\t" << pat_elc[i].vy << "\t" << pat_elc[i].vz << endl;
	}
}
void data_save(Paticle *pat_elc)
{

	/*************************************************/
	time_inter = 1000;
	char s_ele[1000];
	sprintf_s(s_ele, "D:\\PIC\\输出e\\wuchafenxi0.txt");
	ofstream on_ele(s_ele, ofstream::app);
	for (int i = 0; i<30; i++)
	{
		on_ele << pat_elc[i].pr << "\t" << pat_elc[i].py << "\t" << pat_elc[i].pz << "\t" << pat_elc[i].vr << "\t" << pat_elc[i].vy << "\t" << pat_elc[i].vz << endl;
	}

}
void data_save(Grid* G_GPU)
{

	ofstream out_chang("D:\\PIC\\输出e\\chang_gpu.txt");
	if (out_chang)
	{
		for (int i = 0; i < nxx*nzz; i++)
			out_chang << G_GPU[i].ex << " " << G_GPU[i].ey << " " << G_GPU[i].ez << " "
			<< G_GPU[i].hx << " " << G_GPU[i].hy << " " << G_GPU[i].hz << endl;
	}
}
void current_save(Grid* G_GPU)
{
	ofstream out_chang("D:\\PIC\\输出e\\current_gpu.txt");
	if (out_chang)
	{
		for (int i = 0; i < nxx*nzz; i++)
		out_chang << G_GPU[i].jr << " " << G_GPU[i].jy << " " << G_GPU[i].jz <<endl;
	}
}
//void curandnumber(int n)
//{
//	cudaMalloc((void**)&DS_number,n*sizeof(float));
//    HS_number=one_array_malloc(n);
//	srand(1);
//	/*for (int i = 0; i < n; i++)
//	{
//		HS_number[i] = (float)rand() / ((float)RAND_MAX);
//	}*/
//	
//float HS_number[30] = { 0.00125126,0.193304,0.585009,
//		              0.350291,0.82284,0.174108,0.710501,0.303995,0.0914029,
//		              0.147313,0.988525,0.119083,0.0089114,0.531663,0.601764,
//		              0.166234,0.450789,0.0570391,0.783319,0.519883,0.875973,
//					  0.955901,0.539354,0.462081,0.862239,0.779656,0.996796,0.611499,
//					  0.266213,0.840144 };
//cudaMemcpy(DS_number,HS_number,n*sizeof(float),cudaMemcpyHostToDevice);    
//}
//5000个随机数
void curandnumber(int n)
{
	cudaMalloc((void**)&DS_number, n*sizeof(float));
	HS_number = one_array_malloc(n);
	srand(1);
	for (int i = 0; i < n; i++)
	{
		HS_number[i] = (float)rand() / ((float)RAND_MAX);
	}
	cudaMemcpy(DS_number, HS_number, n*sizeof(float), cudaMemcpyHostToDevice);
}

