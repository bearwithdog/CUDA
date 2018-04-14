#include<iostream>
#include<fstream>
#include<string>
#include "allocat_space.cuh"
using namespace std;

void static_magnetic()
{
	
	int n=nx,m=nz,i=0;
	float ** B;

	float ** temp1;
	float ** temp2;
	B=two_array_malloc((nx+1)*(nz+1),6);
	for(int i=0;i<(nx+1)*(nz+1);i++)
	{
		for(int j=0;j<6;j++)
		B[i][j]=0;
	}
	stac_Bx=one_array_malloc(nxx*nzz);
	stac_Bz=one_array_malloc(nxx*nzz);
	//stac_Bz=two_array_malloc(nxx,nzz);

	temp1=two_array_malloc(nx+1,nz+1);
	temp2=two_array_malloc(nx+1,nz+1);

	for(int i=0;i<nxx;i++)		
		for(int j=0;j<nzz;j++)
		{
			int tid=i*nzz+j;
			*(stac_Bx+tid)=0;
			*(stac_Bz+tid)=0;
			temp1[i][j]=0;
			temp1[i][j]=0;
		}
	
	ifstream infile;
	infile.open("D:\\PIC\\static_end30.txt",ios::in);
		if(!infile)
		{
			cout<<"静磁场读入失败"<<endl;
		}
		else
		{
			for(int i=0;i<(nx+1)*(nz+1);i=i+1)
			{
				for(int j=0;j<6;j++)
				{
					infile>>B[i][j];
				//	cout<<B[i][j]<<"\t";
				}
			     //   cout<<endl;
			}

		}
		//ofstream on1("E:lfa1.txt");

		int t=0;
	for(int i=0;i<nx+1;i++)		
		for(int j=0;j<nz+1;j++)
		{
			temp1[i][j]=B[t++][3];
		}
	int t1=0;
	for(int i=0;i<nx+1;i++)		
		for(int j=0;j<nz+1;j++)
		{
			temp2[i][j]=B[t1++][5];
		}
	for(int i=1;i<nx+1;i++)		
		for(int j=1;j<nz+1;j++)
		{
			int tid=(i-1)*nzz+(j-1);
			*(stac_Bx+tid)=temp1[i][j]*0.1;
			*(stac_Bz+tid)=temp2[i][j]*0.1;
		}
		//cout<<temp1[10][61]<<stac_Bx[1869]<<endl;
	////cout<<"电位分布如下：\n";
	ofstream on("D:\\de5\\hhzz.txt");
	for(int i=0;i<n+1;i++)
	 {
		 for(int j=0;j<m+1;j++)
			{ int tid=i*nzz+j;
		 on<<*(stac_Bz+tid)<<"\t";}
		     on<<endl;                            //输出结果
	   }
	ofstream on3("D:\\de5\\hhxx.txt");
	for(int i=0;i<n+1;i++)
	 {
		 for(int j=0;j<m+1;j++)
			{ int tid=i*nzz+j;
		 on3<<*(stac_Bx+tid)<<"\t";}
		     on3<<endl;                            //输出结果
	   }
	//cout<<"迭代次数为：\n";
	//cout<<"n="<<k<<endl;

}

