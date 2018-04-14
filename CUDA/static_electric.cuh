//考虑左下角
#include<iostream>
#include<fstream>
#include<math.h>
#include "allocat_space.cuh"
using namespace std;

void static_electric()
{
	int n=nx,m=nz,i=0;
	int fl=ca_posz/dz;
	int posr=ca_posr/dr;
	int posz=40;
	//float  temp=0.0,error=1.0;
	float ** A;
	int k=0;
	stac_ex=one_array_malloc(nxx*nzz);
	stac_ez=one_array_malloc(nxx*nzz);
	for(int i=0;i<nx+1;i++)		
		for(int j=0;j<nz+1;j++)
		{
			int tid=i*nzz+j;
			*(stac_ex+tid)=0;
			*(stac_ez+tid)=0;
		}
	A=two_array_malloc(n+2,m+2);

	for(int i=0;i<n+2;i++ )
	  for (int j=0;j<m+2;j++)
	 {
		  A[i][j]=0;
	 }

	  for(int i=0;i<n+2;i++ )
		{
			if(i>=posz)
		  A[i][m+1]=wall_p;
			else if(i<posz)
          A[i][m+1]=scr_p; 
	    }
	  for(int i=0;i<m+2;i++ )
		  A[n+1][i]=wall_p;

      for(int i=posr;i<n+2;i++ )
		  A[i][0]=wall_p;
	   for(int i=0;i<posr;i++ )
		   for(int j=0;j<fl;j++)
		{
		     A[i][j]=ca_p;
	    }
	 //for(int i=0;i<posr;i++ )
	 // for(int j=0;j<fl;j++)
		//{
		//	 A[i][j]=0;
		//}

	  float  error=1,temp=0;
	while(error>(1e-3))
	{
		k=k+1;error=0;
	for(int i=1;i<n+1;i++)
	{	
		for(int j=1;j<m+1;j++)
		{
			if((i>0&&i<=posr&&j>=fl)||(i>posr&&j>0))
			{
			temp=(A[i-1][j]+A[i+1][j]+A[i][j+1]+A[i][j-1])/4.0;					
			error=max(error,abs(temp-A[i][j]));
			A[i][j]=temp;	
			}	
		}	
	 }
	 for(int j=1;j<m+1;j++)
	{	
		int i=0;
	   A[0][j]=(4*A[i+1][j]+A[i][j+1]+A[i][j-1])/6.0;//轴线上的点处理
	 }
	// cout<<error<<endl;
	}

	for(int i=1;i<n+1;i++)		
		for(int j=1;j<m+1;j++)
		{ 	
			if((i>0&&i<=posr&&j>=fl)||(i>posr&&j>0))
			{int tid=(i-1)*nzz+j-1;
			*(stac_ex+tid)=-(A[i+1][j]-A[i-1][j])/2/dr;}
			
		}
	for(int i=1;i<n+1;i++)		
		for(int j=1;j<m+1;j++)
		{ 	
			if((i>0&&i<=posr&&j>=fl)||(i>posr&&j>0))
			{
			int tid=(i-1)*nzz+j-1;
			*(stac_ez+tid)=-(A[i][j+1]-A[i][j-1])/2/dz;
			}
		}
		 ofstream onn("D:\\de5\\esh.txt");
	  for(int i=0;i<n+2;i++)
	 {
		 for(int j=0;j<m+2;j++)
			 onn<<A[i][j]<<"\t";
		     onn<<endl;                            //输出结果
	   }
	  	ofstream on1("D:\\de5\\exx.txt");
	for(int i=0;i<n+1;i++)
	 {
		 for(int j=0;j<m+1;j++)
		 { 
			 int tid=i*nzz+j;
			 on1<<*(stac_ex+tid)<<"\t";
		 }
		     on1<<endl;                            //输出结果
	   }
	ofstream on2("D:\\de5\\ezz.txt");
	for(int i=0;i<n;i++)
	 {
		 for(int j=0;j<m;j++)
			 {
			int tid=i*nzz+j;
			 on2<<*(stac_ez+tid)<<"\t";
			}
		     on2<<endl;                            //输出结果
	   }
}




