
#include <iostream>
//#include "grid.h"
#ifndef _ALLOCAT_CUH
#define _ALLOCAT_CUH
float  *one_array_malloc(int n)   //一维数组分配
{
    float  *a;
    a=new float [n];
    return a;
}
float  **two_array_malloc(int n,int m)   //二维数组分配
{
    float ** a=new float *[n];
	for(int i=0;i<n;i++)
		{
			a[i]=new float [m];
		}
	return a;
}

#endif