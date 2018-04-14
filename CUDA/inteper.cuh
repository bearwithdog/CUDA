#include "iostream"
#include "math.h"
#include "ostream"
#include "vector"
#include "fstream"
#include "sstream"
#include "string"
#include <iomanip>
#include<conio.h>
#include<algorithm>
#include"allocat_space.cuh"
using namespace std;


void read_cross(string filename)
{   
	arr=one_array_malloc(1000);
	arr1=one_array_malloc(1000);
	int i=0;
	string  line;
	float tmp,tmp1;
	ifstream infile(filename.c_str(),ios::in);
	if(!infile.is_open())
	{
		cout<<"error openfile"<<endl;
	} 
	while(!infile.eof())
	{
		getline(infile,line);
		if(line.empty())
			continue; 
        istringstream input(line);
		input>>tmp;
		arr[i]=tmp;
		input>>tmp1;
		arr1[i]=tmp1;
		i++;
	}

}

float cross_section(float t)  
{  

	float *x=arr;
	float *y=arr1;
	int n=1000;
    int i,k;  
    float u,f;  
    for(i=0;i<=(n-2);i=i+1)  
    {  
        if(t<=x[i+1])  
        {  
            k=i;  
            break;  
        }  
        else  
            k=n-2;  
    }  
    u=(t-x[k])/(x[k+1]-x[k]);  
    f=y[k]+u*(y[k+1]-y[k]);  
    return (f);  
}  