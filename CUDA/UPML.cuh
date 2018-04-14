#ifndef _UPML_
#define _UPML_
#include "allocat_space.cuh"
#include "common.cuh"
float * sigmaz1;
float * sigmaz;

float   v=1/sqrt(mur*epsl);
float  a6=(v*dt-dz)/(v*dt+dz);
void initialsigma()
{
	float  m=4;
	float  sigmaMax=-1.5*(m+1)*log(1.0e-10)/(2*pml*dz*sqrt(mur/epsl))*0.5;
	for (int i=0; i<nzz; i++) 
	{ 
		sigmaz1[i]=0;
		sigmaz[i]=0;
	} 
	for (int j=0; j<pml; j++) 
	{ 
		float  thick=(float )(j)/(float )pml;   //初始化sigmaz
	    sigmaz1[nz-pml+j]=sigmaMax*pow(thick,m); 
	} 

	for (int j=0; j<pml; j++) 
	{ 
		float  thick=(float )(j+0.5)/(float )pml;   //初始化sigmaz
	    sigmaz[nz-pml+j]=sigmaMax*pow(thick,m); 
	} 

	for (int j=1; j<pml; j++) 
	{ 
		float  thick=(float )(pml-j)/(float )pml;   //初始化sigmaz
		sigmaz1[j]=sigmaMax*pow(thick,m); 
	} 

	for (int j=0; j<pml; j++) 
	{ 
		float  thick=(float )(pml-j-0.5)/(float )pml;   //初始化sigmaz
		sigmaz[j]=sigmaMax*pow(thick,m); 
	} 
}

#endif