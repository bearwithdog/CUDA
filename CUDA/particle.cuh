
#ifndef PARTICLE_CUH
#define PARTICLE_CUH

#include <vector>
class Paticle
{
public:
	float vr,vy,vz;
	float pr,py,pz;
	int species;
	int blei,blek;	
};
class Pre_Paticle
{
public:
	float pr, pz;
};
#endif