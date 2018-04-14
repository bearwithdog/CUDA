#ifndef GRID_CUH
#define GRID_CUH  
class Grid	
{
public:
    float ex,ey,ez,hx,hy,hz;
	float ave_ex,ave_ey,ave_ez,ave_hx,ave_hy,ave_hz;
	float jr,jy,jz,Q,den,jc,jc1,jd1,jd,midu;
	float jr_ion,jy_ion,jz_ion,Q_ion,den_ion;
	int ne[3];
	float Pmax;
};
#endif