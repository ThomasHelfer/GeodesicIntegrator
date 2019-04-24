/*
 * compiled with gcc 5.4:
 * g++-mp-5 -std=c++14 -o rk4 rk4.cc
 *
 */
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm> 
#include <ctime>

#include "tensor.hpp"
#include "schwarzschild.hpp"
#include "rk4.hpp"
#include "render.hpp"
#include "geodesic_shooter.hpp"

using namespace std;

#define FOR1(IDX) for (int IDX = 0; IDX < 4; ++IDX)
#define FOR2(IDX1,IDX2) FOR1(IDX1) FOR1(IDX2)
#define FOR3(IDX1,IDX2,IDX3) FOR2(IDX1,IDX2) FOR1(IDX3)
#define FOR4(IDX1,IDX2,IDX3,IDX4) FOR2(IDX1,IDX2) FOR2(IDX3,IDX4)

int main(void)
{
	// Setting up inital data
		
	double center_x = 100;
	double center_y = 2;
	double center_z = 0.0;
	double start_time = 0.0;
	double velocity_x = -1.0;
	double velocity_y = 0.0;
	double velocity_z = 0.0;
	double lapse  = -1.0;

	const Vec3 initial_data(center_x, center_y, center_z,
		                start_time,
		                velocity_x, velocity_y, velocity_z, 
			         lapse);

	geodesic_shooter<Black_Hole> pewpew;
		
	pewpew.shoot(initial_data,0.3,50);


 	string name_render = "out.ppm";	
	render_black_hole<Black_Hole> rend;
	const double size_x = 20;
	const double size_y = 20;
	int pic[H*H] ; 

	const Vec3 center(100.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0);

	rend.picture(pic,center,size_x,size_y);	

	rend.render(pic,name_render);


	 return 0;
}
