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

using namespace std;

#define FOR1(IDX) for (int IDX = 0; IDX < 4; ++IDX)
#define FOR2(IDX1,IDX2) FOR1(IDX1) FOR1(IDX2)
#define FOR3(IDX1,IDX2,IDX3) FOR2(IDX1,IDX2) FOR1(IDX3)
#define FOR4(IDX1,IDX2,IDX3,IDX4) FOR2(IDX1,IDX2) FOR2(IDX3,IDX4)

template <typename data_t> 
class geodesic_shooter{
	public:
	void shoot(Vec3 center, double shift = 1/5., int shoot = 10 ,  const double TIME_MAXIMUM = 150.0, const double T_START = 0.0, const double DT = 0.1){

		data_t metric;
		auto dy = rk4( metric.eval_diff_eqn ) ;

        	for(int i = 0; i<shoot ; i++){			

			const Vec3 Y_START(0.0, 7+i/5., 0.0, 0.0, -1.0, 0.00, 0.0, -1.0);
        		Vec3 y = Y_START + center;
			double t = T_START;

			// ========= Preparing output ==========
 			string imgname="xpos";
                	char cbuff[20];
                	sprintf (cbuff, "%03d", i);
                	imgname.append(cbuff);
                	imgname.append(".csv");
			ofstream myfile(imgname);

			// ========== Integration and output ===
        		while(t <= TIME_MAXIMUM) {
  	  			y = y + dy(t,y,DT) ; t += DT;			

				myfile << y.x <<"	" << y.y << "	" << y.z << "	" << y.t << "	" << metric.calculate_norm(y) <<   endl;
			}

			// ========== clean up ==================
			myfile.close();

		}
	
	}


};

int main(void)
{
	// Setting up inital data
		
	double center_x = 100;
	double center_y = 4.6;
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
		
	pewpew.shoot(initial_data);


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
