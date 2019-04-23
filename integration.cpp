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


int main(void)
{

        const double TIME_MAXIMUM = 150.0, WHOLE_TOLERANCE = 1e-12 ;
        const double T_START = 0.0, DT = 0.1;
	Black_Hole BH;

	auto dy = rk4( BH.eval_diff_eqn ) ;

        for(int i = 0;i<5;i++){			

		const Vec3 Y_START(100.0, 2.6+i/5., 0.0, 0.0, -1.0, 0.00, 0.0, -1.0);
        	Vec3 y = Y_START;
		double t = T_START;

		cout << "Norm of vector "<< BH.calculate_norm(Y_START) << endl; 
		ofstream myfile("xpos"+std::to_string(i)+".csv");

        	while(t <= TIME_MAXIMUM) {
  	  		y = y + dy(t,y,DT) ; t += DT;			

			myfile << y.x <<"	" << y.y << "	" << y.z << "	" << y.t << "	" << BH.calculate_norm(y) <<   endl;
		}

		myfile.close();

	}



 	string name_render = "out.ppm";	
	render_black_hole<Black_Hole> rend;
	const double size_x = 20;
	const double size_y = 20;
	int pic[H*H] ; 

	const Vec3 center(100.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0);

	rend.picture(pic,center, size_x,size_y);	

	rend.render(pic,name_render);


	 return 0;
}
