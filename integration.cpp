/*
 * compiled with gcc 5.4:
 * g++-mp-5 -std=c++14 -o rk4 rk4.cc
 *
 */
# include <iostream>
#include <fstream>
# include <math.h>
# include <vector>
# include <algorithm> 
#include "tensor.hpp"
#include "schwarzschild.hpp"
using namespace std;


#define FOR1(IDX) for (int IDX = 0; IDX < 4; ++IDX)
#define FOR2(IDX1,IDX2) FOR1(IDX1) FOR1(IDX2)
#define FOR3(IDX1,IDX2,IDX3) FOR2(IDX1,IDX2) FOR1(IDX3)
#define FOR4(IDX1,IDX2,IDX3,IDX4) FOR2(IDX1,IDX2) FOR2(IDX3,IDX4)


auto rk4(Vec3 f(double, Vec3))
{
        return
        [       f            ](double t, Vec3 y, double dt ) -> Vec3{ return
        [t,y,dt,f            ](                    Vec3 dy1) -> Vec3{ return
        [t,y,dt,f,dy1        ](                    Vec3 dy2) -> Vec3{ return
        [t,y,dt,f,dy1,dy2    ](                    Vec3 dy3) -> Vec3{ return
        [t,y,dt,f,dy1,dy2,dy3](                    Vec3 dy4) -> Vec3{ return
        ( dy1 + dy2*2 + dy3*2 + dy4 ) / 6   ;} (
        f( t+dt  , y+dy3   )*dt          );} (
        f( t+dt/2, y+dy2/2 )*dt          );} (
        f( t+dt/2, y+dy1/2 )*dt          );} (
        f( t     , y       )*dt          );} ;
}




int main(void)
{

        const double TIME_MAXIMUM = 200.0, WHOLE_TOLERANCE = 1e-12 ;
        const double T_START = 0.0, DT = 0.001;
	Black_Hole BH;

	auto dy = rk4( BH.eval_diff_eqn ) ;

        for(int i = 0;i<5;i++){			

		const Vec3 Y_START(20.0, 7.0+i/2., 0.0, 0.0, -0.3, 0.00, 0.0, -1.0);
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

  	return 0;
}
