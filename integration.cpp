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
#include <ctime>

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


  	const int H = 2000;
  	const int W = 2000;
	double max_x = 20;
	double max_y = 20;

	std::clock_t start;
	double duration;

  	std::ofstream out("out.ppm");
  	out << "P3\n" << W << ' ' << H << ' ' << "255\n";

	 for (int y1 = 0; y1 < H; ++y1) {
		 start = std::clock();
		 for (int x1 = 0; x1 < W; ++x1) {
        		bool draw = true;
			double x_shot = (double)x1/W*max_x-max_x/2.;
			double y_shot = (double)y1/H*max_y-max_y/2.;
			

			const Vec3 Y_START(100.0, x_shot, y_shot, 0.0, -1.0, 0.00, 0.0, -1.0);
        		Vec3 y = Y_START;
			double t = T_START;
	       		
        		while((t <= TIME_MAXIMUM) && draw) {
  	  			y = y + dy(t,y,DT) ; t += DT;
				double rr = sqrt((y.x)*(y.x)+(y.y)*(y.y))-6.5;
		
				if( rr < 1 && abs(y.z) < 0.2){
					  out << 0 << ' '
		                              << 200 << ' '
	                                      << 0 << '\n';				      
					  draw = false;
				}
			}
			if(draw){
	   		    out << 0 << ' '
        		    << 0 << ' '
        		    << 0 << '\n';
			}
     	     }

	     duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	     cout << " Completed = " <<  (double)y1/(double)H*100 << " % " << endl;	
	     cout << " duration " << duration << " sec" <<endl;
	     cout << " Estmated time to finish " << (double)duration*((double)H-(double)y1)/60.0 << " min " << endl;

	 }

	 return 0;
}
