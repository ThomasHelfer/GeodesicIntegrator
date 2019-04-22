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

using namespace std;

#define FOR1(IDX) for (int IDX = 0; IDX < 4; ++IDX)
#define FOR2(IDX1,IDX2) FOR1(IDX1) FOR1(IDX2)
#define FOR3(IDX1,IDX2,IDX3) FOR2(IDX1,IDX2) FOR1(IDX3)
#define FOR4(IDX1,IDX2,IDX3,IDX4) FOR2(IDX1,IDX2) FOR2(IDX3,IDX4)


#define H 50


template <typename data_t>
class render_black_hole{
	public: 
	bool* picture(bool* picture, double  max_x, double  max_y, 
		      const double TIME_MAXIMUM = 150.0, const double DT = 0.1,
		      const double T_START = 0){

		data_t metric;
		std::clock_t start;
		double duration;
		

		auto dy = rk4( metric.eval_diff_eqn ) ;

	 	for (int y1 = 0; y1 < H; ++y1) {
			 start = std::clock();
		  	 for (int x1 = 0; x1 < H; ++x1) {
        			bool draw = true;
				double x_shot = (double)x1/H*max_x-max_x/2.;
				double y_shot = (double)y1/H*max_y-max_y/2.;
			

				const Vec3 Y_START(100.0, x_shot, y_shot, 0.0, -1.0, 0.00, 0.0, -1.0);
        			Vec3 y = Y_START;
				double t = T_START;
	       		
        			while((t <= TIME_MAXIMUM) && draw) {
  	  				y = y + dy(t,y,DT) ; t += DT;
					double rr = sqrt((y.x)*(y.x)+(y.y)*(y.y))-6.5;
		
					if( rr < 1 && abs(y.z) < 0.2){
						picture[y1*H+x1] = true;
						draw = false;
					}	
				}
				if(draw){
			    		picture[y1*H+x1] = false;
				}
     	     		}

	     		duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	     		cout << " Completed = " <<  (double)y1/(double)H*100 << " % " << endl;	
	     		cout << " duration " << duration << " sec" <<endl;
	     		cout << " Estmated time to finish " << (double)duration*((double)H-(double)y1)/60.0 << " min " << endl;

	 	}
	 	return picture;
	}
	
	void render(bool *picture, char file_name[]){
  	  std::ofstream out("out.ppm");
  	  out << "P3\n" << H << ' ' << H << ' ' << "255\n";

	  for(int i = 0; i < H*H; i++){
		if(picture[i]){
			    out << 0 << ' '
                            << 255 << ' '
                            << 0 << '\n';
		}else{
			 out << 0 << ' '
                            << 0 << ' '
                            << 0 << '\n';
		}
	  }
	}
};


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



 	char name_render[] = "out.ppm";	
	render_black_hole<Black_Hole> rend;
	const double size_x = 20;
	const double size_y = 20;
	bool pic[H*H] ; 

	rend.picture(pic, size_x,size_y);	

	rend.render(pic,name_render);


	 return 0;
}
