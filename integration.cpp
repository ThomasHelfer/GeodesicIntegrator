/*
 * compiled with gcc 5.4:
 * g++-mp-5 -std=c++14 -o rk4 rk4.cc
 *
 */
# include <iostream>
# include <math.h>
# include <vector>
# include<algorithm> 
using namespace std;

struct Vec3 {
  double x,y,z;
  Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
  Vec3 operator + (const Vec3& v) const { return Vec3(x+v.x, y+v.y, z+v.z); }
  Vec3 operator - (const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z); }
  Vec3 operator * (double d) const { return Vec3(x*d, y*d, z*d); }
  Vec3 operator / (double d) const { return Vec3(x/d, y/d, z/d); }
  Vec3 sqrt3(Vec3& v) const { return Vec3(sqrt(v.x),sqrt(v.y),sqrt(v.z)); };
  Vec3 print() const {
	  std::cout <<"x " << x << "\n";
	  std::cout <<"y " << y << "\n";
	  std::cout <<"z " << z << "\n";
  };
  Vec3 normalize() const {
    double mg = sqrt(x*x + y*y + z*z);
    return Vec3(x/mg,y/mg,z/mg);
  }
};

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

        const double TIME_MAXIMUM = 10.0, WHOLE_TOLERANCE = 1e-12 ;
        const double T_START = 0.0, DT = 0.01;
	const Vec3 Y_START(1.0, 1.0, 1.0);

        auto eval_diff_eqn = [               ](double t, Vec3 y)->Vec3{ return y.sqrt3(y)*t                        ; } ;
        auto eval_solution = [               ](double t          )->double{ return pow(t*t+4,2)/16                   ; } ;
        auto find_error    = [eval_solution  ](double t, double y)->double{ return fabs(y-eval_solution(t))          ; } ;
        auto is_whole      = [WHOLE_TOLERANCE](double t          )->bool  { return fabs(t-round(t)) < WHOLE_TOLERANCE; } ;

	auto dy = rk4( eval_diff_eqn ) ;
 
        Vec3 y = Y_START;
	double t = T_START;
        while(t <= TIME_MAXIMUM) {
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.x,find_error(t,y.x)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.y,find_error(t,y.y)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.z,find_error(t,y.z)); }
          y = y + dy(t,y,DT) ; t += DT;
        }

  	return 0;
}
