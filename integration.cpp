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

auto rk4(double f(double, double))
{
        return
        [       f            ](double t, double y, double dt ) -> double{ return
        [t,y,dt,f            ](                    double dy1) -> double{ return
        [t,y,dt,f,dy1        ](                    double dy2) -> double{ return
        [t,y,dt,f,dy1,dy2    ](                    double dy3) -> double{ return
        [t,y,dt,f,dy1,dy2,dy3](                    double dy4) -> double{ return
        ( dy1 + 2*dy2 + 2*dy3 + dy4 ) / 6   ;} (
        dt * f( t+dt  , y+dy3   )          );} (
        dt * f( t+dt/2, y+dy2/2 )          );} (
        dt * f( t+dt/2, y+dy1/2 )          );} (
        dt * f( t     , y       )          );} ;
}

void printVector(vector<double> v) 
{ 
    // lambda expression to print vector 
    for_each(v.begin(), v.end(), [](int i) 
    { 
        std::cout << i << " "; 
    }); 
    cout << endl; 
} 

int main(void)
{
	Vec3 arr1( 4, 5, 8);
	arr1.print();

        const double TIME_MAXIMUM = 10.0, WHOLE_TOLERANCE = 1e-12 ;
        const double T_START = 0.0, Y_START = 1.0, DT = 0.10;

//	auto test 	   = [		     ](double t, double y)->Vec3{ return Vec3(t,y,y);};
        auto eval_diff_eqn = [               ](double t, double y)->double{ return t*sqrt(y)                         ; } ;
        auto eval_solution = [               ](double t          )->double{ return pow(t*t+4,2)/16                   ; } ;
        auto find_error    = [eval_solution  ](double t, double y)->double{ return fabs(y-eval_solution(t))          ; } ;
        auto is_whole      = [WHOLE_TOLERANCE](double t          )->bool  { return fabs(t-round(t)) < WHOLE_TOLERANCE; } ;


        auto dy = rk4( eval_diff_eqn ) ;
 
        double y = Y_START, t = T_START ;
 
        while(t <= TIME_MAXIMUM) {
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y,find_error(t,y)); }
          y += dy(t,y,DT) ; t += DT;
        }
        return 0;
}
