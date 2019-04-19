/*
 * compiled with gcc 5.4:
 * g++-mp-5 -std=c++14 -o rk4 rk4.cc
 *
 */
# include <iostream>
# include <math.h>
# include <vector>
# include <algorithm> 
using namespace std;


#define FOR1(IDX) for (int IDX = 0; IDX < 4; ++IDX)
#define FOR2(IDX1,IDX2) FOR1(IDX1) FOR1(IDX2)
#define FOR3(IDX1,IDX2,IDX3) FOR2(IDX1,IDX2) FOR1(IDX3)
#define FOR4(IDX1,IDX2,IDX3,IDX4) FOR2(IDX1,IDX2) FOR2(IDX3,IDX4)

struct Vec3 {
  double x,y,z,t, vx,vy,vz,vt;
  Vec3(double x, double y, double z, double t,
       double vx, double vy, double vz, double vt ) : 
	  x(x), y(y), z(z), t(t), 
	  vx(vx), vy(vy), vz(vz), vt(vt) {}

  Vec3 operator + (const Vec3& v) const { return Vec3(x+v.x, y+v.y, z+v.z, t+v.t,
		  				      vx+v.vx, vy+v.vy, vz+v.vz, vt+v.vt); }
  Vec3 operator - (const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z, t-v.t,
		  				      vx-v.vx, vy-v.vy, vz-v.vz, vt-v.vt); }
  Vec3 operator * (double d) const { return Vec3(x*d, y*d, z*d, t*d,
		  				 vx*d, vy*d, vz*d, vt*d ); }
  Vec3 operator / (double d) const { return Vec3(x/d, y/d, z/d, t/d,
		  				 vx/d, vy/d, vz/d, vt/d ); }
  Vec3 sqrt3(Vec3& v) const { return Vec3(sqrt(v.x),sqrt(v.y),sqrt(v.z),sqrt(v.t),
						sqrt(v.vx),sqrt(v.vy),sqrt(v.vz),sqrt(v.vt) );};
  Vec3 print() const {
	  std::cout <<"x " << x << "\n";
	  std::cout <<"y " << y << "\n";
	  std::cout <<"z " << z << "\n";
	  std::cout <<"t " << z << "\n";
  };
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

class Black_Hole{
	public : 
	static Vec3 eval_diff_eqn(double t, Vec3 v){

		double mass = 1;
		double g[4][4] = {}; // Metix Index low low 
		double g_spher[4][4] = {}; // Metix Index low low 
		double dg[4][4][4]  = {}; //
		double chris[4][4][4] = {}; // Christoffel index high low low 
		double jacobian[4][4] = {};

		double x[4] = {v.x,v.y,v.z,v.t};
		double dx[4] = {v.vx,v.vy,v.vz,v.vt};
		
		double rr2 =   pow(v.x,2)+ pow(v.y,2) + pow(v.z,2);
  		double rho2 =  pow(v.x,2)+ pow(v.y,2);

		double rr = sqrt(rr2);
  		double rho = sqrt(rho2);
 		 // sinus(theta)
  		double sintheta = rho/rr;
  		// cos(phi)
  		double cosphi = v.x/rho ; 
 	 	// sin(phi)
 		double sinphi = v.y/rho ; 

		// ==================== 

		g_spher[0][0] = 1./(1.0 - 2.0*mass/rr);
  		// Define theta theta component 
  		g_spher[1][1] = rr2;
  		// Define phi phi component 
  		g_spher[2][2] = rr2*pow(sintheta,2);
		// Devine tt component 
		g_spher[3][3] = (1. - 2.0*mass/rr);


		jacobian[0][0] = v.x/rr ;  
  		jacobian[1][0] = cosphi*v.z/rr2 ; 
  		jacobian[2][0] = -v.y/rho2 ; 
  		jacobian[0][1] = v.y/rr ; 
  		jacobian[1][1] = sinphi*v.z/rr2 ; 
  		jacobian[2][1] = v.x/rho2 ; 
  		jacobian[0][2] = v.z/rr ; 
  		jacobian[1][2] = -rho/rr2 ;  
  		jacobian[2][2] = 0.0 ; 	
  		jacobian[3][3] = 1.0 ; 	

		// ====================

 		FOR2(i,j)
  		{  
   		 FOR2(k,l){
			g[i][j] += g_spher[k][l]*jacobian[k][i]*jacobian[l][j]; 
    		    } 			
  		}


       

		return v.sqrt3(v)*t;
	};	


};


int main(void)
{

        const double TIME_MAXIMUM = 10.0, WHOLE_TOLERANCE = 1e-12 ;
        const double T_START = 0.0, DT = 0.001;
	const Vec3 Y_START(1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0, 1.0);
	Black_Hole BH;

        auto eval_solution = [               ](double t          )->double{ return pow(t*t+4,2)/16                   	; } ;
        auto find_error    = [eval_solution  ](double t, double y)->double{ return fabs(y-eval_solution(t))          	; } ;
        auto is_whole      = [WHOLE_TOLERANCE](double t          )->bool  { return fabs(t-round(t)) < WHOLE_TOLERANCE	; } ;

	auto dy = rk4( BH.eval_diff_eqn ) ;
 
        Vec3 y = Y_START;
	double t = T_START;
        while(t <= TIME_MAXIMUM) {
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.x,find_error(t,y.x)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.y,find_error(t,y.y)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.z,find_error(t,y.z)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.t,find_error(t,y.t)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.vx,find_error(t,y.vx)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.vy,find_error(t,y.vy)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.vz,find_error(t,y.vz)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.vt,find_error(t,y.vt)); }
          y = y + dy(t,y,DT) ; t += DT;
        }

  	return 0;
}
