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

	private:
	static void printArr(double a[][4]) {
		FOR1(i){
		   cout << a[i][0] << " " << a[i][1] << " " << a[i][2] << " " << a[i][3]<< endl;
		}
		cout << "-----------" <<endl;
	}

   	static tensor<2, double> get_metric( double M,  double x, double y, double z){

                tensor<2,double> jacobian   = {};
                tensor<2,double> g = {};
                tensor<2,double> g_spher    = {};

		FOR2(i,j){
			g[i][j] = 0;
			g_spher[i][j] = 0;
			jacobian[i][j] = 0;
		}

                double rr2 =   pow(x,2)+ pow(y,2) + pow(z,2);
                double rho2 =  pow(x,2)+ pow(y,2);

                double rr = sqrt(rr2);
                double rho = sqrt(rho2);
                 // sinus(theta)
                double sintheta = rho/rr;
                double costheta = z/rr;
                // cos(phi)
                double cosphi = x/rho ;
                // sin(phi)
                double sinphi = y/rho ;

                g_spher[0][0] = 1./(1.0 - 2.0*M/rr);
                // Define theta theta component
                g_spher[1][1] = rr2;
                // Define phi phi component
                g_spher[2][2] = rr2*pow(sintheta,2);
                // Devine tt component
                g_spher[3][3] = -(1. - 2.0*M/rr);

                jacobian[0][0] = x/rr ;
                jacobian[1][0] = cosphi*z/rr2 ;
                jacobian[2][0] = -y/rho2 ;
                jacobian[0][1] = y/rr ;
                jacobian[1][1] = sinphi*z/rr2 ;
                jacobian[2][1] = x/rho2 ;
                jacobian[0][2] = z/rr ;
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

                return g;
        }


	public:	
	static Vec3 eval_diff_eqn(double t, Vec3 v){

		double M = 1;		
 	        tensor<2,double> g; // Metix Index low low
                tensor<2,double> g_UU;
                tensor<2,double> g_spher;
                tensor<2,double> g_spher_UU; // Metix Index low low
                tensor<3,double> dg; //
                tensor<3,double> dg_spher; //
                tensor<3,double> chris_LLL; // Christoffel index low low low
                tensor<3,double> chris_ULL; // Christoffel index high low low
                tensor<2,double> jacobian = {};

                tensor<1,double> dx = {v.vx,v.vy,v.vz,v.vt};
                tensor<1,double> ddx;

		FOR1(i){ddx[i] = 0;}
		FOR2(i,j){
			g[i][j] = 0;
			g_UU[i][j] = 0;
			g_spher_UU[i][j] = 0;
			jacobian[i][j] = 0;
		}
		FOR3(i,j,k){
			dg[i][j][k] = 0;
			dg_spher[i][j][k] = 0;
			chris_LLL[i][j][k] = 0;
			chris_ULL[i][j][k] = 0;
		}

		double rr2 =   pow(v.x,2) + pow(v.y,2) + pow(v.z,2);
  		double rho2 =  pow(v.x,2) + pow(v.y,2);

		double rr = sqrt(rr2);
  		double rho = sqrt(rho2);
 		 // sinus(theta)
  		double sintheta = rho/rr;
  		double costheta = v.z/rr;
  		// cos(phi)
  		double cosphi = v.x/rho ; 
 	 	// sin(phi)
 		double sinphi = v.y/rho ; 

		// ==================== 

		dg_spher[0][0][0] = (-2.*M)/pow((-2.*M + rr),2.);
		dg_spher[1][1][0] = 2.*rr;
		dg_spher[2][2][0] = 2.*rr*pow(sintheta,2);
		dg_spher[3][3][0] = -2.*M/rr2;

		dg_spher[2][2][1] = rr2*2.*sintheta*costheta;

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


		g  = get_metric(M,v.x,v.y,v.z);

		
		double unity[4][4] = {};
	
		double det = g[0][0]*(g[1][1]*g[2][2]-g[1][2]*g[2][1])-
                      g[0][1]*(g[2][2]*g[1][0]-g[1][2]*g[2][0])+
                      g[0][2]*(g[1][0]*g[2][1]-g[1][1]*g[2][0]);	

		double det_inverse = 1.0/det;

        	g_UU[0][0] = (g[1][1]*g[2][2] - g[1][2]*g[1][2]) * det_inverse;
        	g_UU[0][1] = (g[0][2]*g[1][2] - g[0][1]*g[2][2]) * det_inverse;
        	g_UU[0][2] = (g[0][1]*g[1][2] - g[0][2]*g[1][1]) * det_inverse;
        	g_UU[1][1] = (g[0][0]*g[2][2] - g[0][2]*g[0][2]) * det_inverse;
        	g_UU[1][2] = (g[0][1]*g[0][2] - g[0][0]*g[1][2]) * det_inverse;
        	g_UU[2][2] = (g[0][0]*g[1][1] - g[0][1]*g[0][1]) * det_inverse;
        	g_UU[1][0] = g_UU[0][1];
       		g_UU[2][0] = g_UU[0][2];
        	g_UU[2][1] = g_UU[1][2];

		FOR2(i,j)
		{
			FOR1(k)
			{
				unity[i][j] += g[i][k]*g_UU[k][j];
			}
		}


 		FOR3(i,j,m)
  		{  
   		 FOR3(k,l,n){
			dg[i][j][m] += dg_spher[k][l][n]*jacobian[k][i]*jacobian[l][j]*jacobian[n][m]; 
    		    } 			
  		}

		//=========================

      		FOR3(i,j,k)
        	{
            		chris_LLL[i][j][k] = 0.5*(dg[j][i][k] + dg[k][i][j] - dg[j][k][i]);
        	}
        	FOR3(i,j,k)
        	{
            		chris_ULL[i][j][k] = 0;
            		FOR1(l)
            		{
                 		chris_ULL[i][j][k] += g_UU[i][l]*chris_LLL[l][j][k];
            		}
        	}
		
		//=========================

		FOR1(i){
			FOR2(k,l){
				ddx[i] += - chris_ULL[i][k][l]*dx[k]*dx[l];	
			}
		}
		
		Vec3 out(dx[0],dx[1],dx[2],dx[3],ddx[0],ddx[1],ddx[2],ddx[3]);		

		return out;
		//return v.sqrt3(v)*t;
	};	


};


int main(void)
{

        const double TIME_MAXIMUM = 1000.0, WHOLE_TOLERANCE = 1e-12 ;
        const double T_START = 0.0, DT = 0.01;
	const Vec3 Y_START(10.0, 0.0, 0.0, 0.0, 0.0, 0.34, 0.0, 1.0);
	Black_Hole BH;

        auto eval_solution = [               ](double t          )->double{ return pow(t*t+4,2)/16                   	; } ;
        auto find_error    = [eval_solution  ](double t, double y)->double{ return fabs(y-eval_solution(t))          	; } ;
        auto is_whole      = [WHOLE_TOLERANCE](double t          )->bool  { return fabs(t-round(t)) < WHOLE_TOLERANCE	; } ;

	auto dy = rk4( BH.eval_diff_eqn ) ;
 
        Vec3 y = Y_START;
	double t = T_START;

	ofstream myfile("xpos.csv");

        while(t <= TIME_MAXIMUM) {
/*          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.x,find_error(t,y.x)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.y,find_error(t,y.y)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.z,find_error(t,y.z)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.t,find_error(t,y.t)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.vx,find_error(t,y.vx)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.vy,find_error(t,y.vy)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.vz,find_error(t,y.vz)); }
          if (is_whole(t)) { printf("y(%4.1f)\t=%12.6f \t error: %12.6e\n",t,y.vt,find_error(t,y.vt)); }
*/
  	  y = y + dy(t,y,DT) ; t += DT;			
	  myfile << y.x <<"	" << y.y << "	" << y.z << "	" << y.t << "" <<  endl;
	}

	myfile.close();

  	return 0;
}
