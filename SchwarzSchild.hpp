#include "tensor.hpp"
#include "Rk4vec.hpp"

using namespace std;


#define FOR1(IDX) for (int IDX = 0; IDX < 4; ++IDX)
#define FOR2(IDX1,IDX2) FOR1(IDX1) FOR1(IDX2)
#define FOR3(IDX1,IDX2,IDX3) FOR2(IDX1,IDX2) FOR1(IDX3)
#define FOR4(IDX1,IDX2,IDX3,IDX4) FOR2(IDX1,IDX2) FOR2(IDX3,IDX4)

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

	static tensor<3, double> get_metric_deriv(double M, double x,double y, double z){
		
		tensor<3,double> dg;
		tensor<2,double> g ;
		tensor<2,double> g_dx ;
		tensor<2,double> g_dy ;
		tensor<2,double> g_dz ;
		double h = 0.00001;

		g = get_metric(M,x,y,z);
		g_dx = get_metric(M,x-h,y,z);
		g_dy = get_metric(M,x,y-h,z);
		g_dz = get_metric(M,x,y,z-h);

		FOR2(i,j){
			dg[i][j][0] = (g[i][j]-g_dx[i][j])/h;
			dg[i][j][1] = (g[i][j]-g_dy[i][j])/h;
			dg[i][j][2] = (g[i][j]-g_dz[i][j])/h;
			dg[i][j][3] = 0;
		}
		return dg;
		
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
		double eps = 0.001;
		if(rr < 2*M+eps){
			FOR1(i){ddx[i]=0;dx[i]=0;}
			Vec3 out(dx[0],dx[1],dx[2],dx[3],ddx[0],ddx[1],ddx[2],ddx[3]);		
			return out;
		}
		// ==================== 

		g  = get_metric(M,v.x,v.y,v.z);

	
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


		dg = get_metric_deriv(M,v.x,v.y,v.z);
				
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
	}
};
