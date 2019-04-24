
# define H 10

template <typename data_t>
class render_black_hole{
        public:
        int* picture(int* picture, Vec3 center, double  max_x, double  max_y,
		      const double alpha = 0,
                      const int start_ind = 0, int end_ind = (H*H), 
    		      const double TIME_MAXIMUM = 150.0, const double DT = 0.1,
                      const double T_START = 0){

                data_t metric;
                std::clock_t start;
                double duration;


                auto dy = rk4( metric.eval_diff_eqn ) ;
		if(end_ind > H*H) { end_ind = H*H ;} 

		for (int i = start_ind; i < end_ind; i++){
                                bool draw = true;
				int x1 = i % H;
				int y1 = (i - x1)/H;
 				int i_local = i - start_ind;
				double x_shot = (double)x1/H*max_x-max_x/2.;
                                double y_shot = (double)y1/H*max_y-max_y/2.;


                                Vec3 Y_START(-y_shot*sin(alpha), x_shot, y_shot*cos(alpha), 0.0, 0.0, 0.00, 0.0, 0.0);
				Y_START = Y_START + center ; 	
                                Vec3 y = Y_START;
                		y = metric.set_norm(y,0);
       				double t = T_START;

                                while((t <= TIME_MAXIMUM) && draw) {
                                        y = y + dy(t,y,DT) ; t += DT;

                                        double rr = abs(sqrt((y.x)*(y.x)+(y.y)*(y.y))-6.5);

                                        if( rr < 1 && abs(y.z) < 0.2){
                                                picture[i_local] = (int)(rr*255);
                                                draw = false;
                                        }
                                }
                                if(draw){
                                        picture[i_local] = 0;
                                }


                }

/*
                        start = std::clock();
   			duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
                        cout << " Completed = " <<  (double)y1/(double)H*100 << " % " << endl;
                        cout << " duration " << duration << " sec" <<endl;
                        cout << " Estmated time to finish " << (double)duration*((double)H-(double)y1)/60.0 << " min " << endl;
*/
			return picture;
        }

       void render(int *picture, string file_name){
          std::ofstream out(file_name);
          out << "P3\n" << H << ' ' << H << ' ' << "255\n";

          for(int i = 0; i < H*H; i++){
		  	   if(picture[i] != 0){
			   int red =  picture[i];
			   int blue = abs(255-picture[i]);
			   //cout << red << endl;
			   out << red << ' '
                            << 0 << ' '
                           << blue << '\n';
			   } else {

			   out << 0 << ' '
                               << 0 << ' '
                               << 0 << '\n';
			   }
          }
        }
};
