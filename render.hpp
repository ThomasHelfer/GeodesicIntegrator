
# define H 50

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
