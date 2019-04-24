

template <typename data_t>
class geodesic_shooter{
        public:
        void shoot(Vec3 center, double shift = 1/5., int shoot = 10 , 
		   const double TIME_MAXIMUM = 150.0, const double T_START = 0.0, const double DT = 0.1){

                data_t metric;
                auto dy = rk4( metric.eval_diff_eqn ) ;

                for(int i = 0; i<shoot ; i++){

                        const Vec3 Y_START(0.0, i*shift, 0.0, 0.0, -1.0, 0.00, 0.0, -1.0);
                        Vec3 y = Y_START + center;
                        double t = T_START;

                        // ========= Preparing output ==========
                        string imgname="xpos";
                        char cbuff[20];
                        sprintf (cbuff, "%03d", i);
                        imgname.append(cbuff);
                        imgname.append(".csv");
                        ofstream myfile(imgname);

                        // ========== Integration and output ===
                        while(t <= TIME_MAXIMUM) {
                                y = y + dy(t,y,DT) ; t += DT;

                                myfile << y.x <<"       " << y.y << "   " << y.z << "   " << y.t << "   " << metric.calculate_norm(y) <<   endl;
                        }

                        // ========== clean up ==================
                        myfile.close();

                }

        }


};
