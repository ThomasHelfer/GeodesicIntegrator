#if !defined(GEODESIC_SHOOTER_HPP)
#error "This file should only be included through geodesic_shooter.hpp"
#endif

#ifndef GEODESIC_SHOOTER_IMPL_HPP
#define GEODESIC_SHOOTER_IMPL_HPP


template <typename data_t>
void geodesic_shooter<data_t>::shoot(Vec3 center, double shift, int numberofgeodesics,
                                     bool set_geodesic_null,
                                     const double time_end,
                                     const double time_start, const double dt)
{
    data_t metric;
    gsl_odeiv2_system sys = {metric.eval_diff_eqn,
                            nullptr, 8};

    const int NumberOutputs = 1;
    const double epsabs = 1e-6;
    const double epsrel = 1e-6;
    const double hstart = 1e-6;
    const int nmax = 1000;

    for (int i = 0; i < numberofgeodesics; i++)
    {
        // === Setting up initial data ========
        const Vec3 Y_START(0.0, i * shift, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        Vec3 y_temp = Y_START + center;
        if (set_geodesic_null)
        {
            y_temp = data_t::set_norm(y_temp, 0);
        }

        double y[8];
        y_temp.write_to_array(y);

        double t = time_start;

        // ========= Preparing output ==========
        std::string imgname = "xpos";
        char cbuff[20];
        std::sprintf(cbuff, "%03d", i);
        imgname.append(cbuff);
        imgname.append(".csv");
        std::ofstream myfile(imgname);
        int status = 0;

        // ========== Integration and output ===
        while (t <= time_end && status == 0  )
        {
            status = ODE_Solver(sys, y, t, t+dt, NumberOutputs,hstart,  epsabs, epsrel,nmax);
            myfile << y[0] << "       " << y[1] << "   " << y[2] << "   " << y[3]
                   << "   " << t << "\n";

            t+= dt;

        }

        myfile << std::flush;

        // ========== clean up ==================
        myfile.close();
    }
};

#endif
