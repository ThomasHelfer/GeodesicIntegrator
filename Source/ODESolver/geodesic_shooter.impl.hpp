#if !defined(GEODESIC_SHOOTER_HPP)
#error "This file should only be included through geodesic_shooter.hpp"
#endif

#ifndef GEODESIC_SHOOTER_IMPL_HPP
#define GEODESIC_SHOOTER_IMPL_HPP

template <typename data_t>
void geodesic_shooter<data_t>::shoot(
    Vec3 center, double shift, int numberofgeodesics, bool set_geodesic_null,
    const double time_end, const double time_start, const double dt,
    const double epsabs, const double epsrel,
    const double hstart, const int nmax)
{
    data_t metric;
    gsl_odeiv2_system sys = {metric.eval_diff_eqn, nullptr, 8};

    for (int i = 0; i < numberofgeodesics; i++)
    {
        // === Setting up initial data ========
        const Vec3 Y_START(0.0, i * shift, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        Vec3 y_temp = Y_START + center;
        if (set_geodesic_null)
        {
            y_temp = TensorAlgebra::set_norm<data_t>(y_temp, 0);
        }

        double y[8];
        y_temp.write_to_array(y);

        single_shot(y, i, time_end, time_start, dt, epsabs, epsrel, hstart,
                    nmax);
    }
};

template <typename data_t>
void geodesic_shooter<data_t>::single_shot(double y[], const int index,
                                           const double time_end,
                                           const double time_start,
                                           const double dt, const double epsabs,
                                           const double epsrel,
                                           const double hstart, const int nmax)
{
    data_t metric;
    gsl_odeiv2_system sys = {metric.eval_diff_eqn, nullptr, 8};

    double t = time_start;

    // ========= Preparing output ==========
    std::string imgname = "xpos";
    char cbuff[20];
    std::sprintf(cbuff, "%03d", index);
    imgname.append(cbuff);
    imgname.append(".csv");
    std::ofstream myfile(imgname);
    int status = 0;

    // ========== Integration and output ===
    while (t <= time_end && status == 0)
    {
        tensor<2, double> g = metric.get_metric(1, y[0], y[1], y[2], y[3]);
        const auto norm =
            TensorAlgebra::calculate_norm(y[4], y[5], y[6], y[7], g);
        myfile << y[0] << "       " << y[1] << "   " << y[2] << "   " << y[3]
               << "   " << y[4] << "       " << y[5] << "   " << y[6] << "   "
               << y[7] << "   " << norm << "\n";

        status = ODE_Solver(sys, y, t, t + dt, hstart, epsabs, epsrel, nmax);
        t += dt;
    }

    myfile << std::flush;

    // ========== clean up ==================
    myfile.close();
};
#endif
