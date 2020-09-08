#if !defined(GEODESIC_SHOOTER_HPP)
#error "This file should only be included through geodesic_shooter.hpp"
#endif

#ifndef GEODESIC_SHOOTER_IMPL_HPP
#define GEODESIC_SHOOTER_IMPL_HPP

template <typename data_t>
int geodesic_shooter<data_t>::func(double t, const double y[],
                                          double f[], void *params)
{
    double parameter = *(double *)params;

    Vec3 v(y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7]);

    v = data_t::eval_diff_eqn(t,v);

    v.write_to_array(f);

    return GSL_SUCCESS;
}

template <typename data_t>
int geodesic_shooter<data_t>::jac(double t, const double y[],
                                         double *dfdy, double dfdt[],
                                         void *params)
{
    // EXAMPLE CODE FOR JACOBIAN ( ONLY NEEDED FOR SOME ODE METHODS )
    // HERE WE IGNORE IT FOR THE MOMENT
    /*
    (void)(chi); //avoid unused parameter warning
    double mu = *(double *)params;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix *m = &dfdy_mat.matrix;
    gsl_matrix_set(m, 0, 0, 0.0);
    gsl_matrix_set(m, 0, 1, 1.0);
    gsl_matrix_set(m, 1, 0, -2.0 * mu * y[0] * y[1] - 1.0);
    gsl_matrix_set(m, 1, 1, -mu * (y[0] * y[0] - 1.0));
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
    */
}

template <typename data_t>
void geodesic_shooter<data_t>::shoot(Vec3 center, double shift, int shoot,
                                     bool set_geodesic_null,
                                     const double time_end,
                                     const double time_start, const double dt)
{

    gsl_odeiv2_system sys = {func,
                            nullptr, 8};

    const double chi_start = 1e-10;
    const double chi_end = M_PI/2;
    const int NumberOutputs = 0;
    const double eta_init = -0.6;


    for (int i = 0; i < shoot; i++)
    {

        const Vec3 Y_START(0.0, i * shift, 0.0, 0.0, 0.0, 0.00, 0.0, 0.0);
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

        // ========== Integration and output ===
        while (t <= time_end)
        {

            ODE_Solver(sys, y, t, t+dt, NumberOutputs, 1e-13,1e-13);
            t+= dt;

            myfile << y[0] << "       " << y[1] << "   " << y[2] << "   " << y[3]
                   << "   " << /*data_t::calculate_norm(y)*/42 << "\n";
        }

        myfile << std::flush;

        // ========== clean up ==================
        myfile.close();
    }
};

#endif
