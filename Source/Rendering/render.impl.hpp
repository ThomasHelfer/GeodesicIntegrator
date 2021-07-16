#if !defined(RENDER_HPP)
#error "This file should only be included through render.hpp"
#endif

#ifndef RENDER_IMPL_HPP
#define RENDER_IMPL_HPP

template <typename data_t>
void render_black_hole<data_t>::picture(
    int *red, int *green, int *blue, Vec3 center, double max_x, double max_y,
    const int resolution, const double alpha, const int start_ind, int end_ind,
    const double time_end, const double dt, const double time_start)
{

    data_t metric;
    std::clock_t start;
    double duration;

    gsl_odeiv2_system sys = {metric.eval_diff_eqn, nullptr, 8};

    const double epsabs = 1e-6;
    const double epsrel = 1e-6;
    const double hstart = 1e-6;
    const int nmax = 1000;

    for (int i = start_ind; i < end_ind; i++)
    {
        bool draw = true;
        int x1 = i % resolution;
        int y1 = (i - x1) / resolution;
        int i_local = i - start_ind;
        double x_shot = (double)x1 / resolution * max_x - max_x / 2.;
        double y_shot = (double)y1 / resolution * max_y - max_y / 2.;

        Vec3 Y_START(-y_shot * sin(alpha), x_shot, y_shot * cos(alpha), 0.0,
                     0.0, 0.00, 0.0, 0.0);
        Y_START = Y_START + center;

        // Making it a null - geodesic
        Vec3 y_temp = Y_START;
        y_temp = TensorAlgebra::set_norm<data_t>(y_temp, 0);

        double t = time_start;
        int status = 0;

        double y[8];
        y_temp.write_to_array(y);

        while ((t <= time_end) && draw && status == 0)
        {
            status =
                ODE_Solver(sys, y, t, t + dt, hstart, epsabs, epsrel, nmax);
            t += dt;

            double rr = abs(sqrt((y[0]) * (y[0]) + (y[1]) * (y[1])) - 6.5);

            if (rr < 1 && abs(y[2]) < 0.2)
            {
                double blue_tmp = 1 - t / y[3];
                double red_tmp = (rr);
                red[i_local] = (int)(255 * red_tmp);
                green[i_local] = 0;
                blue[i_local] = (int)(255 * blue_tmp);
                draw = false;
            }
        }
        if (draw)
        {
            red[i_local] = 0;
            green[i_local] = 0;
            blue[i_local] = 0;
        }
    }
}

template <typename data_t>
void render_black_hole<data_t>::render_circle(int *red, int *green, int *blue,
                                              double max_x, double max_y,
                                              const int resolution,
                                              const int start_ind, int end_ind,
                                              double Radius, double thickness)
{

    data_t metric;
    double duration;

    if (end_ind > resolution * resolution)
    {
        end_ind = resolution * resolution;
    }

    for (int i = start_ind; i < end_ind; i++)
    {
        int x1 = i % resolution;
        int y1 = (i - x1) / resolution;
        int i_local = i - start_ind;
        double x_shot = (double)x1 / resolution * max_x - max_x / 2.;
        double y_shot = (double)y1 / resolution * max_y - max_y / 2.;
        double rr = sqrt(x_shot * x_shot + y_shot * y_shot);

        if ((rr > Radius - thickness) && (rr < Radius + thickness))
        {
            red[i_local] = 0;
            green[i_local] = 255;
            blue[i_local] = 0;
        }
    }
}

template <typename data_t>
void render_black_hole<data_t>::render(int *red, int *green, int *blue,
                                       const int resolution,
                                       std::string file_name)
{
    std::ofstream out(file_name);
    out << "P3\n" << resolution << ' ' << resolution << ' ' << "255\n";

    for (int i = 0; i < resolution * resolution; i++)
    {
        out << red[i] << ' ' << green[i] << ' ' << blue[i] << '\n';
    }
}

#endif
