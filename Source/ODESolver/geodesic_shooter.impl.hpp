#if !defined(GEODESIC_SHOOTER_HPP)
#error "This file should only be included through geodesic_shooter.hpp"
#endif

#ifndef GEODESIC_SHOOTER_IMPL_HPP
#define GEODESIC_SHOOTER_IMPL_HPP

template <typename data_t>
void geodesic_shooter<data_t>::shoot(Vec3 center, double shift, int shoot,
                                     bool set_geodesic_null,
                                     const double time_end,
                                     const double time_start, const double dt)
{

    data_t metric;
    auto dy = rk4(metric.eval_diff_eqn);

    for (int i = 0; i < shoot; i++)
    {

        const Vec3 Y_START(0.0, i * shift, 0.0, 0.0, 0.0, 0.00, 0.0, 0.0);
        Vec3 y = Y_START + center;

        if (set_geodesic_null)
        {
            y = metric.set_norm(y, 0);
        }

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
            y = y + dy(t, y, dt);
            t += dt;

            myfile << y.x << "       " << y.y << "   " << y.z << "   " << y.t
                   << "   " << metric.calculate_norm(y) << std::endl;
        }

        // ========== clean up ==================
        myfile.close();
    }
};

#endif
