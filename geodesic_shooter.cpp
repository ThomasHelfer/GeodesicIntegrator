
#include "geodesic_shooter.hpp"

template <typename data_t>
void geodesic_shooter<data_t>::shoot(Vec3 center, double shift, int shoot,
                                     bool set_null, const double TIME_MAXIMUM,
                                     const double T_START, const double DT)
{

    data_t metric;
    auto dy = rk4(metric.eval_diff_eqn);

    for (int i = 0; i < shoot; i++)
    {

        const Vec3 Y_START(0.0, i * shift, 0.0, 0.0, 0.0, 0.00, 0.0, 0.0);
        Vec3 y = Y_START + center;

        if (set_null)
        {
            y = metric.set_norm(y, 0);
        }

        double t = T_START;

        // ========= Preparing output ==========
        std::string imgname = "xpos";
        char cbuff[20];
        std::sprintf(cbuff, "%03d", i);
        imgname.append(cbuff);
        imgname.append(".csv");
        std::ofstream myfile(imgname);

        // ========== Integration and output ===
        while (t <= TIME_MAXIMUM)
        {
            y = y + dy(t, y, DT);
            t += DT;

            myfile << y.x << "       " << y.y << "   " << y.z << "   " << y.t
                   << "   " << metric.calculate_norm(y) << std::endl;
        }

        // ========== clean up ==================
        myfile.close();
    }
};
