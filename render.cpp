
#include "render.hpp"

template <typename data_t>
void render_black_hole::picture(int *red, int *green, int *blue, Vec3 center,
                                double max_x, double max_y,
                                const double alpha = 0, const int start_ind = 0,
                                int end_ind = (H * H),
                                const double TIME_MAXIMUM = 150.0,
                                const double DT = 0.1,
                                const double T_START = 0) {

  data_t metric;
  std::clock_t start;
  double duration;

  auto dy = rk4(metric.eval_diff_eqn);
  if (end_ind > H * H) {
    end_ind = H * H;
  }

  for (int i = start_ind; i < end_ind; i++) {
    bool draw = true;
    int x1 = i % H;
    int y1 = (i - x1) / H;
    int i_local = i - start_ind;
    double x_shot = (double)x1 / H * max_x - max_x / 2.;
    double y_shot = (double)y1 / H * max_y - max_y / 2.;

    Vec3 Y_START(-y_shot * sin(alpha), x_shot, y_shot * cos(alpha), 0.0, 0.0,
                 0.00, 0.0, 0.0);
    Y_START = Y_START + center;
    Vec3 y = Y_START;
    y = metric.set_norm(y, 0);
    double t = T_START;

    while ((t <= TIME_MAXIMUM) && draw) {
      y = y + dy(t, y, DT);
      t += DT;

      double rr = abs(sqrt((y.x) * (y.x) + (y.y) * (y.y)) - 6.5);

      if (rr < 1 && abs(y.z) < 0.2) {
        double blue_tmp = 1 - t / y.t;
        double red_tmp = (rr);
        red[i_local] = (int)(255 * red_tmp);
        green[i_local] = 0;
        blue[i_local] = (int)(255 * blue_tmp);
        draw = false;
      }
    }
    if (draw) {
      red[i_local] = 0;
      green[i_local] = 0;
      blue[i_local] = 0;
    }
  }
}

void render_black_hole::render_circle(int *red, int *green, int *blue,
                                      double max_x, double max_y,
                                      const int start_ind = 0,
                                      int end_ind = (H * H), double Radius = 2,
                                      double thickness = 0.1) {

  data_t metric;
  double duration;

  if (end_ind > H * H) {
    end_ind = H * H;
  }

  for (int i = start_ind; i < end_ind; i++) {
    int x1 = i % H;
    int y1 = (i - x1) / H;
    int i_local = i - start_ind;
    double x_shot = (double)x1 / H * max_x - max_x / 2.;
    double y_shot = (double)y1 / H * max_y - max_y / 2.;
    double rr = sqrt(x_shot * x_shot + y_shot * y_shot);

    if ((rr > Radius - thickness) && (rr < Radius + thickness)) {
      red[i_local] = 0;
      green[i_local] = 255;
      blue[i_local] = 0;
    }
  }
}

void render_black_hole::render(int *red, int *green, int *blue,
                               std::string file_name) {
  std::ofstream out(file_name);
  out << "P3\n" << H << ' ' << H << ' ' << "255\n";

  for (int i = 0; i < H * H; i++) {
    out << red[i] << ' ' << green[i] << ' ' << blue[i] << '\n';
  }
}
