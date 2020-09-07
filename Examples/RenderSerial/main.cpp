#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>


#include "geodesic_shooter.hpp"
#include "render.hpp"
#include "rk4.hpp"
#include "schwarzschild.hpp"
#include "tensor.hpp"
#include "DimensionDefinitions.hpp"

using namespace std;


int main(void)
{

    // ==========================================
    // ========== Render a picture ==============
    // ==========================================
    const int resolution = 10;
    const double alpha = 0;
    const double start_ind = 0;
    const double end_ind = resolution*resolution;
    string name_render = "out.ppm";
    render_black_hole<Black_Hole> rend;
    const double size_x = 20;
    const double size_y = 20;
    int red[resolution*resolution];
    int green[resolution*resolution];
    int blue[resolution*resolution];

    const Vec3 center(15.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0);

    rend.picture(red, green, blue, center, size_x, size_y, resolution, alpha, start_ind, end_ind);

    rend.render(red, green, blue,resolution, name_render);

    return 0;
}
