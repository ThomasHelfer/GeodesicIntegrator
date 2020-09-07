//#include "Rk4vec.hpp"

auto rk4(Vec3 f(double, Vec3))
{
    return [f](double t, Vec3 y, double dt) -> Vec3 {
        return [t, y, dt, f](Vec3 dy1) -> Vec3 {
            return [t, y, dt, f, dy1](Vec3 dy2) -> Vec3 {
                return [t, y, dt, f, dy1, dy2](Vec3 dy3) -> Vec3 {
                    return [t, y, dt, f, dy1, dy2, dy3](Vec3 dy4) -> Vec3 {
                        return (dy1 + dy2 * 2 + dy3 * 2 + dy4) / 6;
                    }(f(t + dt, y + dy3) * dt);
                }(f(t + dt / 2, y + dy2 / 2) * dt);
            }(f(t + dt / 2, y + dy1 / 2) * dt);
        }(f(t, y) * dt);
    };
}
