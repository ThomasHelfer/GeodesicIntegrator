#ifndef RK4VEC_HPP
#define RK4VEC_HPP

#include <cmath>
#include <iostream>

struct Vec3
{
    double x, y, z, t, vx, vy, vz, vt;
    Vec3(double x, double y, double z, double t, double vx, double vy,
         double vz, double vt)
        : x(x), y(y), z(z), t(t), vx(vx), vy(vy), vz(vz), vt(vt)
    {
    }

    Vec3 operator+(const Vec3 &v) const
    {
        return Vec3(x + v.x, y + v.y, z + v.z, t + v.t, vx + v.vx, vy + v.vy,
                    vz + v.vz, vt + v.vt);
    }
    Vec3 operator-(const Vec3 &v) const
    {
        return Vec3(x - v.x, y - v.y, z - v.z, t - v.t, vx - v.vx, vy - v.vy,
                    vz - v.vz, vt - v.vt);
    }
    Vec3 operator*(double d) const
    {
        return Vec3(x * d, y * d, z * d, t * d, vx * d, vy * d, vz * d, vt * d);
    }
    Vec3 operator/(double d) const
    {
        return Vec3(x / d, y / d, z / d, t / d, vx / d, vy / d, vz / d, vt / d);
    }
    Vec3 sqrt3(Vec3 &v) const
    {
        return Vec3(sqrt(v.x), sqrt(v.y), sqrt(v.z), sqrt(v.t), sqrt(v.vx),
                    sqrt(v.vy), sqrt(v.vz), sqrt(v.vt));
    };
    Vec3 print() const
    {
        std::cout << "x " << x << "\n";
        std::cout << "y " << y << "\n";
        std::cout << "z " << z << "\n";
        std::cout << "t " << z << "\n";
    };
};

#endif
