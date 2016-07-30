
/*! \file
\brief Vec object class.


*/

#include "Vec.hpp"
#include <math.h>

// to use std input/output
#include <iostream>
#include <stdexcept>

void Vec::normalise()
{

    double length = pow(u*u + v*v + w*w, 0.5);
    if (length > 0.0) {
        u /= length;
        v /= length;
        w /= length;
    }
}

Vec::Vec() : u(0.0), v (0.0), w (0.0) {}

Vec::Vec(double uu, double vv, double ww)
                : u(uu), v (vv), w (ww)
{
    normalise();
}

Vec::Vec(const Vec& other)
                : u(other.u), v (other.v), w (other.w) {}

Vec& Vec::operator=(const Vec& rhs)
{
    u = rhs.u;
    v = rhs.v;
    w = rhs.w;
    normalise();

    return *this;
}

Vec& Vec::operator+=(const Vec& rhs)
{
    u += rhs.u;
    v += rhs.v;
    w += rhs.w;
    normalise();

    return *this;
}

const Vec Vec::operator+(const Vec& rhs) const
{
    return Vec(*this) += rhs;

}

Vec& Vec::operator-=(const Vec& rhs)
{
    u -= rhs.u;
    v -= rhs.v;
    w -= rhs.w;
    normalise();

    return *this;
}

const Vec Vec::operator-(const Vec& rhs) const
{
    return Vec(*this) -= rhs;

}

Vec& Vec::operator*=(const Vec& rhs)
{
    u = v * rhs.w - w * rhs.v;
    v = w * rhs.u - u * rhs.w;
    w = u * rhs.v - v * rhs.u;
    normalise();

    return *this;
}

const Vec Vec::operator*(const Vec& rhs) const
{
    return Vec(*this) *= rhs;

}

Vec& Vec::operator/=(const Vec& rhs)
{
    throw std::runtime_error("vector division attempted");

    return *this;
}

const Vec Vec::operator/(const Vec& rhs) const
{
    throw std::runtime_error("vector division attempted");
	return *this;

}

bool Vec::operator==(const Vec& rhs) const
{
    return !(*this < rhs) && !(*this > rhs);
}

bool Vec::operator!=(const Vec& rhs) const
{
    return (*this < rhs) || !(*this > rhs);
}

bool Vec::operator>(const Vec& rhs) const
{
    return ((u + v + w) < (rhs.u + rhs.v + rhs.w));

}

bool Vec::operator<(const Vec& rhs) const
{
     return ((u + v + w) > (rhs.u + rhs.v + rhs.w));
}

void Vec::print(std::ostream& os) const
{

    os << u << '\t' << v << '\t' << w;

}


//friend function
std::ostream& operator<<(std::ostream& os, const Vec& c)
{

    c.print(os);

    return os;
}
