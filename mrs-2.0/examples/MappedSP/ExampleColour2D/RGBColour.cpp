
/*! \file
\brief RGBColour object class.


*/

#include "RGBColour.hpp"
//#include "cxsc.hpp"

// to use std input/output
#include <iostream>


void RGBColour::normalise()
{

    if (red > 1.0 || green > 1.0 || blue > 1.0) {
        double total = red + green + blue;
        red /= total;
        green /= total;
        blue /= total;
    }
}

RGBColour::RGBColour() : red(0.0), green (0.0), blue (0.0) {}

RGBColour::RGBColour(double r, double g, double b)
                : red(r), green (g), blue (b)
{
    normalise();
}

RGBColour::RGBColour(const RGBColour& other)
                : red(other.red), green (other.green), blue (other.blue) {}

RGBColour& RGBColour::operator=(const RGBColour& rhs)
{
    red = rhs.red;
    green = rhs.green;
    blue = rhs.blue;
    normalise();

    return *this;
}

RGBColour& RGBColour::operator+=(const RGBColour& rhs)
{
    red += rhs.red;
    green += rhs.green;
    blue += rhs.blue;
    normalise();

    return *this;
}

const RGBColour RGBColour::operator+(const RGBColour& rhs) const
{
    return RGBColour(*this) += rhs;

}

RGBColour& RGBColour::operator-=(const RGBColour& rhs)
{
    red > rhs.red ? red -= rhs.red : 0.0;
    green > rhs.green ? green -= rhs.green : 0.0;
    blue > rhs.blue ? blue -= rhs.blue : 0.0;
    normalise();

    return *this;
}

const RGBColour RGBColour::operator-(const RGBColour& rhs) const
{
    return RGBColour(*this) -= rhs;

}

RGBColour& RGBColour::operator*=(const RGBColour& rhs)
{
    red = red * rhs.red;
    green = green * rhs.green;
    blue = blue * rhs.blue;
    normalise();

    return *this;
}

const RGBColour RGBColour::operator*(const RGBColour& rhs) const
{
    return RGBColour(*this) *= rhs;

}

RGBColour& RGBColour::operator/=(const RGBColour& rhs)
{
    red = (rhs.red > 0.0 ? red / rhs.red : 1.0);
    green = (rhs.green > 0.0 ? red / rhs.green : 1.0);
    blue = (rhs.blue > 0.0 ? red / rhs.blue : 1.0);
    normalise();

    return *this;
}

const RGBColour RGBColour::operator/(const RGBColour& rhs) const
{
    return RGBColour(*this) /= rhs;

}

bool RGBColour::operator==(const RGBColour& rhs) const
{
    return !(*this < rhs) && !(*this > rhs);
}

bool RGBColour::operator!=(const RGBColour& rhs) const
{
    return (*this < rhs) || !(*this > rhs);
}

bool RGBColour::operator>(const RGBColour& rhs) const
{
    return ((red + green + blue) < (rhs.red + rhs.green + rhs.blue));

}

bool RGBColour::operator<(const RGBColour& rhs) const
{
     return ((red + green + blue) > (rhs.red + rhs.green + rhs.blue));
}

RGBColour RGBColour::operator|(const RGBColour& rhs) const
{
    RGBColour newColour(*this);
    if (rhs.red > newColour.red) newColour.red = rhs.red;
    if (rhs.green > newColour.green) newColour.green = rhs.green;
    if (rhs.blue > newColour.blue) newColour.blue = rhs.blue;
    newColour.normalise();

    return newColour;

}

void RGBColour::print(std::ostream& os) const
{

    os << red << '\t' << green << '\t' << blue;

}


//friend function
std::ostream& operator<<(std::ostream& os, const RGBColour& c)
{

    c.print(os);

    return os;
}
