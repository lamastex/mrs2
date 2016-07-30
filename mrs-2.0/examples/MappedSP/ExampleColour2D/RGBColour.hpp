
/*! \file
\brief Declarations for RGBColour class
*/

#ifndef __RGBCOLOUR_HPP__
#define __RGBCOLOUR_HPP__

#include <iosfwd>
//#include <string>

class RGBColour {

    friend std::ostream& operator<<(std::ostream& output, const RGBColour& c);


    private:
        double red;
        double green;
        double blue;

        void normalise();

    public:

        //! no argument constructor
        RGBColour();

        //! parameterised constructor
        RGBColour(double r, double g, double b);

        //! copy constructor
        RGBColour(const RGBColour& other);

        RGBColour& operator=(const RGBColour& rhs);

        RGBColour& operator+=(const RGBColour& rhs);

        const RGBColour operator+(const RGBColour& rhs) const;

        RGBColour& operator-=(const RGBColour& rhs);

        const RGBColour operator-(const RGBColour& rhs) const;

        RGBColour& operator*=(const RGBColour& rhs);

        const RGBColour operator*(const RGBColour& rhs) const;

        RGBColour& operator/=(const RGBColour& rhs);

        const RGBColour operator/(const RGBColour& rhs) const;

		bool operator==(const RGBColour& rhs) const;

        bool operator!=(const RGBColour& rhs) const;

        bool operator>(const RGBColour& rhs) const;

        bool operator<(const RGBColour& rhs) const;

        RGBColour operator|(const RGBColour& rhs) const;

        void print(std::ostream& os) const;
};


#endif
