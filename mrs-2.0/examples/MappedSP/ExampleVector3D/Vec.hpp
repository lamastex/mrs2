
/*! \file
\brief Declarations for Vec class
*/

#ifndef __VEC_HPP__
#define __VEC_HPP__

#include <iosfwd>
//#include <string>

class Vec {

    friend std::ostream& operator<<(std::ostream& output, const Vec& c);


    private:
        double u;
        double v;
        double w;

        void normalise();

    public:

        //! no argument constructor
        Vec();

        //! parameterised constructor
        Vec(double uu, double vv, double ww);

        //! copy constructor
        Vec(const Vec& other);

        Vec& operator=(const Vec& rhs);

        Vec& operator+=(const Vec& rhs);

        const Vec operator+(const Vec& rhs) const;

        Vec& operator-=(const Vec& rhs);

        const Vec operator-(const Vec& rhs) const;

        Vec& operator*=(const Vec& rhs);

        const Vec operator*(const Vec& rhs) const;
		
		Vec& operator/=(const Vec& rhs);

        const Vec operator/(const Vec& rhs) const;

        bool operator==(const Vec& rhs) const;

        bool operator!=(const Vec& rhs) const;

        bool operator>(const Vec& rhs) const;

        bool operator<(const Vec& rhs) const;

        void print(std::ostream& os) const;
};


#endif
