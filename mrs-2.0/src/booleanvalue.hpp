/*
* Copyright (C) 2012 Jennifer Harlow
*
* This file is part of mrs, a C++ class library for statistical set processing.
*
* mrs is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or (at
* your option) any later version.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

/*!/ \file
\brief BooleanValue declarations.
*/

#ifndef __BOOLEANVALUE_HPP__
#define __BOOLEANVALUE_HPP__

#include <algorithm>
#include <iostream>

namespace subpavings {
	
	class BooleanMappedValue {
		public:
		
			/*! \brief No-arguments constructor. */
			BooleanMappedValue();
			
			/*! \brief Constructor.
			 * \param b the boolean value to represent with this. */
			BooleanMappedValue(bool b);
			
			/*! \brief Get the boolean value represented by this. */
			bool getValue() const;
			
			/*! @name Arithmetic operations.*/
			
			//@{
			
			BooleanMappedValue& operator+=(bool b);
			const BooleanMappedValue operator+(bool b) const;
			BooleanMappedValue& operator+=(const BooleanMappedValue& bv);
			const BooleanMappedValue operator+(const BooleanMappedValue& bv) const;
			BooleanMappedValue& operator-=(bool b);
			const BooleanMappedValue operator-(bool b) const;
			BooleanMappedValue& operator-=(const BooleanMappedValue& bv);
			const BooleanMappedValue operator-(const BooleanMappedValue& bv) const;
			BooleanMappedValue& operator*=(bool b);
			const BooleanMappedValue operator*(bool b) const;
			BooleanMappedValue& operator*=(const BooleanMappedValue& bv);
			const BooleanMappedValue operator*(const BooleanMappedValue& bv) const;
			BooleanMappedValue& operator/=(bool b);
			const BooleanMappedValue operator/(bool b) const;
			BooleanMappedValue& operator/=(const BooleanMappedValue& bv);
			const BooleanMappedValue operator/(const BooleanMappedValue& bv) const;
			
			//@}
				
			/*! @name Equality and inequality operations.*/
			
			//@{	
			bool operator==(bool b) const;
			bool operator==(const BooleanMappedValue& bv) const;
			bool operator!=(bool b) const;
			bool operator!=(const BooleanMappedValue& bv) const;
			//@}
			
			/*! \brief Swap contents of this and another %BooleanMappedValue.
			 * \param bv the %BooleanMappedValue to swap contents with.*/
			void swap(BooleanMappedValue& bv);
		
		private :
		
			bool value;
		
		
	};
	
	/*! \brief Output the contents of an BooleanValue object.*/
	std::ostream & operator<<(std::ostream &os, const BooleanMappedValue& bv);


} // end namespace subpavings

// Full specializations of the templates in std namespace can be added in std namespace.
namespace std
{
	template <>
	void swap(subpavings::BooleanMappedValue & v1, 
			subpavings::BooleanMappedValue & v2); // throw ()
	
}


#endif
