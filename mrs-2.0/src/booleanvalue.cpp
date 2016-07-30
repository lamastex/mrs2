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
\brief BooleanValue definitions.
*/

/*header for controlling debugging*/
#include "debug_control.hpp"

#include "booleanvalue.hpp"

using namespace subpavings;
using namespace std;

				
// ------------------------ public member functions ----------------
// ------------BooleanMappedValue public member functions ------------
BooleanMappedValue::BooleanMappedValue()
		: value(false) {}

BooleanMappedValue::BooleanMappedValue(bool b)
		: value(b) {}

bool BooleanMappedValue::getValue() const 
{
	return value;
}


// add := union
BooleanMappedValue& 
			BooleanMappedValue::operator+=(bool b)
{
	value = (value || b);
	return *this;
}
				
const BooleanMappedValue BooleanMappedValue::operator+(bool b) const
{
	BooleanMappedValue tmp(*this);
	tmp+=b;
	return tmp;
}

BooleanMappedValue& BooleanMappedValue::operator+=(
								const BooleanMappedValue& bv)
{
	value = (value || bv.value);
	return *this;
}
				
const BooleanMappedValue BooleanMappedValue::operator+(
								const BooleanMappedValue& bv) const
{
	BooleanMappedValue tmp(*this);
	tmp+=bv;
	return tmp;
}
	
// subtract := symmetric set difference, XOR

BooleanMappedValue& BooleanMappedValue::operator-=(bool b)
{
	value = (value && b ? false : (value || b));
	return *this;
}
				
const BooleanMappedValue BooleanMappedValue::operator-(bool b) const
{
	BooleanMappedValue tmp(*this);
	tmp-=b;
	return tmp;
}

BooleanMappedValue& BooleanMappedValue::operator-=(
									const BooleanMappedValue& bv)
{
	value = (value && bv.value ? false : (value || bv.value));
	return *this;
}
				
const BooleanMappedValue BooleanMappedValue::operator-(
								const BooleanMappedValue& bv) const
{
	BooleanMappedValue tmp(*this);
	tmp-=bv;
	return tmp;
}
		
// multiply := intersection
BooleanMappedValue& BooleanMappedValue::operator*=(bool b)
{
	value = (value && b);
	return *this;
}
				
const BooleanMappedValue BooleanMappedValue::operator*(bool b) const
{
	BooleanMappedValue tmp(*this);
	tmp*=b;
	return tmp;
}

BooleanMappedValue& BooleanMappedValue::operator*=(
								const BooleanMappedValue& bv)
{
	value = (value && bv.value);
	return *this;
}
				
const BooleanMappedValue BooleanMappedValue::operator*(
								const BooleanMappedValue& bv) const
{
	BooleanMappedValue tmp(*this);
	tmp*=bv;
	return tmp;
}
	
// divide := set difference, /

BooleanMappedValue& BooleanMappedValue::operator/=(bool b)
{
	value = (value && !b ? true : false);
	return *this;
}
				
const BooleanMappedValue BooleanMappedValue::operator/(bool b) const
{
	BooleanMappedValue tmp(*this);
	tmp/=b;
	return tmp;
}

BooleanMappedValue& BooleanMappedValue::operator/=(
						const BooleanMappedValue& bv)
{
	value = (value && !bv.value ? true : false);
	return *this;
}
				
const BooleanMappedValue BooleanMappedValue::operator/(
						const BooleanMappedValue& bv) const
{
	BooleanMappedValue tmp(*this);
	tmp/=bv;
	return tmp;
}

bool BooleanMappedValue::operator==(bool b) const
{
	return ((value && b) || (!value && !b));
}

bool BooleanMappedValue::operator==(const BooleanMappedValue& bv) const
{
	return ((value && bv.value) || (!value && !bv.value));
}

bool BooleanMappedValue::operator!=(bool b) const
{
	return ((value && !b) || (!value && b));
}

bool BooleanMappedValue::operator!=(const BooleanMappedValue& bv) const
{
	return ((value && !bv.value) || (!value && bv.value));
}

void BooleanMappedValue::swap(BooleanMappedValue& bv)
{
	if (value && !bv.value) {
		value = false;
		bv.value = true;
	}
	if (!value && bv.value) {
		value = true;
		bv.value = false;
	}
}

// Full specializations of the templates in std namespace can be added in std namespace.
template <>
void std::swap(subpavings::BooleanMappedValue & v1, 
			subpavings::BooleanMappedValue & v2) // throw ()
	{
		v1.swap(v2);
	}

std::ostream & subpavings::operator<<(std::ostream &os, const subpavings::BooleanMappedValue& bv)
{
	os << bv.getValue();
	return os;
}
