/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009, 2010, 2011 Jennifer Harlow
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
\brief RangeCollectionHist definitions
*/

#include "rangecollection_hist.hpp"

#include "subpaving_exception.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

#include <numeric>
#include <algorithm>
#include <functional>
#include <iterator>


using namespace subpavings;
using namespace std;

//public    

RangeCollectionHist::RangeCollectionHist() {}


RangeCollectionHist::RangeCollectionHist(const RealVec& range)
{
	// easy to do when containers are just RVecDatas
	container.assign(range.begin(), range.end());
	
	//for (RVecData::const_iterator it = range.begin(); it < range.end(), ++it) {
	//	container.push_back( RangeCollectionHistValueType(r) );
	//}
}

RangeCollectionHist::RangeCollectionHist(const cxsc::real& r)
{
	container.push_back( RangeCollectionHistPrivateValueType(r) );
}

RangeCollectionHist::RangeCollectionHist(
			const RangeCollectionHist& other)
				: container(other.container) {}


RangeCollectionHist& RangeCollectionHist::operator=(
							RangeCollectionHist rhs)
{
	rhs.swap(*this);
	return *this;
}

// return iterator to start of collection
RangeCollectionHist::RangeCollectionHistItr RangeCollectionHist::begin() const
{
	return container.begin();
}

// return iterator to end of collection
RangeCollectionHist::RangeCollectionHistItr RangeCollectionHist::end() const
{
	return container.end();
}


std::size_t RangeCollectionHist::getSize() const
{
	return container.size();
}

bool RangeCollectionHist::isEmpty() const
{
	return container.empty();
}

RealVec RangeCollectionHist::getRangeCollectionValues() const
{
	//because our value type is in fact a real, we can just return a copy of the container
	return container;
}

cxsc::real RangeCollectionHist::getTotal() const
{
	real init(0.0);
	
	return std::accumulate(container.begin(), container.end(), init);
}

cxsc::real RangeCollectionHist::getAbsTotal() const
{
	real result(0.0);
	
	for (RangeCollectionHistPrivateTypeConstItr cit = container.begin();
			cit < container.end(); ++cit) {
			
			result += cxsc::abs(*cit);
	}
	return result;
}

cxsc::real RangeCollectionHist::getAverageValue() const
{
	if ( isEmpty() ) {
		throw UnfulfillableRequest_Error("RangeCollectionHist::getAverageValue()");
	}
	return getTotal()/(1.0*getSize());
}


RangeCollectionHist RangeCollectionHist::makeAverageRangeCollection() const
{
	return RangeCollectionHist(getAverageValue());
}

RealVec RangeCollectionHist::getAbsDiffsToAverageRangeCollection() const
{
	RealVec result;
	result.reserve( getSize() );
	
	cxsc::real av = getAverageValue();
	for (RangeCollectionHistPrivateTypeConstItr cit = container.begin();
			cit < container.end(); ++cit) {
		result.push_back( cxsc::abs( (*cit) - av ) );
	}
	return result;
}


// parallel adding of two collections of the same size
//ith element of result is addition of ith element of lhs to ith element of this
void RangeCollectionHist::parallelAdd(const RangeCollectionHist& rhs)
{
	if (getSize() != rhs.getSize() ) {
		throw std::logic_error(
			"RangeCollectionHist::parallelAdd(const RangeCollectionHist&) : collections of different sizes");
	}
	
	std::transform(container.begin(), container.end(), 
			rhs.begin(), container.begin(),
			std::plus<RangeCollectionHistPrivateValueType>() );
}


void RangeCollectionHist::reduce()
{
	if (getSize() > 1) {

		RangeCollectionHistPrivateValueType sum(0.0);
		sum = std::accumulate(container.begin(),
										container.end(),
										sum);
										
		RangeCollectionHistPrivateType temp(1,sum);
		temp.swap(container);
	}
}

RangeCollectionHist& RangeCollectionHist::operator+=(
							const RangeCollectionHist& rhs)
{
	container.insert(container.end(),
									rhs.container.begin(),
									rhs.container.end());
	return *this;
}
			
const RangeCollectionHist RangeCollectionHist::operator+(
							const RangeCollectionHist& rhs) const
{
	RangeCollectionHist result = *this;     
	result += rhs;            
	return result; 
}
			
RangeCollectionHist& RangeCollectionHist::operator*=(const cxsc::real& mult)
{
	cxsc::real realMult = cxsc::real(mult);
	
	this->scalarMult(realMult);
	
	//transform (container.begin(), container.end(), container.begin(),
	//			bind1st(multiplies<double>(),mult) );

	/*
	RangeCollectionHistPrivateTypeItr rit = container.begin();
	for (RangeCollectionHistPrivateTypeItr rit = container.begin();
					rit < container.end(); rit++) {
				*rit = mult * (*rit);
	}
	*/
	return *this;
}
			
const RangeCollectionHist RangeCollectionHist::operator*(const cxsc::real& mult) const
{
	RangeCollectionHist result = *this;     
	result *= mult;            
	return result; 
}

RangeCollectionHist& RangeCollectionHist::operator/=(const cxsc::real& div)
{
	if (div == 0.0) {
		throw std::invalid_argument("RangeCollectionHist::operator/= : divide by zero");
	}
	
	cxsc::real realMult = 1.0/cxsc::real(div);
	
	this->scalarMult(realMult);
		
	return *this;
}
			
const RangeCollectionHist RangeCollectionHist::operator/(const cxsc::real& div) const
{
	RangeCollectionHist result = *this;     
	result /= div;            
	return result; 
}




bool RangeCollectionHist::checkRangeCollectionAllPositive() const
{
	RangeCollectionHistPrivateTypeConstItr it;
	it = std::find_if(container.begin(), container.end(), 
			std::bind2nd(std::less<RangeCollectionHistPrivateValueType>(),0) );
	return ( !(it < container.end()) );
}

/*
std::string RangeCollectionHist::toString(int prec) const
{
	try {
		std::string result = "";
		
		if (!isEmpty()) {
			std::ostringstream oss;
			
			oss << std::setprecision(prec);

			std::string delim = "\t";
			
			std::ostream_iterator<RangeCollectionHistValueType> out_it (oss, delim.c_str());
			copy ( container.begin(), container.end(), out_it );
			
			result = oss.str();
			std::size_t found = result.find_last_of(delim);
			if (found < std::string::npos && found > 1) {
				result = result.substr(0, found-1);
			}
		}
		return result;
	}
	catch (exception& e) {
		string msg( e.what() );
		throw RangeCollectionException("Error in toString :\n"
				+ msg);
	}
}	
*/		


std::ostream& RangeCollectionHist::outputTabs(std::ostream &os,
										int prec) const
{
	
	os << cxsc::SaveOpt;
	os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
	
	outputTabs(os);
	
	os << cxsc::RestoreOpt;
	
	return os;
}

std::ostream& RangeCollectionHist::outputTabs(std::ostream &os) const
{
	
	std::string delim = "\t";
	
	if (getSize() > 1 ) {
		RangeCollectionHistPrivateTypeConstItr eit = container.end();
		eit--;

		std::ostream_iterator<RangeCollectionHistPrivateValueType> out_it (os, delim.c_str());
		copy ( container.begin(), eit, out_it );
	}
	if (!isEmpty() ) {
		os << container.back();
	}
	return os;
}

std::ostream& RangeCollectionHist::outputAverage(std::ostream &os,
										int prec) const
{
	os << cxsc::SaveOpt;
	os << cxsc::Variable << cxsc::SetPrecision(prec+2,prec);
	
	outputAverage(os);
	
	os << cxsc::RestoreOpt;
	
	return os;
}

std::ostream& RangeCollectionHist::outputAverage(std::ostream &os) const
{
	os << getAverageValue();
	return os;
}

// non-throwing swap idiom
void RangeCollectionHist::swap(RangeCollectionHist& other)
{
   std::swap(container, other.container);
}


//private


// multiplies each element of the range collection by mult
void RangeCollectionHist::scalarMult(const cxsc::real& mult)
{
	RangeCollectionHistPrivateType tmp_container;
	
	// don't assume the value types are real
	// could do faster with boost ...
	
	for (RangeCollectionHistPrivateTypeConstItr cit = container.begin(); 
					cit < container.end(); ++cit) {
		tmp_container.push_back( (*cit) * mult);
	}
		
	std::swap(this->container, tmp_container);
		
	//transform (container.begin(), container.end(), container.begin(),
	//			bind1st(multiplies<double>(),mult) );

	/*
	RangeCollectionHistPrivateTypeItr rit = container.begin();
	for (RangeCollectionHistPrivateTypeItr rit = container.begin();
					rit < container.end(); rit++) {
				*rit = mult * (*rit);
	}
	*/
}


// Full specializations of the templates in std namespace can be added in std namespace.

template <>
void std::swap<subpavings::RangeCollectionHist>(subpavings::RangeCollectionHist & rh1, 
		subpavings::RangeCollectionHist & rh2) // throw ()
{
	rh1.swap (rh2);
}
