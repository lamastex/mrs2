/*
* Copyright (C) 2010 Jennifer Harlow
*
*/

/*! \file
\brief MappedFobj10D definitions.
*/


#include "mappedFobj10D.hpp"

using namespace subpavings;

cxsc::interval MappedFobj10D::operator()(const cxsc::ivector& ivec) const
{
	return this->operator()(
		ivec[1],
		ivec[2],
		ivec[3],
		ivec[4],
		ivec[5],
		ivec[6],
		ivec[7],
		ivec[8],
		ivec[9],
		ivec[10]);
}

cxsc::real MappedFobj10D::operator()(const cxsc::rvector& rv) const
{
	return this->operator()(
		rv[1],
		rv[2],
		rv[3],
		rv[4],
		rv[5],
		rv[6],
		rv[7],
		rv[8],
		rv[9],
		rv[10]);
}

cxsc::real MappedFobj10D::imageMid(const cxsc::ivector& ivec) const
{
	return this->operator()(
		cxsc::mid(ivec[1]),
		cxsc::mid(ivec[2]),
		cxsc::mid(ivec[3]),
		cxsc::mid(ivec[4]),
		cxsc::mid(ivec[5]),
		cxsc::mid(ivec[6]),
		cxsc::mid(ivec[7]),
		cxsc::mid(ivec[8]),
		cxsc::mid(ivec[9]),
		cxsc::mid(ivec[10]));
}

cxsc::real MappedFobj10D::imageMid(
	const cxsc::interval& ival1,
	const cxsc::interval& ival2,
	const cxsc::interval& ival3,
	const cxsc::interval& ival4,
	const cxsc::interval& ival5,
	const cxsc::interval& ival6,
	const cxsc::interval& ival7,
	const cxsc::interval& ival8,
	const cxsc::interval& ival9,
	const cxsc::interval& ival10) const
{
	return this->operator()(
		cxsc::mid(ival1),
		cxsc::mid(ival2),
		cxsc::mid(ival3),
		cxsc::mid(ival4),
		cxsc::mid(ival5),
		cxsc::mid(ival6),
		cxsc::mid(ival7),
		cxsc::mid(ival8),
		cxsc::mid(ival9),
		cxsc::mid(ival10));
	/*int ivL = Lb(iv);
	int ivU = Ub(iv);
	rvector rv(ivU - ivL +1 );
	for (int i = ivL; i<=ivU; i++) {
		rv[i] = mid(iv[i]);
	}
	return this->operator()(rv);
	*/
}

