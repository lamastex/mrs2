/*! \file GaussianFobj100D.cpp
\brief MappedSPnode example 2-d function object class.

This example is for the standard bivariate gaussian.

*/

#include "GaussianFobj100D.hpp"
#include <cmath> //to use M_PI
//#include "cxsc.hpp"

using namespace cxsc;
using namespace std;
using namespace subpavings;

interval GaussianFobj100D::operator()(
			const cxsc::interval& ival1,
			const cxsc::interval& ival2,
			const cxsc::interval& ival3,
			const cxsc::interval& ival4,
			const cxsc::interval& ival5,
			const cxsc::interval& ival6,
			const cxsc::interval& ival7,
			const cxsc::interval& ival8,
			const cxsc::interval& ival9,
			const cxsc::interval& ival10,
			const cxsc::interval& ival11,
			const cxsc::interval& ival12,
			const cxsc::interval& ival13,
			const cxsc::interval& ival14,
			const cxsc::interval& ival15,
			const cxsc::interval& ival16,
			const cxsc::interval& ival17,
			const cxsc::interval& ival18,
			const cxsc::interval& ival19,
			const cxsc::interval& ival20,
			const cxsc::interval& ival21,
			const cxsc::interval& ival22,
			const cxsc::interval& ival23,
			const cxsc::interval& ival24,
			const cxsc::interval& ival25,
			const cxsc::interval& ival26,
			const cxsc::interval& ival27,
			const cxsc::interval& ival28,
			const cxsc::interval& ival29,
			const cxsc::interval& ival30,
			const cxsc::interval& ival31,
			const cxsc::interval& ival32,
			const cxsc::interval& ival33,
			const cxsc::interval& ival34,
			const cxsc::interval& ival35,
			const cxsc::interval& ival36,
			const cxsc::interval& ival37,
			const cxsc::interval& ival38,
			const cxsc::interval& ival39,
			const cxsc::interval& ival40,
			const cxsc::interval& ival41,
			const cxsc::interval& ival42,
			const cxsc::interval& ival43,
			const cxsc::interval& ival44,
			const cxsc::interval& ival45,
			const cxsc::interval& ival46,
			const cxsc::interval& ival47,
			const cxsc::interval& ival48,
			const cxsc::interval& ival49,
			const cxsc::interval& ival50,
			const cxsc::interval& ival51,
			const cxsc::interval& ival52,
			const cxsc::interval& ival53,
			const cxsc::interval& ival54,
			const cxsc::interval& ival55,
			const cxsc::interval& ival56,
			const cxsc::interval& ival57,
			const cxsc::interval& ival58,
			const cxsc::interval& ival59,
			const cxsc::interval& ival60,
			const cxsc::interval& ival61,
			const cxsc::interval& ival62,
			const cxsc::interval& ival63,
			const cxsc::interval& ival64,
			const cxsc::interval& ival65,
			const cxsc::interval& ival66,
			const cxsc::interval& ival67,
			const cxsc::interval& ival68,
			const cxsc::interval& ival69,
			const cxsc::interval& ival70,
			const cxsc::interval& ival71,
			const cxsc::interval& ival72,
			const cxsc::interval& ival73,
			const cxsc::interval& ival74,
			const cxsc::interval& ival75,
			const cxsc::interval& ival76,
			const cxsc::interval& ival77,
			const cxsc::interval& ival78,
			const cxsc::interval& ival79,
			const cxsc::interval& ival80,
			const cxsc::interval& ival81,
			const cxsc::interval& ival82,
			const cxsc::interval& ival83,
			const cxsc::interval& ival84,
			const cxsc::interval& ival85,
			const cxsc::interval& ival86,
			const cxsc::interval& ival87,
			const cxsc::interval& ival88,
			const cxsc::interval& ival89,
			const cxsc::interval& ival90,
			const cxsc::interval& ival91,
			const cxsc::interval& ival92,
			const cxsc::interval& ival93,
			const cxsc::interval& ival94,
			const cxsc::interval& ival95,
			const cxsc::interval& ival96,
			const cxsc::interval& ival97,
			const cxsc::interval& ival98,
			const cxsc::interval& ival99,
			const cxsc::interval& ival100
			) const
{
	real a = power(2*M_PI, 10.0/2.0);
   interval b = -0.5 * (power(ival1,2)+
power(ival2,2)+
power(ival3,2)+
power(ival4,2)+
power(ival5,2)+
power(ival6,2)+
power(ival7,2)+
power(ival8,2)+
power(ival9,2)+
power(ival10,2)+
power(ival11,2)+
power(ival12,2)+
power(ival13,2)+
power(ival14,2)+
power(ival15,2)+
power(ival16,2)+
power(ival17,2)+
power(ival18,2)+
power(ival19,2)+
power(ival20,2)+
power(ival21,2)+
power(ival22,2)+
power(ival23,2)+
power(ival24,2)+
power(ival25,2)+
power(ival26,2)+
power(ival27,2)+
power(ival28,2)+
power(ival29,2)+
power(ival30,2)+
power(ival31,2)+
power(ival32,2)+
power(ival33,2)+
power(ival34,2)+
power(ival35,2)+
power(ival36,2)+
power(ival37,2)+
power(ival38,2)+
power(ival39,2)+
power(ival40,2)+
power(ival41,2)+
power(ival42,2)+
power(ival43,2)+
power(ival44,2)+
power(ival45,2)+
power(ival46,2)+
power(ival47,2)+
power(ival48,2)+
power(ival49,2)+
power(ival50,2)+
power(ival51,2)+
power(ival52,2)+
power(ival53,2)+
power(ival54,2)+
power(ival55,2)+
power(ival56,2)+
power(ival57,2)+
power(ival58,2)+
power(ival59,2)+
power(ival60,2)+
power(ival61,2)+
power(ival62,2)+
power(ival63,2)+
power(ival64,2)+
power(ival65,2)+
power(ival66,2)+
power(ival67,2)+
power(ival68,2)+
power(ival69,2)+
power(ival70,2)+
power(ival71,2)+
power(ival72,2)+
power(ival73,2)+
power(ival74,2)+
power(ival75,2)+
power(ival76,2)+
power(ival77,2)+
power(ival78,2)+
power(ival79,2)+
power(ival80,2)+
power(ival81,2)+
power(ival82,2)+
power(ival83,2)+
power(ival84,2)+
power(ival85,2)+
power(ival86,2)+
power(ival87,2)+
power(ival88,2)+
power(ival89,2)+
power(ival90,2)+
power(ival91,2)+
power(ival92,2)+
power(ival93,2)+
power(ival94,2)+
power(ival95,2)+
power(ival96,2)+
power(ival97,2)+
power(ival98,2)+
power(ival99,2)+
power(ival100,2)
);
   interval IntPDF = 1.0/a * exp(b);

	return IntPDF;
}

real GaussianFobj100D::operator()(
			const cxsc::real& r1,
			const cxsc::real& r2,
			const cxsc::real& r3,
			const cxsc::real& r4,
			const cxsc::real& r5,
			const cxsc::real& r6,
			const cxsc::real& r7,
			const cxsc::real& r8,
			const cxsc::real& r9,
			const cxsc::real& r10,
			const cxsc::real& r11,
			const cxsc::real& r12,
			const cxsc::real& r13,
			const cxsc::real& r14,
			const cxsc::real& r15,
			const cxsc::real& r16,
			const cxsc::real& r17,
			const cxsc::real& r18,
			const cxsc::real& r19,
			const cxsc::real& r20,
			const cxsc::real& r21,
			const cxsc::real& r22,
			const cxsc::real& r23,
			const cxsc::real& r24,
			const cxsc::real& r25,
			const cxsc::real& r26,
			const cxsc::real& r27,
			const cxsc::real& r28,
			const cxsc::real& r29,
			const cxsc::real& r30,
			const cxsc::real& r31,
			const cxsc::real& r32,
			const cxsc::real& r33,
			const cxsc::real& r34,
			const cxsc::real& r35,
			const cxsc::real& r36,
			const cxsc::real& r37,
			const cxsc::real& r38,
			const cxsc::real& r39,
			const cxsc::real& r40,
			const cxsc::real& r41,
			const cxsc::real& r42,
			const cxsc::real& r43,
			const cxsc::real& r44,
			const cxsc::real& r45,
			const cxsc::real& r46,
			const cxsc::real& r47,
			const cxsc::real& r48,
			const cxsc::real& r49,
			const cxsc::real& r50,
			const cxsc::real& r51,
			const cxsc::real& r52,
			const cxsc::real& r53,
			const cxsc::real& r54,
			const cxsc::real& r55,
			const cxsc::real& r56,
			const cxsc::real& r57,
			const cxsc::real& r58,
			const cxsc::real& r59,
			const cxsc::real& r60,
			const cxsc::real& r61,
			const cxsc::real& r62,
			const cxsc::real& r63,
			const cxsc::real& r64,
			const cxsc::real& r65,
			const cxsc::real& r66,
			const cxsc::real& r67,
			const cxsc::real& r68,
			const cxsc::real& r69,
			const cxsc::real& r70,
			const cxsc::real& r71,
			const cxsc::real& r72,
			const cxsc::real& r73,
			const cxsc::real& r74,
			const cxsc::real& r75,
			const cxsc::real& r76,
			const cxsc::real& r77,
			const cxsc::real& r78,
			const cxsc::real& r79,
			const cxsc::real& r80,
			const cxsc::real& r81,
			const cxsc::real& r82,
			const cxsc::real& r83,
			const cxsc::real& r84,
			const cxsc::real& r85,
			const cxsc::real& r86,
			const cxsc::real& r87,
			const cxsc::real& r88,
			const cxsc::real& r89,
			const cxsc::real& r90,
			const cxsc::real& r91,
			const cxsc::real& r92,
			const cxsc::real& r93,
			const cxsc::real& r94,
			const cxsc::real& r95,
			const cxsc::real& r96,
			const cxsc::real& r97,
			const cxsc::real& r98,
			const cxsc::real& r99,
			const cxsc::real& r100
			) const
{
    real a = power(2*M_PI, 10.0/2.0);
    real b = -0.5 * (		power(r1,2)+
		power(r2,2)+
		power(r3,2)+
		power(r4,2)+
		power(r5,2)+
		power(r6,2)+
		power(r7,2)+
		power(r8,2)+
		power(r9,2)+
		power(r10,2)+
		power(r11,2)+
		power(r12,2)+
		power(r13,2)+
		power(r14,2)+
		power(r15,2)+
		power(r16,2)+
		power(r17,2)+
		power(r18,2)+
		power(r19,2)+
		power(r20,2)+
		power(r21,2)+
		power(r22,2)+
		power(r23,2)+
		power(r24,2)+
		power(r25,2)+
		power(r26,2)+
		power(r27,2)+
		power(r28,2)+
		power(r29,2)+
		power(r30,2)+
		power(r31,2)+
		power(r32,2)+
		power(r33,2)+
		power(r34,2)+
		power(r35,2)+
		power(r36,2)+
		power(r37,2)+
		power(r38,2)+
		power(r39,2)+
		power(r40,2)+
		power(r41,2)+
		power(r42,2)+
		power(r43,2)+
		power(r44,2)+
		power(r45,2)+
		power(r46,2)+
		power(r47,2)+
		power(r48,2)+
		power(r49,2)+
		power(r50,2)+
		power(r51,2)+
		power(r52,2)+
		power(r53,2)+
		power(r54,2)+
		power(r55,2)+
		power(r56,2)+
		power(r57,2)+
		power(r58,2)+
		power(r59,2)+
		power(r60,2)+
		power(r61,2)+
		power(r62,2)+
		power(r63,2)+
		power(r64,2)+
		power(r65,2)+
		power(r66,2)+
		power(r67,2)+
		power(r68,2)+
		power(r69,2)+
		power(r70,2)+
		power(r71,2)+
		power(r72,2)+
		power(r73,2)+
		power(r74,2)+
		power(r75,2)+
		power(r76,2)+
		power(r77,2)+
		power(r78,2)+
		power(r79,2)+
		power(r80,2)+
		power(r81,2)+
		power(r82,2)+
		power(r83,2)+
		power(r84,2)+
		power(r85,2)+
		power(r86,2)+
		power(r87,2)+
		power(r88,2)+
		power(r89,2)+
		power(r90,2)+
		power(r91,2)+
		power(r92,2)+
		power(r93,2)+
		power(r94,2)+
		power(r95,2)+
		power(r96,2)+
		power(r97,2)+
		power(r98,2)+
		power(r99,2)+
		power(r100,2)
		);
    real RePDF = 1.0/a * exp(b);
    
    return RePDF;
}

