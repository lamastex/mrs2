/*! \file RosenFobj100D.cpp
\brief MappedSPnode example 2-d function object class.

This example is for the standard bivariate Rosen.

*/

#include "RosenFobj100D.hpp"
#include <cmath> //to use M_PI
#include <vector>
//#include "cxsc.hpp"

using namespace cxsc;
using namespace std;
using namespace subpavings;


interval RosenFobj100D::operator()(const cxsc::interval& ival1,
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
	 std::vector<interval> R;
   R.push_back(ival1);
R.push_back(ival2);
R.push_back(ival3);
R.push_back(ival4);
R.push_back(ival5);
R.push_back(ival6);
R.push_back(ival7);
R.push_back(ival8);
R.push_back(ival9);
R.push_back(ival10);
R.push_back(ival11);
R.push_back(ival12);
R.push_back(ival13);
R.push_back(ival14);
R.push_back(ival15);
R.push_back(ival16);
R.push_back(ival17);
R.push_back(ival18);
R.push_back(ival19);
R.push_back(ival20);
R.push_back(ival21);
R.push_back(ival22);
R.push_back(ival23);
R.push_back(ival24);
R.push_back(ival25);
R.push_back(ival26);
R.push_back(ival27);
R.push_back(ival28);
R.push_back(ival29);
R.push_back(ival30);
R.push_back(ival31);
R.push_back(ival32);
R.push_back(ival33);
R.push_back(ival34);
R.push_back(ival35);
R.push_back(ival36);
R.push_back(ival37);
R.push_back(ival38);
R.push_back(ival39);
R.push_back(ival40);
R.push_back(ival41);
R.push_back(ival42);
R.push_back(ival43);
R.push_back(ival44);
R.push_back(ival45);
R.push_back(ival46);
R.push_back(ival47);
R.push_back(ival48);
R.push_back(ival49);
R.push_back(ival50);
R.push_back(ival51);
R.push_back(ival52);
R.push_back(ival53);
R.push_back(ival54);
R.push_back(ival55);
R.push_back(ival56);
R.push_back(ival57);
R.push_back(ival58);
R.push_back(ival59);
R.push_back(ival60);
R.push_back(ival61);
R.push_back(ival62);
R.push_back(ival63);
R.push_back(ival64);
R.push_back(ival65);
R.push_back(ival66);
R.push_back(ival67);
R.push_back(ival68);
R.push_back(ival69);
R.push_back(ival70);
R.push_back(ival71);
R.push_back(ival72);
R.push_back(ival73);
R.push_back(ival74);
R.push_back(ival75);
R.push_back(ival76);
R.push_back(ival77);
R.push_back(ival78);
R.push_back(ival79);
R.push_back(ival80);
R.push_back(ival81);
R.push_back(ival82);
R.push_back(ival83);
R.push_back(ival84);
R.push_back(ival85);
R.push_back(ival86);
R.push_back(ival87);
R.push_back(ival88);
R.push_back(ival89);
R.push_back(ival90);
R.push_back(ival91);
R.push_back(ival92);
R.push_back(ival93);
R.push_back(ival94);
R.push_back(ival95);
R.push_back(ival96);
R.push_back(ival97);
R.push_back(ival98);
R.push_back(ival99);
R.push_back(ival100);

    interval result(0,0);
    
    for (size_t i = 0; i < (R.size()-1); i++) 
    {
      result += 100.0 * (power((R[i+1] - power(R[i],2)),2)) +
									power((R[i] - 1), 2);
    }
  
	result = exp (-(1.0 * result));
	
	return result;
	
}

real RosenFobj100D::operator()(const cxsc::real& r1,
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
const cxsc::real& r100) const
{
    std::vector<real> R;
R.push_back(r1);
R.push_back(r2);
R.push_back(r3);
R.push_back(r4);
R.push_back(r5);
R.push_back(r6);
R.push_back(r7);
R.push_back(r8);
R.push_back(r9);
R.push_back(r10);
R.push_back(r11);
R.push_back(r12);
R.push_back(r13);
R.push_back(r14);
R.push_back(r15);
R.push_back(r16);
R.push_back(r17);
R.push_back(r18);
R.push_back(r19);
R.push_back(r20);
R.push_back(r21);
R.push_back(r22);
R.push_back(r23);
R.push_back(r24);
R.push_back(r25);
R.push_back(r26);
R.push_back(r27);
R.push_back(r28);
R.push_back(r29);
R.push_back(r30);
R.push_back(r31);
R.push_back(r32);
R.push_back(r33);
R.push_back(r34);
R.push_back(r35);
R.push_back(r36);
R.push_back(r37);
R.push_back(r38);
R.push_back(r39);
R.push_back(r40);
R.push_back(r41);
R.push_back(r42);
R.push_back(r43);
R.push_back(r44);
R.push_back(r45);
R.push_back(r46);
R.push_back(r47);
R.push_back(r48);
R.push_back(r49);
R.push_back(r50);
R.push_back(r51);
R.push_back(r52);
R.push_back(r53);
R.push_back(r54);
R.push_back(r55);
R.push_back(r56);
R.push_back(r57);
R.push_back(r58);
R.push_back(r59);
R.push_back(r60);
R.push_back(r61);
R.push_back(r62);
R.push_back(r63);
R.push_back(r64);
R.push_back(r65);
R.push_back(r66);
R.push_back(r67);
R.push_back(r68);
R.push_back(r69);
R.push_back(r70);
R.push_back(r71);
R.push_back(r72);
R.push_back(r73);
R.push_back(r74);
R.push_back(r75);
R.push_back(r76);
R.push_back(r77);
R.push_back(r78);
R.push_back(r79);
R.push_back(r80);
R.push_back(r81);
R.push_back(r82);
R.push_back(r83);
R.push_back(r84);
R.push_back(r85);
R.push_back(r86);
R.push_back(r87);
R.push_back(r88);
R.push_back(r89);
R.push_back(r90);
R.push_back(r91);
R.push_back(r92);
R.push_back(r93);
R.push_back(r94);
R.push_back(r95);
R.push_back(r96);
R.push_back(r97);
R.push_back(r98);
R.push_back(r99);
R.push_back(r100);
    
    real result = 0;
    
    for (size_t i = 0; i < (R.size()-1); i++) 
    {
      result += 100.0 * (power((R[i+1] - power(R[i],2)),2)) +
									power((R[i] - 1), 2);
    }

	result = exp (-(1.0 * result));

   return result;
}
