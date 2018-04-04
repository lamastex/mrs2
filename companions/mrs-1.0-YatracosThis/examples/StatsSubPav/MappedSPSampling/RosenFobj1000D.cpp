/*! \file RosenFobj1000D.cpp
\brief MappedSPnode example 2-d function object class.

This example is for the standard bivariate Rosen.

*/

#include "RosenFobj1000D.hpp"
#include <cmath> //to use M_PI
#include <vector>
//#include "cxsc.hpp"

using namespace cxsc;
using namespace std;
using namespace subpavings;


interval RosenFobj1000D::operator()(const cxsc::interval& ival1,
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
const cxsc::interval& ival100,
const cxsc::interval& ival101,
const cxsc::interval& ival102,
const cxsc::interval& ival103,
const cxsc::interval& ival104,
const cxsc::interval& ival105,
const cxsc::interval& ival106,
const cxsc::interval& ival107,
const cxsc::interval& ival108,
const cxsc::interval& ival109,
const cxsc::interval& ival110,
const cxsc::interval& ival111,
const cxsc::interval& ival112,
const cxsc::interval& ival113,
const cxsc::interval& ival114,
const cxsc::interval& ival115,
const cxsc::interval& ival116,
const cxsc::interval& ival117,
const cxsc::interval& ival118,
const cxsc::interval& ival119,
const cxsc::interval& ival120,
const cxsc::interval& ival121,
const cxsc::interval& ival122,
const cxsc::interval& ival123,
const cxsc::interval& ival124,
const cxsc::interval& ival125,
const cxsc::interval& ival126,
const cxsc::interval& ival127,
const cxsc::interval& ival128,
const cxsc::interval& ival129,
const cxsc::interval& ival130,
const cxsc::interval& ival131,
const cxsc::interval& ival132,
const cxsc::interval& ival133,
const cxsc::interval& ival134,
const cxsc::interval& ival135,
const cxsc::interval& ival136,
const cxsc::interval& ival137,
const cxsc::interval& ival138,
const cxsc::interval& ival139,
const cxsc::interval& ival140,
const cxsc::interval& ival141,
const cxsc::interval& ival142,
const cxsc::interval& ival143,
const cxsc::interval& ival144,
const cxsc::interval& ival145,
const cxsc::interval& ival146,
const cxsc::interval& ival147,
const cxsc::interval& ival148,
const cxsc::interval& ival149,
const cxsc::interval& ival150,
const cxsc::interval& ival151,
const cxsc::interval& ival152,
const cxsc::interval& ival153,
const cxsc::interval& ival154,
const cxsc::interval& ival155,
const cxsc::interval& ival156,
const cxsc::interval& ival157,
const cxsc::interval& ival158,
const cxsc::interval& ival159,
const cxsc::interval& ival160,
const cxsc::interval& ival161,
const cxsc::interval& ival162,
const cxsc::interval& ival163,
const cxsc::interval& ival164,
const cxsc::interval& ival165,
const cxsc::interval& ival166,
const cxsc::interval& ival167,
const cxsc::interval& ival168,
const cxsc::interval& ival169,
const cxsc::interval& ival170,
const cxsc::interval& ival171,
const cxsc::interval& ival172,
const cxsc::interval& ival173,
const cxsc::interval& ival174,
const cxsc::interval& ival175,
const cxsc::interval& ival176,
const cxsc::interval& ival177,
const cxsc::interval& ival178,
const cxsc::interval& ival179,
const cxsc::interval& ival180,
const cxsc::interval& ival181,
const cxsc::interval& ival182,
const cxsc::interval& ival183,
const cxsc::interval& ival184,
const cxsc::interval& ival185,
const cxsc::interval& ival186,
const cxsc::interval& ival187,
const cxsc::interval& ival188,
const cxsc::interval& ival189,
const cxsc::interval& ival190,
const cxsc::interval& ival191,
const cxsc::interval& ival192,
const cxsc::interval& ival193,
const cxsc::interval& ival194,
const cxsc::interval& ival195,
const cxsc::interval& ival196,
const cxsc::interval& ival197,
const cxsc::interval& ival198,
const cxsc::interval& ival199,
const cxsc::interval& ival200,
const cxsc::interval& ival201,
const cxsc::interval& ival202,
const cxsc::interval& ival203,
const cxsc::interval& ival204,
const cxsc::interval& ival205,
const cxsc::interval& ival206,
const cxsc::interval& ival207,
const cxsc::interval& ival208,
const cxsc::interval& ival209,
const cxsc::interval& ival210,
const cxsc::interval& ival211,
const cxsc::interval& ival212,
const cxsc::interval& ival213,
const cxsc::interval& ival214,
const cxsc::interval& ival215,
const cxsc::interval& ival216,
const cxsc::interval& ival217,
const cxsc::interval& ival218,
const cxsc::interval& ival219,
const cxsc::interval& ival220,
const cxsc::interval& ival221,
const cxsc::interval& ival222,
const cxsc::interval& ival223,
const cxsc::interval& ival224,
const cxsc::interval& ival225,
const cxsc::interval& ival226,
const cxsc::interval& ival227,
const cxsc::interval& ival228,
const cxsc::interval& ival229,
const cxsc::interval& ival230,
const cxsc::interval& ival231,
const cxsc::interval& ival232,
const cxsc::interval& ival233,
const cxsc::interval& ival234,
const cxsc::interval& ival235,
const cxsc::interval& ival236,
const cxsc::interval& ival237,
const cxsc::interval& ival238,
const cxsc::interval& ival239,
const cxsc::interval& ival240,
const cxsc::interval& ival241,
const cxsc::interval& ival242,
const cxsc::interval& ival243,
const cxsc::interval& ival244,
const cxsc::interval& ival245,
const cxsc::interval& ival246,
const cxsc::interval& ival247,
const cxsc::interval& ival248,
const cxsc::interval& ival249,
const cxsc::interval& ival250,
const cxsc::interval& ival251,
const cxsc::interval& ival252,
const cxsc::interval& ival253,
const cxsc::interval& ival254,
const cxsc::interval& ival255,
const cxsc::interval& ival256,
const cxsc::interval& ival257,
const cxsc::interval& ival258,
const cxsc::interval& ival259,
const cxsc::interval& ival260,
const cxsc::interval& ival261,
const cxsc::interval& ival262,
const cxsc::interval& ival263,
const cxsc::interval& ival264,
const cxsc::interval& ival265,
const cxsc::interval& ival266,
const cxsc::interval& ival267,
const cxsc::interval& ival268,
const cxsc::interval& ival269,
const cxsc::interval& ival270,
const cxsc::interval& ival271,
const cxsc::interval& ival272,
const cxsc::interval& ival273,
const cxsc::interval& ival274,
const cxsc::interval& ival275,
const cxsc::interval& ival276,
const cxsc::interval& ival277,
const cxsc::interval& ival278,
const cxsc::interval& ival279,
const cxsc::interval& ival280,
const cxsc::interval& ival281,
const cxsc::interval& ival282,
const cxsc::interval& ival283,
const cxsc::interval& ival284,
const cxsc::interval& ival285,
const cxsc::interval& ival286,
const cxsc::interval& ival287,
const cxsc::interval& ival288,
const cxsc::interval& ival289,
const cxsc::interval& ival290,
const cxsc::interval& ival291,
const cxsc::interval& ival292,
const cxsc::interval& ival293,
const cxsc::interval& ival294,
const cxsc::interval& ival295,
const cxsc::interval& ival296,
const cxsc::interval& ival297,
const cxsc::interval& ival298,
const cxsc::interval& ival299,
const cxsc::interval& ival300,
const cxsc::interval& ival301,
const cxsc::interval& ival302,
const cxsc::interval& ival303,
const cxsc::interval& ival304,
const cxsc::interval& ival305,
const cxsc::interval& ival306,
const cxsc::interval& ival307,
const cxsc::interval& ival308,
const cxsc::interval& ival309,
const cxsc::interval& ival310,
const cxsc::interval& ival311,
const cxsc::interval& ival312,
const cxsc::interval& ival313,
const cxsc::interval& ival314,
const cxsc::interval& ival315,
const cxsc::interval& ival316,
const cxsc::interval& ival317,
const cxsc::interval& ival318,
const cxsc::interval& ival319,
const cxsc::interval& ival320,
const cxsc::interval& ival321,
const cxsc::interval& ival322,
const cxsc::interval& ival323,
const cxsc::interval& ival324,
const cxsc::interval& ival325,
const cxsc::interval& ival326,
const cxsc::interval& ival327,
const cxsc::interval& ival328,
const cxsc::interval& ival329,
const cxsc::interval& ival330,
const cxsc::interval& ival331,
const cxsc::interval& ival332,
const cxsc::interval& ival333,
const cxsc::interval& ival334,
const cxsc::interval& ival335,
const cxsc::interval& ival336,
const cxsc::interval& ival337,
const cxsc::interval& ival338,
const cxsc::interval& ival339,
const cxsc::interval& ival340,
const cxsc::interval& ival341,
const cxsc::interval& ival342,
const cxsc::interval& ival343,
const cxsc::interval& ival344,
const cxsc::interval& ival345,
const cxsc::interval& ival346,
const cxsc::interval& ival347,
const cxsc::interval& ival348,
const cxsc::interval& ival349,
const cxsc::interval& ival350,
const cxsc::interval& ival351,
const cxsc::interval& ival352,
const cxsc::interval& ival353,
const cxsc::interval& ival354,
const cxsc::interval& ival355,
const cxsc::interval& ival356,
const cxsc::interval& ival357,
const cxsc::interval& ival358,
const cxsc::interval& ival359,
const cxsc::interval& ival360,
const cxsc::interval& ival361,
const cxsc::interval& ival362,
const cxsc::interval& ival363,
const cxsc::interval& ival364,
const cxsc::interval& ival365,
const cxsc::interval& ival366,
const cxsc::interval& ival367,
const cxsc::interval& ival368,
const cxsc::interval& ival369,
const cxsc::interval& ival370,
const cxsc::interval& ival371,
const cxsc::interval& ival372,
const cxsc::interval& ival373,
const cxsc::interval& ival374,
const cxsc::interval& ival375,
const cxsc::interval& ival376,
const cxsc::interval& ival377,
const cxsc::interval& ival378,
const cxsc::interval& ival379,
const cxsc::interval& ival380,
const cxsc::interval& ival381,
const cxsc::interval& ival382,
const cxsc::interval& ival383,
const cxsc::interval& ival384,
const cxsc::interval& ival385,
const cxsc::interval& ival386,
const cxsc::interval& ival387,
const cxsc::interval& ival388,
const cxsc::interval& ival389,
const cxsc::interval& ival390,
const cxsc::interval& ival391,
const cxsc::interval& ival392,
const cxsc::interval& ival393,
const cxsc::interval& ival394,
const cxsc::interval& ival395,
const cxsc::interval& ival396,
const cxsc::interval& ival397,
const cxsc::interval& ival398,
const cxsc::interval& ival399,
const cxsc::interval& ival400,
const cxsc::interval& ival401,
const cxsc::interval& ival402,
const cxsc::interval& ival403,
const cxsc::interval& ival404,
const cxsc::interval& ival405,
const cxsc::interval& ival406,
const cxsc::interval& ival407,
const cxsc::interval& ival408,
const cxsc::interval& ival409,
const cxsc::interval& ival410,
const cxsc::interval& ival411,
const cxsc::interval& ival412,
const cxsc::interval& ival413,
const cxsc::interval& ival414,
const cxsc::interval& ival415,
const cxsc::interval& ival416,
const cxsc::interval& ival417,
const cxsc::interval& ival418,
const cxsc::interval& ival419,
const cxsc::interval& ival420,
const cxsc::interval& ival421,
const cxsc::interval& ival422,
const cxsc::interval& ival423,
const cxsc::interval& ival424,
const cxsc::interval& ival425,
const cxsc::interval& ival426,
const cxsc::interval& ival427,
const cxsc::interval& ival428,
const cxsc::interval& ival429,
const cxsc::interval& ival430,
const cxsc::interval& ival431,
const cxsc::interval& ival432,
const cxsc::interval& ival433,
const cxsc::interval& ival434,
const cxsc::interval& ival435,
const cxsc::interval& ival436,
const cxsc::interval& ival437,
const cxsc::interval& ival438,
const cxsc::interval& ival439,
const cxsc::interval& ival440,
const cxsc::interval& ival441,
const cxsc::interval& ival442,
const cxsc::interval& ival443,
const cxsc::interval& ival444,
const cxsc::interval& ival445,
const cxsc::interval& ival446,
const cxsc::interval& ival447,
const cxsc::interval& ival448,
const cxsc::interval& ival449,
const cxsc::interval& ival450,
const cxsc::interval& ival451,
const cxsc::interval& ival452,
const cxsc::interval& ival453,
const cxsc::interval& ival454,
const cxsc::interval& ival455,
const cxsc::interval& ival456,
const cxsc::interval& ival457,
const cxsc::interval& ival458,
const cxsc::interval& ival459,
const cxsc::interval& ival460,
const cxsc::interval& ival461,
const cxsc::interval& ival462,
const cxsc::interval& ival463,
const cxsc::interval& ival464,
const cxsc::interval& ival465,
const cxsc::interval& ival466,
const cxsc::interval& ival467,
const cxsc::interval& ival468,
const cxsc::interval& ival469,
const cxsc::interval& ival470,
const cxsc::interval& ival471,
const cxsc::interval& ival472,
const cxsc::interval& ival473,
const cxsc::interval& ival474,
const cxsc::interval& ival475,
const cxsc::interval& ival476,
const cxsc::interval& ival477,
const cxsc::interval& ival478,
const cxsc::interval& ival479,
const cxsc::interval& ival480,
const cxsc::interval& ival481,
const cxsc::interval& ival482,
const cxsc::interval& ival483,
const cxsc::interval& ival484,
const cxsc::interval& ival485,
const cxsc::interval& ival486,
const cxsc::interval& ival487,
const cxsc::interval& ival488,
const cxsc::interval& ival489,
const cxsc::interval& ival490,
const cxsc::interval& ival491,
const cxsc::interval& ival492,
const cxsc::interval& ival493,
const cxsc::interval& ival494,
const cxsc::interval& ival495,
const cxsc::interval& ival496,
const cxsc::interval& ival497,
const cxsc::interval& ival498,
const cxsc::interval& ival499,
const cxsc::interval& ival500,
const cxsc::interval& ival501,
const cxsc::interval& ival502,
const cxsc::interval& ival503,
const cxsc::interval& ival504,
const cxsc::interval& ival505,
const cxsc::interval& ival506,
const cxsc::interval& ival507,
const cxsc::interval& ival508,
const cxsc::interval& ival509,
const cxsc::interval& ival510,
const cxsc::interval& ival511,
const cxsc::interval& ival512,
const cxsc::interval& ival513,
const cxsc::interval& ival514,
const cxsc::interval& ival515,
const cxsc::interval& ival516,
const cxsc::interval& ival517,
const cxsc::interval& ival518,
const cxsc::interval& ival519,
const cxsc::interval& ival520,
const cxsc::interval& ival521,
const cxsc::interval& ival522,
const cxsc::interval& ival523,
const cxsc::interval& ival524,
const cxsc::interval& ival525,
const cxsc::interval& ival526,
const cxsc::interval& ival527,
const cxsc::interval& ival528,
const cxsc::interval& ival529,
const cxsc::interval& ival530,
const cxsc::interval& ival531,
const cxsc::interval& ival532,
const cxsc::interval& ival533,
const cxsc::interval& ival534,
const cxsc::interval& ival535,
const cxsc::interval& ival536,
const cxsc::interval& ival537,
const cxsc::interval& ival538,
const cxsc::interval& ival539,
const cxsc::interval& ival540,
const cxsc::interval& ival541,
const cxsc::interval& ival542,
const cxsc::interval& ival543,
const cxsc::interval& ival544,
const cxsc::interval& ival545,
const cxsc::interval& ival546,
const cxsc::interval& ival547,
const cxsc::interval& ival548,
const cxsc::interval& ival549,
const cxsc::interval& ival550,
const cxsc::interval& ival551,
const cxsc::interval& ival552,
const cxsc::interval& ival553,
const cxsc::interval& ival554,
const cxsc::interval& ival555,
const cxsc::interval& ival556,
const cxsc::interval& ival557,
const cxsc::interval& ival558,
const cxsc::interval& ival559,
const cxsc::interval& ival560,
const cxsc::interval& ival561,
const cxsc::interval& ival562,
const cxsc::interval& ival563,
const cxsc::interval& ival564,
const cxsc::interval& ival565,
const cxsc::interval& ival566,
const cxsc::interval& ival567,
const cxsc::interval& ival568,
const cxsc::interval& ival569,
const cxsc::interval& ival570,
const cxsc::interval& ival571,
const cxsc::interval& ival572,
const cxsc::interval& ival573,
const cxsc::interval& ival574,
const cxsc::interval& ival575,
const cxsc::interval& ival576,
const cxsc::interval& ival577,
const cxsc::interval& ival578,
const cxsc::interval& ival579,
const cxsc::interval& ival580,
const cxsc::interval& ival581,
const cxsc::interval& ival582,
const cxsc::interval& ival583,
const cxsc::interval& ival584,
const cxsc::interval& ival585,
const cxsc::interval& ival586,
const cxsc::interval& ival587,
const cxsc::interval& ival588,
const cxsc::interval& ival589,
const cxsc::interval& ival590,
const cxsc::interval& ival591,
const cxsc::interval& ival592,
const cxsc::interval& ival593,
const cxsc::interval& ival594,
const cxsc::interval& ival595,
const cxsc::interval& ival596,
const cxsc::interval& ival597,
const cxsc::interval& ival598,
const cxsc::interval& ival599,
const cxsc::interval& ival600,
const cxsc::interval& ival601,
const cxsc::interval& ival602,
const cxsc::interval& ival603,
const cxsc::interval& ival604,
const cxsc::interval& ival605,
const cxsc::interval& ival606,
const cxsc::interval& ival607,
const cxsc::interval& ival608,
const cxsc::interval& ival609,
const cxsc::interval& ival610,
const cxsc::interval& ival611,
const cxsc::interval& ival612,
const cxsc::interval& ival613,
const cxsc::interval& ival614,
const cxsc::interval& ival615,
const cxsc::interval& ival616,
const cxsc::interval& ival617,
const cxsc::interval& ival618,
const cxsc::interval& ival619,
const cxsc::interval& ival620,
const cxsc::interval& ival621,
const cxsc::interval& ival622,
const cxsc::interval& ival623,
const cxsc::interval& ival624,
const cxsc::interval& ival625,
const cxsc::interval& ival626,
const cxsc::interval& ival627,
const cxsc::interval& ival628,
const cxsc::interval& ival629,
const cxsc::interval& ival630,
const cxsc::interval& ival631,
const cxsc::interval& ival632,
const cxsc::interval& ival633,
const cxsc::interval& ival634,
const cxsc::interval& ival635,
const cxsc::interval& ival636,
const cxsc::interval& ival637,
const cxsc::interval& ival638,
const cxsc::interval& ival639,
const cxsc::interval& ival640,
const cxsc::interval& ival641,
const cxsc::interval& ival642,
const cxsc::interval& ival643,
const cxsc::interval& ival644,
const cxsc::interval& ival645,
const cxsc::interval& ival646,
const cxsc::interval& ival647,
const cxsc::interval& ival648,
const cxsc::interval& ival649,
const cxsc::interval& ival650,
const cxsc::interval& ival651,
const cxsc::interval& ival652,
const cxsc::interval& ival653,
const cxsc::interval& ival654,
const cxsc::interval& ival655,
const cxsc::interval& ival656,
const cxsc::interval& ival657,
const cxsc::interval& ival658,
const cxsc::interval& ival659,
const cxsc::interval& ival660,
const cxsc::interval& ival661,
const cxsc::interval& ival662,
const cxsc::interval& ival663,
const cxsc::interval& ival664,
const cxsc::interval& ival665,
const cxsc::interval& ival666,
const cxsc::interval& ival667,
const cxsc::interval& ival668,
const cxsc::interval& ival669,
const cxsc::interval& ival670,
const cxsc::interval& ival671,
const cxsc::interval& ival672,
const cxsc::interval& ival673,
const cxsc::interval& ival674,
const cxsc::interval& ival675,
const cxsc::interval& ival676,
const cxsc::interval& ival677,
const cxsc::interval& ival678,
const cxsc::interval& ival679,
const cxsc::interval& ival680,
const cxsc::interval& ival681,
const cxsc::interval& ival682,
const cxsc::interval& ival683,
const cxsc::interval& ival684,
const cxsc::interval& ival685,
const cxsc::interval& ival686,
const cxsc::interval& ival687,
const cxsc::interval& ival688,
const cxsc::interval& ival689,
const cxsc::interval& ival690,
const cxsc::interval& ival691,
const cxsc::interval& ival692,
const cxsc::interval& ival693,
const cxsc::interval& ival694,
const cxsc::interval& ival695,
const cxsc::interval& ival696,
const cxsc::interval& ival697,
const cxsc::interval& ival698,
const cxsc::interval& ival699,
const cxsc::interval& ival700,
const cxsc::interval& ival701,
const cxsc::interval& ival702,
const cxsc::interval& ival703,
const cxsc::interval& ival704,
const cxsc::interval& ival705,
const cxsc::interval& ival706,
const cxsc::interval& ival707,
const cxsc::interval& ival708,
const cxsc::interval& ival709,
const cxsc::interval& ival710,
const cxsc::interval& ival711,
const cxsc::interval& ival712,
const cxsc::interval& ival713,
const cxsc::interval& ival714,
const cxsc::interval& ival715,
const cxsc::interval& ival716,
const cxsc::interval& ival717,
const cxsc::interval& ival718,
const cxsc::interval& ival719,
const cxsc::interval& ival720,
const cxsc::interval& ival721,
const cxsc::interval& ival722,
const cxsc::interval& ival723,
const cxsc::interval& ival724,
const cxsc::interval& ival725,
const cxsc::interval& ival726,
const cxsc::interval& ival727,
const cxsc::interval& ival728,
const cxsc::interval& ival729,
const cxsc::interval& ival730,
const cxsc::interval& ival731,
const cxsc::interval& ival732,
const cxsc::interval& ival733,
const cxsc::interval& ival734,
const cxsc::interval& ival735,
const cxsc::interval& ival736,
const cxsc::interval& ival737,
const cxsc::interval& ival738,
const cxsc::interval& ival739,
const cxsc::interval& ival740,
const cxsc::interval& ival741,
const cxsc::interval& ival742,
const cxsc::interval& ival743,
const cxsc::interval& ival744,
const cxsc::interval& ival745,
const cxsc::interval& ival746,
const cxsc::interval& ival747,
const cxsc::interval& ival748,
const cxsc::interval& ival749,
const cxsc::interval& ival750,
const cxsc::interval& ival751,
const cxsc::interval& ival752,
const cxsc::interval& ival753,
const cxsc::interval& ival754,
const cxsc::interval& ival755,
const cxsc::interval& ival756,
const cxsc::interval& ival757,
const cxsc::interval& ival758,
const cxsc::interval& ival759,
const cxsc::interval& ival760,
const cxsc::interval& ival761,
const cxsc::interval& ival762,
const cxsc::interval& ival763,
const cxsc::interval& ival764,
const cxsc::interval& ival765,
const cxsc::interval& ival766,
const cxsc::interval& ival767,
const cxsc::interval& ival768,
const cxsc::interval& ival769,
const cxsc::interval& ival770,
const cxsc::interval& ival771,
const cxsc::interval& ival772,
const cxsc::interval& ival773,
const cxsc::interval& ival774,
const cxsc::interval& ival775,
const cxsc::interval& ival776,
const cxsc::interval& ival777,
const cxsc::interval& ival778,
const cxsc::interval& ival779,
const cxsc::interval& ival780,
const cxsc::interval& ival781,
const cxsc::interval& ival782,
const cxsc::interval& ival783,
const cxsc::interval& ival784,
const cxsc::interval& ival785,
const cxsc::interval& ival786,
const cxsc::interval& ival787,
const cxsc::interval& ival788,
const cxsc::interval& ival789,
const cxsc::interval& ival790,
const cxsc::interval& ival791,
const cxsc::interval& ival792,
const cxsc::interval& ival793,
const cxsc::interval& ival794,
const cxsc::interval& ival795,
const cxsc::interval& ival796,
const cxsc::interval& ival797,
const cxsc::interval& ival798,
const cxsc::interval& ival799,
const cxsc::interval& ival800,
const cxsc::interval& ival801,
const cxsc::interval& ival802,
const cxsc::interval& ival803,
const cxsc::interval& ival804,
const cxsc::interval& ival805,
const cxsc::interval& ival806,
const cxsc::interval& ival807,
const cxsc::interval& ival808,
const cxsc::interval& ival809,
const cxsc::interval& ival810,
const cxsc::interval& ival811,
const cxsc::interval& ival812,
const cxsc::interval& ival813,
const cxsc::interval& ival814,
const cxsc::interval& ival815,
const cxsc::interval& ival816,
const cxsc::interval& ival817,
const cxsc::interval& ival818,
const cxsc::interval& ival819,
const cxsc::interval& ival820,
const cxsc::interval& ival821,
const cxsc::interval& ival822,
const cxsc::interval& ival823,
const cxsc::interval& ival824,
const cxsc::interval& ival825,
const cxsc::interval& ival826,
const cxsc::interval& ival827,
const cxsc::interval& ival828,
const cxsc::interval& ival829,
const cxsc::interval& ival830,
const cxsc::interval& ival831,
const cxsc::interval& ival832,
const cxsc::interval& ival833,
const cxsc::interval& ival834,
const cxsc::interval& ival835,
const cxsc::interval& ival836,
const cxsc::interval& ival837,
const cxsc::interval& ival838,
const cxsc::interval& ival839,
const cxsc::interval& ival840,
const cxsc::interval& ival841,
const cxsc::interval& ival842,
const cxsc::interval& ival843,
const cxsc::interval& ival844,
const cxsc::interval& ival845,
const cxsc::interval& ival846,
const cxsc::interval& ival847,
const cxsc::interval& ival848,
const cxsc::interval& ival849,
const cxsc::interval& ival850,
const cxsc::interval& ival851,
const cxsc::interval& ival852,
const cxsc::interval& ival853,
const cxsc::interval& ival854,
const cxsc::interval& ival855,
const cxsc::interval& ival856,
const cxsc::interval& ival857,
const cxsc::interval& ival858,
const cxsc::interval& ival859,
const cxsc::interval& ival860,
const cxsc::interval& ival861,
const cxsc::interval& ival862,
const cxsc::interval& ival863,
const cxsc::interval& ival864,
const cxsc::interval& ival865,
const cxsc::interval& ival866,
const cxsc::interval& ival867,
const cxsc::interval& ival868,
const cxsc::interval& ival869,
const cxsc::interval& ival870,
const cxsc::interval& ival871,
const cxsc::interval& ival872,
const cxsc::interval& ival873,
const cxsc::interval& ival874,
const cxsc::interval& ival875,
const cxsc::interval& ival876,
const cxsc::interval& ival877,
const cxsc::interval& ival878,
const cxsc::interval& ival879,
const cxsc::interval& ival880,
const cxsc::interval& ival881,
const cxsc::interval& ival882,
const cxsc::interval& ival883,
const cxsc::interval& ival884,
const cxsc::interval& ival885,
const cxsc::interval& ival886,
const cxsc::interval& ival887,
const cxsc::interval& ival888,
const cxsc::interval& ival889,
const cxsc::interval& ival890,
const cxsc::interval& ival891,
const cxsc::interval& ival892,
const cxsc::interval& ival893,
const cxsc::interval& ival894,
const cxsc::interval& ival895,
const cxsc::interval& ival896,
const cxsc::interval& ival897,
const cxsc::interval& ival898,
const cxsc::interval& ival899,
const cxsc::interval& ival900,
const cxsc::interval& ival901,
const cxsc::interval& ival902,
const cxsc::interval& ival903,
const cxsc::interval& ival904,
const cxsc::interval& ival905,
const cxsc::interval& ival906,
const cxsc::interval& ival907,
const cxsc::interval& ival908,
const cxsc::interval& ival909,
const cxsc::interval& ival910,
const cxsc::interval& ival911,
const cxsc::interval& ival912,
const cxsc::interval& ival913,
const cxsc::interval& ival914,
const cxsc::interval& ival915,
const cxsc::interval& ival916,
const cxsc::interval& ival917,
const cxsc::interval& ival918,
const cxsc::interval& ival919,
const cxsc::interval& ival920,
const cxsc::interval& ival921,
const cxsc::interval& ival922,
const cxsc::interval& ival923,
const cxsc::interval& ival924,
const cxsc::interval& ival925,
const cxsc::interval& ival926,
const cxsc::interval& ival927,
const cxsc::interval& ival928,
const cxsc::interval& ival929,
const cxsc::interval& ival930,
const cxsc::interval& ival931,
const cxsc::interval& ival932,
const cxsc::interval& ival933,
const cxsc::interval& ival934,
const cxsc::interval& ival935,
const cxsc::interval& ival936,
const cxsc::interval& ival937,
const cxsc::interval& ival938,
const cxsc::interval& ival939,
const cxsc::interval& ival940,
const cxsc::interval& ival941,
const cxsc::interval& ival942,
const cxsc::interval& ival943,
const cxsc::interval& ival944,
const cxsc::interval& ival945,
const cxsc::interval& ival946,
const cxsc::interval& ival947,
const cxsc::interval& ival948,
const cxsc::interval& ival949,
const cxsc::interval& ival950,
const cxsc::interval& ival951,
const cxsc::interval& ival952,
const cxsc::interval& ival953,
const cxsc::interval& ival954,
const cxsc::interval& ival955,
const cxsc::interval& ival956,
const cxsc::interval& ival957,
const cxsc::interval& ival958,
const cxsc::interval& ival959,
const cxsc::interval& ival960,
const cxsc::interval& ival961,
const cxsc::interval& ival962,
const cxsc::interval& ival963,
const cxsc::interval& ival964,
const cxsc::interval& ival965,
const cxsc::interval& ival966,
const cxsc::interval& ival967,
const cxsc::interval& ival968,
const cxsc::interval& ival969,
const cxsc::interval& ival970,
const cxsc::interval& ival971,
const cxsc::interval& ival972,
const cxsc::interval& ival973,
const cxsc::interval& ival974,
const cxsc::interval& ival975,
const cxsc::interval& ival976,
const cxsc::interval& ival977,
const cxsc::interval& ival978,
const cxsc::interval& ival979,
const cxsc::interval& ival980,
const cxsc::interval& ival981,
const cxsc::interval& ival982,
const cxsc::interval& ival983,
const cxsc::interval& ival984,
const cxsc::interval& ival985,
const cxsc::interval& ival986,
const cxsc::interval& ival987,
const cxsc::interval& ival988,
const cxsc::interval& ival989,
const cxsc::interval& ival990,
const cxsc::interval& ival991,
const cxsc::interval& ival992,
const cxsc::interval& ival993,
const cxsc::interval& ival994,
const cxsc::interval& ival995,
const cxsc::interval& ival996,
const cxsc::interval& ival997,
const cxsc::interval& ival998,
const cxsc::interval& ival999,
const cxsc::interval& ival1000) const
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
R.push_back(ival101);
R.push_back(ival102);
R.push_back(ival103);
R.push_back(ival104);
R.push_back(ival105);
R.push_back(ival106);
R.push_back(ival107);
R.push_back(ival108);
R.push_back(ival109);
R.push_back(ival110);
R.push_back(ival111);
R.push_back(ival112);
R.push_back(ival113);
R.push_back(ival114);
R.push_back(ival115);
R.push_back(ival116);
R.push_back(ival117);
R.push_back(ival118);
R.push_back(ival119);
R.push_back(ival120);
R.push_back(ival121);
R.push_back(ival122);
R.push_back(ival123);
R.push_back(ival124);
R.push_back(ival125);
R.push_back(ival126);
R.push_back(ival127);
R.push_back(ival128);
R.push_back(ival129);
R.push_back(ival130);
R.push_back(ival131);
R.push_back(ival132);
R.push_back(ival133);
R.push_back(ival134);
R.push_back(ival135);
R.push_back(ival136);
R.push_back(ival137);
R.push_back(ival138);
R.push_back(ival139);
R.push_back(ival140);
R.push_back(ival141);
R.push_back(ival142);
R.push_back(ival143);
R.push_back(ival144);
R.push_back(ival145);
R.push_back(ival146);
R.push_back(ival147);
R.push_back(ival148);
R.push_back(ival149);
R.push_back(ival150);
R.push_back(ival151);
R.push_back(ival152);
R.push_back(ival153);
R.push_back(ival154);
R.push_back(ival155);
R.push_back(ival156);
R.push_back(ival157);
R.push_back(ival158);
R.push_back(ival159);
R.push_back(ival160);
R.push_back(ival161);
R.push_back(ival162);
R.push_back(ival163);
R.push_back(ival164);
R.push_back(ival165);
R.push_back(ival166);
R.push_back(ival167);
R.push_back(ival168);
R.push_back(ival169);
R.push_back(ival170);
R.push_back(ival171);
R.push_back(ival172);
R.push_back(ival173);
R.push_back(ival174);
R.push_back(ival175);
R.push_back(ival176);
R.push_back(ival177);
R.push_back(ival178);
R.push_back(ival179);
R.push_back(ival180);
R.push_back(ival181);
R.push_back(ival182);
R.push_back(ival183);
R.push_back(ival184);
R.push_back(ival185);
R.push_back(ival186);
R.push_back(ival187);
R.push_back(ival188);
R.push_back(ival189);
R.push_back(ival190);
R.push_back(ival191);
R.push_back(ival192);
R.push_back(ival193);
R.push_back(ival194);
R.push_back(ival195);
R.push_back(ival196);
R.push_back(ival197);
R.push_back(ival198);
R.push_back(ival199);
R.push_back(ival200);
R.push_back(ival201);
R.push_back(ival202);
R.push_back(ival203);
R.push_back(ival204);
R.push_back(ival205);
R.push_back(ival206);
R.push_back(ival207);
R.push_back(ival208);
R.push_back(ival209);
R.push_back(ival210);
R.push_back(ival211);
R.push_back(ival212);
R.push_back(ival213);
R.push_back(ival214);
R.push_back(ival215);
R.push_back(ival216);
R.push_back(ival217);
R.push_back(ival218);
R.push_back(ival219);
R.push_back(ival220);
R.push_back(ival221);
R.push_back(ival222);
R.push_back(ival223);
R.push_back(ival224);
R.push_back(ival225);
R.push_back(ival226);
R.push_back(ival227);
R.push_back(ival228);
R.push_back(ival229);
R.push_back(ival230);
R.push_back(ival231);
R.push_back(ival232);
R.push_back(ival233);
R.push_back(ival234);
R.push_back(ival235);
R.push_back(ival236);
R.push_back(ival237);
R.push_back(ival238);
R.push_back(ival239);
R.push_back(ival240);
R.push_back(ival241);
R.push_back(ival242);
R.push_back(ival243);
R.push_back(ival244);
R.push_back(ival245);
R.push_back(ival246);
R.push_back(ival247);
R.push_back(ival248);
R.push_back(ival249);
R.push_back(ival250);
R.push_back(ival251);
R.push_back(ival252);
R.push_back(ival253);
R.push_back(ival254);
R.push_back(ival255);
R.push_back(ival256);
R.push_back(ival257);
R.push_back(ival258);
R.push_back(ival259);
R.push_back(ival260);
R.push_back(ival261);
R.push_back(ival262);
R.push_back(ival263);
R.push_back(ival264);
R.push_back(ival265);
R.push_back(ival266);
R.push_back(ival267);
R.push_back(ival268);
R.push_back(ival269);
R.push_back(ival270);
R.push_back(ival271);
R.push_back(ival272);
R.push_back(ival273);
R.push_back(ival274);
R.push_back(ival275);
R.push_back(ival276);
R.push_back(ival277);
R.push_back(ival278);
R.push_back(ival279);
R.push_back(ival280);
R.push_back(ival281);
R.push_back(ival282);
R.push_back(ival283);
R.push_back(ival284);
R.push_back(ival285);
R.push_back(ival286);
R.push_back(ival287);
R.push_back(ival288);
R.push_back(ival289);
R.push_back(ival290);
R.push_back(ival291);
R.push_back(ival292);
R.push_back(ival293);
R.push_back(ival294);
R.push_back(ival295);
R.push_back(ival296);
R.push_back(ival297);
R.push_back(ival298);
R.push_back(ival299);
R.push_back(ival300);
R.push_back(ival301);
R.push_back(ival302);
R.push_back(ival303);
R.push_back(ival304);
R.push_back(ival305);
R.push_back(ival306);
R.push_back(ival307);
R.push_back(ival308);
R.push_back(ival309);
R.push_back(ival310);
R.push_back(ival311);
R.push_back(ival312);
R.push_back(ival313);
R.push_back(ival314);
R.push_back(ival315);
R.push_back(ival316);
R.push_back(ival317);
R.push_back(ival318);
R.push_back(ival319);
R.push_back(ival320);
R.push_back(ival321);
R.push_back(ival322);
R.push_back(ival323);
R.push_back(ival324);
R.push_back(ival325);
R.push_back(ival326);
R.push_back(ival327);
R.push_back(ival328);
R.push_back(ival329);
R.push_back(ival330);
R.push_back(ival331);
R.push_back(ival332);
R.push_back(ival333);
R.push_back(ival334);
R.push_back(ival335);
R.push_back(ival336);
R.push_back(ival337);
R.push_back(ival338);
R.push_back(ival339);
R.push_back(ival340);
R.push_back(ival341);
R.push_back(ival342);
R.push_back(ival343);
R.push_back(ival344);
R.push_back(ival345);
R.push_back(ival346);
R.push_back(ival347);
R.push_back(ival348);
R.push_back(ival349);
R.push_back(ival350);
R.push_back(ival351);
R.push_back(ival352);
R.push_back(ival353);
R.push_back(ival354);
R.push_back(ival355);
R.push_back(ival356);
R.push_back(ival357);
R.push_back(ival358);
R.push_back(ival359);
R.push_back(ival360);
R.push_back(ival361);
R.push_back(ival362);
R.push_back(ival363);
R.push_back(ival364);
R.push_back(ival365);
R.push_back(ival366);
R.push_back(ival367);
R.push_back(ival368);
R.push_back(ival369);
R.push_back(ival370);
R.push_back(ival371);
R.push_back(ival372);
R.push_back(ival373);
R.push_back(ival374);
R.push_back(ival375);
R.push_back(ival376);
R.push_back(ival377);
R.push_back(ival378);
R.push_back(ival379);
R.push_back(ival380);
R.push_back(ival381);
R.push_back(ival382);
R.push_back(ival383);
R.push_back(ival384);
R.push_back(ival385);
R.push_back(ival386);
R.push_back(ival387);
R.push_back(ival388);
R.push_back(ival389);
R.push_back(ival390);
R.push_back(ival391);
R.push_back(ival392);
R.push_back(ival393);
R.push_back(ival394);
R.push_back(ival395);
R.push_back(ival396);
R.push_back(ival397);
R.push_back(ival398);
R.push_back(ival399);
R.push_back(ival400);
R.push_back(ival401);
R.push_back(ival402);
R.push_back(ival403);
R.push_back(ival404);
R.push_back(ival405);
R.push_back(ival406);
R.push_back(ival407);
R.push_back(ival408);
R.push_back(ival409);
R.push_back(ival410);
R.push_back(ival411);
R.push_back(ival412);
R.push_back(ival413);
R.push_back(ival414);
R.push_back(ival415);
R.push_back(ival416);
R.push_back(ival417);
R.push_back(ival418);
R.push_back(ival419);
R.push_back(ival420);
R.push_back(ival421);
R.push_back(ival422);
R.push_back(ival423);
R.push_back(ival424);
R.push_back(ival425);
R.push_back(ival426);
R.push_back(ival427);
R.push_back(ival428);
R.push_back(ival429);
R.push_back(ival430);
R.push_back(ival431);
R.push_back(ival432);
R.push_back(ival433);
R.push_back(ival434);
R.push_back(ival435);
R.push_back(ival436);
R.push_back(ival437);
R.push_back(ival438);
R.push_back(ival439);
R.push_back(ival440);
R.push_back(ival441);
R.push_back(ival442);
R.push_back(ival443);
R.push_back(ival444);
R.push_back(ival445);
R.push_back(ival446);
R.push_back(ival447);
R.push_back(ival448);
R.push_back(ival449);
R.push_back(ival450);
R.push_back(ival451);
R.push_back(ival452);
R.push_back(ival453);
R.push_back(ival454);
R.push_back(ival455);
R.push_back(ival456);
R.push_back(ival457);
R.push_back(ival458);
R.push_back(ival459);
R.push_back(ival460);
R.push_back(ival461);
R.push_back(ival462);
R.push_back(ival463);
R.push_back(ival464);
R.push_back(ival465);
R.push_back(ival466);
R.push_back(ival467);
R.push_back(ival468);
R.push_back(ival469);
R.push_back(ival470);
R.push_back(ival471);
R.push_back(ival472);
R.push_back(ival473);
R.push_back(ival474);
R.push_back(ival475);
R.push_back(ival476);
R.push_back(ival477);
R.push_back(ival478);
R.push_back(ival479);
R.push_back(ival480);
R.push_back(ival481);
R.push_back(ival482);
R.push_back(ival483);
R.push_back(ival484);
R.push_back(ival485);
R.push_back(ival486);
R.push_back(ival487);
R.push_back(ival488);
R.push_back(ival489);
R.push_back(ival490);
R.push_back(ival491);
R.push_back(ival492);
R.push_back(ival493);
R.push_back(ival494);
R.push_back(ival495);
R.push_back(ival496);
R.push_back(ival497);
R.push_back(ival498);
R.push_back(ival499);
R.push_back(ival500);
R.push_back(ival501);
R.push_back(ival502);
R.push_back(ival503);
R.push_back(ival504);
R.push_back(ival505);
R.push_back(ival506);
R.push_back(ival507);
R.push_back(ival508);
R.push_back(ival509);
R.push_back(ival510);
R.push_back(ival511);
R.push_back(ival512);
R.push_back(ival513);
R.push_back(ival514);
R.push_back(ival515);
R.push_back(ival516);
R.push_back(ival517);
R.push_back(ival518);
R.push_back(ival519);
R.push_back(ival520);
R.push_back(ival521);
R.push_back(ival522);
R.push_back(ival523);
R.push_back(ival524);
R.push_back(ival525);
R.push_back(ival526);
R.push_back(ival527);
R.push_back(ival528);
R.push_back(ival529);
R.push_back(ival530);
R.push_back(ival531);
R.push_back(ival532);
R.push_back(ival533);
R.push_back(ival534);
R.push_back(ival535);
R.push_back(ival536);
R.push_back(ival537);
R.push_back(ival538);
R.push_back(ival539);
R.push_back(ival540);
R.push_back(ival541);
R.push_back(ival542);
R.push_back(ival543);
R.push_back(ival544);
R.push_back(ival545);
R.push_back(ival546);
R.push_back(ival547);
R.push_back(ival548);
R.push_back(ival549);
R.push_back(ival550);
R.push_back(ival551);
R.push_back(ival552);
R.push_back(ival553);
R.push_back(ival554);
R.push_back(ival555);
R.push_back(ival556);
R.push_back(ival557);
R.push_back(ival558);
R.push_back(ival559);
R.push_back(ival560);
R.push_back(ival561);
R.push_back(ival562);
R.push_back(ival563);
R.push_back(ival564);
R.push_back(ival565);
R.push_back(ival566);
R.push_back(ival567);
R.push_back(ival568);
R.push_back(ival569);
R.push_back(ival570);
R.push_back(ival571);
R.push_back(ival572);
R.push_back(ival573);
R.push_back(ival574);
R.push_back(ival575);
R.push_back(ival576);
R.push_back(ival577);
R.push_back(ival578);
R.push_back(ival579);
R.push_back(ival580);
R.push_back(ival581);
R.push_back(ival582);
R.push_back(ival583);
R.push_back(ival584);
R.push_back(ival585);
R.push_back(ival586);
R.push_back(ival587);
R.push_back(ival588);
R.push_back(ival589);
R.push_back(ival590);
R.push_back(ival591);
R.push_back(ival592);
R.push_back(ival593);
R.push_back(ival594);
R.push_back(ival595);
R.push_back(ival596);
R.push_back(ival597);
R.push_back(ival598);
R.push_back(ival599);
R.push_back(ival600);
R.push_back(ival601);
R.push_back(ival602);
R.push_back(ival603);
R.push_back(ival604);
R.push_back(ival605);
R.push_back(ival606);
R.push_back(ival607);
R.push_back(ival608);
R.push_back(ival609);
R.push_back(ival610);
R.push_back(ival611);
R.push_back(ival612);
R.push_back(ival613);
R.push_back(ival614);
R.push_back(ival615);
R.push_back(ival616);
R.push_back(ival617);
R.push_back(ival618);
R.push_back(ival619);
R.push_back(ival620);
R.push_back(ival621);
R.push_back(ival622);
R.push_back(ival623);
R.push_back(ival624);
R.push_back(ival625);
R.push_back(ival626);
R.push_back(ival627);
R.push_back(ival628);
R.push_back(ival629);
R.push_back(ival630);
R.push_back(ival631);
R.push_back(ival632);
R.push_back(ival633);
R.push_back(ival634);
R.push_back(ival635);
R.push_back(ival636);
R.push_back(ival637);
R.push_back(ival638);
R.push_back(ival639);
R.push_back(ival640);
R.push_back(ival641);
R.push_back(ival642);
R.push_back(ival643);
R.push_back(ival644);
R.push_back(ival645);
R.push_back(ival646);
R.push_back(ival647);
R.push_back(ival648);
R.push_back(ival649);
R.push_back(ival650);
R.push_back(ival651);
R.push_back(ival652);
R.push_back(ival653);
R.push_back(ival654);
R.push_back(ival655);
R.push_back(ival656);
R.push_back(ival657);
R.push_back(ival658);
R.push_back(ival659);
R.push_back(ival660);
R.push_back(ival661);
R.push_back(ival662);
R.push_back(ival663);
R.push_back(ival664);
R.push_back(ival665);
R.push_back(ival666);
R.push_back(ival667);
R.push_back(ival668);
R.push_back(ival669);
R.push_back(ival670);
R.push_back(ival671);
R.push_back(ival672);
R.push_back(ival673);
R.push_back(ival674);
R.push_back(ival675);
R.push_back(ival676);
R.push_back(ival677);
R.push_back(ival678);
R.push_back(ival679);
R.push_back(ival680);
R.push_back(ival681);
R.push_back(ival682);
R.push_back(ival683);
R.push_back(ival684);
R.push_back(ival685);
R.push_back(ival686);
R.push_back(ival687);
R.push_back(ival688);
R.push_back(ival689);
R.push_back(ival690);
R.push_back(ival691);
R.push_back(ival692);
R.push_back(ival693);
R.push_back(ival694);
R.push_back(ival695);
R.push_back(ival696);
R.push_back(ival697);
R.push_back(ival698);
R.push_back(ival699);
R.push_back(ival700);
R.push_back(ival701);
R.push_back(ival702);
R.push_back(ival703);
R.push_back(ival704);
R.push_back(ival705);
R.push_back(ival706);
R.push_back(ival707);
R.push_back(ival708);
R.push_back(ival709);
R.push_back(ival710);
R.push_back(ival711);
R.push_back(ival712);
R.push_back(ival713);
R.push_back(ival714);
R.push_back(ival715);
R.push_back(ival716);
R.push_back(ival717);
R.push_back(ival718);
R.push_back(ival719);
R.push_back(ival720);
R.push_back(ival721);
R.push_back(ival722);
R.push_back(ival723);
R.push_back(ival724);
R.push_back(ival725);
R.push_back(ival726);
R.push_back(ival727);
R.push_back(ival728);
R.push_back(ival729);
R.push_back(ival730);
R.push_back(ival731);
R.push_back(ival732);
R.push_back(ival733);
R.push_back(ival734);
R.push_back(ival735);
R.push_back(ival736);
R.push_back(ival737);
R.push_back(ival738);
R.push_back(ival739);
R.push_back(ival740);
R.push_back(ival741);
R.push_back(ival742);
R.push_back(ival743);
R.push_back(ival744);
R.push_back(ival745);
R.push_back(ival746);
R.push_back(ival747);
R.push_back(ival748);
R.push_back(ival749);
R.push_back(ival750);
R.push_back(ival751);
R.push_back(ival752);
R.push_back(ival753);
R.push_back(ival754);
R.push_back(ival755);
R.push_back(ival756);
R.push_back(ival757);
R.push_back(ival758);
R.push_back(ival759);
R.push_back(ival760);
R.push_back(ival761);
R.push_back(ival762);
R.push_back(ival763);
R.push_back(ival764);
R.push_back(ival765);
R.push_back(ival766);
R.push_back(ival767);
R.push_back(ival768);
R.push_back(ival769);
R.push_back(ival770);
R.push_back(ival771);
R.push_back(ival772);
R.push_back(ival773);
R.push_back(ival774);
R.push_back(ival775);
R.push_back(ival776);
R.push_back(ival777);
R.push_back(ival778);
R.push_back(ival779);
R.push_back(ival780);
R.push_back(ival781);
R.push_back(ival782);
R.push_back(ival783);
R.push_back(ival784);
R.push_back(ival785);
R.push_back(ival786);
R.push_back(ival787);
R.push_back(ival788);
R.push_back(ival789);
R.push_back(ival790);
R.push_back(ival791);
R.push_back(ival792);
R.push_back(ival793);
R.push_back(ival794);
R.push_back(ival795);
R.push_back(ival796);
R.push_back(ival797);
R.push_back(ival798);
R.push_back(ival799);
R.push_back(ival800);
R.push_back(ival801);
R.push_back(ival802);
R.push_back(ival803);
R.push_back(ival804);
R.push_back(ival805);
R.push_back(ival806);
R.push_back(ival807);
R.push_back(ival808);
R.push_back(ival809);
R.push_back(ival810);
R.push_back(ival811);
R.push_back(ival812);
R.push_back(ival813);
R.push_back(ival814);
R.push_back(ival815);
R.push_back(ival816);
R.push_back(ival817);
R.push_back(ival818);
R.push_back(ival819);
R.push_back(ival820);
R.push_back(ival821);
R.push_back(ival822);
R.push_back(ival823);
R.push_back(ival824);
R.push_back(ival825);
R.push_back(ival826);
R.push_back(ival827);
R.push_back(ival828);
R.push_back(ival829);
R.push_back(ival830);
R.push_back(ival831);
R.push_back(ival832);
R.push_back(ival833);
R.push_back(ival834);
R.push_back(ival835);
R.push_back(ival836);
R.push_back(ival837);
R.push_back(ival838);
R.push_back(ival839);
R.push_back(ival840);
R.push_back(ival841);
R.push_back(ival842);
R.push_back(ival843);
R.push_back(ival844);
R.push_back(ival845);
R.push_back(ival846);
R.push_back(ival847);
R.push_back(ival848);
R.push_back(ival849);
R.push_back(ival850);
R.push_back(ival851);
R.push_back(ival852);
R.push_back(ival853);
R.push_back(ival854);
R.push_back(ival855);
R.push_back(ival856);
R.push_back(ival857);
R.push_back(ival858);
R.push_back(ival859);
R.push_back(ival860);
R.push_back(ival861);
R.push_back(ival862);
R.push_back(ival863);
R.push_back(ival864);
R.push_back(ival865);
R.push_back(ival866);
R.push_back(ival867);
R.push_back(ival868);
R.push_back(ival869);
R.push_back(ival870);
R.push_back(ival871);
R.push_back(ival872);
R.push_back(ival873);
R.push_back(ival874);
R.push_back(ival875);
R.push_back(ival876);
R.push_back(ival877);
R.push_back(ival878);
R.push_back(ival879);
R.push_back(ival880);
R.push_back(ival881);
R.push_back(ival882);
R.push_back(ival883);
R.push_back(ival884);
R.push_back(ival885);
R.push_back(ival886);
R.push_back(ival887);
R.push_back(ival888);
R.push_back(ival889);
R.push_back(ival890);
R.push_back(ival891);
R.push_back(ival892);
R.push_back(ival893);
R.push_back(ival894);
R.push_back(ival895);
R.push_back(ival896);
R.push_back(ival897);
R.push_back(ival898);
R.push_back(ival899);
R.push_back(ival900);
R.push_back(ival901);
R.push_back(ival902);
R.push_back(ival903);
R.push_back(ival904);
R.push_back(ival905);
R.push_back(ival906);
R.push_back(ival907);
R.push_back(ival908);
R.push_back(ival909);
R.push_back(ival910);
R.push_back(ival911);
R.push_back(ival912);
R.push_back(ival913);
R.push_back(ival914);
R.push_back(ival915);
R.push_back(ival916);
R.push_back(ival917);
R.push_back(ival918);
R.push_back(ival919);
R.push_back(ival920);
R.push_back(ival921);
R.push_back(ival922);
R.push_back(ival923);
R.push_back(ival924);
R.push_back(ival925);
R.push_back(ival926);
R.push_back(ival927);
R.push_back(ival928);
R.push_back(ival929);
R.push_back(ival930);
R.push_back(ival931);
R.push_back(ival932);
R.push_back(ival933);
R.push_back(ival934);
R.push_back(ival935);
R.push_back(ival936);
R.push_back(ival937);
R.push_back(ival938);
R.push_back(ival939);
R.push_back(ival940);
R.push_back(ival941);
R.push_back(ival942);
R.push_back(ival943);
R.push_back(ival944);
R.push_back(ival945);
R.push_back(ival946);
R.push_back(ival947);
R.push_back(ival948);
R.push_back(ival949);
R.push_back(ival950);
R.push_back(ival951);
R.push_back(ival952);
R.push_back(ival953);
R.push_back(ival954);
R.push_back(ival955);
R.push_back(ival956);
R.push_back(ival957);
R.push_back(ival958);
R.push_back(ival959);
R.push_back(ival960);
R.push_back(ival961);
R.push_back(ival962);
R.push_back(ival963);
R.push_back(ival964);
R.push_back(ival965);
R.push_back(ival966);
R.push_back(ival967);
R.push_back(ival968);
R.push_back(ival969);
R.push_back(ival970);
R.push_back(ival971);
R.push_back(ival972);
R.push_back(ival973);
R.push_back(ival974);
R.push_back(ival975);
R.push_back(ival976);
R.push_back(ival977);
R.push_back(ival978);
R.push_back(ival979);
R.push_back(ival980);
R.push_back(ival981);
R.push_back(ival982);
R.push_back(ival983);
R.push_back(ival984);
R.push_back(ival985);
R.push_back(ival986);
R.push_back(ival987);
R.push_back(ival988);
R.push_back(ival989);
R.push_back(ival990);
R.push_back(ival991);
R.push_back(ival992);
R.push_back(ival993);
R.push_back(ival994);
R.push_back(ival995);
R.push_back(ival996);
R.push_back(ival997);
R.push_back(ival998);
R.push_back(ival999);
R.push_back(ival1000);

    interval result(0,0);
    
    for (size_t i = 0; i < (R.size()-1); i++) 
    {
      result += 100.0 * (power((R[i+1] - power(R[i],2)),2)) +
									power((R[i] - 1), 2);
    }
  
	result = exp (-(1.0 * result));
	
	return result;
	
}

real RosenFobj1000D::operator()(const cxsc::real& r1,
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
const cxsc::real& r100,
const cxsc::real& r101,
const cxsc::real& r102,
const cxsc::real& r103,
const cxsc::real& r104,
const cxsc::real& r105,
const cxsc::real& r106,
const cxsc::real& r107,
const cxsc::real& r108,
const cxsc::real& r109,
const cxsc::real& r110,
const cxsc::real& r111,
const cxsc::real& r112,
const cxsc::real& r113,
const cxsc::real& r114,
const cxsc::real& r115,
const cxsc::real& r116,
const cxsc::real& r117,
const cxsc::real& r118,
const cxsc::real& r119,
const cxsc::real& r120,
const cxsc::real& r121,
const cxsc::real& r122,
const cxsc::real& r123,
const cxsc::real& r124,
const cxsc::real& r125,
const cxsc::real& r126,
const cxsc::real& r127,
const cxsc::real& r128,
const cxsc::real& r129,
const cxsc::real& r130,
const cxsc::real& r131,
const cxsc::real& r132,
const cxsc::real& r133,
const cxsc::real& r134,
const cxsc::real& r135,
const cxsc::real& r136,
const cxsc::real& r137,
const cxsc::real& r138,
const cxsc::real& r139,
const cxsc::real& r140,
const cxsc::real& r141,
const cxsc::real& r142,
const cxsc::real& r143,
const cxsc::real& r144,
const cxsc::real& r145,
const cxsc::real& r146,
const cxsc::real& r147,
const cxsc::real& r148,
const cxsc::real& r149,
const cxsc::real& r150,
const cxsc::real& r151,
const cxsc::real& r152,
const cxsc::real& r153,
const cxsc::real& r154,
const cxsc::real& r155,
const cxsc::real& r156,
const cxsc::real& r157,
const cxsc::real& r158,
const cxsc::real& r159,
const cxsc::real& r160,
const cxsc::real& r161,
const cxsc::real& r162,
const cxsc::real& r163,
const cxsc::real& r164,
const cxsc::real& r165,
const cxsc::real& r166,
const cxsc::real& r167,
const cxsc::real& r168,
const cxsc::real& r169,
const cxsc::real& r170,
const cxsc::real& r171,
const cxsc::real& r172,
const cxsc::real& r173,
const cxsc::real& r174,
const cxsc::real& r175,
const cxsc::real& r176,
const cxsc::real& r177,
const cxsc::real& r178,
const cxsc::real& r179,
const cxsc::real& r180,
const cxsc::real& r181,
const cxsc::real& r182,
const cxsc::real& r183,
const cxsc::real& r184,
const cxsc::real& r185,
const cxsc::real& r186,
const cxsc::real& r187,
const cxsc::real& r188,
const cxsc::real& r189,
const cxsc::real& r190,
const cxsc::real& r191,
const cxsc::real& r192,
const cxsc::real& r193,
const cxsc::real& r194,
const cxsc::real& r195,
const cxsc::real& r196,
const cxsc::real& r197,
const cxsc::real& r198,
const cxsc::real& r199,
const cxsc::real& r200,
const cxsc::real& r201,
const cxsc::real& r202,
const cxsc::real& r203,
const cxsc::real& r204,
const cxsc::real& r205,
const cxsc::real& r206,
const cxsc::real& r207,
const cxsc::real& r208,
const cxsc::real& r209,
const cxsc::real& r210,
const cxsc::real& r211,
const cxsc::real& r212,
const cxsc::real& r213,
const cxsc::real& r214,
const cxsc::real& r215,
const cxsc::real& r216,
const cxsc::real& r217,
const cxsc::real& r218,
const cxsc::real& r219,
const cxsc::real& r220,
const cxsc::real& r221,
const cxsc::real& r222,
const cxsc::real& r223,
const cxsc::real& r224,
const cxsc::real& r225,
const cxsc::real& r226,
const cxsc::real& r227,
const cxsc::real& r228,
const cxsc::real& r229,
const cxsc::real& r230,
const cxsc::real& r231,
const cxsc::real& r232,
const cxsc::real& r233,
const cxsc::real& r234,
const cxsc::real& r235,
const cxsc::real& r236,
const cxsc::real& r237,
const cxsc::real& r238,
const cxsc::real& r239,
const cxsc::real& r240,
const cxsc::real& r241,
const cxsc::real& r242,
const cxsc::real& r243,
const cxsc::real& r244,
const cxsc::real& r245,
const cxsc::real& r246,
const cxsc::real& r247,
const cxsc::real& r248,
const cxsc::real& r249,
const cxsc::real& r250,
const cxsc::real& r251,
const cxsc::real& r252,
const cxsc::real& r253,
const cxsc::real& r254,
const cxsc::real& r255,
const cxsc::real& r256,
const cxsc::real& r257,
const cxsc::real& r258,
const cxsc::real& r259,
const cxsc::real& r260,
const cxsc::real& r261,
const cxsc::real& r262,
const cxsc::real& r263,
const cxsc::real& r264,
const cxsc::real& r265,
const cxsc::real& r266,
const cxsc::real& r267,
const cxsc::real& r268,
const cxsc::real& r269,
const cxsc::real& r270,
const cxsc::real& r271,
const cxsc::real& r272,
const cxsc::real& r273,
const cxsc::real& r274,
const cxsc::real& r275,
const cxsc::real& r276,
const cxsc::real& r277,
const cxsc::real& r278,
const cxsc::real& r279,
const cxsc::real& r280,
const cxsc::real& r281,
const cxsc::real& r282,
const cxsc::real& r283,
const cxsc::real& r284,
const cxsc::real& r285,
const cxsc::real& r286,
const cxsc::real& r287,
const cxsc::real& r288,
const cxsc::real& r289,
const cxsc::real& r290,
const cxsc::real& r291,
const cxsc::real& r292,
const cxsc::real& r293,
const cxsc::real& r294,
const cxsc::real& r295,
const cxsc::real& r296,
const cxsc::real& r297,
const cxsc::real& r298,
const cxsc::real& r299,
const cxsc::real& r300,
const cxsc::real& r301,
const cxsc::real& r302,
const cxsc::real& r303,
const cxsc::real& r304,
const cxsc::real& r305,
const cxsc::real& r306,
const cxsc::real& r307,
const cxsc::real& r308,
const cxsc::real& r309,
const cxsc::real& r310,
const cxsc::real& r311,
const cxsc::real& r312,
const cxsc::real& r313,
const cxsc::real& r314,
const cxsc::real& r315,
const cxsc::real& r316,
const cxsc::real& r317,
const cxsc::real& r318,
const cxsc::real& r319,
const cxsc::real& r320,
const cxsc::real& r321,
const cxsc::real& r322,
const cxsc::real& r323,
const cxsc::real& r324,
const cxsc::real& r325,
const cxsc::real& r326,
const cxsc::real& r327,
const cxsc::real& r328,
const cxsc::real& r329,
const cxsc::real& r330,
const cxsc::real& r331,
const cxsc::real& r332,
const cxsc::real& r333,
const cxsc::real& r334,
const cxsc::real& r335,
const cxsc::real& r336,
const cxsc::real& r337,
const cxsc::real& r338,
const cxsc::real& r339,
const cxsc::real& r340,
const cxsc::real& r341,
const cxsc::real& r342,
const cxsc::real& r343,
const cxsc::real& r344,
const cxsc::real& r345,
const cxsc::real& r346,
const cxsc::real& r347,
const cxsc::real& r348,
const cxsc::real& r349,
const cxsc::real& r350,
const cxsc::real& r351,
const cxsc::real& r352,
const cxsc::real& r353,
const cxsc::real& r354,
const cxsc::real& r355,
const cxsc::real& r356,
const cxsc::real& r357,
const cxsc::real& r358,
const cxsc::real& r359,
const cxsc::real& r360,
const cxsc::real& r361,
const cxsc::real& r362,
const cxsc::real& r363,
const cxsc::real& r364,
const cxsc::real& r365,
const cxsc::real& r366,
const cxsc::real& r367,
const cxsc::real& r368,
const cxsc::real& r369,
const cxsc::real& r370,
const cxsc::real& r371,
const cxsc::real& r372,
const cxsc::real& r373,
const cxsc::real& r374,
const cxsc::real& r375,
const cxsc::real& r376,
const cxsc::real& r377,
const cxsc::real& r378,
const cxsc::real& r379,
const cxsc::real& r380,
const cxsc::real& r381,
const cxsc::real& r382,
const cxsc::real& r383,
const cxsc::real& r384,
const cxsc::real& r385,
const cxsc::real& r386,
const cxsc::real& r387,
const cxsc::real& r388,
const cxsc::real& r389,
const cxsc::real& r390,
const cxsc::real& r391,
const cxsc::real& r392,
const cxsc::real& r393,
const cxsc::real& r394,
const cxsc::real& r395,
const cxsc::real& r396,
const cxsc::real& r397,
const cxsc::real& r398,
const cxsc::real& r399,
const cxsc::real& r400,
const cxsc::real& r401,
const cxsc::real& r402,
const cxsc::real& r403,
const cxsc::real& r404,
const cxsc::real& r405,
const cxsc::real& r406,
const cxsc::real& r407,
const cxsc::real& r408,
const cxsc::real& r409,
const cxsc::real& r410,
const cxsc::real& r411,
const cxsc::real& r412,
const cxsc::real& r413,
const cxsc::real& r414,
const cxsc::real& r415,
const cxsc::real& r416,
const cxsc::real& r417,
const cxsc::real& r418,
const cxsc::real& r419,
const cxsc::real& r420,
const cxsc::real& r421,
const cxsc::real& r422,
const cxsc::real& r423,
const cxsc::real& r424,
const cxsc::real& r425,
const cxsc::real& r426,
const cxsc::real& r427,
const cxsc::real& r428,
const cxsc::real& r429,
const cxsc::real& r430,
const cxsc::real& r431,
const cxsc::real& r432,
const cxsc::real& r433,
const cxsc::real& r434,
const cxsc::real& r435,
const cxsc::real& r436,
const cxsc::real& r437,
const cxsc::real& r438,
const cxsc::real& r439,
const cxsc::real& r440,
const cxsc::real& r441,
const cxsc::real& r442,
const cxsc::real& r443,
const cxsc::real& r444,
const cxsc::real& r445,
const cxsc::real& r446,
const cxsc::real& r447,
const cxsc::real& r448,
const cxsc::real& r449,
const cxsc::real& r450,
const cxsc::real& r451,
const cxsc::real& r452,
const cxsc::real& r453,
const cxsc::real& r454,
const cxsc::real& r455,
const cxsc::real& r456,
const cxsc::real& r457,
const cxsc::real& r458,
const cxsc::real& r459,
const cxsc::real& r460,
const cxsc::real& r461,
const cxsc::real& r462,
const cxsc::real& r463,
const cxsc::real& r464,
const cxsc::real& r465,
const cxsc::real& r466,
const cxsc::real& r467,
const cxsc::real& r468,
const cxsc::real& r469,
const cxsc::real& r470,
const cxsc::real& r471,
const cxsc::real& r472,
const cxsc::real& r473,
const cxsc::real& r474,
const cxsc::real& r475,
const cxsc::real& r476,
const cxsc::real& r477,
const cxsc::real& r478,
const cxsc::real& r479,
const cxsc::real& r480,
const cxsc::real& r481,
const cxsc::real& r482,
const cxsc::real& r483,
const cxsc::real& r484,
const cxsc::real& r485,
const cxsc::real& r486,
const cxsc::real& r487,
const cxsc::real& r488,
const cxsc::real& r489,
const cxsc::real& r490,
const cxsc::real& r491,
const cxsc::real& r492,
const cxsc::real& r493,
const cxsc::real& r494,
const cxsc::real& r495,
const cxsc::real& r496,
const cxsc::real& r497,
const cxsc::real& r498,
const cxsc::real& r499,
const cxsc::real& r500,
const cxsc::real& r501,
const cxsc::real& r502,
const cxsc::real& r503,
const cxsc::real& r504,
const cxsc::real& r505,
const cxsc::real& r506,
const cxsc::real& r507,
const cxsc::real& r508,
const cxsc::real& r509,
const cxsc::real& r510,
const cxsc::real& r511,
const cxsc::real& r512,
const cxsc::real& r513,
const cxsc::real& r514,
const cxsc::real& r515,
const cxsc::real& r516,
const cxsc::real& r517,
const cxsc::real& r518,
const cxsc::real& r519,
const cxsc::real& r520,
const cxsc::real& r521,
const cxsc::real& r522,
const cxsc::real& r523,
const cxsc::real& r524,
const cxsc::real& r525,
const cxsc::real& r526,
const cxsc::real& r527,
const cxsc::real& r528,
const cxsc::real& r529,
const cxsc::real& r530,
const cxsc::real& r531,
const cxsc::real& r532,
const cxsc::real& r533,
const cxsc::real& r534,
const cxsc::real& r535,
const cxsc::real& r536,
const cxsc::real& r537,
const cxsc::real& r538,
const cxsc::real& r539,
const cxsc::real& r540,
const cxsc::real& r541,
const cxsc::real& r542,
const cxsc::real& r543,
const cxsc::real& r544,
const cxsc::real& r545,
const cxsc::real& r546,
const cxsc::real& r547,
const cxsc::real& r548,
const cxsc::real& r549,
const cxsc::real& r550,
const cxsc::real& r551,
const cxsc::real& r552,
const cxsc::real& r553,
const cxsc::real& r554,
const cxsc::real& r555,
const cxsc::real& r556,
const cxsc::real& r557,
const cxsc::real& r558,
const cxsc::real& r559,
const cxsc::real& r560,
const cxsc::real& r561,
const cxsc::real& r562,
const cxsc::real& r563,
const cxsc::real& r564,
const cxsc::real& r565,
const cxsc::real& r566,
const cxsc::real& r567,
const cxsc::real& r568,
const cxsc::real& r569,
const cxsc::real& r570,
const cxsc::real& r571,
const cxsc::real& r572,
const cxsc::real& r573,
const cxsc::real& r574,
const cxsc::real& r575,
const cxsc::real& r576,
const cxsc::real& r577,
const cxsc::real& r578,
const cxsc::real& r579,
const cxsc::real& r580,
const cxsc::real& r581,
const cxsc::real& r582,
const cxsc::real& r583,
const cxsc::real& r584,
const cxsc::real& r585,
const cxsc::real& r586,
const cxsc::real& r587,
const cxsc::real& r588,
const cxsc::real& r589,
const cxsc::real& r590,
const cxsc::real& r591,
const cxsc::real& r592,
const cxsc::real& r593,
const cxsc::real& r594,
const cxsc::real& r595,
const cxsc::real& r596,
const cxsc::real& r597,
const cxsc::real& r598,
const cxsc::real& r599,
const cxsc::real& r600,
const cxsc::real& r601,
const cxsc::real& r602,
const cxsc::real& r603,
const cxsc::real& r604,
const cxsc::real& r605,
const cxsc::real& r606,
const cxsc::real& r607,
const cxsc::real& r608,
const cxsc::real& r609,
const cxsc::real& r610,
const cxsc::real& r611,
const cxsc::real& r612,
const cxsc::real& r613,
const cxsc::real& r614,
const cxsc::real& r615,
const cxsc::real& r616,
const cxsc::real& r617,
const cxsc::real& r618,
const cxsc::real& r619,
const cxsc::real& r620,
const cxsc::real& r621,
const cxsc::real& r622,
const cxsc::real& r623,
const cxsc::real& r624,
const cxsc::real& r625,
const cxsc::real& r626,
const cxsc::real& r627,
const cxsc::real& r628,
const cxsc::real& r629,
const cxsc::real& r630,
const cxsc::real& r631,
const cxsc::real& r632,
const cxsc::real& r633,
const cxsc::real& r634,
const cxsc::real& r635,
const cxsc::real& r636,
const cxsc::real& r637,
const cxsc::real& r638,
const cxsc::real& r639,
const cxsc::real& r640,
const cxsc::real& r641,
const cxsc::real& r642,
const cxsc::real& r643,
const cxsc::real& r644,
const cxsc::real& r645,
const cxsc::real& r646,
const cxsc::real& r647,
const cxsc::real& r648,
const cxsc::real& r649,
const cxsc::real& r650,
const cxsc::real& r651,
const cxsc::real& r652,
const cxsc::real& r653,
const cxsc::real& r654,
const cxsc::real& r655,
const cxsc::real& r656,
const cxsc::real& r657,
const cxsc::real& r658,
const cxsc::real& r659,
const cxsc::real& r660,
const cxsc::real& r661,
const cxsc::real& r662,
const cxsc::real& r663,
const cxsc::real& r664,
const cxsc::real& r665,
const cxsc::real& r666,
const cxsc::real& r667,
const cxsc::real& r668,
const cxsc::real& r669,
const cxsc::real& r670,
const cxsc::real& r671,
const cxsc::real& r672,
const cxsc::real& r673,
const cxsc::real& r674,
const cxsc::real& r675,
const cxsc::real& r676,
const cxsc::real& r677,
const cxsc::real& r678,
const cxsc::real& r679,
const cxsc::real& r680,
const cxsc::real& r681,
const cxsc::real& r682,
const cxsc::real& r683,
const cxsc::real& r684,
const cxsc::real& r685,
const cxsc::real& r686,
const cxsc::real& r687,
const cxsc::real& r688,
const cxsc::real& r689,
const cxsc::real& r690,
const cxsc::real& r691,
const cxsc::real& r692,
const cxsc::real& r693,
const cxsc::real& r694,
const cxsc::real& r695,
const cxsc::real& r696,
const cxsc::real& r697,
const cxsc::real& r698,
const cxsc::real& r699,
const cxsc::real& r700,
const cxsc::real& r701,
const cxsc::real& r702,
const cxsc::real& r703,
const cxsc::real& r704,
const cxsc::real& r705,
const cxsc::real& r706,
const cxsc::real& r707,
const cxsc::real& r708,
const cxsc::real& r709,
const cxsc::real& r710,
const cxsc::real& r711,
const cxsc::real& r712,
const cxsc::real& r713,
const cxsc::real& r714,
const cxsc::real& r715,
const cxsc::real& r716,
const cxsc::real& r717,
const cxsc::real& r718,
const cxsc::real& r719,
const cxsc::real& r720,
const cxsc::real& r721,
const cxsc::real& r722,
const cxsc::real& r723,
const cxsc::real& r724,
const cxsc::real& r725,
const cxsc::real& r726,
const cxsc::real& r727,
const cxsc::real& r728,
const cxsc::real& r729,
const cxsc::real& r730,
const cxsc::real& r731,
const cxsc::real& r732,
const cxsc::real& r733,
const cxsc::real& r734,
const cxsc::real& r735,
const cxsc::real& r736,
const cxsc::real& r737,
const cxsc::real& r738,
const cxsc::real& r739,
const cxsc::real& r740,
const cxsc::real& r741,
const cxsc::real& r742,
const cxsc::real& r743,
const cxsc::real& r744,
const cxsc::real& r745,
const cxsc::real& r746,
const cxsc::real& r747,
const cxsc::real& r748,
const cxsc::real& r749,
const cxsc::real& r750,
const cxsc::real& r751,
const cxsc::real& r752,
const cxsc::real& r753,
const cxsc::real& r754,
const cxsc::real& r755,
const cxsc::real& r756,
const cxsc::real& r757,
const cxsc::real& r758,
const cxsc::real& r759,
const cxsc::real& r760,
const cxsc::real& r761,
const cxsc::real& r762,
const cxsc::real& r763,
const cxsc::real& r764,
const cxsc::real& r765,
const cxsc::real& r766,
const cxsc::real& r767,
const cxsc::real& r768,
const cxsc::real& r769,
const cxsc::real& r770,
const cxsc::real& r771,
const cxsc::real& r772,
const cxsc::real& r773,
const cxsc::real& r774,
const cxsc::real& r775,
const cxsc::real& r776,
const cxsc::real& r777,
const cxsc::real& r778,
const cxsc::real& r779,
const cxsc::real& r780,
const cxsc::real& r781,
const cxsc::real& r782,
const cxsc::real& r783,
const cxsc::real& r784,
const cxsc::real& r785,
const cxsc::real& r786,
const cxsc::real& r787,
const cxsc::real& r788,
const cxsc::real& r789,
const cxsc::real& r790,
const cxsc::real& r791,
const cxsc::real& r792,
const cxsc::real& r793,
const cxsc::real& r794,
const cxsc::real& r795,
const cxsc::real& r796,
const cxsc::real& r797,
const cxsc::real& r798,
const cxsc::real& r799,
const cxsc::real& r800,
const cxsc::real& r801,
const cxsc::real& r802,
const cxsc::real& r803,
const cxsc::real& r804,
const cxsc::real& r805,
const cxsc::real& r806,
const cxsc::real& r807,
const cxsc::real& r808,
const cxsc::real& r809,
const cxsc::real& r810,
const cxsc::real& r811,
const cxsc::real& r812,
const cxsc::real& r813,
const cxsc::real& r814,
const cxsc::real& r815,
const cxsc::real& r816,
const cxsc::real& r817,
const cxsc::real& r818,
const cxsc::real& r819,
const cxsc::real& r820,
const cxsc::real& r821,
const cxsc::real& r822,
const cxsc::real& r823,
const cxsc::real& r824,
const cxsc::real& r825,
const cxsc::real& r826,
const cxsc::real& r827,
const cxsc::real& r828,
const cxsc::real& r829,
const cxsc::real& r830,
const cxsc::real& r831,
const cxsc::real& r832,
const cxsc::real& r833,
const cxsc::real& r834,
const cxsc::real& r835,
const cxsc::real& r836,
const cxsc::real& r837,
const cxsc::real& r838,
const cxsc::real& r839,
const cxsc::real& r840,
const cxsc::real& r841,
const cxsc::real& r842,
const cxsc::real& r843,
const cxsc::real& r844,
const cxsc::real& r845,
const cxsc::real& r846,
const cxsc::real& r847,
const cxsc::real& r848,
const cxsc::real& r849,
const cxsc::real& r850,
const cxsc::real& r851,
const cxsc::real& r852,
const cxsc::real& r853,
const cxsc::real& r854,
const cxsc::real& r855,
const cxsc::real& r856,
const cxsc::real& r857,
const cxsc::real& r858,
const cxsc::real& r859,
const cxsc::real& r860,
const cxsc::real& r861,
const cxsc::real& r862,
const cxsc::real& r863,
const cxsc::real& r864,
const cxsc::real& r865,
const cxsc::real& r866,
const cxsc::real& r867,
const cxsc::real& r868,
const cxsc::real& r869,
const cxsc::real& r870,
const cxsc::real& r871,
const cxsc::real& r872,
const cxsc::real& r873,
const cxsc::real& r874,
const cxsc::real& r875,
const cxsc::real& r876,
const cxsc::real& r877,
const cxsc::real& r878,
const cxsc::real& r879,
const cxsc::real& r880,
const cxsc::real& r881,
const cxsc::real& r882,
const cxsc::real& r883,
const cxsc::real& r884,
const cxsc::real& r885,
const cxsc::real& r886,
const cxsc::real& r887,
const cxsc::real& r888,
const cxsc::real& r889,
const cxsc::real& r890,
const cxsc::real& r891,
const cxsc::real& r892,
const cxsc::real& r893,
const cxsc::real& r894,
const cxsc::real& r895,
const cxsc::real& r896,
const cxsc::real& r897,
const cxsc::real& r898,
const cxsc::real& r899,
const cxsc::real& r900,
const cxsc::real& r901,
const cxsc::real& r902,
const cxsc::real& r903,
const cxsc::real& r904,
const cxsc::real& r905,
const cxsc::real& r906,
const cxsc::real& r907,
const cxsc::real& r908,
const cxsc::real& r909,
const cxsc::real& r910,
const cxsc::real& r911,
const cxsc::real& r912,
const cxsc::real& r913,
const cxsc::real& r914,
const cxsc::real& r915,
const cxsc::real& r916,
const cxsc::real& r917,
const cxsc::real& r918,
const cxsc::real& r919,
const cxsc::real& r920,
const cxsc::real& r921,
const cxsc::real& r922,
const cxsc::real& r923,
const cxsc::real& r924,
const cxsc::real& r925,
const cxsc::real& r926,
const cxsc::real& r927,
const cxsc::real& r928,
const cxsc::real& r929,
const cxsc::real& r930,
const cxsc::real& r931,
const cxsc::real& r932,
const cxsc::real& r933,
const cxsc::real& r934,
const cxsc::real& r935,
const cxsc::real& r936,
const cxsc::real& r937,
const cxsc::real& r938,
const cxsc::real& r939,
const cxsc::real& r940,
const cxsc::real& r941,
const cxsc::real& r942,
const cxsc::real& r943,
const cxsc::real& r944,
const cxsc::real& r945,
const cxsc::real& r946,
const cxsc::real& r947,
const cxsc::real& r948,
const cxsc::real& r949,
const cxsc::real& r950,
const cxsc::real& r951,
const cxsc::real& r952,
const cxsc::real& r953,
const cxsc::real& r954,
const cxsc::real& r955,
const cxsc::real& r956,
const cxsc::real& r957,
const cxsc::real& r958,
const cxsc::real& r959,
const cxsc::real& r960,
const cxsc::real& r961,
const cxsc::real& r962,
const cxsc::real& r963,
const cxsc::real& r964,
const cxsc::real& r965,
const cxsc::real& r966,
const cxsc::real& r967,
const cxsc::real& r968,
const cxsc::real& r969,
const cxsc::real& r970,
const cxsc::real& r971,
const cxsc::real& r972,
const cxsc::real& r973,
const cxsc::real& r974,
const cxsc::real& r975,
const cxsc::real& r976,
const cxsc::real& r977,
const cxsc::real& r978,
const cxsc::real& r979,
const cxsc::real& r980,
const cxsc::real& r981,
const cxsc::real& r982,
const cxsc::real& r983,
const cxsc::real& r984,
const cxsc::real& r985,
const cxsc::real& r986,
const cxsc::real& r987,
const cxsc::real& r988,
const cxsc::real& r989,
const cxsc::real& r990,
const cxsc::real& r991,
const cxsc::real& r992,
const cxsc::real& r993,
const cxsc::real& r994,
const cxsc::real& r995,
const cxsc::real& r996,
const cxsc::real& r997,
const cxsc::real& r998,
const cxsc::real& r999,
const cxsc::real& r1000) const
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
R.push_back(r101);
R.push_back(r102);
R.push_back(r103);
R.push_back(r104);
R.push_back(r105);
R.push_back(r106);
R.push_back(r107);
R.push_back(r108);
R.push_back(r109);
R.push_back(r110);
R.push_back(r111);
R.push_back(r112);
R.push_back(r113);
R.push_back(r114);
R.push_back(r115);
R.push_back(r116);
R.push_back(r117);
R.push_back(r118);
R.push_back(r119);
R.push_back(r120);
R.push_back(r121);
R.push_back(r122);
R.push_back(r123);
R.push_back(r124);
R.push_back(r125);
R.push_back(r126);
R.push_back(r127);
R.push_back(r128);
R.push_back(r129);
R.push_back(r130);
R.push_back(r131);
R.push_back(r132);
R.push_back(r133);
R.push_back(r134);
R.push_back(r135);
R.push_back(r136);
R.push_back(r137);
R.push_back(r138);
R.push_back(r139);
R.push_back(r140);
R.push_back(r141);
R.push_back(r142);
R.push_back(r143);
R.push_back(r144);
R.push_back(r145);
R.push_back(r146);
R.push_back(r147);
R.push_back(r148);
R.push_back(r149);
R.push_back(r150);
R.push_back(r151);
R.push_back(r152);
R.push_back(r153);
R.push_back(r154);
R.push_back(r155);
R.push_back(r156);
R.push_back(r157);
R.push_back(r158);
R.push_back(r159);
R.push_back(r160);
R.push_back(r161);
R.push_back(r162);
R.push_back(r163);
R.push_back(r164);
R.push_back(r165);
R.push_back(r166);
R.push_back(r167);
R.push_back(r168);
R.push_back(r169);
R.push_back(r170);
R.push_back(r171);
R.push_back(r172);
R.push_back(r173);
R.push_back(r174);
R.push_back(r175);
R.push_back(r176);
R.push_back(r177);
R.push_back(r178);
R.push_back(r179);
R.push_back(r180);
R.push_back(r181);
R.push_back(r182);
R.push_back(r183);
R.push_back(r184);
R.push_back(r185);
R.push_back(r186);
R.push_back(r187);
R.push_back(r188);
R.push_back(r189);
R.push_back(r190);
R.push_back(r191);
R.push_back(r192);
R.push_back(r193);
R.push_back(r194);
R.push_back(r195);
R.push_back(r196);
R.push_back(r197);
R.push_back(r198);
R.push_back(r199);
R.push_back(r200);
R.push_back(r201);
R.push_back(r202);
R.push_back(r203);
R.push_back(r204);
R.push_back(r205);
R.push_back(r206);
R.push_back(r207);
R.push_back(r208);
R.push_back(r209);
R.push_back(r210);
R.push_back(r211);
R.push_back(r212);
R.push_back(r213);
R.push_back(r214);
R.push_back(r215);
R.push_back(r216);
R.push_back(r217);
R.push_back(r218);
R.push_back(r219);
R.push_back(r220);
R.push_back(r221);
R.push_back(r222);
R.push_back(r223);
R.push_back(r224);
R.push_back(r225);
R.push_back(r226);
R.push_back(r227);
R.push_back(r228);
R.push_back(r229);
R.push_back(r230);
R.push_back(r231);
R.push_back(r232);
R.push_back(r233);
R.push_back(r234);
R.push_back(r235);
R.push_back(r236);
R.push_back(r237);
R.push_back(r238);
R.push_back(r239);
R.push_back(r240);
R.push_back(r241);
R.push_back(r242);
R.push_back(r243);
R.push_back(r244);
R.push_back(r245);
R.push_back(r246);
R.push_back(r247);
R.push_back(r248);
R.push_back(r249);
R.push_back(r250);
R.push_back(r251);
R.push_back(r252);
R.push_back(r253);
R.push_back(r254);
R.push_back(r255);
R.push_back(r256);
R.push_back(r257);
R.push_back(r258);
R.push_back(r259);
R.push_back(r260);
R.push_back(r261);
R.push_back(r262);
R.push_back(r263);
R.push_back(r264);
R.push_back(r265);
R.push_back(r266);
R.push_back(r267);
R.push_back(r268);
R.push_back(r269);
R.push_back(r270);
R.push_back(r271);
R.push_back(r272);
R.push_back(r273);
R.push_back(r274);
R.push_back(r275);
R.push_back(r276);
R.push_back(r277);
R.push_back(r278);
R.push_back(r279);
R.push_back(r280);
R.push_back(r281);
R.push_back(r282);
R.push_back(r283);
R.push_back(r284);
R.push_back(r285);
R.push_back(r286);
R.push_back(r287);
R.push_back(r288);
R.push_back(r289);
R.push_back(r290);
R.push_back(r291);
R.push_back(r292);
R.push_back(r293);
R.push_back(r294);
R.push_back(r295);
R.push_back(r296);
R.push_back(r297);
R.push_back(r298);
R.push_back(r299);
R.push_back(r300);
R.push_back(r301);
R.push_back(r302);
R.push_back(r303);
R.push_back(r304);
R.push_back(r305);
R.push_back(r306);
R.push_back(r307);
R.push_back(r308);
R.push_back(r309);
R.push_back(r310);
R.push_back(r311);
R.push_back(r312);
R.push_back(r313);
R.push_back(r314);
R.push_back(r315);
R.push_back(r316);
R.push_back(r317);
R.push_back(r318);
R.push_back(r319);
R.push_back(r320);
R.push_back(r321);
R.push_back(r322);
R.push_back(r323);
R.push_back(r324);
R.push_back(r325);
R.push_back(r326);
R.push_back(r327);
R.push_back(r328);
R.push_back(r329);
R.push_back(r330);
R.push_back(r331);
R.push_back(r332);
R.push_back(r333);
R.push_back(r334);
R.push_back(r335);
R.push_back(r336);
R.push_back(r337);
R.push_back(r338);
R.push_back(r339);
R.push_back(r340);
R.push_back(r341);
R.push_back(r342);
R.push_back(r343);
R.push_back(r344);
R.push_back(r345);
R.push_back(r346);
R.push_back(r347);
R.push_back(r348);
R.push_back(r349);
R.push_back(r350);
R.push_back(r351);
R.push_back(r352);
R.push_back(r353);
R.push_back(r354);
R.push_back(r355);
R.push_back(r356);
R.push_back(r357);
R.push_back(r358);
R.push_back(r359);
R.push_back(r360);
R.push_back(r361);
R.push_back(r362);
R.push_back(r363);
R.push_back(r364);
R.push_back(r365);
R.push_back(r366);
R.push_back(r367);
R.push_back(r368);
R.push_back(r369);
R.push_back(r370);
R.push_back(r371);
R.push_back(r372);
R.push_back(r373);
R.push_back(r374);
R.push_back(r375);
R.push_back(r376);
R.push_back(r377);
R.push_back(r378);
R.push_back(r379);
R.push_back(r380);
R.push_back(r381);
R.push_back(r382);
R.push_back(r383);
R.push_back(r384);
R.push_back(r385);
R.push_back(r386);
R.push_back(r387);
R.push_back(r388);
R.push_back(r389);
R.push_back(r390);
R.push_back(r391);
R.push_back(r392);
R.push_back(r393);
R.push_back(r394);
R.push_back(r395);
R.push_back(r396);
R.push_back(r397);
R.push_back(r398);
R.push_back(r399);
R.push_back(r400);
R.push_back(r401);
R.push_back(r402);
R.push_back(r403);
R.push_back(r404);
R.push_back(r405);
R.push_back(r406);
R.push_back(r407);
R.push_back(r408);
R.push_back(r409);
R.push_back(r410);
R.push_back(r411);
R.push_back(r412);
R.push_back(r413);
R.push_back(r414);
R.push_back(r415);
R.push_back(r416);
R.push_back(r417);
R.push_back(r418);
R.push_back(r419);
R.push_back(r420);
R.push_back(r421);
R.push_back(r422);
R.push_back(r423);
R.push_back(r424);
R.push_back(r425);
R.push_back(r426);
R.push_back(r427);
R.push_back(r428);
R.push_back(r429);
R.push_back(r430);
R.push_back(r431);
R.push_back(r432);
R.push_back(r433);
R.push_back(r434);
R.push_back(r435);
R.push_back(r436);
R.push_back(r437);
R.push_back(r438);
R.push_back(r439);
R.push_back(r440);
R.push_back(r441);
R.push_back(r442);
R.push_back(r443);
R.push_back(r444);
R.push_back(r445);
R.push_back(r446);
R.push_back(r447);
R.push_back(r448);
R.push_back(r449);
R.push_back(r450);
R.push_back(r451);
R.push_back(r452);
R.push_back(r453);
R.push_back(r454);
R.push_back(r455);
R.push_back(r456);
R.push_back(r457);
R.push_back(r458);
R.push_back(r459);
R.push_back(r460);
R.push_back(r461);
R.push_back(r462);
R.push_back(r463);
R.push_back(r464);
R.push_back(r465);
R.push_back(r466);
R.push_back(r467);
R.push_back(r468);
R.push_back(r469);
R.push_back(r470);
R.push_back(r471);
R.push_back(r472);
R.push_back(r473);
R.push_back(r474);
R.push_back(r475);
R.push_back(r476);
R.push_back(r477);
R.push_back(r478);
R.push_back(r479);
R.push_back(r480);
R.push_back(r481);
R.push_back(r482);
R.push_back(r483);
R.push_back(r484);
R.push_back(r485);
R.push_back(r486);
R.push_back(r487);
R.push_back(r488);
R.push_back(r489);
R.push_back(r490);
R.push_back(r491);
R.push_back(r492);
R.push_back(r493);
R.push_back(r494);
R.push_back(r495);
R.push_back(r496);
R.push_back(r497);
R.push_back(r498);
R.push_back(r499);
R.push_back(r500);
R.push_back(r501);
R.push_back(r502);
R.push_back(r503);
R.push_back(r504);
R.push_back(r505);
R.push_back(r506);
R.push_back(r507);
R.push_back(r508);
R.push_back(r509);
R.push_back(r510);
R.push_back(r511);
R.push_back(r512);
R.push_back(r513);
R.push_back(r514);
R.push_back(r515);
R.push_back(r516);
R.push_back(r517);
R.push_back(r518);
R.push_back(r519);
R.push_back(r520);
R.push_back(r521);
R.push_back(r522);
R.push_back(r523);
R.push_back(r524);
R.push_back(r525);
R.push_back(r526);
R.push_back(r527);
R.push_back(r528);
R.push_back(r529);
R.push_back(r530);
R.push_back(r531);
R.push_back(r532);
R.push_back(r533);
R.push_back(r534);
R.push_back(r535);
R.push_back(r536);
R.push_back(r537);
R.push_back(r538);
R.push_back(r539);
R.push_back(r540);
R.push_back(r541);
R.push_back(r542);
R.push_back(r543);
R.push_back(r544);
R.push_back(r545);
R.push_back(r546);
R.push_back(r547);
R.push_back(r548);
R.push_back(r549);
R.push_back(r550);
R.push_back(r551);
R.push_back(r552);
R.push_back(r553);
R.push_back(r554);
R.push_back(r555);
R.push_back(r556);
R.push_back(r557);
R.push_back(r558);
R.push_back(r559);
R.push_back(r560);
R.push_back(r561);
R.push_back(r562);
R.push_back(r563);
R.push_back(r564);
R.push_back(r565);
R.push_back(r566);
R.push_back(r567);
R.push_back(r568);
R.push_back(r569);
R.push_back(r570);
R.push_back(r571);
R.push_back(r572);
R.push_back(r573);
R.push_back(r574);
R.push_back(r575);
R.push_back(r576);
R.push_back(r577);
R.push_back(r578);
R.push_back(r579);
R.push_back(r580);
R.push_back(r581);
R.push_back(r582);
R.push_back(r583);
R.push_back(r584);
R.push_back(r585);
R.push_back(r586);
R.push_back(r587);
R.push_back(r588);
R.push_back(r589);
R.push_back(r590);
R.push_back(r591);
R.push_back(r592);
R.push_back(r593);
R.push_back(r594);
R.push_back(r595);
R.push_back(r596);
R.push_back(r597);
R.push_back(r598);
R.push_back(r599);
R.push_back(r600);
R.push_back(r601);
R.push_back(r602);
R.push_back(r603);
R.push_back(r604);
R.push_back(r605);
R.push_back(r606);
R.push_back(r607);
R.push_back(r608);
R.push_back(r609);
R.push_back(r610);
R.push_back(r611);
R.push_back(r612);
R.push_back(r613);
R.push_back(r614);
R.push_back(r615);
R.push_back(r616);
R.push_back(r617);
R.push_back(r618);
R.push_back(r619);
R.push_back(r620);
R.push_back(r621);
R.push_back(r622);
R.push_back(r623);
R.push_back(r624);
R.push_back(r625);
R.push_back(r626);
R.push_back(r627);
R.push_back(r628);
R.push_back(r629);
R.push_back(r630);
R.push_back(r631);
R.push_back(r632);
R.push_back(r633);
R.push_back(r634);
R.push_back(r635);
R.push_back(r636);
R.push_back(r637);
R.push_back(r638);
R.push_back(r639);
R.push_back(r640);
R.push_back(r641);
R.push_back(r642);
R.push_back(r643);
R.push_back(r644);
R.push_back(r645);
R.push_back(r646);
R.push_back(r647);
R.push_back(r648);
R.push_back(r649);
R.push_back(r650);
R.push_back(r651);
R.push_back(r652);
R.push_back(r653);
R.push_back(r654);
R.push_back(r655);
R.push_back(r656);
R.push_back(r657);
R.push_back(r658);
R.push_back(r659);
R.push_back(r660);
R.push_back(r661);
R.push_back(r662);
R.push_back(r663);
R.push_back(r664);
R.push_back(r665);
R.push_back(r666);
R.push_back(r667);
R.push_back(r668);
R.push_back(r669);
R.push_back(r670);
R.push_back(r671);
R.push_back(r672);
R.push_back(r673);
R.push_back(r674);
R.push_back(r675);
R.push_back(r676);
R.push_back(r677);
R.push_back(r678);
R.push_back(r679);
R.push_back(r680);
R.push_back(r681);
R.push_back(r682);
R.push_back(r683);
R.push_back(r684);
R.push_back(r685);
R.push_back(r686);
R.push_back(r687);
R.push_back(r688);
R.push_back(r689);
R.push_back(r690);
R.push_back(r691);
R.push_back(r692);
R.push_back(r693);
R.push_back(r694);
R.push_back(r695);
R.push_back(r696);
R.push_back(r697);
R.push_back(r698);
R.push_back(r699);
R.push_back(r700);
R.push_back(r701);
R.push_back(r702);
R.push_back(r703);
R.push_back(r704);
R.push_back(r705);
R.push_back(r706);
R.push_back(r707);
R.push_back(r708);
R.push_back(r709);
R.push_back(r710);
R.push_back(r711);
R.push_back(r712);
R.push_back(r713);
R.push_back(r714);
R.push_back(r715);
R.push_back(r716);
R.push_back(r717);
R.push_back(r718);
R.push_back(r719);
R.push_back(r720);
R.push_back(r721);
R.push_back(r722);
R.push_back(r723);
R.push_back(r724);
R.push_back(r725);
R.push_back(r726);
R.push_back(r727);
R.push_back(r728);
R.push_back(r729);
R.push_back(r730);
R.push_back(r731);
R.push_back(r732);
R.push_back(r733);
R.push_back(r734);
R.push_back(r735);
R.push_back(r736);
R.push_back(r737);
R.push_back(r738);
R.push_back(r739);
R.push_back(r740);
R.push_back(r741);
R.push_back(r742);
R.push_back(r743);
R.push_back(r744);
R.push_back(r745);
R.push_back(r746);
R.push_back(r747);
R.push_back(r748);
R.push_back(r749);
R.push_back(r750);
R.push_back(r751);
R.push_back(r752);
R.push_back(r753);
R.push_back(r754);
R.push_back(r755);
R.push_back(r756);
R.push_back(r757);
R.push_back(r758);
R.push_back(r759);
R.push_back(r760);
R.push_back(r761);
R.push_back(r762);
R.push_back(r763);
R.push_back(r764);
R.push_back(r765);
R.push_back(r766);
R.push_back(r767);
R.push_back(r768);
R.push_back(r769);
R.push_back(r770);
R.push_back(r771);
R.push_back(r772);
R.push_back(r773);
R.push_back(r774);
R.push_back(r775);
R.push_back(r776);
R.push_back(r777);
R.push_back(r778);
R.push_back(r779);
R.push_back(r780);
R.push_back(r781);
R.push_back(r782);
R.push_back(r783);
R.push_back(r784);
R.push_back(r785);
R.push_back(r786);
R.push_back(r787);
R.push_back(r788);
R.push_back(r789);
R.push_back(r790);
R.push_back(r791);
R.push_back(r792);
R.push_back(r793);
R.push_back(r794);
R.push_back(r795);
R.push_back(r796);
R.push_back(r797);
R.push_back(r798);
R.push_back(r799);
R.push_back(r800);
R.push_back(r801);
R.push_back(r802);
R.push_back(r803);
R.push_back(r804);
R.push_back(r805);
R.push_back(r806);
R.push_back(r807);
R.push_back(r808);
R.push_back(r809);
R.push_back(r810);
R.push_back(r811);
R.push_back(r812);
R.push_back(r813);
R.push_back(r814);
R.push_back(r815);
R.push_back(r816);
R.push_back(r817);
R.push_back(r818);
R.push_back(r819);
R.push_back(r820);
R.push_back(r821);
R.push_back(r822);
R.push_back(r823);
R.push_back(r824);
R.push_back(r825);
R.push_back(r826);
R.push_back(r827);
R.push_back(r828);
R.push_back(r829);
R.push_back(r830);
R.push_back(r831);
R.push_back(r832);
R.push_back(r833);
R.push_back(r834);
R.push_back(r835);
R.push_back(r836);
R.push_back(r837);
R.push_back(r838);
R.push_back(r839);
R.push_back(r840);
R.push_back(r841);
R.push_back(r842);
R.push_back(r843);
R.push_back(r844);
R.push_back(r845);
R.push_back(r846);
R.push_back(r847);
R.push_back(r848);
R.push_back(r849);
R.push_back(r850);
R.push_back(r851);
R.push_back(r852);
R.push_back(r853);
R.push_back(r854);
R.push_back(r855);
R.push_back(r856);
R.push_back(r857);
R.push_back(r858);
R.push_back(r859);
R.push_back(r860);
R.push_back(r861);
R.push_back(r862);
R.push_back(r863);
R.push_back(r864);
R.push_back(r865);
R.push_back(r866);
R.push_back(r867);
R.push_back(r868);
R.push_back(r869);
R.push_back(r870);
R.push_back(r871);
R.push_back(r872);
R.push_back(r873);
R.push_back(r874);
R.push_back(r875);
R.push_back(r876);
R.push_back(r877);
R.push_back(r878);
R.push_back(r879);
R.push_back(r880);
R.push_back(r881);
R.push_back(r882);
R.push_back(r883);
R.push_back(r884);
R.push_back(r885);
R.push_back(r886);
R.push_back(r887);
R.push_back(r888);
R.push_back(r889);
R.push_back(r890);
R.push_back(r891);
R.push_back(r892);
R.push_back(r893);
R.push_back(r894);
R.push_back(r895);
R.push_back(r896);
R.push_back(r897);
R.push_back(r898);
R.push_back(r899);
R.push_back(r900);
R.push_back(r901);
R.push_back(r902);
R.push_back(r903);
R.push_back(r904);
R.push_back(r905);
R.push_back(r906);
R.push_back(r907);
R.push_back(r908);
R.push_back(r909);
R.push_back(r910);
R.push_back(r911);
R.push_back(r912);
R.push_back(r913);
R.push_back(r914);
R.push_back(r915);
R.push_back(r916);
R.push_back(r917);
R.push_back(r918);
R.push_back(r919);
R.push_back(r920);
R.push_back(r921);
R.push_back(r922);
R.push_back(r923);
R.push_back(r924);
R.push_back(r925);
R.push_back(r926);
R.push_back(r927);
R.push_back(r928);
R.push_back(r929);
R.push_back(r930);
R.push_back(r931);
R.push_back(r932);
R.push_back(r933);
R.push_back(r934);
R.push_back(r935);
R.push_back(r936);
R.push_back(r937);
R.push_back(r938);
R.push_back(r939);
R.push_back(r940);
R.push_back(r941);
R.push_back(r942);
R.push_back(r943);
R.push_back(r944);
R.push_back(r945);
R.push_back(r946);
R.push_back(r947);
R.push_back(r948);
R.push_back(r949);
R.push_back(r950);
R.push_back(r951);
R.push_back(r952);
R.push_back(r953);
R.push_back(r954);
R.push_back(r955);
R.push_back(r956);
R.push_back(r957);
R.push_back(r958);
R.push_back(r959);
R.push_back(r960);
R.push_back(r961);
R.push_back(r962);
R.push_back(r963);
R.push_back(r964);
R.push_back(r965);
R.push_back(r966);
R.push_back(r967);
R.push_back(r968);
R.push_back(r969);
R.push_back(r970);
R.push_back(r971);
R.push_back(r972);
R.push_back(r973);
R.push_back(r974);
R.push_back(r975);
R.push_back(r976);
R.push_back(r977);
R.push_back(r978);
R.push_back(r979);
R.push_back(r980);
R.push_back(r981);
R.push_back(r982);
R.push_back(r983);
R.push_back(r984);
R.push_back(r985);
R.push_back(r986);
R.push_back(r987);
R.push_back(r988);
R.push_back(r989);
R.push_back(r990);
R.push_back(r991);
R.push_back(r992);
R.push_back(r993);
R.push_back(r994);
R.push_back(r995);
R.push_back(r996);
R.push_back(r997);
R.push_back(r998);
R.push_back(r999);
R.push_back(r1000);
    
    real result = 0;
    
    for (size_t i = 0; i < (R.size()-1); i++) 
    {
      result += 100.0 * (power((R[i+1] - power(R[i],2)),2)) +
									power((R[i] - 1), 2);
    }

	result = exp (-(1.0 * result));

   return result;
}
