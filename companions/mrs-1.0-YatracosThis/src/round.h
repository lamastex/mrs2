/*   File: round.h

     A header file that defines the directed  
     rounding modes for Linux and Sparc.

     Author: Warwick Tucker <warwick@math.uu.se>
     Latest edit: Wed Oct 15, 2001
*/


#if defined(__linux)
#include <fenv.h>
#define ROUND_DOWN  FE_DOWNWARD
#define ROUND_UP    FE_UPWARD
#define ROUND_NEAR  FE_TONEAREST
#define ROUND_ZERO  FE_TOWARDZERO
#define setRound    fesetround
#endif

#if defined(__sparc)
#include <ieeefp.h>
#define ROUND_DOWN  FP_RM
#define ROUND_UP    FP_RP
#define ROUND_NEAR  FP_RN
#define ROUND_ZERO  FP_RZ
#define setRound    fpsetround
#endif

void setRoundDown() { setRound(ROUND_DOWN); }
void setRoundUp  () { setRound(ROUND_UP);   }
void setRoundNear() { setRound(ROUND_NEAR); }
void setRoundZero() { setRound(ROUND_ZERO); }

