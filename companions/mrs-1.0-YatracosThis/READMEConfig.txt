To turn off debugging option -g, open file configure.

Change "-g" to "-Wall":

if test "$ac_test_CXXFLAGS" = set; then
  CXXFLAGS=$ac_save_CXXFLAGS
elif test $ac_cv_prog_cxx_g = yes; then
  if test "$GXX" = yes; then
    #CXXFLAGS="-g -O2"
    CXXFLAGS="-Wall -O2"
  else
    CXXFLAGS="-g"
  fi
else
  if test "$GXX" = yes; then
    CXXFLAGS="-O2"
  else
    CXXFLAGS=
  fi
fi
