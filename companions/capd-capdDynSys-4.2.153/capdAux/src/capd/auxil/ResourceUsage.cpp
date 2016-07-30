/////////////////////////////////////////////////////////////////////////////
/// @file ResourceUsage.cpp
///
/// @author Mateusz Juda <mateusz.juda@{ii.uj.edu.pl,gmail.com}>
///
/// @date 2014-04-30
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#include <capd/auxil/ResourceUsage.h>
#include <capd/config-capdAux.h>

#ifdef HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif

#include <iostream>

#include <ctime>
#include <cstring>



using namespace capd::auxil;

namespace
{
#ifdef HAVE_SYS_RESOURCE_H
  rusage GetRUsage()
  {
    rusage result;
    memset(&result, 0, sizeof(result));

    rusage self;
    int res_s = getrusage(RUSAGE_SELF, &self);

    if (res_s == 0) {
      result = self;
    }
    return result;
  }

  double timevalToSec(const timeval& t)
  {
    return t.tv_sec + t.tv_usec/1000000.0;
  }
#endif
}

const double ResourceUsage::START_TIME = ResourceUsage::time();


ResourceUsage::ResourceUsage()
{
  _startStat = (*this)();
  _prevStat = _startStat;
}


size_t ResourceUsage::maxMemory() const
{
#ifdef HAVE_SYS_RESOURCE_H
  return GetRUsage().ru_maxrss;
#else
  return 0;
#endif
}

size_t ResourceUsage::currentMemory() const
{
#ifdef HAVE_SYS_RESOURCE_H
  rusage ru = GetRUsage();
  return ru.ru_idrss + ru.ru_isrss;
#else
  return 0;
#endif
}


double ResourceUsage::userCPU() const
{
#ifdef HAVE_SYS_RESOURCE_H
  return timevalToSec(GetRUsage().ru_utime);
#else
  return 0;
#endif
}

double ResourceUsage::sysCPU() const
{
#ifdef HAVE_SYS_RESOURCE_H
  return timevalToSec(GetRUsage().ru_stime);
#else
  return 0;
#endif
}

double ResourceUsage::elapsedTime() const
{
  return ResourceUsage::time() - START_TIME;
}

double ResourceUsage::time()
{
  return difftime(::time(NULL), time_t(0));
}


ResourceUsage::Stat ResourceUsage::operator()() const
{
  ResourceUsage::Stat stat;

  stat.maxMemory = maxMemory();
  stat.currentMemory = currentMemory();

  stat.userCPU = userCPU();
  stat.sysCPU = sysCPU();
  stat.elapsedTime = elapsedTime();

  return stat;
}

ResourceUsage::Stat ResourceUsage::usage() const
{
  return (*this)() - _startStat;
}

ResourceUsage::Stat ResourceUsage::delta()
{
  Stat curr = (*this)();
  Stat res =  curr - _prevStat;
  _prevStat = curr;

  return res;
}

std::ostream& operator<<(std::ostream& out, const capd::auxil::ResourceUsage::Stat& stat)
{
  out << "user CPU: " << stat.userCPU << "\n"
      << "sys CPU: " << stat.sysCPU << "\n"
      << "elapsedTime: " << stat.elapsedTime <<  "\n"
      << "memory: " << stat.currentMemory << " (" << stat.maxMemory << ")";


  return out;
}

  // out << "user CPU: " << stat.userCPU << " (" << stat.totalUserCPU << ")\n"
  //     << "sys CPU: " << stat.sysCPU << " (" << stat.totalSysCPU << ")\n"
  //     << "elapsedTime: " << stat.elapsedTime << " (" << stat.totalElapsedTime << ")\n"
  //     << "memory: " << stat.currentMemory << " (" << stat.maxMemory << ")";
