/////////////////////////////////////////////////////////////////////////////
/// @file ResourceUsage.h
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

#ifndef CAPD_FILE_RESOURCEUSAGE_H
#define CAPD_FILE_RESOURCEUSAGE_H

#include <iosfwd>
#include <ctime>

namespace capd
{

  namespace auxil
  {
    class ResourceUsage
    {
    public:

      struct Stat
      {
	Stat():
	  userCPU(0),
	  sysCPU(0),
	  elapsedTime(0),
	  maxMemory(0),
	  currentMemory(0)
	{}

	double userCPU;
	double sysCPU;
	double elapsedTime;
	size_t maxMemory;
	size_t currentMemory;

	Stat operator-(const Stat& rhs) const
	{
	  Stat res = *this;
	  res.userCPU -= rhs.userCPU;
	  res.sysCPU -= rhs.sysCPU;
	  res.elapsedTime -= rhs.elapsedTime;
	  res.maxMemory -= rhs.maxMemory;
	  res.currentMemory -= rhs.currentMemory;

	  return res;
	}
      };

      static const double START_TIME;

      ResourceUsage();

      size_t maxMemory() const;
      size_t currentMemory() const;

      double userCPU() const;
      double sysCPU() const;
      double elapsedTime() const;

      static double time();

      Stat operator()() const;
      Stat usage() const;
      Stat delta();

    private:
      Stat _startStat;
      Stat _prevStat;
    };
  }

}

std::ostream& operator<<(std::ostream& out, const capd::auxil::ResourceUsage::Stat& stat);
std::ostream& operator<<(std::ostream& out, const capd::auxil::ResourceUsage& ru);

#endif // CAPD_FILE_RESOURCEUSAGE_H
