mrs 2.0, a C++ class library for statistical set processing and computer-aided proofs in statistics.
====================================================================================================

created:	Sat Jul 30 18:02:29 CEST 2016

This software is copy-lefted and distributed under the terms of the 
GNU General Public License (GPL).
see http://www.gnu.org/copyleft/gpl.html

See specific source files for author lists by year.

The authors make no warranties with respect to the adequacies of
this program/code for any particular purpose or with respect to
any particular result.  They make no warranties, express or implied,
that the programs are free of error, or are consistent with any
standard of merchantibility, or that they will meet your requirements
for any particulr purpose or application.  In particular, they should 
not be relied on for solving a problem whose incorrect solution could
result in injury to a person or loss of property.  If you do use the
programs in such a manner, it is entirely at your own risk.  In no
event shall the authors be liable for special, direct, indirect or
consequential damages, losses, costs, charges, claims, demands or
claim for lost profits, fees or expenses of any nature or kind.
Also see the disclaimers in the C-XSC library that this program 
builds upon. 

Dependencies (see ./companions):
================================

1.
You are assumed to have C-XSC2.5* already installed
http://www2.math.uni-wuppertal.de/~xsc/
C-XSC is distributed under the terms of the GNU General Public License (GPL).

2.
You are also expected to have GSL already installed
http://www.gnu.org/software/gsl/
GSL is distributed under the terms of the GNU General Public License (GPL).
You may need libgsl0-dev also.

3.
Computer-aided Proofs in Dynamics:
http://capd.ii.uj.edu.pl/

4.
The Computational Geometry Algorithms Library:
http://www.cgal.org/

Other dependencies include boost++ and standard GNU packages and collections.

Distribution:
https://github.com/raazesh-sainudiin/mrs2

Citation:
"mrs 2.0, a C++ class library for statistical set processing and computer-aided proofs in statistics", Raazesh Sainudiin and Thomas York,

mrs-2.0 configure and make instructions:
========================================

Modify paths to companion libraries in `custom_config.sh` and then:

        cd mrs-2.0
        ./custom_config.sh
        make

For a basic installation, run
=============================

        make install

In addition, the (highly incomplete!) documentation can be compiled with

        make doxygen-doc

If the required dependencies are located in places other than the default,
will you will need to specify their locations as arguments to the configure
script. For example:

	./configure CPPFLAGS='-I /usr/local/cxsc/include' \
                    LDFLAGS='-L /usr/local/cxsc/lib' \
	            --prefix='/home/user_name/mrs'

The file `custom_config.sh' can also be used, which contains a variation on the 
above command.  Make sure you export the environmet variable LD_LIBRARY_PATH 
properly (for example see last few lines of custom_config.sh):
export LD_LIBRARY_PATH=/usr/local/cxsc/lib:$LD_LIBRARY_PATH

Supported by
============

* 2003 by NSF grant DGE-9870631 (United States)
* 2004, 2005 by joint NSF/NIGMS grant DMS-02-01037 (United States)
* 2006, 2007 by Research Fellowship of the Royal Commission for the Exhibition of 1851 (United Kingdom)
* 2007, 2008 by internal research grants from Mathematics and Statistics Department, University of Canterbury (New Zealand)
* 2009 by statistical consulting revenues of Raazesh Sainudiin
* 2012 University of Canterbury Postgraduate Scholarship 
* 2015, 2016 by the project CORCON: Correctness by Construction, Seventh Framework Programme of the European Union, Marie Curie Actions-People, International Research Staff Exchange Scheme with counter-part funding by The Royal Society of New Zealand
