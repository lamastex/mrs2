/* 
 * Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009 Raazesh Sainudiin 
 * Copyright (C) 2009 Jennifer Harlow
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
/*! \file examples/GOpt/PhyloPOT/PhyloPOT.cpp
\brief Example to use FPhyloPOT and GOpt to do global optimisation 
for phylogenetic problem.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include <string>

#include <math.h>
#include <time.h>

//CSXC includes
#include "interval.hpp"		// Include interval arithmetic package
#include "imath.hpp"		// Include interval standard functions
#include "rmath.hpp"		// Include real standard functions
#include "intvector.hpp"
#include "ivector.hpp"
#include "rvector.hpp"
#include "imatrix.hpp"

#include <stdio.h>
#include <stdlib.h>

#include "toolz.hpp"
#include "SmallClasses.hpp"
#include "Fobj.hpp"
#include "FPhyloPOT.hpp"
#include "PhyloTree.hpp"
#include "GOpt.hpp"


using namespace std;
using namespace cxsc;



int main (int argc, char **argv)
{
	ios::sync_with_stdio ();	// call this function so iostream works with stdio
	cout << SetPrecision (20, 15);	// Number of mantissa digits in I/O
	
	// default value
	double Alb = 1.0;// partition until lower bound on Acceptance Prob. is > Alb
	
	// default value
	int n_boxes = 500;
	
	// default value
	int n_samples = 100;
	
	// default value
	unsigned theSeed = 1234;
	
	// default value
	bool CenteredForm = true;
        int CenFrm=1;

	// default value
	bool UseLogPi = true;
	
	// default value
	bool use_f_scale = true;
	
	// default value
	int prior_type = 0;
	
	// parameters for FPhyloPOT
	// default value
	//interval DomainCFN3(0.0000000001,1);
	interval DomainCFN3(0.0000000001,0.5);
	
	// default value
	int tree_space = 3;	// number of taxa
	int character_space = 4;	// number of DNA states 4=JC69, 2=CFN model
	
	// check the input 
	if (argc >= 2)
	{
		sscanf (argv[1], "%i", &tree_space);
                if (argc>=3)
		  	sscanf (argv[2], "%i", &character_space);
		if (argc >= 4)
			sscanf (argv[3], "%i", &CenFrm);
		cout << "# Usage: PhyloPOT <tree_space> " << tree_space  
		     << " <character_space> " << character_space 
		     << " <Centered_Form> " << CenteredForm << endl
                     << "extra arguments ignored.\n";
	}
	
        CenteredForm=CenFrm;

	cout << "# tree_space: " << tree_space << " # characters: " << character_space 
		<< " Centered Form " << CenteredForm << endl;
	
	
	FPhyloPOT FPhylo(tree_space, character_space, DomainCFN3, CenteredForm, UseLogPi, prior_type);
	
	// ***************** global optimisation ****************
	
		
	// set up a search box
  // get the dimensions of the domain for this first tree
	int domain_dim = FPhylo.getLabeledDomainDim(0); 
	ivector search (1, domain_dim);
	for (int i = 1; i <= domain_dim; i++)
	{
		search[i] = DomainCFN3;			// note smaller domain used here for gopt
	}
	
	real tolerance;	// set a tolerance
	tolerance = 1e-8;
	
	// cycle through the trees and do global optimisation for minimums and maximums
	int noTrees = FPhylo.getNoTrees();
	std::cout << "Note there are " << noTrees 
            << " trees for this function object." << std::endl << std::endl;
	for(int i = 0; i < noTrees; i++) {
		
		std::cout << "Performing global optimisation only on tree number " 
              << i << std::endl;
				
		// Minimums
    // call GoptMin with pointer to FPhylo for tree i
		//GOptMin(&FPhylo, search, tolerance, i);
		
		// Maximums
    // call GoptMax with pointer to FPhylo for tree i
		GOptMax(&FPhylo, search, tolerance, i); 
    
    /*! \todo Need to do simultaneous glabal opoptimisation over all labeles
      see Raaz's 2004 code with hacked older C-XSC libs*/
			
	}
	
	
	return 0;			
}
