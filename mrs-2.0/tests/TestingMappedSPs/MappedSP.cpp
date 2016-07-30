/*
* Copyright (C) 2009, 2010, 2011, 2012 Jennifer Harlow
*/

/*! \file
\brief Testing MappedSPnode class
*/


#include "MappedSPTesting.hpp"

// to use std input/output
#include <iostream>

// to use exceptions
#include <stdexcept>

#include <signal.h> // signal handling


using namespace std;

// see http://publications.gbdirect.co.uk/c_book/chapter9/signal_handling.html
// also see http://www.ibm.com/developerworks/linux/library/l-cppexcep/index.html

void div_by_zero_handler(int sig) {
	cerr << "Div by zero exception about to be thrown" << endl;
	throw std::runtime_error("Division by zero");

}
	
    
int main()
{
	
	// div by zero handler
	(void) signal(SIGFPE, div_by_zero_handler);

	
    try {
        
		// testing general arithmetic
        testingArithmetic();
		
		// testing for ints
        testingInts();

        // testing for reals
        testingReals();

		testingRealsMarginaliseFail();
		
		testingRealsMarginalise();
		

        // testing for intervals
        testingIntervals();

        // testing for intervals
        testingRvectors();
		
		// testing for bools
        testingBools();
        
        cout << "End of test" << endl;
		
		return 0;

   }
   catch (exception& e) {
            cerr << "std::exception error in testing:\n" << e.what() << std::endl;
   }
   catch (...) {
			cout << "some other exception" << endl;
	}

 
} // end of test

