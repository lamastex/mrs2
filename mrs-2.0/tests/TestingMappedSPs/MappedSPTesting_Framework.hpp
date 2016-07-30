/*
* Copyright (C) 2009, 2010, 2011, 2012 Jennifer Harlow
*/

/*! \file
\brief Templatised procedures for testing MappedSPnode class
*/
#ifndef __MAPPEDSP_TESTING_FRAMEWORK_HPP__
#define __MAPPEDSP_TESTING_FRAMEWORK_HPP__

#include "mappedspnode.hpp"

namespace subpavings{
	/* This is really bad - supply a dummy implementation for these
	 * operators because they are not defined:  REALLY should
	 * not have an example for rvectors since it does not
	 * implement all of {+, - , *, /}*/
	cxsc::rvector& operator /=(cxsc::rvector& a, const cxsc::rvector& b)
	{ throw std::runtime_error("Rvector division attempted");}

	cxsc::rvector operator /(const cxsc::rvector& a, const cxsc::rvector& b)
	{ throw std::runtime_error("Rvector division attempted");}
}

using namespace cxsc;
using namespace std;
using namespace subpavings;


template<typename T>
void printTestTrees(const MappedSPnode<T>& lhs, const MappedSPnode<T>& rhs)
{
    cout << "lhs is " << endl;
    lhs.nodesAllOutput(cout, 0);
    cout << endl;

    cout << "rhs is " << endl;
    rhs.nodesAllOutput(cout, 0);
    cout << endl;

}


template<typename T>
void printTestTree(const MappedSPnode<T>& node)
{
    cout << "tree is " << endl;
    node.nodesAllOutput(cout, 0);
    cout << endl;

}




template<typename T>
void mappedAdditionTest(const MappedSPnode<T>& additiveIdentity,
						const MappedSPnode<T>& lhs,
						const MappedSPnode<T>& rhs)
{
   try {

		cout << "\n** addition test **\n" << endl;
     
		printTestTrees(lhs, rhs);

        cout << "Trying + operator, creating add = rhs + lhs" << endl;

        MappedSPnode<T> add = lhs + rhs;

		cout << "add is " << endl;
        add.nodesAllOutput(cout, 0);
        cout << endl;
		
		cout << "Trying addMinimised = minimiseLeaves on a copy of add" << endl;

        MappedSPnode<T> addMinimised(add);
		addMinimised.minimiseLeaves();

        cout << "addMinimised is " << endl;
        addMinimised.nodesAllOutput(cout, 0);
        cout << endl;

        
    }
    catch (exception& e) {
		std::cerr << "Error in mappedAdditionTest:\n" << e.what() << std::endl;
            throw;
    }


} // end of mapped addition test


template<typename T>
void mappedSubtractTest(const MappedSPnode<T>& additiveIdentity,
						const MappedSPnode<T>& lhs,
						const MappedSPnode<T>& rhs)
{
    try {

		cout << "\n** subtraction test **\n" << endl;
       
        printTestTrees(lhs, rhs);

        cout << "creating sub = rhs - lhs" << endl;
        cout << endl;

        MappedSPnode<T> sub = lhs - rhs;

        cout << "sub is " << endl;
        sub.nodesAllOutput(cout, 0);
        cout << endl;

		cout << "Trying subMinimised = minimiseLeaves on a copy of sub" << endl;

        MappedSPnode<T> subMinimised(sub);
		subMinimised.minimiseLeaves();

        cout << "subMinimised is " << endl;
        subMinimised.nodesAllOutput(cout, 0);
        cout << endl;

    }
    catch (exception& e) {
		std::cerr << "Error in mappedSubtractionTest:\n" << e.what() << std::endl;
            throw;
    }
}

template<typename T>
void mappedProductTest(const MappedSPnode<T>& multiplicativeIdentity,
						const MappedSPnode<T>& lhs,
						const MappedSPnode<T>& rhs)
{
    try {
		cout << "\n** multiply test **\n" << endl;
       
		printTestTrees(lhs, rhs);

        cout << "creating prod = rhs * lhs" << endl;
        cout << endl;

        MappedSPnode<T> prod = lhs * rhs;

        cout << "prod is " << endl;
        prod.nodesAllOutput(cout, 0);
        cout << endl;
		
		cout << "Trying prodMinimised = minimiseLeaves on a copy of prod" << endl;

        MappedSPnode<T> prodMinimised(prod);
		prodMinimised.minimiseLeaves();

        cout << "prodMinimised is " << endl;
        prodMinimised.nodesAllOutput(cout, 0);
        cout << endl;

    }
    catch (exception& e) {
		std::cerr << "Error in mappedProductTest:\n" << e.what() << std::endl;
            throw;
    }
}

template<typename T>
void mappedDivisionTest(const MappedSPnode<T>& multiplicativeIdentity,
						const MappedSPnode<T>& lhs,
						const MappedSPnode<T>& rhs)
{
    try {
		cout << "\n** division test **\n" << endl;
       
		printTestTrees(lhs, rhs);

        cout << "creating div = rhs / lhs" << endl;
        cout << endl;

        MappedSPnode<T> div = lhs / rhs;

        cout << "div is " << endl;
        div.nodesAllOutput(cout, 0);
        cout << endl;

		cout << "Trying divMinimised = minimiseLeaves on a copy of div" << endl;

        MappedSPnode<T> divMinimised(div);
		divMinimised.minimiseLeaves();

        cout << "divMinimised is " << endl;
        divMinimised.nodesAllOutput(cout, 0);
        cout << endl;

    }
    catch (exception& e) {
		std::cerr << "Error in mappedDivisionTest:\n" << e.what() << std::endl;
            throw;
    }
}

template<typename T>
void mappedScalarMultTest1(const MappedSPnode<T>& lhs,
                            const MappedSPnode<T>& rhs,
                            const T multiplier)
{
    try {

		cout << "\n** scalar multiplication test **\n" << endl;
       
        printTestTrees(lhs, rhs);

		cout << "creating lhs2 = lhs" << endl;
        cout << endl;
		MappedSPnode<T> lhs2(lhs);
        
		cout << "Trying own type scalar multiplication multLhs = lhs * " << multiplier << endl;
        cout << endl;

        MappedSPnode<T> multLhs = lhs * multiplier;

        cout << "multLhs now is " << endl;
        multLhs.nodesAllOutput(cout, 0);
        cout << endl;
		
		cout << "Trying own type scalar multiplication on lhs2*=" << multiplier << endl;
        cout << endl;

        lhs2*=multiplier;

        cout << "lhs2 now is " << endl;
        lhs2.nodesAllOutput(cout, 0);
        cout << endl;
		
		cout << "Trying lhs2Minimised = minimiseLeaves on a copy of lhs2" << endl;

        MappedSPnode<T> lhs2Minimised(lhs2);
		lhs2Minimised.minimiseLeaves();

        cout << "lhs2Minimised is " << endl;
        lhs2Minimised.nodesAllOutput(cout, 0);
        cout << endl;


        
    }
    catch (exception& e) {
		std::cerr << "Error in mappedScalarMultTest1:\n" << e.what() << std::endl;
            throw;
    }
}
/*
template<typename T, typename M>
void mappedScalarMultTest2(const MappedSPnode<T>& lhs,
                            const MappedSPnode<T>& rhs,
                            const M multiplier)
{
    try {


        cout << "creating addedLazy1 rhs lazy(+) lhs" << endl;
        cout << endl;

        MappedSPnode<T> addedLazy1 = MappedSPnode<T>::lazyCollationNonMinimalUnion(lhs, rhs);

        cout << "addedLazy1 is " << endl;
        addedLazy1.nodesAllOutput(cout, 0);
        cout << endl;

		cout << "creating addedLazy2 = addedLazy1" << endl;
        cout << endl;
		MappedSPnode<T> addedLazy2 = addedLazy1;

	
		cout << "creating lhs2 = lhs" << endl;
        cout << endl;
		MappedSPnode<T> lhs2(lhs);
        
		cout << "Trying different type scalar multiplication multLhs = lhs * " << multiplier << endl;
        cout << endl;

        MappedSPnode<T> multLhs = lhs * multiplier;

        cout << "multLhs now is " << endl;
        multLhs.nodesAllOutput(cout, 0);
        cout << endl;
		
		cout << "Trying different type scalar multiplication on lhs2 scalarMult(lhs2, " << multiplier << ")" << endl;
        cout << endl;

        scalarMult(lhs2, multiplier);

        cout << "lhs2 now is " << endl;
        lhs2.nodesAllOutput(cout, 0);
        cout << endl;

        cout << "Trying different type scalar multiplication scalarMult(adddedLazy1, " << multiplier << ")" << endl;
        cout << endl;

        scalarMult(addedLazy1, multiplier);

        cout << "addedLazy1 now is " << endl;
        addedLazy1.nodesAllOutput(cout, 0);
        cout << endl;

		cout << "Trying different type scalar multiplication multAdded2 = addedLazy2 * " << multiplier << endl;
        cout << endl;

        MappedSPnode<T> multAdded2 = addedLazy2 * multiplier;

        cout << "multAdded2 now is " << endl;
        multAdded2.nodesAllOutput(cout, 0);
        cout << endl;

    }
    catch (exception& e) {
		std::cerr << "Error in mappedScalarMultTest2:\n" << e.what() << std::endl;
            throw;
    }
}
*/
template<typename T>
void mappedScalarDivTest1(const MappedSPnode<T>& lhs,
                            const MappedSPnode<T>& rhs,
                            const T divisor)
{
    try {

		cout << "\n** scalar division test **\n" << endl;
       
		printTestTrees(lhs, rhs);

		cout << "creating lhs2 = lhs" << endl;
        cout << endl;
		MappedSPnode<T> lhs2(lhs);
       
		cout << "Trying own type scalar division on lhs2 lhs2/=" << divisor << endl;
        cout << endl;

        lhs2/=divisor;

        cout << "lhs2 now is " << endl;
        lhs2.nodesAllOutput(cout, 0);
        cout << endl;
		
		cout << "Trying lhs2Minimised = minimiseLeaves on a copy of lhs2" << endl;

        MappedSPnode<T> lhs2Minimised(lhs2);
		lhs2Minimised.minimiseLeaves();

        cout << "lhs2Minimised is " << endl;
        lhs2Minimised.nodesAllOutput(cout, 0);
        cout << endl;


		cout << "Trying own type scalar divison divLhs = lhs / " << divisor << endl;
        cout << endl;

        MappedSPnode<T> divLhs = lhs/divisor;

        cout << "divLhs now is " << endl;
        divLhs.nodesAllOutput(cout, 0);
        cout << endl;

        
    }
    catch (exception& e) {
		std::cerr << "Error in mappedScalarDivTest1:\n" << e.what() << std::endl;
            throw;
    }
}

/*
template<typename T, typename M>
void mappedScalarDivTest2(const MappedSPnode<T>& lhs,
                            const MappedSPnode<T>& rhs,
                            const M divisor)
{
    try {


		cout << "creating addedLazy1 rhs lazy(+) lhs" << endl;
        cout << endl;

        MappedSPnode<T> addedLazy1 = MappedSPnode<T>::lazyCollationNonMinimalUnion(lhs, rhs);

        cout << "addedLazy1 is " << endl;
        addedLazy1.nodesAllOutput(cout, 0);
        cout << endl;

		cout << "creating addedLazy2 = addedLazy1" << endl;
        cout << endl;
		MappedSPnode<T> addedLazy2 = addedLazy1;

        cout << "creating lhs2 = lhs" << endl;
        cout << endl;
		MappedSPnode<T> lhs2(lhs);
       
		cout << "Trying different type scalar division on lhs2 scalarDiv(lhs2, " << divisor << ")" << endl;
        cout << endl;

        scalarDiv(lhs2, divisor);

        cout << "lhs2 now is " << endl;
        lhs2.nodesAllOutput(cout, 0);
        cout << endl;

		cout << "Trying differeent type scalar divison divLhs = lhs / " << divisor << endl;
         cout << endl;

        MappedSPnode<T> divLhs = lhs/divisor;

        cout << "divLhs now is " << endl;
        divLhs.nodesAllOutput(cout, 0);
        cout << endl;

        cout << "Trying different type scalar division scalarDiv(addedLazy1, " << divisor << ")" << endl;
        cout << endl;

        scalarDiv(addedLazy1, divisor);

        cout << "addedLazy1 now is " << endl;
        addedLazy1.nodesAllOutput(cout, 0);
        cout << endl;

		cout << "Trying different type scalar division divAdded2 = addedLazy2 / " << divisor << endl;
        cout << endl;

        MappedSPnode<T> divAdded2 = addedLazy2 / divisor;

        cout << "divAdded2 now is " << endl;
        divAdded2.nodesAllOutput(cout, 0);
        cout << endl;

    }
    catch (exception& e) {
		std::cerr << "Error in mappedScalarDivTest2:\n" << e.what() << std::endl;
            throw;
    }
}
*/

template<typename T>
void mappedMinimiseTest(const MappedSPnode<T>& lhs, const MappedSPnode<T>& rhs)
{
   try {

		cout << "\n** minimise test **\n" << endl;
        printTestTrees(lhs, rhs);

		cout << "Creating add = rhs + lhs" << endl;

        MappedSPnode<T> add = rhs + lhs;
		
        cout << "addis " << endl;
        add.nodesAllOutput(cout, 0);
        cout << endl;

        cout << "Trying minimiseLeaves() on add" << endl;
        add.minimiseLeaves();

        cout << "add now is " << endl;
        add.nodesAllOutput(cout, 0);
        cout << endl;

    }
    catch (exception& e) {
		std::cerr << "Error in mappedMinimiseTest:\n" << e.what() << std::endl;
            throw;
    }


} // end of mapped minimise test

template<typename T>
void mappedArithmeticTest(const MappedSPnode<T>& empty,
						const MappedSPnode<T>& nonEmpty)
{
   try {

		cout << "\n** general arithmetic test **\n" << endl;
       
	    try {
			cout << "\nTrying addition empty + empty" << endl;
			
			MappedSPnode<T> tmp = empty + empty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		
		try {
			cout << "\nTrying addition empty + nonEmpty" << endl;
			
			MappedSPnode<T> tmp = empty + nonEmpty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		try {
			cout << "\nTrying addition nonEmpty + empty" << endl;
			
			MappedSPnode<T> tmp = nonEmpty + empty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		try {
			cout << "\nTrying addition nonEmpty + nonEmpty" << endl;
			
			MappedSPnode<T> tmp = nonEmpty + nonEmpty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}

try {
			cout << "\nTrying subtraction empty - empty" << endl;
			
			MappedSPnode<T> tmp = empty - empty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		
		try {
			cout << "\nTrying subtraction empty - nonEmpty" << endl;
			
			MappedSPnode<T> tmp = empty - nonEmpty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		try {
			cout << "\nTrying subtraction nonEmpty - empty" << endl;
			
			MappedSPnode<T> tmp = nonEmpty - empty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		try {
			cout << "\nTrying subtraction nonEmpty - nonEmpty" << endl;
			
			MappedSPnode<T> tmp = nonEmpty - nonEmpty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		
		
		try {
			cout << "\nTrying multiplication empty * empty" << endl;
			
			MappedSPnode<T> tmp = empty * empty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		
		try {
			cout << "\nTrying multiplication empty * nonEmpty" << endl;
			
			MappedSPnode<T> tmp = empty * nonEmpty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		try {
			cout << "\nTrying multiplication nonEmpty * empty" << endl;
			
			MappedSPnode<T> tmp = nonEmpty * empty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		try {
			cout << "\nTrying multiplication nonEmpty * nonEmpty" << endl;
			
			MappedSPnode<T> tmp = nonEmpty * nonEmpty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		
			
		try {
			cout << "\nTrying division empty / empty" << endl;
			
			MappedSPnode<T> tmp = empty / empty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		
		try {
			cout << "\nTrying division empty / nonEmpty" << endl;
			
			MappedSPnode<T> tmp = empty / nonEmpty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		try {
			cout << "\nTrying division nonEmpty / empty" << endl;
			
			MappedSPnode<T> tmp = nonEmpty / empty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}
		try {
			cout << "\nTrying division nonEmpty / nonEmpty" << endl;
			
			MappedSPnode<T> tmp = nonEmpty / nonEmpty;
			cout << "...sucess" << endl;
			tmp.nodesAllOutput(cout, 0);
		}
		catch (subpavings::IncompatibleDimensions_Error& e) {
			cout << "IncompatibleDimensions_Error:\n" << e.what() << endl;
		}
		catch (std::exception& ee) {
			cout << "Exception:\n" << ee.what() << endl;
			throw;
		}
		catch (...) {
			cout << "Exception" << endl;
			throw;
		}


    }
    catch (exception& e) {
		std::cerr << "Error in mappedArithmeticTest:\n" << e.what() << std::endl;
            throw;
    }


} // end of mapped arithmetic test

#endif
