/*
* Copyright (C) 2009, 2010, 2011, 2012 Jennifer Harlow
*/
/*! \file
\brief Testing MappedSPnode class
*/

#include "MappedSPTesting_Framework.hpp"
#include "realmappedspnode.hpp"

#include <iterator>
#include <stdexcept>

using namespace cxsc;
using namespace std;
using namespace subpavings;


void doMarginals(const subpavings::RealMappedSPnode& node);
void doMarginal(const subpavings::RealMappedSPnode& node, 
				const std::vector < int >& reqDims);

void doMarginalsFail(const subpavings::RealMappedSPnode& node);
void doMarginalFail(const subpavings::RealMappedSPnode& node, 
				const std::vector < int >& reqDims);





void testingArithmetic()
{
    try {
        // testing arithmetic with reals
		
		cout << "\nTesting arithmetic (with reals)\n" << endl;

        int dims = 1;
        ivector pavingBox(dims);
        interval pavingInterval(-4,4);
        for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;

		RealMappedSPnode empty;
        RealMappedSPnode nonEmpty(pavingBox); // make another MappedSPnode object

        // allocate ranges
		cxsc::real range3(3.0);
       
        std::vector<real> ranges;
        ranges.push_back(range3);
        
        nonEmpty.allocateRanges(ranges);

		mappedArithmeticTest(empty, nonEmpty);
        
    }
    catch (exception& e) {
		std::cerr << "Error in testingArithmetic:\n" << e.what() << std::endl;
            throw;
    }



} // end of test for ints


void testingInts()
{
    try {
        // testing for ints
		
		cout << "\nTesting for ranges as ints\n" << endl;

        int dims = 1;
        ivector pavingBox(dims);
        interval pavingInterval(-4,4);
        for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;

		MappedSPnode<int> additiveIdentity(pavingBox); 
		MappedSPnode<int> multiplicativeIdentity(pavingBox); 
        MappedSPnode<int> nodeOne(pavingBox); // make a MappedSPnode object
        MappedSPnode<int> nodeTwo(pavingBox); // make another MappedSPnode object

        // split the nodes to particular shapes
        nodeOne.splitToShape("1,1"); // split it
        nodeTwo.splitToShape("2,2,1"); // split it

        // allocate ranges
		int range0 = 0;
        int range1 = 1;
		int range2 = 2;
		int range3 = 3;
		int range4 = 4;
		int range5 = 5;
        
        std::vector<int> rangesOne;
        rangesOne.push_back(range2);
        rangesOne.push_back(range1);
        rangesOne.push_back(range3);

        nodeOne.allocateRanges(rangesOne, 0);


        // allocate ranges for nodeTwo
        //std::vector<int>  rangeTwo2(1, 3);
        //std::vector<int>  rangeTwo3(1, 1);
        //std::vector<int>  rangeTwo5(, 4);

        vector<int> rangesTwo;
        rangesTwo.push_back(range4);
        rangesTwo.push_back(range3);
        rangesTwo.push_back(range4);
        rangesTwo.push_back(range2);
        rangesTwo.push_back(range5);

        nodeTwo.allocateRanges(rangesTwo, 0);

		std::vector<int> rangesAddId;
		rangesAddId.push_back(range0);
		
		additiveIdentity.allocateRanges(rangesAddId);

		std::vector<int> rangesMultId;
		rangesMultId.push_back(range1);
		
		multiplicativeIdentity.allocateRanges(rangesMultId);

		mappedAdditionTest(additiveIdentity, nodeOne, nodeTwo);
        mappedSubtractTest(additiveIdentity, nodeOne, nodeTwo);
        mappedProductTest(multiplicativeIdentity, nodeOne, nodeTwo);
		mappedDivisionTest(multiplicativeIdentity, nodeOne, nodeTwo);

        {
			int multiplier = 4; // supply a multiplier of the right type

			mappedScalarMultTest1(nodeOne, nodeTwo, multiplier);
			mappedScalarDivTest1(nodeOne, nodeTwo, multiplier);
		}

		MappedSPnode<int> nodeThree(pavingBox); 
        MappedSPnode<int> nodeFour(pavingBox); 
		MappedSPnode<int> nodeFive(pavingBox);
		
		// split the nodes to particular shapes
        nodeThree.splitToShape("2,3,3,1"); // split it
		nodeFour.splitToShape("2,3,3,1"); // split it
        nodeFive.splitToShape("2,3,3,2,3,3"); // split it

        std::vector<int> rangesThree;
        rangesThree.push_back(range3);
        rangesThree.push_back(range3);
        rangesThree.push_back(range3);
		rangesThree.push_back(range3);
		rangesThree.push_back(range2);
		rangesThree.push_back(range4);
		rangesThree.push_back(range3);

        nodeThree.allocateRanges(rangesThree, 0);

		std::vector<int> rangesFour;
        rangesFour.push_back(range3);
        rangesFour.push_back(range3);
        rangesFour.push_back(range3);
		rangesFour.push_back(range3);
		rangesFour.push_back(range4);
		rangesFour.push_back(range2);
		rangesFour.push_back(range3);

        nodeFour.allocateRanges(rangesFour, 0);
		
		std::vector<int> rangesFive;
        rangesFive.push_back(range3);
        rangesFive.push_back(range3);
        rangesFive.push_back(range3);
		rangesFive.push_back(range3);
		rangesFive.push_back(range4);
		rangesFive.push_back(range2);
		rangesFive.push_back(range3);
        rangesFive.push_back(range3);
		rangesFive.push_back(range3);
		rangesFive.push_back(range2);
		rangesFive.push_back(range4);
		

        nodeFive.allocateRanges(rangesFive, 0);
		
		mappedMinimiseTest(nodeThree, nodeFour);
		mappedMinimiseTest(nodeThree, nodeFive);
		

    }
    catch (exception& e) {
		std::cerr << "Error in testingInts:\n" << e.what() << std::endl;
            throw;
    }



} // end of test for ints


void testingReals()
{
    try {

        // testing for reals

        cout << "\nTesting for ranges as csxc::reals\n" << endl;

        int dims = 1;
        ivector pavingBox(dims);
        interval pavingInterval(-4,4);
        for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;

		RealMappedSPnode additiveIdentity(pavingBox); 
		RealMappedSPnode multiplicativeIdentity(pavingBox); 
       
        RealMappedSPnode nodeOne(pavingBox); // make a MappedSPnode object
        RealMappedSPnode nodeTwo(pavingBox); // make another MappedSPnode object

        // split the nodes to particular shapes
        nodeOne.splitToShape("1,1"); // split it
        nodeTwo.splitToShape("2,2,1"); // split it

        // allocate ranges
		cxsc::real range0(0.0);
        cxsc::real range1(1.0);
        cxsc::real range2(2.0);
        cxsc::real range3(3.0);
		cxsc::real range4(4.0);
        cxsc::real range5(5.0);
        
        std::vector<real> rangesOne;
        rangesOne.push_back(range2);
        rangesOne.push_back(range1);
        rangesOne.push_back(range3);

        nodeOne.allocateRanges(rangesOne, 0);


        // allocate ranges for nodeTwo
        //std::vector<int>  rangeTwo2(3);
        //std::vector<int>  rangeTwo3(1);
        //std::vector<int>  rangeTwo5(4);

        std::vector<real> rangesTwo;
        rangesTwo.push_back(range4);
        rangesTwo.push_back(range3);
        rangesTwo.push_back(range4);
        rangesTwo.push_back(range2);
        rangesTwo.push_back(range5);

        nodeTwo.allocateRanges(rangesTwo, 0);
		
		std::vector<real> rangesAddId;
		rangesAddId.push_back(range0);
		
		additiveIdentity.allocateRanges(rangesAddId);

		std::vector<real> rangesMultId;
		rangesMultId.push_back(range1);
		
		multiplicativeIdentity.allocateRanges(rangesMultId);


        printTestTrees(nodeOne, nodeTwo);

		mappedAdditionTest(additiveIdentity, nodeOne, nodeTwo);
        mappedSubtractTest(additiveIdentity, nodeOne, nodeTwo);
        mappedProductTest(multiplicativeIdentity, nodeOne, nodeTwo);
		mappedDivisionTest(multiplicativeIdentity, nodeOne, nodeTwo);
		
        real multiplier(4.0); // supply a multiplier of the right type

        mappedScalarMultTest1(nodeOne, nodeTwo, multiplier);
		mappedScalarDivTest1(nodeOne, nodeTwo, multiplier);
		
        
    }
    catch (exception& e) {
		std::cerr << "Error in testingReals:\n" << e.what() << std::endl;
            throw;
    }


} // end of test for reals


void testingRealsMarginalise()
{
    try {

        // testing for reals with marginalising

        cout << "\nTesting for marginalising as csxc::reals\n" << endl;

		int maxdims = 3;
		
		for (int dims = 2; dims <= maxdims; ++dims) {
			
			ivector pavingBox(dims);
			interval pavingInterval(0,1);
			for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;

			//RealMappedSPnode additiveIdentity(pavingBox); 
			//RealMappedSPnode multiplicativeIdentity(pavingBox); 
		   
			if (dims > 0) {
				RealMappedSPnode node(pavingBox); // make a MappedSPnode object
			
				node.splitToShape("1,1"); // split it
				
				std::vector<real> ranges;
				
				// allocate ranges
				ranges.push_back( real(1.0));
				ranges.push_back( real(1.0/3*2));
				ranges.push_back( real(2.0/3*2));
				node.allocateRanges(ranges, 0);
				
				doMarginals(node);
				
			}
			
			if (dims > 1) {
				RealMappedSPnode node(pavingBox); // make a MappedSPnode object
			
				node.splitToShape("2,2,1"); // split it
				
				std::vector<real> ranges;
				
				// allocate ranges
				ranges.push_back( real(1.0));
				ranges.push_back( real(3.0/7*2));
				ranges.push_back( real(1.0/7*4));
				ranges.push_back( real(2.0/7*4));
				ranges.push_back( real(4.0/7*2));
				
				node.allocateRanges(ranges, 0);
				
				doMarginals(node);
			}
			
			if (dims > 2) {
				RealMappedSPnode node(pavingBox); // make a MappedSPnode object
			
				node.splitToShape("2,3,3,1"); // split it
				
				std::vector<real> ranges;
				
				// allocate ranges
				ranges.push_back( real(1.0));
				ranges.push_back( real(6.0/10*2));
				ranges.push_back( real(1.0/10*4));
				ranges.push_back( real(5.0/10*4));
				ranges.push_back( real(2.0/10*8));
				ranges.push_back( real(3.0/10*8));
				ranges.push_back( real(4.0/10*2));
				
				node.allocateRanges(ranges, 0);
				
				doMarginals(node);

			}
		}
			
    }
    catch (exception& e) {
		std::cerr << "Error in testingRealsMarginalise:\n" << e.what() << std::endl;
            throw;
    }


} // end of test for reals marginalising

void testingRealsMarginaliseFail()
{
    try {

        // testing for reals with marginalising - cases that should fail

        cout << "\nTesting for marginalising as csxc::reals\n" << endl;

		int maxdims = 3;
		
		for (int dims = 2; dims <= maxdims; ++dims) {
			
			ivector pavingBox(dims);
			interval pavingInterval(0,1);
			for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;

			//RealMappedSPnode additiveIdentity(pavingBox); 
			//RealMappedSPnode multiplicativeIdentity(pavingBox); 
		   
			if (dims == 1 ) {
				RealMappedSPnode node(pavingBox); // make a MappedSPnode object
			
				node.splitToShape("1,1"); // split it
				
				std::vector<real> ranges;
				
				// allocate ranges
				ranges.push_back( real(1.0));
				ranges.push_back( real(1.0/3*2));
				ranges.push_back( real(2.0/3*2));
				node.allocateRanges(ranges, 0);
				
				doMarginalsFail(node);
				
			}
			
			if (dims == 2) {
				RealMappedSPnode node(pavingBox); // make a MappedSPnode object
			
				node.splitToShape("2,2,1"); // split it
				
				std::vector<real> ranges;
				
				// allocate ranges
				ranges.push_back( real(1.0));
				ranges.push_back( real(3.0/7*2));
				ranges.push_back( real(1.0/7*4));
				ranges.push_back( real(2.0/7*4));
				ranges.push_back( real(4.0/7*2));
				
				node.allocateRanges(ranges, 0);
				
				doMarginalsFail(node);
			}
			
			if (dims > 2) {
				RealMappedSPnode node(pavingBox); // make a MappedSPnode object
			
				node.splitToShape("2,3,3,1"); // split it
				
				std::vector<real> ranges;
				
				// allocate ranges
				ranges.push_back( real(1.0));
				ranges.push_back( real(6.0/10*2));
				ranges.push_back( real(1.0/10*4));
				ranges.push_back( real(5.0/10*4));
				ranges.push_back( real(2.0/10*8));
				ranges.push_back( real(3.0/10*8));
				ranges.push_back( real(4.0/10*2));
				
				node.allocateRanges(ranges, 0);
				
				doMarginalsFail(node);

			}
		}
			
    }
    catch (exception& e) {
		std::cerr << "Error in testingRealsMarginaliseFail:\n" << e.what() << std::endl;
            throw;
    }


} // end of test for reals


void doMarginals(const RealMappedSPnode& node) {

	int dims = node.getDimension();
	cout << "\n\nTest node with box dimensions = " << dims << endl;
	printTestTree(node);
	cout << "Marginals :\n" << endl;
				
	for (int i = 1; i <= dims; ++i) {
		
		// do the 1-d ones
		std::vector <int> reqDims1;
		reqDims1.push_back(i);
		doMarginal(node, reqDims1);
		
		for (int j = i+1; j <= dims; ++j) {
			
			//do the 2-d ones
			std::vector <int> reqDims2 = reqDims1;
			reqDims2.push_back(j);
			doMarginal(node, reqDims2);
			
			for (int k = j+1; k <= dims; ++k) {
				
				// do the 3-d ones
				std::vector <int> reqDims3 = reqDims2;
				reqDims3.push_back(k);
				doMarginal(node, reqDims3);
				
			}
			
		}
		
	}
}

// cases that should fail
void doMarginalsFail(const RealMappedSPnode& node) {

	int dims = node.getDimension();
	cout << "\n\nTest marginalisations that should fail with node with box dimensions = " << dims << endl;
	
	std::vector <int> reqDims;
	reqDims.push_back(-1);
	
	doMarginalFail(node, reqDims);
				
	for (int i = dims; i <= dims+1; ++i) {
		
		std::vector <int> reqDims1;
		reqDims1.push_back(i);
			
		// do the 1-d ones
		if (i > dims) {
			doMarginalFail(node, reqDims1);
		}
		
		for (int j = dims; j <= dims+1; ++j) {
				
			//do the 2-d ones
			std::vector <int> reqDims2 = reqDims1;
			reqDims2.push_back(j);
			
			if (i > dims || j > dims) {
				doMarginalFail(node, reqDims2);
			}
				
			for (int k = dims; k <= dims+1; ++k) {
					
					// do the 3-d ones
					std::vector <int> reqDims3 = reqDims2;
					reqDims3.push_back(k);
					
					if (i > dims || j > dims || k > dims) {
						doMarginalFail(node, reqDims3);
					}
			}
				
		}
	}
	
}

void doMarginal(const RealMappedSPnode& node, 
				const std::vector < int >& reqDims) {
		
		RealMappedSPnode marg = node.makeMarginalised(reqDims);
		cout << "After marginalising on ( ";
		ostream_iterator< int > out_it (cout," ");
		copy ( reqDims.begin(), reqDims.end(), out_it );
		cout << ") ";

		printTestTree(marg);
}

void doMarginalFail(const RealMappedSPnode& node, 
				const std::vector < int >& reqDims) {
	try {
		cout << "Try to marginalise on ( ";
		ostream_iterator< int > out_it (cout," ");
		copy ( reqDims.begin(), reqDims.end(), out_it );
		cout << ") " << endl;

		RealMappedSPnode marg = node.makeMarginalised(reqDims);
		throw std::logic_error("Should not be able to do that");
	}
	catch (std::invalid_argument& e) {
		cout << "invalid_argument exception :\n" << e.what() << endl;
	}
}		
		
void testingIntervals()
{

    try {
        // testing for intervals

        cout << "\nTesting for ranges as intervals\n" << endl;

        int dims = 1;
        ivector pavingBox(dims);
        interval pavingInterval(-4,4);
        for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;

		MappedSPnode<interval> additiveIdentity(pavingBox); 
		MappedSPnode<interval> multiplicativeIdentity(pavingBox); 
       
        MappedSPnode<interval> nodeOne(pavingBox); // make a MappedSPnode object
        MappedSPnode<interval> nodeTwo(pavingBox); // make another MappedSPnode object

        // split the nodes to particular shapes
        nodeOne.splitToShape("1,1"); // split it
        nodeTwo.splitToShape("2,2,1"); // split it

        // allocate ranges for nodeOne
        cxsc::interval  rangeOne1(1.0,5.0);
        cxsc::interval rangeOne2(1.5,3.0);
        cxsc::interval rangeOne3(2.5,4.0);

        std::vector<interval> rangesOne;
        rangesOne.push_back(rangeOne1);
        rangesOne.push_back(rangeOne2);
        rangesOne.push_back(rangeOne3);

        nodeOne.allocateRanges(rangesOne, 0);


        // allocate ranges for nodeTwo
        cxsc::interval rangeTwo1(1.0,6.0);
        cxsc::interval rangeTwo2(2.0,5.0);
        cxsc::interval rangeTwo3(2.5,3.5);
        cxsc::interval rangeTwo4(3.0,4.5);
        cxsc::interval rangeTwo5(1.5,5.5);

        std::vector<interval> rangesTwo;
        rangesTwo.push_back(rangeTwo1);
        rangesTwo.push_back(rangeTwo2);
        rangesTwo.push_back(rangeTwo3);
        rangesTwo.push_back(rangeTwo4);
        rangesTwo.push_back(rangeTwo5);

        nodeTwo.allocateRanges(rangesTwo, 0);

        cxsc::interval  rangeAddId(0.0,0.0);
         std::vector<interval> rangesAddId;
        rangesAddId.push_back(rangeAddId);
		additiveIdentity.allocateRanges(rangesAddId, 0);

		cxsc::interval  rangeMultId(1.0,1.0);
        std::vector<interval> rangesMultId;
        rangesMultId.push_back(rangeMultId);
		multiplicativeIdentity.allocateRanges(rangesMultId, 0);

        mappedAdditionTest(additiveIdentity, nodeOne, nodeTwo);
        mappedSubtractTest(additiveIdentity, nodeOne, nodeTwo);
        mappedProductTest(multiplicativeIdentity, nodeOne, nodeTwo);
		mappedDivisionTest(multiplicativeIdentity, nodeOne, nodeTwo);



        interval multiplier(1.0,4.0); // supply a multiplier of the right type

        mappedScalarMultTest1(nodeOne, nodeTwo, multiplier);
		mappedScalarDivTest1(nodeOne, nodeTwo, multiplier);


    }
    catch (exception& e) {
		std::cerr << "Error in testingIntervals:\n" << e.what() << std::endl;
            throw;
    }


} // end of test for intervals

rvector make2Drvector(real r1, real r2)
{
    rvector rv(2);
    rv[1] = r1;
    rv[2] = r2;
    return rv;

}



void testingRvectors()
{

   try {

        // testing for rvectors

        cout << "\nTesting for ranges as rvectors\n" << endl;

        int dims = 2;
        ivector pavingBox(dims);
        interval pavingInterval(-4,4);
        for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;


        MappedSPnode<rvector> additiveIdentity(pavingBox); 
		MappedSPnode<rvector> nodeOne(pavingBox); // make a MappedSPnode object
        MappedSPnode<rvector> nodeTwo(pavingBox); // make another MappedSPnode object

        // split the nodes to particular shapes
        nodeOne.splitToShape("1,1"); // split it
        nodeTwo.splitToShape("2,2,1"); // split it

        // allocate ranges for nodeOne
        cxsc::rvector  rangeOne1= make2Drvector(2.0,3.5);
        cxsc::rvector rangeOne2 = make2Drvector(1.5,3.0);
        cxsc::rvector rangeOne3= make2Drvector(2.5,4.0);

        std::vector<rvector> rangesOne;
        rangesOne.push_back(rangeOne1);
        rangesOne.push_back(rangeOne2);
        rangesOne.push_back(rangeOne3);

        nodeOne.allocateRanges(rangesOne, 0);


        // allocate ranges for nodeTwo
        cxsc::rvector rangeTwo1 = make2Drvector(2.0,4.5);
        cxsc::rvector rangeTwo2 = make2Drvector(2.5,4.0);
        cxsc::rvector rangeTwo3 = make2Drvector(2.0,3.5);
        cxsc::rvector rangeTwo4 = make2Drvector(3.0,4.5);
        cxsc::rvector rangeTwo5 = make2Drvector(1.5,5.0);

        std::vector<rvector> rangesTwo;
        rangesTwo.push_back(rangeTwo1);
        rangesTwo.push_back(rangeTwo2);
        rangesTwo.push_back(rangeTwo3);
        rangesTwo.push_back(rangeTwo4);
        rangesTwo.push_back(rangeTwo5);

        nodeTwo.allocateRanges(rangesTwo, 0);

        cxsc::rvector rangeAddId = make2Drvector(0.0,0.0);
        std::vector<rvector> rangesAddId;
        rangesAddId.push_back(rangeAddId);
		additiveIdentity.allocateRanges(rangesAddId, 0);


        mappedAdditionTest(additiveIdentity, nodeOne, nodeTwo);
        mappedSubtractTest(additiveIdentity, nodeOne, nodeTwo);
        //mappedProductTest(multiplicativeIdentity, nodeOne, nodeTwo);
		cout << "Multiplication with rvectors gives dot product as a real so no multiplication test" << endl;
		//mappedDivisionTest(multiplicativeIdentity, nodeOne, nodeTwo);
		cout << "No division operation defined for rvectors so no division test" << endl;
       
        rvector multiplier = make2Drvector(1.0,4.0); // supply a multiplier of the right type

        //mappedScalarMultTest1(nodeOne, nodeTwo, multiplier);
		//mappedScalarDivTest1(nodeOne, nodeTwo, multiplier);
		cout << "No division operation defined for rvectors so no division test" << endl;
	   

	}
    catch (exception& e) {
		std::cerr << "Error in testingRvectors:\n" << e.what() << std::endl;
            throw;
    }



} // end of test for rvectors






void testingBools()
{
	
	try {
        // testing for bools
		
		cout << "\nTesting for ranges as bools\n" << endl;
		
		
		bool bf = false;
		bool bt = true;
		
		#if 0
			cout << "bf = " << bf << endl;
			cout << "bt = " << bt << endl;
			cout << "bf + ft = " << (bf + bf) << endl;
			cout << "bf + bt = " << (bf + bt) << endl;
			cout << "bt + bt = " << (bt + bt) << endl;
			cout << "bf - bf = " << (bf - bf) << endl;
			cout << "bf - bt = " << (bf - bt) << endl;
			cout << "bt - bf = " << (bt - bf) << endl;
			cout << "bt - bt = " << (bt - bt) << endl;
			cout << "bf * bf = " << (bf * bf) << endl;
			cout << "bf * bt = " << (bf * bt) << endl;
			cout << "bt * bt = " << (bt * bt) << endl;
			cout << "bf / bf = " << (bf / bf) << endl;
			cout << "bf / bt = " << (bf / bt) << endl;
			//cout << "bt / bf = " << (bt / bf) << endl;
			cout << "bt / bt = " << (bt / bt) << endl;
		#endif
		
        bool dims = 2;
        ivector pavingBox(dims);
        interval pavingInterval(-4,4);
        for(int k=1; k <= dims; k++) pavingBox[k] = pavingInterval;

		MappedSPnode<bool> additiveIdentity(pavingBox); 
		MappedSPnode<bool> multiplicativeIdentity(pavingBox); 
        MappedSPnode<bool> nodeOne(pavingBox); // make a MappedSPnode object
        MappedSPnode<bool> nodeTwo(pavingBox); // make another MappedSPnode object
		MappedSPnode<bool> nodeThree(pavingBox); // make another MappedSPnode object


        // split the nodes to particular shapes
        nodeOne.splitToShape("1,2,2"); // split it
        nodeTwo.splitToShape("2,2,1"); // split it
		nodeThree.splitToShape("2,2,1"); // split it

		
        // allocate ranges
		bool range0 = false;
        bool range1 = true;
        
        std::vector<bool> rangesOne;
        rangesOne.push_back(range1);
        rangesOne.push_back(range0);
        rangesOne.push_back(range1);
		rangesOne.push_back(range1);
		rangesOne.push_back(range0);

        nodeOne.allocateRanges(rangesOne, 0);

		std::vector<bool> rangesTwo;
        rangesTwo.push_back(range1);
        rangesTwo.push_back(range1);
        rangesTwo.push_back(range1);
        rangesTwo.push_back(range0);
        rangesTwo.push_back(range0);

        nodeTwo.allocateRanges(rangesTwo, 0);
		
		std::vector<bool> rangesThree;
        rangesThree.push_back(range1);
        rangesThree.push_back(range1);
        rangesThree.push_back(range1);
        rangesThree.push_back(range0);
        rangesThree.push_back(range1);
		
        nodeThree.allocateRanges(rangesThree, 0);

		std::vector<bool> rangesAddId;
		rangesAddId.push_back(range0);
		
		additiveIdentity.allocateRanges(rangesAddId);

		std::vector<bool> rangesMultId;
		rangesMultId.push_back(range1);
		
		multiplicativeIdentity.allocateRanges(rangesMultId);

		//mappedAdditionTest(additiveIdentity, nodeOne, nodeTwo);
		mappedAdditionTest(additiveIdentity, nodeOne, nodeThree);
        //mappedSubtractTest(additiveIdentity, nodeOne, nodeTwo);
		mappedSubtractTest(additiveIdentity, nodeOne, nodeThree);
        //mappedProductTest(multiplicativeIdentity, nodeOne, nodeTwo);
		mappedProductTest(multiplicativeIdentity, nodeOne, nodeThree);
		//mappedDivisionTest(multiplicativeIdentity, nodeOne, nodeTwo);
		cout << "Don't do division tests - will get floating point exceptions if we divide by false" << endl;


        {
			bool multiplier = true; // supply a multiplier of the right type

			mappedScalarMultTest1(nodeOne, nodeTwo, multiplier);
			//mappedScalarDivTest1(nodeOne, nodeTwo, multiplier);
		}

		{
			bool multiplier = false; // supply a multiplier of the right type

			mappedScalarMultTest1(nodeOne, nodeTwo, multiplier);
			//mappedScalarDivTest1(nodeOne, nodeTwo, multiplier);
		}
		
		/*
		MappedSPnode<bool> nodeThree(pavingBox); 
        MappedSPnode<bool> nodeFour(pavingBox); 
		MappedSPnode<bool> nodeFive(pavingBox);
		
		// split the nodes to particular shapes
        nodeThree.splitToShape("2,3,3,1"); // split it
		nodeFour.splitToShape("2,3,3,1"); // split it
        nodeFive.splitToShape("2,3,3,2,3,3"); // split it

        std::vector<bool> rangesThree;
        rangesThree.push_back(range3);
        rangesThree.push_back(range3);
        rangesThree.push_back(range3);
		rangesThree.push_back(range3);
		rangesThree.push_back(range2);
		rangesThree.push_back(range4);
		rangesThree.push_back(range3);

        nodeThree.allocateRanges(rangesThree, 0);

		std::vector<bool> rangesFour;
        rangesFour.push_back(range3);
        rangesFour.push_back(range3);
        rangesFour.push_back(range3);
		rangesFour.push_back(range3);
		rangesFour.push_back(range4);
		rangesFour.push_back(range2);
		rangesFour.push_back(range3);

        nodeFour.allocateRanges(rangesFour, 0);
		
		std::vector<bool> rangesFive;
        rangesFive.push_back(range3);
        rangesFive.push_back(range3);
        rangesFive.push_back(range3);
		rangesFive.push_back(range3);
		rangesFive.push_back(range4);
		rangesFive.push_back(range2);
		rangesFive.push_back(range3);
        rangesFive.push_back(range3);
		rangesFive.push_back(range3);
		rangesFive.push_back(range2);
		rangesFive.push_back(range4);
		

        nodeFive.allocateRanges(rangesFive, 0);
		
		mappedMinimiseTest(nodeThree, nodeFour);
		mappedMinimiseTest(nodeThree, nodeFive);
		*/

    }
    catch (exception& e) {
		std::cerr << "Error in testingBools:\n" << e.what() << std::endl;
            throw;
    }



} // end of test for bools


