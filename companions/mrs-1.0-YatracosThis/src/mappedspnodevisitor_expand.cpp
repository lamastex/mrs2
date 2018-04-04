/*
* Copyright (C) 2010 Jenny Harlow
*/


/*! \file
\brief definitions for MappedSPnodeVisitorExpander

*/

#include "mappedspnodevisitor_expand.hpp"
#include "nodecompobjmapped.hpp"
#include <assert.h>


using namespace cxsc;
using namespace std;
using namespace subpavings;


//gat41
// a class for comparison
class MyCompare
{
    const NodeCompObjMapped& myNC;

    public:
    MyCompare(const NodeCompObjMapped& nc) : myNC(nc) {}

    bool operator()   (const SPnode * const lhs,
                      const SPnode * const rhs) const

    { return myNC(lhs, rhs); }
};




MappedSPnodeVisitorExpand::MappedSPnodeVisitorExpand(MappedFobj& f, real tol)
			: SPnodeVisitor(), fobj(f), tolerance(tol) {}


void MappedSPnodeVisitorExpand::visit(SPnode * mspn)
{

	//std::cout << "in visit, for " << mspn->getNodeName() << std::endl;

	// check if we need to split
	ivector box = mspn->getBox();
	//std::cout << "this box is " << box << " and has volume " << mspn->nodeVolume() << std::endl;
	interval thisRange = fobj(box);
	//std::cout << "this range is " << thisRange << std::endl;
	real thisMidImage = fobj.imageMid(box);
	//std::cout << "this midImage is " << thisMidImage << std::endl;
	
	// split if so and then visit children
	
cout << diam(thisRange) << '\t' << thisRange << endl;
	// wants the simple function approx using the mid image to be within epsilon of the worst-case scenario
	if (max(Sup(thisRange) - thisMidImage, thisMidImage - Inf(thisRange)) > tolerance) {
		//std::cout << "expanding" << std::endl;
		mspn->nodeExpand();

		//visit the children
		//std::cout << "check children" << std::endl;
		(mspn->getLeftChild())->accept(*this);
		(mspn->getRightChild())->accept(*this);
	}
}

real MappedSPnodeVisitorExpand::tellMe(SPnode * mspn)
{
	return fobj.imageMid(mspn->getBox());
}

/*
real MappedSPnodeVisitorExpand::getSPArea(SPnode * mspn)
{
	std::cout << "in visit, for " << mspn->getNodeName() << std::endl;

	ivector box = mspn->getBox();
	interval thisRange = fobj(box);
	real area = Volume(box) * diam(thisRange);
	
	return area;
}
*/

//gat41
bool MappedSPnodeVisitorExpand::priorityVisit(SPnode * mspn, size_t critLeaves, std::vector<real>& eps)
{
	 bool retValue = false;
	 gsl_rng * rgsl = NULL;
	  // set up a random number generator for uniform rvs
	  const gsl_rng_type * tgsl;
	  // set the library variables *gsl_rng_default and
	  // gsl_rng_default_seed to default environmental vars
	  gsl_rng_env_setup();
	  tgsl = gsl_rng_default; // make tgsl the default type
	  rgsl = gsl_rng_alloc (tgsl); // set up with default seed
	  retValue = priorityVisit(mspn, critLeaves, rgsl, eps);
	  gsl_rng_free (rgsl);
	 return retValue;
}

//gat41
bool MappedSPnodeVisitorExpand::priorityVisit(SPnode * mspn, size_t critLeaves, gsl_rng * rgsl, std::vector<real>& eps)
{    
   //cout << "Calling priority visit: " << endl;

   bool cancontinue = false;
   size_t numNodes = 0;
   //real normConst = 0.0;
   
   //comparison function
   CompSPArea compTest(fobj);
   // a multiset for the queue (key values are not necessarily unique)
   multiset<SPnode*, MyCompare> pq((MyCompare(compTest)));

   //cout << "get fobj of box " << endl;
   ivector box = mspn->getBox();
   interval thisRange = fobj(box);
   
   //cout << "get fobj of mid point " << endl;
   //real thisMidImage = fobj.imageMid(box);
   //normConst = (mspn->nodeVolume() * thisMidImage);

   interval RiemannDiff = (mspn->nodeVolume()) * (thisRange);
   real MaxEps = diam(RiemannDiff);
   //normalize MaxEps by normConst
   interval TotEps = interval(MaxEps); 
   real midTotEps = mid(TotEps);
   pq.insert(mspn);
   mspn->collectRange(*this);
   numNodes++;
	// this is optional (collecting midTotEps)
   eps.push_back(midTotEps);
   
   cancontinue = (!pq.empty());
   if(!cancontinue) {
			std::cout << "No splittable leaves to split - aborting" << std::endl;
    } 

   // split until have desired number of leaf nodes
   // we only put splittable nodes into the set, so we don't have to check
   // that they are splittable when we take them out
   while (cancontinue && (numNodes < critLeaves) && (midTotEps > tolerance))
   {
		SPnode* largest = *(pq.rbegin ()); // the last largest in the set
		SPnode* chosenLargest;
		
		// find if there are any more equal to largest around
		multiset<SPnode*, MyCompare>::iterator mit;
		pair<multiset<SPnode*, MyCompare>::iterator,
			 multiset<SPnode*, MyCompare>::iterator> equalLargest;

		equalLargest = pq.equal_range(largest); // everything that = largest
		size_t numberLargest = pq.count(largest); // number of =largest

		if (numberLargest > 1) {
			 // draw a random number in [0,1)
			 double rand = gsl_rng_uniform(rgsl);
			 real sum = 0.0;

			 // random selection of the =largest node to chose
			 for (mit=equalLargest.first; mit!=equalLargest.second; ++mit) {
				  sum += 1.0/(1.0*numberLargest);
				  if (rand < sum) {
						break;
				  }
			 }
			 chosenLargest = *(mit); // the chosen largest in the set
			 pq.erase(mit);// take the iterator to chosen largest out of the set
			 numNodes--; //check if pq.size() == numNodes;
			 assert(numNodes == pq.size());
		}
		else {
			 chosenLargest = *(pq.rbegin ()); // the only largest
			 multiset<SPnode*, MyCompare>::iterator it = pq.end();
			 it--;
			 pq.erase(it);// take this largest out of the set
			 numNodes--; //check if pq.size() == numNodes
			 assert(numNodes == pq.size());
		}

		// split the biggest one
		//cout << "---------------splitting " << chosenLargest->getNodeName() << "-----" << endl;
		
		//cout << "get fobj of box " << endl;
		box = chosenLargest->getBox();
		thisRange = fobj(box);
		//thisRange = thisRange/normConst; //normalize the heights
		RiemannDiff = (chosenLargest->nodeVolume()) * (thisRange);
		MaxEps = diam(RiemannDiff);
		
		// now update TotEps
		TotEps = TotEps - interval(MaxEps);
		
		 // now update the normalizing constant since a node is removed
		 //thisMidImage = fobj.imageMid(box);
		 //real remNormConst = thisMidImage * (chosenLargest->nodeVolume());
		 //cout << "to remove: " << chosenLargest->getNodeName() << "\t" << remNormConst << endl;
		 //assert(normConst >= remNormConst);
		 //normConst = normConst - remNormConst;
		 
		 //now split
		chosenLargest->nodeExpand();
		
		//name the new children
		(chosenLargest->getLeftChild())->recursiveRename();
		(chosenLargest->getLeftChild())->collectRange(*this);

		(chosenLargest->getRightChild())->recursiveRename();
		(chosenLargest->getRightChild())->collectRange(*this);
		
		// insert these nodes into the priority queue
		pq.insert(chosenLargest->getLeftChild());
		pq.insert(chosenLargest->getRightChild());
		numNodes = numNodes+2;
		assert(pq.size() == numNodes);
		
		// update normConst with the addition of the left and right child nodes
		//normConst = normConst 
		//		+ (((chosenLargest->getLeftChild())->nodeVolume())* 
		//		   (fobj.imageMid((chosenLargest->getLeftChild())->getBox()))
		//		  );

		//normConst = normConst 
		//		+ (((chosenLargest->getRightChild())->nodeVolume())* 
		//		   (fobj.imageMid((chosenLargest->getRightChild())->getBox()))
		//		  );

		//normalize TotEps to the current normConst

		
		// get the epsilon for the left child node
		//cout << "get fobj of box and mid point " << endl;
		MaxEps = diam(((chosenLargest->getLeftChild())->nodeVolume()) * (fobj((chosenLargest->getLeftChild())->getBox())));
		TotEps = TotEps + interval(MaxEps);

		// get the epsilon for the right child node
		MaxEps = diam(((chosenLargest->getRightChild())->nodeVolume()) * (fobj((chosenLargest->getRightChild())->getBox())));
		TotEps = TotEps + interval(MaxEps);

		//get the mid point
		midTotEps = mid(TotEps);
		//optional
		//if (normConst != 0.0) 
		eps.push_back(midTotEps); 

		cancontinue = (!pq.empty());
		if (!cancontinue)
			 std::cout << "Terminated splitting: no splittable nodes left"
				  << std::endl;
	}
    return (cancontinue);
}

// end of MappedSPnodeVisitorExpand class
