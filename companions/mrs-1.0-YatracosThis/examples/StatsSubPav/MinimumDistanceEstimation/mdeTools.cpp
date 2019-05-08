/*! \file
\brief Common routines for MDE

 */
 
#include "mdeTools.hpp"
#include <vector> //std::vector
#include <map>	//std::multimap
#include <algorithm>    // std::sort
#include <iostream>

using namespace std;

// return a vector of the top k indices
void topk(std::vector<double> a, std::vector<int> & indtop, size_t k)
{
	std::multimap<double, size_t> m; // mapping from value to its index
	std::vector<double>::iterator it;

	for (it = a.begin(); it != a.end(); ++it)
		m.insert(std::make_pair(*it, it - a.begin()));

	std::multimap<double, size_t>::iterator itm; // mapping from value to its index
	size_t indx=0;
	double val =0;
	for (itm = m.begin(); itm != m.end(); ++itm){
		//std::cout << itm->first <<" , "<< itm->second << std::endl;
		//if (itm->first != val) { 
			indtop.push_back(itm->second);
			indx++;
		//}
		//val = itm->first;
		if ( indx == k) break; 
	}	
	
	std::sort(indtop.begin(), indtop.end());

} // end of topk

// return a vector of the sequence of leaf nodes for collation
void getSequence(std::vector<int> & seq, int temp, int critLeaves, int increment)
{
	while ( temp < critLeaves) {
		temp += increment;
	 	seq.push_back(temp); 
	 }
	std::sort(seq.begin(), seq.end());
	seq.erase( unique( seq.begin(), seq.end() ), seq.end() );
} // end of getSequence
