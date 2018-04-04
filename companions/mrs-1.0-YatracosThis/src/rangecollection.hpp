/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
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

/*!/ \file
\brief RangeCollection definitions.
*/

#ifndef __RANGECOLLECTION_HPP__
#define __RANGECOLLECTION_HPP__

// put it all in the header for the moment and sort out the template issues later

#include "rangecollectionexception.hpp"

#include "cxsc.hpp"

// to use std input/output
#include <iostream>

// to use exceptions
//#include <exception>

//to use accumulate
#include <numeric>
#include <algorithm>
#include <functional>

#include <vector>



//using namespace std;


class Vec;

namespace subpavings {

    // templatised hull operator function, for most types
    template<typename T>
    T hullOperator(const T& x, const T& y)
    {
        return x > y ? x : y;
    }

    // specialisation for intervals
    template<>
    cxsc::interval hullOperator<cxsc::interval>(const cxsc::interval& x, const cxsc::interval& y);

    // specialisation for interval vectors
    template<>
    cxsc::ivector hullOperator<cxsc::ivector>(const cxsc::ivector& x, const cxsc::ivector& y);

    // specialisation for vectors
    template<>
    Vec hullOperator<Vec>(const Vec& x, const Vec& y);

/*
    template<typename T>
    struct transformMult : public binary_function<int, T, T > {
      T operator() (const int a, const T b) const {return (a*b);}
    };
*/


/*! \brief A class for range collection objects.

A RangeCollectionClass is an implementation for collections of ranges used by
MappedSPnode objects.  The type of the range object is the template parameter T.

The underlying container should be some kind of sequence container.

This implementation uses an stl::vector as the implementation container.

*/


// implementation of RangeCollection using vector for container
template <typename T>
    class RangeCollectionClass {

        typedef std::vector<T> RangeCollectionType;
        typedef typename std::vector<T>::iterator RangeCollectionItr;
        typedef typename std::vector<T>::const_iterator RangeCollectionConstItr;

        private:
            RangeCollectionType container;

        public:

        // no-argument constructor,
        RangeCollectionClass() {}

        // constructor with initial value
        RangeCollectionClass(const T& range)
        {
            pushBack(range);
        }

        // copy constructor
        RangeCollectionClass(const RangeCollectionClass<T>& other)
        : container(other.container) {}

        // Copy assignment operator
        RangeCollectionClass<T>& operator=(const RangeCollectionClass<T>& rhs)
        {
            container = rhs.container;
            return *this;
        }

        // get the size of the RangeCollection, ie number of elements
        size_t getSize() const
        {
            return container.size();
        }

        // get whether container is empty
        bool isEmpty() const
        {
            return container.empty();
        }

        // add to the back of the container
        void pushBack(const T& toAdd)
        {
            container.push_back(toAdd);
        }

        // reduce the rangeCollection to one value, the sum of the present values
        // Fudge to find an initial value and use std::accumulation
        void addReduce()
        {
            if (!container.empty()) {

                // use std::accumulate for this but from second value in collection on
                // so that we can use first value as intial one
                // avoids having to know the equivalent of '0' for T

                T accumulation = *(container.begin());

                if (container.size() > 1) {
                    accumulation = std::accumulate(++container.begin(),
                                    container.end(), accumulation);
                }

                container.clear();
                container.push_back(accumulation);
            }
        }

        // reduce the rangeCollection to one value
        // as the product of the present values
        void productReduce()
        {
            if (!container.empty()) {

                RangeCollectionItr rit = container.begin();
                T product = *(rit);

                if (container.size() > 1) {
                    rit++;
                    for (;
                            rit < container.end(); rit++) {
                        // don't use *= in case type T does not support it
                        product = product * (*rit);
                    }
               }

                container.clear();
                container.push_back(product);
            }
        }

        // return absolute value of range collection
        T absValue()
        {
            T accumulation; // default constructor

            if (!container.empty()) {

                accumulation = abs(*(container.begin()));


                if (container.size() > 1) {
                    RangeCollectionItr rit = container.begin();

                    rit++;
                    for (;
                            rit < container.end(); rit++) {
                        // don't use *= in case type T does not support it
                        accumulation = accumulation + abs(*rit);
                    }
                }

            }
            return accumulation;
        }

        // incorporate another rangeCollection into this rangeCollection
        // by taking the pair-wise interval hull
        void hullCollection(const RangeCollectionClass<T>& other)
        {

            size_t i = 0;

            size_t n = other.getSize();

            if (isEmpty()) { // nothing here already
                // copy the other's container for this
                container = other.container;
            }
            else { // has something in

                // we have to make sure the rangeCollection for this matches
                // that of the children
                // number of elements in this rangeCollection should = child
                if (getSize() != n) {
                    throw RangeCollectionException("rangeCollections do not match");
                }

                // store current container temporarily
                RangeCollectionType _temp = container;

                container.clear();

                // put into this container the hull operator of
                // the pairs of values from the other's container and the this's container
                for (i=0; i < n; i++) {
                    container.push_back(
                                hullOperator((other.container)[i],
                                _temp[i]));
                }
            }
        }
/*
        // multiply each element in the collection by mult
        void intScalarMult(int mult)
        {
            transform (container.begin(), container.end(), container.begin(),
                        bind1st(transformMult<T>(),mult));
        }
*/
        // multiply each element in the collection by mult when mult is of type T
        void scalarMult(const T& mult)
        {
            //transform (container.begin(), container.end(), container.begin(),
            //            bind1st(transformMult<T>(),mult));

            RangeCollectionItr rit = container.begin();
            for (RangeCollectionItr rit = container.begin();
                            rit < container.end(); rit++) {
                        // don't use *= in case type T does not support it
                        *rit = mult * (*rit);
            }
        }

        // replace contents of this's container with contents of lhs and rhs containers
        void hullCollection(const RangeCollectionClass<T>& lhs,
                                const RangeCollectionClass<T>& rhs)
        {
            container = lhs.container;
            hullCollection(rhs);
        }


        // combines contents of two other collections into this one
        // replaces contents of this
        void combineCollection(const RangeCollectionClass<T>& lhs,
                                const RangeCollectionClass<T>& rhs)
        {
            // now give this node a range collection which contains all the
            // elements in both lhs and rhs range collections, going l->r
            container.clear(); // make sure current range collection empty
            container.insert(container.begin(),
                                            lhs.container.begin(),
                                            lhs.container.end());
            container.insert(container.end(),
                                            rhs.container.begin(),
                                            rhs.container.end());
        }

        // replaces contents of this with the difference of the given collections
        void subtractCollection(const RangeCollectionClass<T>& lhs,
                                const RangeCollectionClass<T>& rhs)
        {
            size_t i = 0;
            size_t n = rhs.getSize();
            size_t m = lhs.getSize();
            // we have to make sure the rangeCollection for this matches
            // that of the other
            // number of elements in this rangeCollection should = child
            if (m < n) {
                throw RangeCollectionException("rangeCollections do not match");
            }
            container.clear();
            // put into this container the difference of
            // the pairs of values from the two containers
            for (i=0; i < n; i++) {
                container.push_back(
                            lhs.container[i] - rhs.container[i]);
            }
            for (i = n; i < m; i++) {
                container.push_back(lhs.container[i]);
            }
        }

        // output contents
        std::ostream& outputTabs(std::ostream &os) const
        {
            RangeCollectionConstItr cit;
            for (cit = container.begin(); cit < container.end(); cit++) {
                os << "\t" << (*cit);
            }
        }

		//gloria added the following functions because she doesn't know how to
		//access the elements of the container from the outside. 
		//get the ranges in the range collection and output to a vector
		//can only do this for rangeCollection.getSize() == 1
		//the weight is volume times height of this
		void getWeights(std::vector<double> & WeightsVector, 
								std::vector<interval> & WeightsInt, double boxVol)
		{
			RangeCollectionItr cit;
			if (container.size() == 1) {
				for (cit = container.begin(); cit < container.end(); cit++) {
						
						interval vol = interval(boxVol);
						interval height = interval(*cit);
						//std::cout << vol << "\t" << height << std::endl;
						
						interval area = vol*height;
						WeightsInt.push_back(area);
						WeightsVector.push_back(_double(mid(area)));
						//i know this is horrible but i don't know how else to do this
						//for now and i don't have time to improve this.. sorry..
						//mainly i need a double for the gsl_ran_discrete_preproc
						//which only takes in doubles, hence the cast.
				}
			}
      }
		
		//gloria added the following functions because she doesn't know how to
		//access the elements of the container from the outside. 
		//get the ranges in the range collection and output to a vector
		//can only do this for rangeCollection.getSize() == 1
		//get the height
		void getHeight(std::vector<real> & HeightsVector)
		{
			RangeCollectionItr cit;
			if (container.size() == 1) {
				for (cit = container.begin(); cit < container.end(); cit++) {
						HeightsVector.push_back((*cit));
				}
			}
      }
		
		//gloria added the following function because she doesn't know how to
		//access the elements of the container from the outside. 
		//get the ranges in the range collection and output to a vector
		//can only do this for rangeCollection.getSize() == 1
		real getNodeIAE(double boxVol)
		{
			RangeCollectionItr cit;
			real IAE;
			if (container.size() == 1) {
				for (cit = container.begin(); cit < container.end(); cit++) {
						real diff = (*cit);
						if (diff < 0) { diff = diff * (-1.0); } //to make it positive
						IAE = (boxVol*diff);
						//std::cout << diff << "\t" << boxVol << std::endl;
				}
			}
			return IAE;
      }
      
      //get the ranges in the range collection and output to a vector
		//can only do this for rangeCollection.getSize() == 1
		real getNodeIAERegHist(double boxVol, double height)
		{
			RangeCollectionItr cit;
			real IAE;
			if (container.size() == 1) {
				for (cit = container.begin(); cit < container.end(); cit++) {
						real diff = (*cit) - height;
						if (diff < 0) { diff = diff * (-1.0); } //to make it positive
						IAE = (boxVol*diff);
				}
			}
			return IAE;
      }
      
      //get the area of this node
      real getNodeArea(double boxVol)
      {
			RangeCollectionItr cit;
			real area = 0.0;
			//std::cout << "get node area: " << std::endl;
			if (container.size() == 1) {
				for (cit = container.begin(); cit < container.end(); cit++) {
						//std::cout << (*cit) << "\t" << boxVol << std::endl;
						area = (*cit)*boxVol;
				}
			}
			return area;
		}
			
		real normNodeHeight(double totalArea)
		{
			RangeCollectionItr cit;
			real newHeight;
			if (container.size() == 1) {
				for (cit = container.begin(); cit < container.end(); cit++) {
						newHeight = (*cit)/totalArea;
				}
			}
			return newHeight;
		} // end of normNodeHeight
				

	//gat41
	//clear the current container
   void clearAll()
   {
		if (!container.empty()) {
			 container.clear();
		}
   }



	}; // end class RangeCollectionClass
} // end namespace subpavings

#endif
