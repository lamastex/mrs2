/*
* Copyright (C) 2007, 2008, 2009 Raazesh Sainudiin
* Copyright (C) 2009, 2010, 2011 Jennifer Harlow
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
\brief RangeCollectionHist declarations
*/

#ifndef __RANGECOLLECTION_HIST_HPP__
#define __RANGECOLLECTION_HIST_HPP__

#include "sptypes.hpp"
#include "cxsc.hpp"

#include <string>
#include <vector>

    

/*! \brief A class for range collections for histograms.

The underlying container should be some kind of sequence container.

This implementation uses an stl::vector as the implementation container


*/


// implementation of RangeCollections for Histograms using vector for container
namespace subpavings {
	
	class RangeCollectionHist {
		
        public:
		
			/*! \brief Definition for the type for the values 
			accessible from a rangeCollection.*/ 
			typedef cxsc::real 	RangeCollectionHistValueType;
			
			/*! \brief Definition for the type for iterating through
			the values in a rangeCollection.*/ 
			typedef std::vector < RangeCollectionHistValueType >::const_iterator 
									RangeCollectionHistItr;

            /*! \brief No-argument constructor.*/
			RangeCollectionHist();

			/*! \brief Constructor initialised with a container
			of values for the rangeCollection.*/
			explicit RangeCollectionHist(const RealVec& range);

			/*! \brief Constructor initialised with a single
			value for the rangeCollection.*/
			explicit RangeCollectionHist(const cxsc::real& r);

			/*! \brief Copy constructor.*/
			RangeCollectionHist(const RangeCollectionHist& other);
			
			
			/*! \brief Copy assignment operator. */
			RangeCollectionHist& operator=(RangeCollectionHist rhs);

			/*! \brief Return an iterator to the start of the collection.
			
			The idiom for iterating through all the valuse in a 
			rangeCollection rc is:
			
			for (RangeCollectionHistItr it = rc.begin();
					it < rc.end(),
					++it) {
				...
			}
			*/
			RangeCollectionHistItr begin() const;
			
			/*! \brief Return an iterator to the end of the collection.*/
			RangeCollectionHistItr end() const;

			/*! \brief Get the size of this.
			 * 
			 Size is number of elements in this.*/
			std::size_t getSize() const;

			/*! \brief Get whether the rangeCollection is empty.
			
			\return true if this has no elements, false otherwise.*/
			bool isEmpty() const;
			
			/*! \brief Get a copy of the elements of this as a 
			collection of reals.*/
			RealVec getRangeCollectionValues() const;
			
			/*! \brief Get the total of all the elements in this.*/
			cxsc::real getTotal() const;
			
			/*! \brief Get the total of the absolute values of all the
			elements in this.*/
			cxsc::real getAbsTotal() const;

			/*! \brief Get the average of all the elements in this.
			
			Throws an UnfulfillableRequest_Error if this is empty.*/
			cxsc::real getAverageValue() const;

			
			/*! \brief Make a rangeCollection containing just the average of
			this.
			
			Throws an UnfulfillableRequest_Error if this is empty.*/
			RangeCollectionHist makeAverageRangeCollection() const;
			
			/*! \brief Make a rangeCollection containing the absolute 
			 differences of each element to their average. 
			
			Throws an UnfulfillableRequest_Error if this is empty.*/
			RealVec getAbsDiffsToAverageRangeCollection() const;
			
			/*! \brief Reduce the rangeCollection to the sum of the
			 present values.*/
			void reduce();
			
			/*! \brief Add another the rangeCollection to this so that 
			each element += parallel element in \a rhs.
			
			Throws a std::logic_error exception if this and \a rhs 
			are different sizes.*/
			void parallelAdd(const RangeCollectionHist& rhs);


			/*! \brief Append the contents of another rangeCollection 
			to this.*/
			RangeCollectionHist& operator+=(const RangeCollectionHist& rhs);
			
			/*! \brief Get the results of appending the contents of another rangeCollection 
			to this.*/
			const RangeCollectionHist operator+(const RangeCollectionHist& rhs) const;
			
			/*! \brief Multiply every element in this by \a mult.*/
			RangeCollectionHist& operator*=(const cxsc::real& mult);
			
			/*! \brief Get the results of multipling every element in this by \a mult.*/
			const RangeCollectionHist operator*(const cxsc::real& mult) const;
			
			/*! \brief Divide every element in this by \a div.
			
			Throws a std::invalid_argument exception if \a div == 0.0.*/
			RangeCollectionHist& operator/=(const cxsc::real& div);
			
			/*! \brief Get the results of dividing every element in this by \a div.
			
			Throws a std::invalid_argument exception if \a div == 0.0.*/
			const RangeCollectionHist operator/(const cxsc::real& div) const;

			/*! \brief Get whether the values in this are all positive.*/ 
			bool checkRangeCollectionAllPositive() const;
			
			/*! @name Output this to an out stream with values 
			separated by tabs.
			
			\param os is the stream to sent the output to.
			\param prec is the precision to use in the output.*/
			//@{
			std::ostream& outputTabs(std::ostream &os,
										int prec) const;
			std::ostream& outputTabs(std::ostream &os) const;
			//@}
			
			/*! @name Output the average of this to an out stream.
			
			\param os is the stream to sent the output to.
			\param prec is the precision to use in the output.*/
			//@{
			std::ostream& outputAverage(std::ostream &os,
										int prec) const;
			std::ostream& outputAverage(std::ostream &os) const;
			//@}
			
			/*! \brief Swap the contents of this and other.*/
			void swap(RangeCollectionHist &other);
			
			
		private:
			
			typedef cxsc::real 
									RangeCollectionHistPrivateValueType;
			typedef std::vector<RangeCollectionHistPrivateValueType> 
									RangeCollectionHistPrivateType;
			typedef RangeCollectionHistPrivateType::iterator 
									RangeCollectionHistPrivateTypeItr;
			typedef RangeCollectionHistPrivateType::const_iterator 
									RangeCollectionHistPrivateTypeConstItr;

            RangeCollectionHistPrivateType container;
			

			void scalarMult(const cxsc::real& mult);
			

    }; // end class RangeCollectionHist


} // end namespace subpavings

/*! A specialisation of std::swap for RangeCollectionHist types.*/
namespace std
{
	template <>
	void swap(subpavings::RangeCollectionHist & rh1, 
			subpavings::RangeCollectionHist & rh2); // throw ()
	
}



#endif
