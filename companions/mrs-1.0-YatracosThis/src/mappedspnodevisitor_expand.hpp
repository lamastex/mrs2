/*
* Copyright (C) 2010 Jenny Harlow
*
*/

#ifndef ___SPNODEVISITOREXPAND_HPP__
#define ___SPNODEVISITOREXPAND_HPP__


/*! \file
\brief declarations for MappedSPnodeVisitorExpander

An interface for a type that visits MappedSPnodes and expands them

*/

#include "spnodevisitor.hpp"
#include "mappedFobj.hpp"
#include "spnode.hpp"
#include <algorithm>


namespace subpavings {

    class MappedFobj;

    //template <typename T>
    //class MappedSPnode;


    class MappedSPnodeVisitorExpand : public SPnodeVisitor {

        private:

            MappedFobj& fobj;
            cxsc::real tolerance;

        public:
            MappedSPnodeVisitorExpand(MappedFobj& f, cxsc::real tol);

            virtual void visit(SPnode * spn);

            virtual cxsc::real tellMe(SPnode * spn);

				//gat41
				virtual bool priorityVisit(SPnode * spn, size_t critLeaves, std::vector<real>& eps);
				virtual bool priorityVisit(SPnode * spn, size_t critLeaves, gsl_rng * rgsl, std::vector<real>& eps);
				//virtual cxsc::real getSPArea(SPnode * mspn);
    };
    // end of MappedSPnodeVisitorExpand class

} // end namespace subpavings

#endif
